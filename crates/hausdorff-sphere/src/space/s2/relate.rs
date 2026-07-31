//! `s2::SpatialPredicate` — the OGC binary spatial predicates evaluated with
//! **great-circle edges on the unit sphere**, the spherical peer of
//! [`hausdorff_planar::space::r2::SpatialPredicate`] (which treats lon/lat as a
//! Cartesian plane). Predicates over the *same* WKB lon/lat geometries the planar
//! path uses, but with edges interpreted as geodesic arcs — so results are correct
//! across the antimeridian, at the poles, and for geometries spanning many
//! degrees, where the planar interpretation bows away from the true edges.
//!
//! **Scope (v1): point-in-areal.** [`SpatialPredicate::eval`] returns `Some(bool)`
//! for a point / multipoint tested against a polygon / multipolygon (the dominant
//! spatial-join case — points against regions — and exactly where the planar
//! lon/lat answer diverges from the spherical one), backed by the proof-checked
//! [`Polygon::contains_point`]. For every other geometry-type combination it
//! returns `None`, the signal for a caller to fall back to the planar predicate —
//! so wiring it in never regresses a case it does not yet cover. Areal×areal and
//! linear predicates are the documented follow-on.

use hausdorff_base::earth::GeoPoint;
use hausdorff_planar::space::r2::Geometry;

use super::cover::{polygon_region, valid_lonlat};
use super::{Polygon as S2Polygon, S2Point};

/// The OGC binary spatial predicates, evaluated on the sphere. Mirrors
/// [`hausdorff_planar::space::r2::SpatialPredicate`]; carry one as data into a join
/// and call [`eval`](SpatialPredicate::eval).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SpatialPredicate {
    /// The geometries share at least one point (`ST_Intersects`).
    Intersects,
    /// `a` contains `b` (`ST_Contains`).
    Contains,
    /// `a` is within `b` (`ST_Within`).
    Within,
    /// `a` covers `b` (`ST_Covers`).
    Covers,
    /// `a` is covered by `b` (`ST_CoveredBy`).
    CoveredBy,
    /// The geometries have no point in common (`ST_Disjoint`).
    Disjoint,
}

impl SpatialPredicate {
    /// Evaluate `a <self> b` on the sphere, or `None` when this predicate over
    /// these geometry types is outside the v1 point-in-areal scope (the caller
    /// then falls back to the planar predicate).
    ///
    /// Supported: one side a point/multipoint, the other a polygon/multipolygon.
    /// Containment uses the shell/hole interior test [`Polygon::contains_point`];
    /// as there (and as in the planar path), a point exactly on a boundary is not
    /// guaranteed either way. `Covers`/`CoveredBy` coincide with
    /// `Contains`/`Within` for a point argument (they differ only on the boundary).
    #[must_use]
    pub fn eval(self, a: &Geometry, b: &Geometry) -> Option<bool> {
        if is_point_like(a) && is_areal(b) {
            return Some(self.eval_points_in_areal(&s2_points(a), &s2_polygons(b)?));
        }
        if is_areal(a) && is_point_like(b) {
            // Reduce the areal×point case to the point×areal one with the
            // arguments and predicate transposed (`a P b` ⟺ `b Pᵀ a`).
            return Some(
                self.transpose()
                    .eval_points_in_areal(&s2_points(b), &s2_polygons(a)?),
            );
        }
        if is_areal(a) && is_areal(b) {
            return self.eval_areal_areal(&s2_polygons(a)?, &s2_polygons(b)?);
        }
        None
    }

    /// Evaluate `areal <self> areal`, exactly on the sphere for the tractable case:
    /// **Intersects / Disjoint over hole-free polygons**, via shell-loop
    /// intersection ([`Loop::intersects`] — boundary crossing or nesting). Returns
    /// `None` (defer to the planar predicate) for containment predicates and for any
    /// polygon carrying holes, which need a full spherical boolean / DE-9IM engine
    /// this does not yet provide.
    fn eval_areal_areal(self, a: &[S2Polygon], b: &[S2Polygon]) -> Option<bool> {
        if !matches!(
            self,
            SpatialPredicate::Intersects | SpatialPredicate::Disjoint
        ) {
            return None;
        }
        if a.iter().chain(b.iter()).any(|p| !p.holes().is_empty()) {
            return None;
        }
        let hit = a
            .iter()
            .any(|x| b.iter().any(|y| x.shell().intersects(y.shell())));
        Some(match self {
            SpatialPredicate::Disjoint => !hit,
            _ => hit,
        })
    }

    /// The predicate with its arguments swapped: `Contains`⇄`Within`,
    /// `Covers`⇄`CoveredBy`; the rest are their own transpose.
    #[must_use]
    fn transpose(self) -> Self {
        match self {
            SpatialPredicate::Contains => SpatialPredicate::Within,
            SpatialPredicate::Within => SpatialPredicate::Contains,
            SpatialPredicate::Covers => SpatialPredicate::CoveredBy,
            SpatialPredicate::CoveredBy => SpatialPredicate::Covers,
            other => other,
        }
    }

    /// Evaluate `points <self> areal`: the predicate with the point side as `a`.
    fn eval_points_in_areal(self, points: &[S2Point], areal: &[S2Polygon]) -> bool {
        // `Intersects`/`Disjoint`/`CoveredBy` include the boundary (a point on
        // an edge shares that point with the region); `Within` is the strict
        // interior (a boundary point is covered by, but not within, the region).
        let inside_closed = |p: S2Point| areal.iter().any(|poly| poly.contains_point_closed(p));
        let inside_open = |p: S2Point| areal.iter().any(|poly| poly.contains_point(p));
        match self {
            // Points and an areal share a point iff some point is on/in the region.
            SpatialPredicate::Intersects => points.iter().any(|&p| inside_closed(p)),
            SpatialPredicate::Disjoint => !points.iter().any(|&p| inside_closed(p)),
            // Covered by the areal iff every point is inside or on its boundary.
            SpatialPredicate::CoveredBy => {
                !points.is_empty() && points.iter().all(|&p| inside_closed(p))
            }
            // Within iff every point is in the strict interior.
            SpatialPredicate::Within => {
                !points.is_empty() && points.iter().all(|&p| inside_open(p))
            }
            // A point set (dimension 0) can never contain/cover an areal region.
            SpatialPredicate::Contains | SpatialPredicate::Covers => false,
        }
    }
}

/// The lon/lat points of a point or multipoint geometry as [`S2Point`]s
/// (invalid coordinates dropped); empty for a non-point geometry.
fn s2_points(geom: &Geometry) -> Vec<S2Point> {
    let point = |x: f64, y: f64| valid_lonlat(x, y).then(|| point_at(x, y));
    match geom {
        Geometry::Point(p) => point(p.x(), p.y()).into_iter().collect(),
        Geometry::MultiPoint(mp) => mp.iter().filter_map(|p| point(p.x(), p.y())).collect(),
        _ => Vec::new(),
    }
}

fn point_at(lon: f64, lat: f64) -> S2Point {
    S2Point::from_geo_point(GeoPoint::from_lon_lat(lon, lat))
}

/// The spherical polygons of an areal geometry; `None` when a ring has an invalid
/// coordinate or the geometry is not areal.
fn s2_polygons(geom: &Geometry) -> Option<Vec<S2Polygon>> {
    match geom {
        Geometry::Polygon(poly) => polygon_region(poly).map(|p| vec![p]),
        Geometry::MultiPolygon(mp) => mp.iter().map(polygon_region).collect(),
        Geometry::Rect(r) => polygon_region(&r.to_polygon()).map(|p| vec![p]),
        Geometry::Triangle(t) => polygon_region(&t.to_polygon()).map(|p| vec![p]),
        _ => None,
    }
}

fn is_point_like(geom: &Geometry) -> bool {
    matches!(geom, Geometry::Point(_) | Geometry::MultiPoint(_))
}

fn is_areal(geom: &Geometry) -> bool {
    matches!(
        geom,
        Geometry::Polygon(_)
            | Geometry::MultiPolygon(_)
            | Geometry::Rect(_)
            | Geometry::Triangle(_)
    )
}

#[cfg(test)]
mod tests;
