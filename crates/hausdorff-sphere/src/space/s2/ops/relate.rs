//! `s2::SpatialPredicate` — the OGC binary spatial predicates evaluated with
//! **great-circle edges on the unit sphere**, the spherical peer of
//! [`hausdorff_planar::space::r2::SpatialPredicate`] (which treats lon/lat as a
//! Cartesian plane). Predicates over the *same* WKB lon/lat geometries the planar
//! path uses, but with edges interpreted as geodesic arcs — so results are correct
//! across the antimeridian, at the poles, and for geometries spanning many
//! degrees, where the planar interpretation bows away from the true edges.
//!
//! **Scope: point-in-areal and areal-in-areal.** [`SpatialPredicate::eval`]
//! returns `Some(bool)` for
//!
//! - a point / multipoint against a polygon / multipolygon, backed by the
//!   proof-checked [`S2Polygon::contains_point`]; and
//! - a polygon / multipolygon against another, **including operands with holes**,
//!   decided through the spherical boolean engine in [`boolean_op`.
//!
//! Areal x areal used to return `None` for the containment predicates and for any
//! holed operand, deferring to the planar engine. That was the wrong place to
//! defer: areal geometry is precisely what spans the antimeridian and the poles,
//! so the deferred cases were the ones the planar reading is most likely to get
//! wrong. With a boolean engine available they are now answered here.
//!
//! What still returns `None` is anything **linear** — a linestring operand — and
//! the four dimension-dependent predicates (`Equals`, `Overlaps`, `Touches`,
//! `Crosses`) over *mixed* dimensions, where deciding them means duplicating the
//! planar DE-9IM engine's dimension bookkeeping for cases in which the spherical
//! and planar answers cannot differ. `None` remains the signal to fall back, so
//! wiring this in never regresses a case it does not cover.

use hausdorff_base::earth::GeoPoint;
use hausdorff_planar::space::r2::Geometry;

use super::super::cell::cover::{polygon_region, valid_lonlat};
use crate::space::s2::{Polygon as S2Polygon, S2Point};

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
    /// The geometries are topologically equal (`ST_Equals`).
    Equals,
    /// The geometries overlap at the same dimension (`ST_Overlaps`).
    Overlaps,
    /// The geometries touch only at their boundaries (`ST_Touches`).
    Touches,
    /// The geometries cross (`ST_Crosses`).
    Crosses,
}

impl SpatialPredicate {
    /// Evaluate `a <self> b` on the sphere, or `None` when this predicate over
    /// these geometry types is outside the v1 point-in-areal scope (the caller
    /// then falls back to the planar predicate).
    ///
    /// Supported: one side a point/multipoint, the other a polygon/multipolygon.
    /// Containment uses the shell/hole interior test [`S2Polygon::contains_point`];
    /// as there (and as in the planar path), a point exactly on a boundary is not
    /// guaranteed either way. `Covers`/`CoveredBy` coincide with
    /// `Contains`/`Within` for a point argument (they differ only on the boundary).
    #[must_use]
    pub fn eval(self, a: &Geometry, b: &Geometry) -> Option<bool> {
        if is_point_like(a) && is_areal(b) {
            return self.eval_points_in_areal(&s2_points(a), &s2_polygons(b)?);
        }
        if is_areal(a) && is_point_like(b) {
            // Reduce the areal×point case to the point×areal one with the
            // arguments and predicate transposed (`a P b` ⟺ `b Pᵀ a`).
            return self
                .transpose()
                .eval_points_in_areal(&s2_points(b), &s2_polygons(a)?);
        }
        if is_areal(a) && is_areal(b) {
            return self.eval_areal_areal(&s2_polygons(a)?, &s2_polygons(b)?);
        }
        None
    }

    /// Evaluate `areal <self> areal` on the sphere, for **every** predicate and for
    /// polygons **with holes**.
    ///
    /// This used to return `None` for the containment predicates and for any
    /// operand carrying a hole, deferring to the planar engine — which is exactly
    /// where the planar answer is most likely to be wrong, since it is the areal
    /// cases that span the antimeridian and the poles. The spherical boolean
    /// engine makes the whole family decidable, by reducing each predicate to
    /// areas of set combinations:
    ///
    /// | predicate | in terms of area |
    /// |---|---|
    /// | `Covers` | `|b \ a| = 0` |
    /// | `CoveredBy` | `|a \ b| = 0` |
    /// | `Equals` | `|a \ b| = 0` and `|b \ a| = 0` |
    /// | `Overlaps` | all three of `|a n b|`, `|a \ b|`, `|b \ a|` positive |
    /// | `Contains` / `Within` | `Covers` / `CoveredBy` plus a positive overlap |
    /// | `Intersects` | positive overlap, or boundary contact |
    /// | `Disjoint` | the negation of `Intersects` |
    /// | `Touches` | boundary contact with no interior overlap |
    ///
    /// **Two limits, both inherent to deciding a predicate from areas.**
    ///
    /// A measure-zero distinction is invisible here. OGC `Contains` differs from
    /// `Covers` only when `b` lies entirely within `a`'s *boundary*, which has no
    /// area, so the two agree for every operand with a non-degenerate interior.
    /// Likewise `Contains` and `Within` are reported through their `Covers`
    /// counterparts with an overlap check rather than by a separate interior test.
    ///
    /// Boundary contact without interior overlap is detected by testing each
    /// operand's vertices against the other's closure. That catches a shared
    /// vertex and a shared edge — the cases that actually occur in tiled data —
    /// and would miss two boundaries meeting tangentially at a point interior to
    /// an edge of both, a configuration with no area and no vertex to find. A
    /// transversal crossing is not affected: it produces real overlap and is
    /// caught by the area test.
    fn eval_areal_areal(self, a: &[S2Polygon], b: &[S2Polygon]) -> Option<bool> {
        use crate::space::s2::{Operation, boolean_op, total_area};

        // Scale the "is this zero" threshold to the operands rather than fixing it
        // absolutely: a sliver left by the arrangement is proportional to the
        // geometry it came from, and an absolute epsilon would be either blind on
        // continent-sized polygons or paranoid on city blocks.
        let (area_a, area_b) = (total_area(a), total_area(b));
        let scale = area_a.max(area_b).max(f64::MIN_POSITIVE);
        let eps = 1e-9 * scale;

        let overlap = total_area(&boolean_op(a, b, Operation::Intersection));
        let only_a = total_area(&boolean_op(a, b, Operation::Difference));
        let only_b = total_area(&boolean_op(b, a, Operation::Difference));

        let covers = only_b <= eps;
        let covered_by = only_a <= eps;
        let overlaps_interior = overlap > eps;
        let touching = !overlaps_interior && boundaries_touch(a, b);

        Some(match self {
            SpatialPredicate::Intersects => overlaps_interior || touching,
            SpatialPredicate::Disjoint => !(overlaps_interior || touching),
            SpatialPredicate::Covers => covers,
            SpatialPredicate::CoveredBy => covered_by,
            SpatialPredicate::Contains => covers && overlaps_interior,
            SpatialPredicate::Within => covered_by && overlaps_interior,
            SpatialPredicate::Equals => covers && covered_by,
            SpatialPredicate::Overlaps => overlaps_interior && only_a > eps && only_b > eps,
            SpatialPredicate::Touches => touching,
            // `Crosses` is defined for operands of *differing* dimension, so it is
            // never true for two areal geometries.
            SpatialPredicate::Crosses => false,
        })
    }

    /// The predicate with its arguments swapped: `Contains`⇄`Within`,
    /// `Covers`⇄`CoveredBy`. The rest — `Intersects`, `Disjoint`, `Equals`,
    /// `Overlaps`, `Touches`, `Crosses` — are symmetric, so each is its own
    /// transpose.
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
    ///
    /// Returns `None` for `Equals`, `Overlaps`, `Touches`, and `Crosses`. Those
    /// are not unsupported so much as *dimension-dependent*: `Overlaps` requires
    /// both operands to share a dimension, `Crosses` requires them to differ, and
    /// both are defined through the dimension of the intersection rather than
    /// through containment. Deciding them here would mean duplicating the planar
    /// DE-9IM engine's dimension bookkeeping for the one case where the spherical
    /// and planar answers cannot differ — a point either is or is not in a region,
    /// and no seam changes that. So they defer, as they always have.
    fn eval_points_in_areal(self, points: &[S2Point], areal: &[S2Polygon]) -> Option<bool> {
        // `Intersects`/`Disjoint`/`CoveredBy` include the boundary (a point on
        // an edge shares that point with the region); `Within` is the strict
        // interior (a boundary point is covered by, but not within, the region).
        let inside_closed = |p: S2Point| areal.iter().any(|poly| poly.contains_point_closed(p));
        let inside_open = |p: S2Point| areal.iter().any(|poly| poly.contains_point(p));
        Some(match self {
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
            SpatialPredicate::Equals
            | SpatialPredicate::Overlaps
            | SpatialPredicate::Touches
            | SpatialPredicate::Crosses => return None,
        })
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
/// Whether the two polygon sets share any boundary point, tested by vertex
/// containment in the other's closure.
///
/// Only consulted when the interiors do **not** overlap, so it is answering a
/// narrow question: do these two otherwise-separate regions touch? A shared
/// vertex or a shared edge puts a vertex of one inside the other's closure, which
/// is the case that arises in tiled or adjacency-built data.
fn boundaries_touch(a: &[S2Polygon], b: &[S2Polygon]) -> bool {
    let touches = |xs: &[S2Polygon], ys: &[S2Polygon]| {
        xs.iter().any(|x| {
            std::iter::once(x.shell())
                .chain(x.holes())
                .flat_map(|ring| ring.vertices().iter().copied())
                .any(|v| ys.iter().any(|y| y.contains_point_closed(v)))
        })
    };
    touches(a, b) || touches(b, a)
}

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
