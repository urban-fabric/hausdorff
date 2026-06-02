//! The OGC binary spatial predicates as data — [`SpatialPredicate`] and its
//! evaluator, which reads the DE-9IM [`relate`](super::relate) matrix and
//! short-circuits the cheap cases.

use super::super::{Contains, Geometry, Intersects};
use super::{dimension, relate};

/// The OGC binary spatial predicates, evaluated natively over [`Geometry`].
/// Carry one as data (e.g. into a spatial join) and call
/// [`eval`](SpatialPredicate::eval); it reads the DE-9IM [`relate`] matrix,
/// short-circuiting the cheap cases.
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
    /// The geometries touch only at their boundaries (`ST_Touches`).
    Touches,
    /// The geometries cross (`ST_Crosses`).
    Crosses,
    /// The geometries are topologically equal (`ST_Equals`).
    Equals,
    /// The geometries overlap at the same dimension (`ST_Overlaps`).
    Overlaps,
}

impl SpatialPredicate {
    /// Evaluate `a <self> b`. `Intersects`/`Disjoint` take the direct
    /// [`Intersects`] path; a point tested against an areal geometry takes the
    /// interior ray-cast fast path (identical to the matrix's contains/within
    /// but without building it — far cheaper for a many-vertex boundary);
    /// everything else reads the DE-9IM [`relate`] matrix.
    #[must_use]
    pub fn eval(self, a: &Geometry, b: &Geometry) -> bool {
        match self {
            SpatialPredicate::Intersects => return a.intersects(b),
            SpatialPredicate::Disjoint => return !a.intersects(b),
            SpatialPredicate::Contains => {
                if let Some(hit) = areal_contains_point(a, b) {
                    return hit;
                }
            }
            SpatialPredicate::Within => {
                if let Some(hit) = areal_contains_point(b, a) {
                    return hit;
                }
            }
            _ => {}
        }
        let matrix = relate(a, b);
        let (da, db) = (dimension(a), dimension(b));
        match self {
            SpatialPredicate::Contains => matrix.is_contains(),
            SpatialPredicate::Within => matrix.is_within(),
            SpatialPredicate::Covers => matrix.is_covers(),
            SpatialPredicate::CoveredBy => matrix.is_coveredby(),
            SpatialPredicate::Touches => matrix.is_touches(),
            SpatialPredicate::Crosses => matrix.is_crosses(da, db),
            SpatialPredicate::Equals => matrix.is_equal_topo(),
            SpatialPredicate::Overlaps => matrix.is_overlaps(da, db),
            SpatialPredicate::Intersects | SpatialPredicate::Disjoint => unreachable!(),
        }
    }

    /// The predicate with its arguments swapped: `a P b` ⟺ `b P.transpose() a`.
    /// `Contains`/`Within` and `Covers`/`CoveredBy` swap; the rest are their own
    /// transpose. Useful when a join indexes one side but evaluates
    /// `predicate(left, right)`.
    #[must_use]
    pub fn transpose(self) -> Self {
        match self {
            SpatialPredicate::Contains => SpatialPredicate::Within,
            SpatialPredicate::Within => SpatialPredicate::Contains,
            SpatialPredicate::Covers => SpatialPredicate::CoveredBy,
            SpatialPredicate::CoveredBy => SpatialPredicate::Covers,
            other => other,
        }
    }
}

/// Interior membership of a point in an areal geometry — the fast path for
/// `Contains`/`Within` of a point. `None` when `area` is not areal or `point` is
/// not a [`Geometry::Point`], so the caller falls back to the full matrix.
fn areal_contains_point(area: &Geometry, point: &Geometry) -> Option<bool> {
    let Geometry::Point(p) = point else {
        return None;
    };
    match area {
        Geometry::Polygon(poly) => Some(poly.contains(p)),
        Geometry::MultiPolygon(mp) => Some(mp.contains(p)),
        Geometry::Rect(r) => Some(r.contains(p)),
        _ => None,
    }
}
