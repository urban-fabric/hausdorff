//! Exact planar predicates over integer coordinates.
//!
//! This is the fixed-grid backend pointed at by `riemann`'s exact-arithmetic
//! proofs: `ExactArith/RobustOrientation.triDet` for orientation and
//! `ExactArith/InCircle.inCircle` for Delaunay in-circle tests. The formulas are
//! the same determinants used by the floating-point geometry kernels, but they
//! run over checked `i128` arithmetic. That gives exact signs for coordinates
//! that callers have already placed on an integer grid, with overflow reported
//! instead of silently wrapping.

/// Result type for checked exact predicates.
pub type ExactResult<T> = Result<T, ExactPredicateError>;

/// Failure modes for exact predicate evaluation.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ExactPredicateError {
    /// An intermediate value exceeded `i128`.
    Overflow,
    /// The defining triangle is collinear, so no triangle interior exists.
    DegenerateTriangle,
    /// The defining triangle is collinear, so no unique circumcircle exists.
    DegenerateCircumcircle,
}

/// A planar coordinate on an exact integer grid.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct ExactCoord {
    pub x: i128,
    pub y: i128,
}

impl ExactCoord {
    /// Construct an exact coordinate.
    #[must_use]
    pub const fn new(x: i128, y: i128) -> Self {
        Self { x, y }
    }
}

impl From<(i128, i128)> for ExactCoord {
    fn from((x, y): (i128, i128)) -> Self {
        Self::new(x, y)
    }
}

impl From<(i64, i64)> for ExactCoord {
    fn from((x, y): (i64, i64)) -> Self {
        Self::new(i128::from(x), i128::from(y))
    }
}

/// Exact orientation sign for a point triple.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum OrientationSign {
    Clockwise,
    Collinear,
    CounterClockwise,
}

impl OrientationSign {
    #[must_use]
    pub const fn is_collinear(self) -> bool {
        matches!(self, Self::Collinear)
    }

    #[must_use]
    pub const fn is_clockwise(self) -> bool {
        matches!(self, Self::Clockwise)
    }

    #[must_use]
    pub const fn is_counter_clockwise(self) -> bool {
        matches!(self, Self::CounterClockwise)
    }

    fn from_determinant(det: i128) -> Self {
        if det > 0 {
            Self::CounterClockwise
        } else if det < 0 {
            Self::Clockwise
        } else {
            Self::Collinear
        }
    }
}

#[cfg(euler_extract)]
macro_rules! euler_targets {
    ($($target:ident),+ $(,)?) => {
        #[allow(dead_code, non_snake_case)]
        mod __euler_targets {
            $(pub fn $target() {})+
        }
    };
}

#[cfg(euler_extract)]
euler_targets!(
    OrientationSign__from_determinant,
    relative,
    det2,
    orientation_determinant,
    orient,
);

/// The position of a query point relative to an oriented triangle's
/// circumcircle, normalized so the result is independent of triangle winding.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum InCirclePosition {
    Inside,
    OnCircle,
    Outside,
}

/// Exact segment-segment intersection classification.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ExactSegmentIntersection {
    Disjoint,
    /// Segments share exactly one endpoint or one endpoint lies on the other
    /// segment.
    Touch,
    /// Non-collinear segments cross in their relative interiors.
    ProperCross,
    /// Collinear segments share a positive-length subsegment.
    CollinearOverlap,
}

impl ExactSegmentIntersection {
    #[must_use]
    pub const fn intersects(self) -> bool {
        !matches!(self, Self::Disjoint)
    }
}

/// Exact point-in-triangle classification.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ExactTrianglePosition {
    Inside,
    OnBoundary,
    Outside,
}

impl ExactTrianglePosition {
    #[must_use]
    pub const fn is_inside_or_boundary(self) -> bool {
        matches!(self, Self::Inside | Self::OnBoundary)
    }
}

#[inline]
fn checked_add(a: i128, b: i128) -> ExactResult<i128> {
    a.checked_add(b).ok_or(ExactPredicateError::Overflow)
}

#[inline]
fn checked_sub(a: i128, b: i128) -> ExactResult<i128> {
    a.checked_sub(b).ok_or(ExactPredicateError::Overflow)
}

#[inline]
fn checked_mul(a: i128, b: i128) -> ExactResult<i128> {
    a.checked_mul(b).ok_or(ExactPredicateError::Overflow)
}

#[inline]
fn checked_square(a: i128) -> ExactResult<i128> {
    checked_mul(a, a)
}

#[inline]
fn det2(ax: i128, ay: i128, bx: i128, by: i128) -> ExactResult<i128> {
    checked_sub(checked_mul(ax, by)?, checked_mul(ay, bx)?)
}

#[inline]
fn relative(p: ExactCoord, origin: ExactCoord) -> ExactResult<(i128, i128)> {
    Ok((checked_sub(p.x, origin.x)?, checked_sub(p.y, origin.y)?))
}

#[inline]
fn sq_norm(x: i128, y: i128) -> ExactResult<i128> {
    checked_add(checked_square(x)?, checked_square(y)?)
}

/// Exact orientation determinant `triDet(a,b,c)`.
///
/// Positive means `a -> b -> c` is counter-clockwise, negative means clockwise,
/// and zero means collinear. This is `riemann`
/// `ExactArith/RobustOrientation.triDet`, evaluated with checked integer
/// arithmetic.
pub fn orientation_determinant(a: ExactCoord, b: ExactCoord, c: ExactCoord) -> ExactResult<i128> {
    let (ax, ay) = relative(a, c)?;
    let (bx, by) = relative(b, c)?;
    det2(ax, ay, bx, by)
}

/// Exact orientation sign for a point triple.
pub fn orient(a: ExactCoord, b: ExactCoord, c: ExactCoord) -> ExactResult<OrientationSign> {
    orientation_determinant(a, b, c).map(OrientationSign::from_determinant)
}

#[inline]
fn strictly_opposite(a: OrientationSign, b: OrientationSign) -> bool {
    matches!(
        (a, b),
        (
            OrientationSign::Clockwise,
            OrientationSign::CounterClockwise
        ) | (
            OrientationSign::CounterClockwise,
            OrientationSign::Clockwise
        )
    )
}

#[inline]
fn between_closed(a: i128, b: i128, p: i128) -> bool {
    a.min(b) <= p && p <= a.max(b)
}

/// True iff `p` lies on the closed segment `a-b`.
pub fn point_on_segment(a: ExactCoord, b: ExactCoord, p: ExactCoord) -> ExactResult<bool> {
    Ok(orient(a, b, p)?.is_collinear()
        && between_closed(a.x, b.x, p.x)
        && between_closed(a.y, b.y, p.y))
}

/// True iff the two segments cross in their relative interiors.
///
/// This is the double-straddle test proved by `riemann`
/// `ExactArith/SegmentCross.segments_properly_cross`.
pub fn segments_properly_cross(
    a: ExactCoord,
    b: ExactCoord,
    c: ExactCoord,
    d: ExactCoord,
) -> ExactResult<bool> {
    let ab_c = orient(a, b, c)?;
    let ab_d = orient(a, b, d)?;
    let cd_a = orient(c, d, a)?;
    let cd_b = orient(c, d, b)?;
    Ok(strictly_opposite(ab_c, ab_d) && strictly_opposite(cd_a, cd_b))
}

fn bbox_intersects(a: ExactCoord, b: ExactCoord, c: ExactCoord, d: ExactCoord) -> bool {
    a.x.min(b.x) <= c.x.max(d.x)
        && c.x.min(d.x) <= a.x.max(b.x)
        && a.y.min(b.y) <= c.y.max(d.y)
        && c.y.min(d.y) <= a.y.max(b.y)
}

fn collinear_overlap_kind(
    a: ExactCoord,
    b: ExactCoord,
    c: ExactCoord,
    d: ExactCoord,
) -> ExactSegmentIntersection {
    if !bbox_intersects(a, b, c, d) {
        return ExactSegmentIntersection::Disjoint;
    }
    let x_positive = a.x.min(b.x).max(c.x.min(d.x)) < a.x.max(b.x).min(c.x.max(d.x));
    let y_positive = a.y.min(b.y).max(c.y.min(d.y)) < a.y.max(b.y).min(c.y.max(d.y));
    if x_positive || y_positive {
        ExactSegmentIntersection::CollinearOverlap
    } else {
        ExactSegmentIntersection::Touch
    }
}

/// Classify the intersection of two closed line segments.
///
/// Proper crossings are the exact double-straddle predicate from `riemann`
/// `ExactArith/SegmentCross`; endpoint and collinear cases add the closed AABB
/// checks used by `ExactArith/AABB`.
pub fn segment_intersection(
    a: ExactCoord,
    b: ExactCoord,
    c: ExactCoord,
    d: ExactCoord,
) -> ExactResult<ExactSegmentIntersection> {
    let ab_c = orient(a, b, c)?;
    let ab_d = orient(a, b, d)?;
    let cd_a = orient(c, d, a)?;
    let cd_b = orient(c, d, b)?;

    if strictly_opposite(ab_c, ab_d) && strictly_opposite(cd_a, cd_b) {
        return Ok(ExactSegmentIntersection::ProperCross);
    }

    if ab_c.is_collinear() && ab_d.is_collinear() && cd_a.is_collinear() && cd_b.is_collinear() {
        return Ok(collinear_overlap_kind(a, b, c, d));
    }

    if (ab_c.is_collinear() && point_on_segment(a, b, c)?)
        || (ab_d.is_collinear() && point_on_segment(a, b, d)?)
        || (cd_a.is_collinear() && point_on_segment(c, d, a)?)
        || (cd_b.is_collinear() && point_on_segment(c, d, b)?)
    {
        return Ok(ExactSegmentIntersection::Touch);
    }

    Ok(ExactSegmentIntersection::Disjoint)
}

/// Classify `p` against the closed triangle `a,b,c`.
///
/// This is the same-sign edge-orientation test backed by `riemann`
/// `ExactArith/OrientationLaws.point_in_triangle_edges`, normalized for either
/// triangle winding.
pub fn point_in_triangle(
    a: ExactCoord,
    b: ExactCoord,
    c: ExactCoord,
    p: ExactCoord,
) -> ExactResult<ExactTrianglePosition> {
    let whole = orient(a, b, c)?;
    if whole.is_collinear() {
        return Err(ExactPredicateError::DegenerateTriangle);
    }

    let edges = [orient(a, b, p)?, orient(b, c, p)?, orient(c, a, p)?];
    let outside = match whole {
        OrientationSign::CounterClockwise => edges.iter().any(|s| s.is_clockwise()),
        OrientationSign::Clockwise => edges.iter().any(|s| s.is_counter_clockwise()),
        OrientationSign::Collinear => unreachable!("handled above"),
    };
    if outside {
        return Ok(ExactTrianglePosition::Outside);
    }
    if edges.iter().any(|s| s.is_collinear()) {
        Ok(ExactTrianglePosition::OnBoundary)
    } else {
        Ok(ExactTrianglePosition::Inside)
    }
}

/// True iff `p` is inside the closed triangle `a,b,c`.
pub fn triangle_contains_point(
    a: ExactCoord,
    b: ExactCoord,
    c: ExactCoord,
    p: ExactCoord,
) -> ExactResult<bool> {
    point_in_triangle(a, b, c, p).map(ExactTrianglePosition::is_inside_or_boundary)
}

/// Exact Delaunay in-circle determinant.
///
/// This is the expanded `riemann` `ExactArith/InCircle.inCircle` determinant:
/// the rows are `(p_i - d, ||p_i - d||^2)` for `p_i` in `{a,b,c}`. Its raw sign
/// depends on the winding of `a,b,c`; use [`in_circle_position`] for the
/// orientation-normalized predicate.
pub fn in_circle_determinant(
    a: ExactCoord,
    b: ExactCoord,
    c: ExactCoord,
    d: ExactCoord,
) -> ExactResult<i128> {
    let (ax, ay) = relative(a, d)?;
    let (bx, by) = relative(b, d)?;
    let (cx, cy) = relative(c, d)?;

    let a2 = sq_norm(ax, ay)?;
    let b2 = sq_norm(bx, by)?;
    let c2 = sq_norm(cx, cy)?;

    let term1 = checked_mul(ax, checked_sub(checked_mul(by, c2)?, checked_mul(b2, cy)?)?)?;
    let term2 = checked_mul(ay, checked_sub(checked_mul(bx, c2)?, checked_mul(b2, cx)?)?)?;
    let term3 = checked_mul(a2, det2(bx, by, cx, cy)?)?;

    checked_add(checked_sub(term1, term2)?, term3)
}

/// Classify `d` against the circumcircle through non-collinear `a,b,c`.
///
/// The raw determinant is positive for inside when `a,b,c` are
/// counter-clockwise and negative for inside when they are clockwise; this
/// function normalizes by the exact orientation sign.
pub fn in_circle_position(
    a: ExactCoord,
    b: ExactCoord,
    c: ExactCoord,
    d: ExactCoord,
) -> ExactResult<InCirclePosition> {
    let orientation = orient(a, b, c)?;
    if orientation.is_collinear() {
        return Err(ExactPredicateError::DegenerateCircumcircle);
    }

    let det = in_circle_determinant(a, b, c, d)?;
    if det == 0 {
        return Ok(InCirclePosition::OnCircle);
    }

    match (orientation, det > 0) {
        (OrientationSign::CounterClockwise, true) | (OrientationSign::Clockwise, false) => {
            Ok(InCirclePosition::Inside)
        }
        (OrientationSign::CounterClockwise, false) | (OrientationSign::Clockwise, true) => {
            Ok(InCirclePosition::Outside)
        }
        (OrientationSign::Collinear, _) => unreachable!("handled above"),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn near_collinear_scaled_point_keeps_exact_sign() {
        let scale = 1_000_000_000_000_000_000_i128;
        let a = ExactCoord::new(0, 0);
        let b = ExactCoord::new(3 * scale, 3 * scale);
        let c = ExactCoord::new(scale, scale + 1);

        assert_eq!(orientation_determinant(a, b, c).unwrap(), 3 * scale);
        assert_eq!(orient(a, b, c).unwrap(), OrientationSign::CounterClockwise);
    }

    #[test]
    fn incircle_classifies_inside_boundary_outside() {
        let a = ExactCoord::new(0, 0);
        let b = ExactCoord::new(2, 0);
        let c = ExactCoord::new(0, 2);

        assert_eq!(
            in_circle_position(a, b, c, ExactCoord::new(1, 1)).unwrap(),
            InCirclePosition::Inside
        );
        assert_eq!(
            in_circle_position(a, b, c, ExactCoord::new(2, 2)).unwrap(),
            InCirclePosition::OnCircle
        );
        assert_eq!(
            in_circle_position(a, b, c, ExactCoord::new(3, 3)).unwrap(),
            InCirclePosition::Outside
        );
    }

    #[test]
    fn incircle_position_is_independent_of_triangle_winding() {
        let a = ExactCoord::new(0, 0);
        let b = ExactCoord::new(2, 0);
        let c = ExactCoord::new(0, 2);
        let query = ExactCoord::new(1, 1);

        assert_eq!(
            in_circle_position(a, b, c, query).unwrap(),
            in_circle_position(c, b, a, query).unwrap()
        );
    }

    #[test]
    fn segments_classify_cross_touch_overlap_and_disjoint() {
        assert_eq!(
            segment_intersection(
                ExactCoord::new(0, 0),
                ExactCoord::new(4, 4),
                ExactCoord::new(0, 4),
                ExactCoord::new(4, 0),
            )
            .unwrap(),
            ExactSegmentIntersection::ProperCross
        );
        assert!(
            segments_properly_cross(
                ExactCoord::new(0, 0),
                ExactCoord::new(4, 4),
                ExactCoord::new(0, 4),
                ExactCoord::new(4, 0),
            )
            .unwrap()
        );
        assert_eq!(
            segment_intersection(
                ExactCoord::new(0, 0),
                ExactCoord::new(2, 0),
                ExactCoord::new(2, 0),
                ExactCoord::new(2, 3),
            )
            .unwrap(),
            ExactSegmentIntersection::Touch
        );
        assert_eq!(
            segment_intersection(
                ExactCoord::new(0, 0),
                ExactCoord::new(4, 0),
                ExactCoord::new(2, 0),
                ExactCoord::new(6, 0),
            )
            .unwrap(),
            ExactSegmentIntersection::CollinearOverlap
        );
        assert_eq!(
            segment_intersection(
                ExactCoord::new(0, 0),
                ExactCoord::new(1, 0),
                ExactCoord::new(2, 0),
                ExactCoord::new(3, 0),
            )
            .unwrap(),
            ExactSegmentIntersection::Disjoint
        );
    }

    #[test]
    fn point_in_triangle_classifies_boundary_and_winding() {
        let a = ExactCoord::new(0, 0);
        let b = ExactCoord::new(4, 0);
        let c = ExactCoord::new(0, 4);

        assert_eq!(
            point_in_triangle(a, b, c, ExactCoord::new(1, 1)).unwrap(),
            ExactTrianglePosition::Inside
        );
        assert_eq!(
            point_in_triangle(a, b, c, ExactCoord::new(2, 0)).unwrap(),
            ExactTrianglePosition::OnBoundary
        );
        assert_eq!(
            point_in_triangle(a, b, c, ExactCoord::new(3, 3)).unwrap(),
            ExactTrianglePosition::Outside
        );
        assert_eq!(
            point_in_triangle(c, b, a, ExactCoord::new(1, 1)).unwrap(),
            ExactTrianglePosition::Inside
        );
    }
}
