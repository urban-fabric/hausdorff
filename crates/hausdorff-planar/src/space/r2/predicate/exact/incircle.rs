use super::*;

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
