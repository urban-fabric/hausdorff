//! Contour tracing: turning the sweep's surviving edges back into rings.
//!
//! The reference implementation walks the result events in sorted order and, at
//! a point shared by several edges, continues along whichever neighbouring event
//! happens to sit next in that order. Where the result touches itself — a pinch,
//! four or more edges meeting at one point — that order is arbitrary, so the
//! walk can leave along an edge belonging to a different part of the boundary.
//! The contour then dead-ends and gets closed with a chord that was never in the
//! result; two such rings overlap and their areas double-count (#44).
//!
//! So the choice is not left to the event order. Every edge is given a direction
//! with the result's interior on its left — from the ring winding it came in
//! with, read against the [`Operation`] — and at a shared point the continuation
//! is the first outgoing edge clockwise from the one just arrived along. That is
//! the usual planar face-traversal rule, and it follows from the labelling
//! rather than guessing at the geometry.
//!
//! **What this does not depend on.** An earlier version required the surviving
//! edges to have an even number meeting at every point, on the reasoning that a
//! closed boundary has no odd-degree vertex and the walk would otherwise have to
//! dead-end. That is true of a boundary that never touches itself, and false in
//! general: where a vertex of one operand lands exactly on a vertex of another,
//! removing the overlapping edge legitimately leaves odd parity, the walk does
//! dead-end there, and [`assemble`] recovers the rings regardless. A segment's
//! offset rectangle meeting a disc at its corner — what `buffer` produces — does
//! all three and returns an answer bit-identical to the same union computed
//! without the coincidence.
//!
//! So what is checked is the output instead: no two parts of the result may
//! share interior area. That is the invariant #44 actually broke, and it has no
//! false positives, because boundary contact between parts is permitted and
//! interior overlap never is.

use super::Operation;
use super::assemble::assemble;
use super::event::{Event, compare_events};
use crate::space::r2::{Coord, Polygon};

pub(super) fn connect_edges(events: &[Event], sorted: &[u32], op: Operation) -> Vec<Polygon> {
    let mut result: Vec<u32> = sorted
        .iter()
        .copied()
        .filter(|&e| {
            let ev = &events[e as usize];
            (ev.left && ev.in_result) || (!ev.left && events[ev.other.unwrap() as usize].in_result)
        })
        .collect();
    result.sort_by(|&a, &b| compare_events(events, a, b));

    let n = result.len();
    if n == 0 {
        return Vec::new();
    }
    let mut index_of = vec![u32::MAX; events.len()];
    for (i, &e) in result.iter().enumerate() {
        index_of[e as usize] = i as u32;
    }
    let other_pos: Vec<usize> = result
        .iter()
        .map(|&e| {
            let other = events[e as usize].other.unwrap();
            let p = index_of[other as usize];
            if p == u32::MAX {
                usize::MAX
            } else {
                p as usize
            }
        })
        .collect();
    if other_pos.contains(&usize::MAX) {
        return Vec::new();
    }

    let is_tail = orient(events, &result, &other_pos, op);
    let pairing = pair_at_shared_points(events, &result, &other_pos, &is_tail);

    let mut visited = vec![false; n];
    let mut rings: Vec<Vec<Coord>> = Vec::new();
    for start in 0..n {
        if visited[start] || !is_tail[start] {
            continue;
        }
        let mut ring: Vec<Coord> = Vec::new();
        let mut pos = start;
        loop {
            visited[pos] = true;
            let head = other_pos[pos];
            visited[head] = true;
            // Only the tail of each edge is recorded: the head is the next
            // edge's tail, the same point.
            ring.push(events[result[pos] as usize].point);
            let next = pairing[head];
            if next == usize::MAX || next == start || visited[next] {
                break;
            }
            pos = next;
        }
        if ring.len() >= 3 {
            rings.push(ring);
        }
    }

    let polygons = assemble(rings);
    debug_assert!(
        interiors_are_disjoint(&polygons),
        "{op:?}: two parts of the result share interior area, so tracing closed a \
         ring with a chord that was never in the result and the overlap is being \
         double counted"
    );
    polygons
}

/// Do no two parts of the result share interior area?
///
/// This is the invariant a multipolygon has to satisfy, and the one the tracing
/// bug in #44 broke: a contour closed with an invented chord produces two rings
/// that overlap, and their areas double-count. Boundary contact is permitted —
/// two parts may meet at a point or along a line — so the test is on the
/// interior/interior cell of the DE-9IM matrix rather than on intersection.
///
/// Checked here rather than on the sweep's input because the obvious
/// precondition does not hold. An even number of result edges at every point is
/// *not* necessary: where a vertex of one operand coincides exactly with a
/// vertex of another, overlapping-edge removal legitimately leaves odd parity,
/// the traversal dead-ends there, and `assemble` still recovers the correct
/// rings. A segment's offset rectangle meeting a disc at its corner does all
/// three and produces an answer identical to the same union without the
/// coincidence, so parity and dead-ends are both false positives for the thing
/// actually feared.
#[cfg(debug_assertions)]
fn interiors_are_disjoint(polygons: &[Polygon]) -> bool {
    use crate::space::r2::Geometry;
    use crate::space::r2::relate;

    polygons.iter().enumerate().all(|(index, polygon)| {
        polygons[index + 1..].iter().all(|other| {
            relate(
                &Geometry::Polygon(polygon.clone()),
                &Geometry::Polygon(other.clone()),
            )
            .interior_interior_dim()
                < 0
        })
    })
}

#[cfg(not(debug_assertions))]
fn interiors_are_disjoint(_polygons: &[Polygon]) -> bool {
    true
}

/// Give every result edge a direction with the result's interior on its left,
/// returning for each event position whether it is that directed edge's tail.
///
/// `in_out` is false when the edge's *own* polygon lies above it, and
/// `other_in_out` is false when the edge lies inside the *other* polygon. What
/// the result's interior is depends on the operation: a clipping edge kept by a
/// difference, for instance, bounds the result on the side its own polygon is
/// *not*. An edge with the interior above runs left-to-right, one with the
/// interior below runs right-to-left, which puts the interior on the left of
/// both.
fn orient(events: &[Event], result: &[u32], other_pos: &[usize], op: Operation) -> Vec<bool> {
    let mut is_tail = vec![false; result.len()];
    for i in 0..result.len() {
        let ev = &events[result[i] as usize];
        if !ev.left {
            continue;
        }
        // The ring runs left-to-right exactly when it leads away from the left
        // event, and a normalised ring keeps its own interior on its left — so
        // this says the edge's own polygon lies above it.
        let own_above = ev.ring_forward;
        let interior_above = match op {
            // The result is bounded on the same side its own polygon is; a
            // coincident pair (`SameTransition`) agrees on which side that is.
            Operation::Union | Operation::Intersection => own_above,
            // Subject edges keep their own side; clipping edges bound the result
            // on the far side from the clip. This also covers the coincident
            // opposed pair (`DifferentTransition`), where the surviving event may
            // belong to either operand.
            Operation::Difference => ev.subject == own_above,
            // Outside the other polygon the edge keeps its own side; inside it,
            // its own interior is covered by both operands and so is excluded,
            // leaving the result on the far side.
            Operation::Xor => ev.other_in_out == own_above,
        };
        is_tail[i] = interior_above;
        is_tail[other_pos[i]] = !interior_above;
    }
    is_tail
}

/// For every event position that a directed edge *arrives* at, the position the
/// contour leaves along. Events sharing a point are adjacent in `result`, so
/// each such run is one vertex of the planar graph.
fn pair_at_shared_points(
    events: &[Event],
    result: &[u32],
    other_pos: &[usize],
    is_tail: &[bool],
) -> Vec<usize> {
    let n = result.len();
    let mut pairing = vec![usize::MAX; n];
    let mut start = 0;
    while start < n {
        let point = events[result[start] as usize].point;
        let mut end = start;
        while end < n && events[result[end] as usize].point == point {
            end += 1;
        }
        pair_group(events, result, other_pos, is_tail, start, end, &mut pairing);
        start = end;
    }
    pairing
}

fn pair_group(
    events: &[Event],
    result: &[u32],
    other_pos: &[usize],
    is_tail: &[bool],
    start: usize,
    end: usize,
    pairing: &mut [usize],
) {
    // The direction away from the shared point along each incident edge. For a
    // tail that is the direction travelled; for a head it is the reverse of the
    // direction arrived along, which is the reference the clockwise rule needs.
    let bearing = |pos: usize| {
        let here = events[result[pos] as usize].point;
        let far = events[result[other_pos[pos]] as usize].point;
        (far.y - here.y).atan2(far.x - here.x)
    };

    let mut order: Vec<usize> = (start..end).collect();
    order.sort_by(|&a, &b| {
        bearing(a)
            .partial_cmp(&bearing(b))
            .unwrap_or(std::cmp::Ordering::Equal)
            .then(a.cmp(&b))
    });

    let m = order.len();
    let mut taken = vec![false; m];

    // The first outgoing edge clockwise from the one arrived along; with the
    // interior on the left of every edge this is the continuation that keeps it
    // there.
    for i in 0..m {
        if is_tail[order[i]] {
            continue;
        }
        for step in 1..=m {
            let k = (i + m - step) % m;
            if is_tail[order[k]] && !taken[k] {
                taken[k] = true;
                pairing[order[i]] = order[k];
                break;
            }
        }
    }

    // An arrival left unpaired here means the vertex has no free edge leaving it
    // — the sweep mislabelled which edges belong to the result, and the walk has
    // to stop. Continuing along some other incident edge instead was tried: it
    // recovers a little area on the worst band pair and loses more on three
    // others, which is what picking arbitrarily always looks like.
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::space::r2::LineString;

    fn square(origin: f64) -> Polygon {
        Polygon::new(
            LineString::from(vec![
                (origin, origin),
                (origin + 2.0, origin),
                (origin + 2.0, origin + 2.0),
                (origin, origin + 2.0),
                (origin, origin),
            ]),
            Vec::new(),
        )
    }

    #[test]
    fn the_output_invariant_rejects_overlapping_parts() {
        // The #44 symptom: a contour closed with an invented chord leaves two
        // rings whose areas double-count. Interior overlap is what that looks
        // like from the outside, so the check has to reject it.
        assert!(!interiors_are_disjoint(&[square(0.0), square(1.0)]));
    }

    /// A coincident vertex must not change the answer.
    ///
    /// A segment's offset rectangle plus a disc at each end, with a disc vertex
    /// landing exactly on the rectangle's corner — what `buffer` produced before
    /// its discs were rotated. This is the configuration that used to abort
    /// tracing: it leaves odd parity at the shared point and the walk dead-ends
    /// there. The answer must still match the same union computed with the discs
    /// rotated so nothing is coincident.
    #[test]
    fn a_coincident_vertex_does_not_change_the_union() {
        use super::super::unary_union;
        use crate::space::r2::{Area as _, MultiPolygon};

        let radius = 0.75;
        let rectangle = Polygon::new(
            LineString::from(vec![
                (0.0, -radius),
                (3.0, -radius),
                (3.0, radius),
                (0.0, radius),
                (0.0, -radius),
            ]),
            Vec::new(),
        );
        let disc = |centre: f64, segments: u32, offset: f64| {
            let mut ring: Vec<(f64, f64)> = (0..segments)
                .map(|step| {
                    let angle =
                        offset + std::f64::consts::TAU * f64::from(step) / f64::from(segments);
                    (centre + radius * angle.cos(), radius * angle.sin())
                })
                .collect();
            ring.push(ring[0]);
            Polygon::new(LineString::from(ring), Vec::new())
        };
        // Four and 1024 segments both place a vertex at exactly 90 degrees,
        // which is exactly the rectangle's corner height.
        for segments in [4u32, 16, 1024] {
            let half_step = std::f64::consts::TAU / (2.0 * f64::from(segments));
            let coincident: Vec<Polygon> = vec![
                rectangle.clone(),
                disc(0.0, segments, 0.0),
                disc(3.0, segments, 0.0),
            ];
            let rotated: Vec<Polygon> = vec![
                rectangle.clone(),
                disc(0.0, segments, half_step),
                disc(3.0, segments, half_step),
            ];
            let merged: MultiPolygon = unary_union(coincident.iter());
            let reference: MultiPolygon = unary_union(rotated.iter());
            assert_eq!(
                merged.0.len(),
                reference.0.len(),
                "{segments} segments: part count moved"
            );
            assert!(
                (merged.unsigned_area() - reference.unsigned_area()).abs() < 1e-9,
                "{segments} segments: area {} against {}",
                merged.unsigned_area(),
                reference.unsigned_area()
            );
        }
    }

    #[test]
    fn the_output_invariant_permits_boundary_contact() {
        // Two parts of a multipolygon may meet at a point or along a line.
        // Rejecting that would make the check fail on correct results — which
        // is the mistake the parity precondition made.
        assert!(interiors_are_disjoint(&[square(0.0), square(2.0)]));
        assert!(interiors_are_disjoint(&[square(0.0), square(10.0)]));
        assert!(interiors_are_disjoint(&[square(0.0)]));
        assert!(interiors_are_disjoint(&[]));
    }
}
