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
//! **What this depends on.** The rule needs the surviving edges to form a closed
//! boundary. Where the sweep mislabels which edges belong to the result it hands
//! over a set with an *odd-degree* vertex — a point where an odd number of result
//! edges meet, which no closed boundary has — and then the walk must dead-end
//! there whatever it does. Six of the thirty-two real band pairs used to arrive
//! that way; see the ordering and re-placement notes in [`super::sweep`] for what
//! that was. Odd parity is the thing to measure if this regresses: count result
//! edges per point, rather than comparing areas, which only says something is
//! wrong and not where.

use super::super::{Coord, Polygon};
use super::Operation;
use super::assemble::assemble;
use super::event::{Event, compare_events};

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

    debug_assert!(
        closed_boundary(events, &result),
        "{op:?}: the sweep left a point with an odd number of result edges at it, \
         which no closed boundary has — the region labelling is wrong, and tracing \
         cannot recover it"
    );

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

    assemble(rings)
}

/// Does an even number of result edges meet at every point? Anything else is a
/// dangling boundary, and the tracer below would have to invent a chord to close
/// it — so this is checked rather than worked around. It is the cheap form of the
/// question areas only answer indirectly: a residual says a result is wrong, this
/// says where. Events sharing a point are adjacent in `result`, so it is one pass.
fn closed_boundary(events: &[Event], result: &[u32]) -> bool {
    let mut start = 0;
    while start < result.len() {
        let point = events[result[start] as usize].point;
        let mut end = start;
        while end < result.len() && events[result[end] as usize].point == point {
            end += 1;
        }
        if (end - start) % 2 == 1 {
            return false;
        }
        start = end;
    }
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
