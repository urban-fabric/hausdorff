//! Turning the surviving arrangement edges into closed rings.
//!
//! The edges handed here bound the odd-even region, so they form closed curves
//! and every vertex has even degree. What is *not* known per edge is which of
//! its sides is inside, and a walk that guesses gets pinches wrong: where four
//! edges meet at a point, pairing the wrong two traces a figure-eight — one
//! self-touching ring instead of two rings meeting at a point.
//!
//! So this traces faces rather than rings. Each edge is walked once in each
//! direction, and at every vertex the walk takes the next edge clockwise from
//! the one it arrived on, which keeps the same face on its left the whole way
//! round. Every directed edge belongs to exactly one face, so the walk is a
//! partition and needs no arbitrary starting choice.
//!
//! Bounded faces come out counter-clockwise and the unbounded face clockwise,
//! so keeping the positively-signed cycles keeps exactly the bounded ones. That
//! includes the inside of a hole, which is bounded — telling a hole from a shell
//! is containment depth's job, not this walk's.

use super::node::Nodes;
use crate::space::r2::Coord;

/// Walk `edges` into the bounded faces of the arrangement they form.
pub(crate) fn rings(edges: &[(usize, usize)], nodes: &Nodes) -> Option<Vec<Vec<Coord>>> {
    if edges.is_empty() {
        return Some(Vec::new());
    }
    let vertex_count = edges.iter().map(|&(a, b)| a.max(b)).max()? + 1;

    // Directed edge 2i is edges[i] forward, 2i+1 is it reversed.
    let head = |d: usize| {
        let (a, b) = edges[d / 2];
        if d.is_multiple_of(2) { b } else { a }
    };
    let tail = |d: usize| {
        let (a, b) = edges[d / 2];
        if d.is_multiple_of(2) { a } else { b }
    };

    // Outgoing directed edges per vertex, sorted by heading so the "next
    // clockwise" choice is a step in a precomputed cyclic order rather than a
    // scan. Sorting also makes the walk independent of the order `edges`
    // arrived in, which a hash map does not fix.
    let mut outgoing: Vec<Vec<usize>> = vec![Vec::new(); vertex_count];
    for d in 0..edges.len() * 2 {
        outgoing[tail(d)].push(d);
    }
    let heading = |d: usize| {
        let (from, to) = (nodes.coord(tail(d)), nodes.coord(head(d)));
        (to.y - from.y).atan2(to.x - from.x)
    };
    for out in &mut outgoing {
        out.sort_by(|&x, &y| heading(x).total_cmp(&heading(y)).then(x.cmp(&y)));
    }
    // Where each directed edge sits in its tail's cyclic order.
    let mut slot = vec![0usize; edges.len() * 2];
    for out in &outgoing {
        for (i, &d) in out.iter().enumerate() {
            slot[d] = i;
        }
    }

    let twin = |d: usize| d ^ 1;
    let mut used = vec![false; edges.len() * 2];
    let mut faces: Vec<Vec<Coord>> = Vec::new();

    for start in 0..edges.len() * 2 {
        if used[start] {
            continue;
        }
        let mut ring: Vec<Coord> = Vec::new();
        let mut d = start;
        loop {
            used[d] = true;
            ring.push(nodes.coord(tail(d)));
            // Arrive at head(d); leave along the edge one step clockwise from
            // the way we came in, which is the predecessor of our twin in the
            // vertex's counter-clockwise ordering.
            let back = twin(d);
            let at = tail(back);
            let ring_of = &outgoing[at];
            let i = slot[back];
            d = ring_of[(i + ring_of.len() - 1) % ring_of.len()];
            if d == start {
                break;
            }
            if used[d] {
                // Only reachable if the arrangement is not a clean planar
                // graph, which would mean the caller's edge set is malformed.
                return None;
            }
        }
        if ring.len() >= 3 && crate::space::r2::util::signed_area2(&ring) > 0.0 {
            faces.push(ring);
        }
    }
    Some(faces)
}
