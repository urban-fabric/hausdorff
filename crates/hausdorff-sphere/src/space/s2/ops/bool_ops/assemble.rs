//! Classify the noded arcs, then walk the survivors into loops and nest them.
//!
//! After the arrangement, every arc of `a` lies entirely inside `b`, entirely
//! outside it, or exactly along its boundary — it cannot be partly in and partly
//! out, because any transition would have been a crossing and crossings are now
//! vertices. That trichotomy is what makes classification a single point test per
//! arc rather than anything geometric.
//!
//! **The test point is the arc's midpoint, and it is not offset.** Offsetting to
//! one side is the usual way to classify an arrangement edge, and it is fragile:
//! the offset has to be small enough not to cross a neighbouring arc and large
//! enough to be numerically meaningful, and there is no bound that is reliably
//! both. Here it is unnecessary. An arc of `a` is already known to be on `a`'s
//! boundary, so the only open question is its relationship to `b`, and its
//! midpoint is generically *not* on `b`'s boundary. The one case where it is —
//! the two boundaries running along each other — is detected structurally as a
//! duplicate arc and handled by the table below, never by a point test.
//!
//! ## The selection rules
//!
//! With `in_other` meaning "this arc's midpoint is inside the other operand":
//!
//! | operation | keep arcs of `a` | keep arcs of `b` |
//! |---|---|---|
//! | union | `!in_b` | `!in_a` |
//! | intersection | `in_b` | `in_a` |
//! | difference `a \ b` | `!in_b` | `in_a`, **reversed** |
//!
//! The reversal in `difference` is the whole content of that operation: the part
//! of `b`'s boundary inside `a` becomes a *hole* boundary, and a hole runs the
//! other way round so that the surviving region stays on the left of every arc.
//!
//! Symmetric difference is not in the table. It is computed as
//! `union(a \ b, b \ a)`, which is its definition, rather than as a fourth
//! orientation rule. That costs an extra pass and buys correctness from code
//! already exercised by the other three — the alternative is a sign convention
//! that is easy to get backwards and whose errors look like plausible geometry.
//!
//! ## Coincident boundaries
//!
//! When the two operands share a stretch of boundary, the arrangement produces two
//! arcs on the same pair of vertices. A point test cannot resolve those (the
//! midpoint is on both boundaries), so they are matched structurally:
//!
//! | | union | intersection | difference |
//! |---|---|---|---|
//! | same direction (regions on the same side) | keep one | keep one | drop both |
//! | opposite direction (regions face each other) | drop both | drop both | keep `a`'s |
//!
//! "Same direction" means the two regions lie on the same side of the shared
//! boundary, so it is interior to the union and to the intersection alike. Opposite
//! means they meet along it, so the union swallows it and the intersection is
//! empty there.

use std::collections::{BTreeMap, HashMap};

use super::Operation;
use super::arrangement::{Arc, Side, VertexPool};
use crate::space::s2::{Loop, Polygon, S2Point, polygon_signed_area_vos};

/// Keep the arcs the operation's boundary is made of, oriented so the result lies
/// to their left.
pub(super) fn classify(
    pool: &VertexPool,
    arcs: &[Arc],
    op: Operation,
    a_contains: &dyn Fn(S2Point) -> bool,
    b_contains: &dyn Fn(S2Point) -> bool,
) -> Vec<(usize, usize)> {
    // Index the arcs by their unordered endpoint pair so a coincident boundary
    // shows up as two entries on one key.
    //
    // **Ordered, not hashed, and that is load-bearing.** This map is *iterated*
    // below, so its order is the order of `kept` — which is the order `trace`
    // seeds its walks in, which fixes each traced ring's starting vertex. `nest`
    // then measures a ring with `polygon_signed_area_vos`, a fan sum whose apex is
    // `ring[0]`, and for a ring enclosing more than a hemisphere that apex decides
    // the *sign*. Through a `HashMap` the whole chain ran off the per-process hash
    // seed: a 200°-wide lune minus an interior quad returned its own complement on
    // some runs and the right answer on others, from identical input. A `BTreeMap`
    // on pool indices makes the order a function of the input alone.
    let mut by_pair: BTreeMap<(usize, usize), Vec<Arc>> = BTreeMap::new();
    for arc in arcs {
        let key = if arc.start <= arc.end {
            (arc.start, arc.end)
        } else {
            (arc.end, arc.start)
        };
        by_pair.entry(key).or_default().push(*arc);
    }

    let mut kept = Vec::new();
    for (_, group) in by_pair {
        let from_a: Vec<Arc> = group
            .iter()
            .copied()
            .filter(|x| x.side == Side::A)
            .collect();
        let from_b: Vec<Arc> = group
            .iter()
            .copied()
            .filter(|x| x.side == Side::B)
            .collect();

        if !from_a.is_empty() && !from_b.is_empty() {
            // Coincident boundary: resolved by direction, not by a point test.
            let a_arc = from_a[0];
            let b_arc = from_b[0];
            let same_direction = a_arc.start == b_arc.start && a_arc.end == b_arc.end;
            match (op, same_direction) {
                (Operation::Union | Operation::Intersection, true) => {
                    kept.push((a_arc.start, a_arc.end));
                }
                (Operation::Union | Operation::Intersection, false) => {}
                (Operation::Difference, true) => {}
                (Operation::Difference, false) => kept.push((a_arc.start, a_arc.end)),
            }
            continue;
        }

        for arc in group {
            let mid = midpoint(pool.point(arc.start), pool.point(arc.end));
            let inside_other = match arc.side {
                Side::A => b_contains(mid),
                Side::B => a_contains(mid),
            };
            let keep = match (op, arc.side) {
                (Operation::Union, _) => !inside_other,
                (Operation::Intersection, _) => inside_other,
                (Operation::Difference, Side::A) => !inside_other,
                (Operation::Difference, Side::B) => inside_other,
            };
            if !keep {
                continue;
            }
            // `b`'s arcs bound a hole in a difference, so they run the other way.
            if op == Operation::Difference && arc.side == Side::B {
                kept.push((arc.end, arc.start));
            } else {
                kept.push((arc.start, arc.end));
            }
        }
    }
    kept
}

/// The midpoint of a great-circle arc: the normalized sum of its endpoints.
///
/// Valid because the arc is shorter than π, so the endpoints are not antipodal and
/// their sum does not vanish.
fn midpoint(a: S2Point, b: S2Point) -> S2Point {
    S2Point::from_direction(a.as_vector().add(b.as_vector()).normalize())
}

/// Walk the kept arcs into closed loops.
///
/// Every vertex of a boundary has equal in- and out-degree — it is the boundary of
/// a region, so a walk that enters must leave — which is what guarantees the walk
/// closes. Where several arcs leave one vertex, the **sharpest right turn** is
/// taken relative to the incoming direction. That is not arbitrary: at a vertex
/// where two components of the result touch, turning as tightly as possible keeps
/// the walk inside the component it is currently tracing instead of hopping to the
/// other one and welding two faces into one self-touching ring.
pub(super) fn trace(pool: &VertexPool, arcs: &[(usize, usize)]) -> Vec<Vec<S2Point>> {
    // Hashed is fine here, unlike `by_pair`: this map is only ever looked up by
    // key, never iterated, so its order cannot reach the output. The order that
    // matters is `arcs`, and each adjacency list keeps that insertion order.
    let mut outgoing: HashMap<usize, Vec<usize>> = HashMap::new();
    for (i, (start, _)) in arcs.iter().enumerate() {
        outgoing.entry(*start).or_default().push(i);
    }
    let mut used = vec![false; arcs.len()];
    let mut loops = Vec::new();

    for seed in 0..arcs.len() {
        if used[seed] {
            continue;
        }
        let mut ring = Vec::new();
        let mut current = seed;
        let start_vertex = arcs[seed].0;

        loop {
            used[current] = true;
            ring.push(pool.point(arcs[current].0));
            let at = arcs[current].1;
            if at == start_vertex {
                break;
            }
            let Some(candidates) = outgoing.get(&at) else {
                // A dangling arc: the walk cannot close. Treated as a dropped
                // ring rather than a panic, because a caller feeding
                // self-intersecting input should get a best effort, not a crash.
                ring.clear();
                break;
            };
            let available: Vec<usize> = candidates.iter().copied().filter(|i| !used[*i]).collect();
            let Some(next) = pick_sharpest_right(pool, arcs, current, &available) else {
                ring.clear();
                break;
            };
            current = next;
        }

        if ring.len() >= 3 {
            loops.push(ring);
        }
    }
    loops
}

/// Of the arcs leaving the current arc's end, the one turning most sharply to the
/// right.
fn pick_sharpest_right(
    pool: &VertexPool,
    arcs: &[(usize, usize)],
    current: usize,
    available: &[usize],
) -> Option<usize> {
    match available.len() {
        0 => None,
        1 => Some(available[0]),
        _ => {
            let pivot = pool.point(arcs[current].1);
            let incoming = pool.point(arcs[current].0);
            // Work in the tangent plane at the pivot: project both the arriving
            // direction and each candidate's departing direction, then compare
            // signed angles. A tangent-plane comparison is the spherical analogue
            // of the planar angular sort and needs no projection of the geometry.
            let basis = pivot.tangent_basis();
            let back = tangent_angle(&basis, pivot, incoming);
            available.iter().copied().min_by(|&x, &y| {
                let ax = turn_from(&basis, pivot, pool.point(arcs[x].1), back);
                let ay = turn_from(&basis, pivot, pool.point(arcs[y].1), back);
                ax.total_cmp(&ay)
            })
        }
    }
}

/// The angle of the direction from `pivot` toward `target`, in `pivot`'s tangent
/// plane.
fn tangent_angle(
    basis: &[hausdorff_base::space::r3::Vector<f64>; 2],
    pivot: S2Point,
    target: S2Point,
) -> f64 {
    let v = pivot.as_vector();
    let d = target.as_vector().sub(v.mul(v.dot(target.as_vector())));
    d.dot(basis[1]).atan2(d.dot(basis[0]))
}

/// How far right the departure toward `target` turns, measured from the direction
/// the walk arrived from. Smaller is sharper.
fn turn_from(
    basis: &[hausdorff_base::space::r3::Vector<f64>; 2],
    pivot: S2Point,
    target: S2Point,
    back: f64,
) -> f64 {
    let out = tangent_angle(basis, pivot, target);
    (out - back).rem_euclid(std::f64::consts::TAU)
}

/// Group traced rings into polygons: counter-clockwise rings are shells,
/// clockwise rings are holes, and each hole joins the smallest shell containing
/// it.
///
/// "Smallest containing shell" rather than "any containing shell" is what makes
/// nested rings come out right — a hole inside an island inside a lake belongs to
/// the island, not the outer shell, and only the area ordering distinguishes them.
///
/// **Known limit on the output side.** Where there is more than one traced ring,
/// shell-versus-hole is decided from the sign of `polygon_signed_area_vos`, and
/// that sign is ambiguous for a ring enclosing more than a hemisphere: the van
/// Oosterom sum returns the spherical excess, so `A > 2π` comes back as `A − 4π`
/// and reads as clockwise.
///
/// Worse than ambiguous, in fact: that function is a *fan sum* apexed at
/// `ring[0]`, so above a hemisphere the sign depends on which vertex the ring
/// happens to start at, and a ring is the same ring whatever its rotation. Two
/// consequences follow. A 250° lune minus an interior quad returns the complement
/// — both rings read as holes, so the `shells.is_empty()` branch below reverses
/// them. And the answer was once not even stable: ring start order came from a
/// `HashMap` (see `classify`), so identical input gave different areas on
/// different runs. The ordering is fixed; the measurement is not.
///
/// Two of the three places this could bite are closed. Inputs are immune, because
/// their orientation is read off the `Polygon` type rather than measured (see
/// `push_ring`). A single-ring result is immune, because a lone ring is taken as a
/// shell by construction rather than measured. What remains unverified is a result
/// that is *both* larger than a hemisphere *and* multi-ring — which needs operands
/// already near-global — and that case is not claimed.
///
/// **Hole rings are reversed back to counter-clockwise before being stored.**
/// There are two different conventions in play and conflating them produces a
/// polygon whose area is zero while looking structurally perfect. The
/// *arrangement* orients every arc so the region lies to its left, which makes a
/// hole run clockwise — that is what lets the classifier be a pure point test.
/// `Polygon`, however, stores holes with the *same* winding as the shell and
/// carries the hole semantics in the struct: `Polygon::area` subtracts
/// `Loop::area` directly, and `Loop::area` is the area to the ring's left. Leave
/// a hole clockwise and `Loop::area` returns the complement — about 4π — so the
/// subtraction goes negative and clamps to zero.
pub(super) fn nest(rings: Vec<Vec<S2Point>>) -> Vec<Polygon> {
    // A single traced ring is **always** a shell, with no measurement involved.
    //
    // Every kept arc carries the result region on its left, so a lone ring's
    // left-region *is* the result — which is exactly what `Loop::area` returns.
    // Deciding this by signed area instead is what broke a region larger than a
    // hemisphere: the van Oosterom sum gives the spherical excess, so `A > 2π`
    // reports `A − 4π` and reads as clockwise. A 200°-wide lune unioned with a lune
    // inside it came back as 5.59 sr — the complement of the correct 6.98.
    //
    // Taking the ring at its word removes the measurement, and with it the
    // ambiguity, for every simply-connected result at any size.
    if rings.len() == 1 {
        return vec![Polygon::from_shell(Loop::new(
            rings.into_iter().next().unwrap(),
        ))];
    }

    let mut shells: Vec<Vec<S2Point>> = Vec::new();
    let mut holes: Vec<Vec<S2Point>> = Vec::new();
    for ring in rings {
        if polygon_signed_area_vos(&ring) >= 0.0 {
            shells.push(ring);
        } else {
            holes.push(ring);
        }
    }

    // A clockwise ring with no shell to belong to is the complement of a region.
    // Reversing it makes it the shell it describes.
    if shells.is_empty() && !holes.is_empty() {
        for hole in &mut holes {
            hole.reverse();
        }
        shells.append(&mut holes);
    }

    let mut built: Vec<(Loop, Vec<Loop>)> = shells
        .into_iter()
        .map(|s| (Loop::new(s), Vec::new()))
        .collect();

    for hole in holes {
        let probe = hole[0];
        let mut best: Option<usize> = None;
        for (i, (shell, _)) in built.iter().enumerate() {
            if !shell.contains_point_closed(probe) {
                continue;
            }
            let smaller = match best {
                None => true,
                Some(j) => shell.area() < built[j].0.area(),
            };
            if smaller {
                best = Some(i);
            }
        }
        if let Some(i) = best {
            // Back to the shell's winding for storage; see the note above.
            let mut ccw = hole;
            ccw.reverse();
            built[i].1.push(Loop::new(ccw));
        }
        // A hole inside nothing is dropped: it bounds no region of the result.
    }

    built
        .into_iter()
        .map(|(shell, holes)| Polygon::new(shell, holes))
        .collect()
}
