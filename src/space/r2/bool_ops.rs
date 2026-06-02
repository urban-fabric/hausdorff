//! `r2::bool_ops` — polygon boolean operations (intersection, union,
//! difference, symmetric difference) via the Martinez–Rueda–Feito sweep-line
//! algorithm. Pure, no external deps.
//!
//! Ported from the reference implementation by F. Martínez et al. (2009/2013);
//! the [`BooleanOps`] trait mirrors `geo`'s so call sites migrate by import.
//! Input winding is normalised (exterior CCW, holes CW) so arbitrarily-wound
//! polygons work. Validated against `geo`'s `BooleanOps` over fixed
//! configurations (disjoint, nested, crossing, concave, touching/collinear,
//! identical, holes, multipolygons) plus ~900 randomized rectangle/triangle
//! pairs; result areas and the set-algebra invariants match.
//!
//! Sweep events live in a flat arena (`Vec<Event>`) and reference each other by
//! `u32` index rather than `Rc<RefCell<…>>`: this removes a heap allocation and
//! reference-count per event, and the runtime borrow-checks the comparators paid
//! on every comparison. The event queue is a small index min-heap whose ordering
//! reads the live arena (the comparison keys can change as segments are divided,
//! exactly as the reference algorithm assumes).
//!
//! **Known limitation:** a *T-vertex* — a vertex of one polygon lying exactly on
//! the interior of the other's edge — can mislabel the sweep's in/out fields and
//! drop part of the result. This is rare for real (float lon/lat) geometry but
//! common on integer grids. [`crate::space::r2::relate`] avoids the clipper for this
//! reason; the overlay kernels do not yet handle it.

use std::cmp::Ordering;

use super::{Coord, CoordsIter, LineString, MultiPolygon, Polygon};

/// Which boolean operation to compute.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Operation {
    Intersection,
    Union,
    Difference,
    Xor,
}

/// Polygon set operations on [`Polygon`]/[`MultiPolygon`].
pub trait BooleanOps {
    /// Region in both operands.
    fn intersection(&self, other: &Self) -> MultiPolygon;
    /// Region in either operand.
    fn union(&self, other: &Self) -> MultiPolygon;
    /// Region in `self` but not `other`.
    fn difference(&self, other: &Self) -> MultiPolygon;
    /// Region in exactly one operand.
    fn xor(&self, other: &Self) -> MultiPolygon;
}

impl BooleanOps for MultiPolygon {
    fn intersection(&self, other: &Self) -> MultiPolygon {
        run(self, other, Operation::Intersection)
    }
    fn union(&self, other: &Self) -> MultiPolygon {
        run(self, other, Operation::Union)
    }
    fn difference(&self, other: &Self) -> MultiPolygon {
        run(self, other, Operation::Difference)
    }
    fn xor(&self, other: &Self) -> MultiPolygon {
        run(self, other, Operation::Xor)
    }
}

impl BooleanOps for Polygon {
    fn intersection(&self, other: &Self) -> MultiPolygon {
        run(&one(self), &one(other), Operation::Intersection)
    }
    fn union(&self, other: &Self) -> MultiPolygon {
        run(&one(self), &one(other), Operation::Union)
    }
    fn difference(&self, other: &Self) -> MultiPolygon {
        run(&one(self), &one(other), Operation::Difference)
    }
    fn xor(&self, other: &Self) -> MultiPolygon {
        run(&one(self), &one(other), Operation::Xor)
    }
}

fn one(p: &Polygon) -> MultiPolygon {
    MultiPolygon::new(vec![p.clone()])
}

// ---- sweep events -----------------------------------------------------------

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum EdgeType {
    Normal,
    NonContributing,
    SameTransition,
    DifferentTransition,
}

/// A sweep event in the arena. `other`/`prev_in_result` are arena indices.
#[derive(Clone, Copy)]
struct Event {
    point: Coord,
    left: bool,
    other: Option<u32>,
    subject: bool,
    edge_type: EdgeType,
    in_out: bool,
    other_in_out: bool,
    in_result: bool,
    prev_in_result: Option<u32>,
}

/// Push a fresh event into the arena and return its index.
fn new_event(events: &mut Vec<Event>, point: Coord, subject: bool) -> u32 {
    let idx = events.len() as u32;
    events.push(Event {
        point,
        left: false,
        other: None,
        subject,
        edge_type: EdgeType::Normal,
        in_out: false,
        other_in_out: false,
        in_result: false,
        prev_in_result: None,
    });
    idx
}

/// The other endpoint of event `e`'s segment.
fn other_point(events: &[Event], e: u32) -> Coord {
    events[events[e as usize].other.unwrap() as usize].point
}

/// Is `p` below the segment `(a → b)` (a is the left endpoint)?  Uses the signed
/// area of the triangle.
fn signed_area(a: Coord, b: Coord, c: Coord) -> f64 {
    (a.x - c.x) * (b.y - c.y) - (b.x - c.x) * (a.y - c.y)
}

/// Order for the event queue: leftmost-lowest first; at a shared point, right
/// events before left events, then the lower segment first.
fn compare_events(events: &[Event], e1: u32, e2: u32) -> Ordering {
    let a = &events[e1 as usize];
    let b = &events[e2 as usize];
    if a.point.x > b.point.x {
        return Ordering::Greater;
    }
    if a.point.x < b.point.x {
        return Ordering::Less;
    }
    if a.point.y != b.point.y {
        return if a.point.y > b.point.y { Ordering::Greater } else { Ordering::Less };
    }
    // same point
    if a.left != b.left {
        return if a.left { Ordering::Greater } else { Ordering::Less };
    }
    // both left or both right at the same point: the lower segment is processed
    // first (matches Martinez `!e1.isBelow(e2.other)`).
    let s = signed_area(a.point, other_point(events, e1), other_point(events, e2));
    if s != 0.0 {
        return if s > 0.0 { Ordering::Less } else { Ordering::Greater };
    }
    // collinear at the same point: subject before clipping
    if a.subject != b.subject {
        return if a.subject { Ordering::Less } else { Ordering::Greater };
    }
    // Fully tied: break by arena index for a strict, stable total order. The
    // index min-heap relies on a total order (it has no implicit insertion-order
    // tiebreak the way `std::BinaryHeap` happened to provide).
    e1.cmp(&e2)
}

// ---- index min-heap event queue --------------------------------------------
//
// `std::BinaryHeap`'s `Ord` can't see the arena, so we keep a small binary
// min-heap over event indices whose comparator is `compare_events` against the
// live arena — preserving the reference algorithm's "compare current state"
// behaviour as segments are divided.

/// Push event index `idx` and sift it up to its place.
fn heap_push(heap: &mut Vec<u32>, events: &[Event], idx: u32) {
    heap.push(idx);
    let mut i = heap.len() - 1;
    while i > 0 {
        let parent = (i - 1) / 2;
        if compare_events(events, heap[i], heap[parent]) == Ordering::Less {
            heap.swap(i, parent);
            i = parent;
        } else {
            break;
        }
    }
}

/// Pop the `compare_events`-minimum event index, or `None` if empty.
fn heap_pop(heap: &mut Vec<u32>, events: &[Event]) -> Option<u32> {
    if heap.is_empty() {
        return None;
    }
    let top = heap.swap_remove(0); // move the last element into the root
    let n = heap.len();
    let mut i = 0;
    loop {
        let (l, r) = (2 * i + 1, 2 * i + 2);
        let mut min = i;
        if l < n && compare_events(events, heap[l], heap[min]) == Ordering::Less {
            min = l;
        }
        if r < n && compare_events(events, heap[r], heap[min]) == Ordering::Less {
            min = r;
        }
        if min == i {
            break;
        }
        heap.swap(i, min);
        i = min;
    }
    Some(top)
}

/// Sweep-line segment order: which left-event's segment is lower at the sweep.
fn compare_segments(events: &[Event], le1: u32, le2: u32) -> Ordering {
    if le1 == le2 {
        return Ordering::Equal;
    }
    let (p1, p1o) = (events[le1 as usize].point, other_point(events, le1));
    let (p2, p2o) = (events[le2 as usize].point, other_point(events, le2));

    if signed_area(p1, p1o, p2) != 0.0 || signed_area(p1, p1o, p2o) != 0.0 {
        // not collinear
        if p1 == p2 {
            // same left endpoint: order by the right endpoints
            return if below(p1, p1o, p2o) { Ordering::Less } else { Ordering::Greater };
        }
        if compare_events(events, le1, le2) == Ordering::Greater {
            // le1 is processed after le2: is le2's segment below le1's left point?
            return if !below(p2, p2o, p1) { Ordering::Less } else { Ordering::Greater };
        }
        return if below(p1, p1o, p2) { Ordering::Less } else { Ordering::Greater };
    }
    // collinear segments
    let (s1, s2) = (events[le1 as usize].subject, events[le2 as usize].subject);
    if s1 == s2 {
        if p1 == p2 {
            // identical: break ties consistently by arena index
            le1.cmp(&le2)
        } else {
            compare_events(events, le1, le2)
        }
    } else if s1 {
        Ordering::Greater
    } else {
        Ordering::Less
    }
}

fn below(a: Coord, b: Coord, p: Coord) -> bool {
    signed_area(a, b, p) > 0.0
}

// ---- the algorithm ----------------------------------------------------------

fn run(subject: &MultiPolygon, clipping: &MultiPolygon, op: Operation) -> MultiPolygon {
    // The sweep's in/out labelling assumes canonical winding (exterior CCW,
    // holes CW); normalise both operands so arbitrarily-wound input works.
    use super::winding::{Direction, Orient};
    let subject = subject.orient(Direction::Default);
    let clipping = clipping.orient(Direction::Default);

    // Two events per edge initially, plus splits; one growing arena instead of a
    // heap allocation per event.
    let cap = 2 * (subject.coords_count() + clipping.coords_count()) + 16;
    let mut events: Vec<Event> = Vec::with_capacity(cap);
    let mut queue: Vec<u32> = Vec::with_capacity(cap);
    for p in subject.iter() {
        fill_polygon(&mut events, &mut queue, p, true);
    }
    for p in clipping.iter() {
        fill_polygon(&mut events, &mut queue, p, false);
    }

    let sorted = subdivide(&mut events, &mut queue, op);
    let contours = connect_edges(&events, &sorted);
    MultiPolygon::new(contours)
}

fn fill_polygon(events: &mut Vec<Event>, queue: &mut Vec<u32>, p: &Polygon, subject: bool) {
    fill_ring(events, queue, &p.exterior().0, subject);
    for hole in p.interiors() {
        fill_ring(events, queue, &hole.0, subject);
    }
}

fn fill_ring(events: &mut Vec<Event>, queue: &mut Vec<u32>, coords: &[Coord], subject: bool) {
    let v = ring_vertices(coords);
    let n = v.len();
    if n < 3 {
        return;
    }
    for i in 0..n {
        let a = v[i];
        let b = v[(i + 1) % n];
        if a == b {
            continue;
        }
        let e1 = new_event(events, a, subject);
        let e2 = new_event(events, b, subject);
        events[e1 as usize].other = Some(e2);
        events[e2 as usize].other = Some(e1);
        // the leftmost endpoint is the "left" event
        if compare_events(events, e1, e2) == Ordering::Less {
            events[e1 as usize].left = true;
            events[e2 as usize].left = false;
        } else {
            events[e1 as usize].left = false;
            events[e2 as usize].left = true;
        }
        heap_push(queue, events, e1);
        heap_push(queue, events, e2);
    }
}

fn ring_vertices(ring: &[Coord]) -> &[Coord] {
    if ring.len() >= 2 && ring[0] == ring[ring.len() - 1] {
        &ring[..ring.len() - 1]
    } else {
        ring
    }
}

fn subdivide(events: &mut Vec<Event>, queue: &mut Vec<u32>, op: Operation) -> Vec<u32> {
    let mut sorted: Vec<u32> = Vec::new();
    // sweep-line status: left-events ordered by compare_segments
    let mut sweep: Vec<u32> = Vec::new();

    while let Some(event) = heap_pop(queue, events) {
        sorted.push(event);

        if events[event as usize].left {
            let pos = sweep_insert(&mut sweep, events, event);
            let prev = if pos > 0 { Some(sweep[pos - 1]) } else { None };
            let next = sweep.get(pos + 1).copied();

            compute_fields(events, event, prev, op);

            if let Some(next) = next
                && possible_intersection(events, queue, event, next) == 2
            {
                compute_fields(events, event, prev, op);
                compute_fields(events, next, Some(event), op);
            }
            if let Some(prev) = prev
                && possible_intersection(events, queue, prev, event) == 2
            {
                let prev_prev = sweep_prev_of(&sweep, prev);
                compute_fields(events, prev, prev_prev, op);
                compute_fields(events, event, Some(prev), op);
            }
        } else {
            // right event: operate on its left counterpart
            let left = events[event as usize].other.unwrap();
            if let Some(pos) = sweep.iter().position(|&e| e == left) {
                let prev = if pos > 0 { Some(sweep[pos - 1]) } else { None };
                let next = sweep.get(pos + 1).copied();
                sweep.remove(pos);
                if let (Some(prev), Some(next)) = (prev, next) {
                    possible_intersection(events, queue, prev, next);
                }
            }
        }
    }
    sorted
}

fn sweep_insert(sweep: &mut Vec<u32>, events: &[Event], event: u32) -> usize {
    let pos = sweep.partition_point(|&e| compare_segments(events, e, event) == Ordering::Less);
    sweep.insert(pos, event);
    pos
}

fn sweep_prev_of(sweep: &[u32], event: u32) -> Option<u32> {
    let pos = sweep.iter().position(|&e| e == event)?;
    if pos > 0 { Some(sweep[pos - 1]) } else { None }
}

fn is_vertical(events: &[Event], e: u32) -> bool {
    events[e as usize].point.x == other_point(events, e).x
}

fn compute_fields(events: &mut [Event], event: u32, prev: Option<u32>, op: Operation) {
    match prev {
        None => {
            let ev = &mut events[event as usize];
            ev.in_out = false;
            ev.other_in_out = true;
        }
        Some(prev) => {
            // Snapshot prev (a `Copy` of the event) and its verticality before
            // taking the mutable borrow of `event`.
            let pr = events[prev as usize];
            let prev_vertical = is_vertical(events, prev);
            let ev_subject = events[event as usize].subject;
            let ev = &mut events[event as usize];
            if ev_subject == pr.subject {
                ev.in_out = !pr.in_out;
                ev.other_in_out = pr.other_in_out;
            } else {
                ev.in_out = !pr.other_in_out;
                ev.other_in_out = if prev_vertical { !pr.in_out } else { pr.in_out };
            }
        }
    }
    let res = in_result(events, event, op);
    events[event as usize].in_result = res;
    // prev_in_result chain
    let pir = match prev {
        Some(prev) => {
            if !events[prev as usize].in_result || is_vertical(events, prev) {
                events[prev as usize].prev_in_result
            } else {
                Some(prev)
            }
        }
        None => None,
    };
    events[event as usize].prev_in_result = pir;
}

fn in_result(events: &[Event], event: u32, op: Operation) -> bool {
    let ev = events[event as usize];
    match ev.edge_type {
        EdgeType::Normal => match op {
            Operation::Intersection => !ev.other_in_out,
            Operation::Union => ev.other_in_out,
            Operation::Difference => {
                (ev.subject && ev.other_in_out) || (!ev.subject && !ev.other_in_out)
            }
            Operation::Xor => true,
        },
        EdgeType::SameTransition => {
            matches!(op, Operation::Intersection | Operation::Union)
        }
        EdgeType::DifferentTransition => op == Operation::Difference,
        EdgeType::NonContributing => false,
    }
}

/// Segment intersection: returns the number of intersection points (0, 1, or 2)
/// and the points.
fn intersect_segments(a1: Coord, a2: Coord, b1: Coord, b2: Coord) -> (u8, Coord, Coord) {
    let va = a2 - a1;
    let vb = b2 - b1;
    let e = b1 - a1;
    let kross = va.x * vb.y - va.y * vb.x;
    let sqr_len_a = va.x * va.x + va.y * va.y;

    if kross != 0.0 {
        let s = (e.x * vb.y - e.y * vb.x) / kross;
        if !(0.0..=1.0).contains(&s) {
            return (0, a1, a1);
        }
        let t = (e.x * va.y - e.y * va.x) / kross;
        if !(0.0..=1.0).contains(&t) {
            return (0, a1, a1);
        }
        return (1, a1 + va * s, a1);
    }

    // parallel
    let kross_e_va = e.x * va.y - e.y * va.x;
    if kross_e_va != 0.0 {
        return (0, a1, a1); // parallel, non-collinear
    }
    // collinear: project b1, b2 onto a
    let sa = (e.x * va.x + e.y * va.y) / sqr_len_a;
    let sb = sa + (vb.x * va.x + vb.y * va.y) / sqr_len_a;
    let smin = sa.min(sb);
    let smax = sa.max(sb);
    if smin > 1.0 || smax < 0.0 {
        return (0, a1, a1);
    }
    let lo = smin.max(0.0);
    let hi = smax.min(1.0);
    if lo == hi {
        (1, a1 + va * lo, a1)
    } else {
        (2, a1 + va * lo, a1 + va * hi)
    }
}

fn possible_intersection(events: &mut Vec<Event>, queue: &mut Vec<u32>, e1: u32, e2: u32) -> u8 {
    let a1 = events[e1 as usize].point;
    let a2 = other_point(events, e1);
    let b1 = events[e2 as usize].point;
    let b2 = other_point(events, e2);

    let (n, ip1, _ip2) = intersect_segments(a1, a2, b1, b2);

    if n == 0 {
        return 0;
    }
    if n == 1 && (a1 == b1 || a2 == b2 || a1 == b2 || a2 == b1) {
        // intersect only at a shared endpoint
        return 0;
    }

    if n == 1 {
        if a1 != ip1 && a2 != ip1 {
            divide_segment(events, queue, e1, ip1);
        }
        if b1 != ip1 && b2 != ip1 {
            divide_segment(events, queue, e2, ip1);
        }
        return 1;
    }

    // n == 2: collinear overlapping segments. Overlap within the same operand is
    // degenerate input — ignore it.
    if events[e1 as usize].subject == events[e2 as usize].subject {
        return 0;
    }

    // Faithful port of Martinez `possibleIntersection`'s overlap branch: order
    // the up-to-four endpoints, then subdivide so the overlapping portions become
    // a coincident pair. The coincident pair is later re-processed via the
    // left-coincide path, which assigns Same/DifferentTransition.
    let e1o = events[e1 as usize].other.unwrap();
    let e2o = events[e2 as usize].other.unwrap();
    let mut order: Vec<u32> = Vec::new();
    let left_coincide = a1 == b1;
    let right_coincide = a2 == b2;

    if left_coincide {
        // share left endpoint: handled below
    } else if compare_events(events, e1, e2) == Ordering::Greater {
        order.push(e2);
        order.push(e1);
    } else {
        order.push(e1);
        order.push(e2);
    }

    if right_coincide {
        // share right endpoint: handled below
    } else if compare_events(events, e1o, e2o) == Ordering::Greater {
        order.push(e2o);
        order.push(e1o);
    } else {
        order.push(e1o);
        order.push(e2o);
    }

    if left_coincide {
        // the two segments share their left endpoint: one is non-contributing,
        // the other carries the overlap transition.
        events[e2 as usize].edge_type = EdgeType::NonContributing;
        let same = events[e2 as usize].in_out == events[e1 as usize].in_out;
        events[e1 as usize].edge_type =
            if same { EdgeType::SameTransition } else { EdgeType::DifferentTransition };
        if !right_coincide {
            // trim the longer segment at the shorter one's right endpoint
            let target = events[order[1] as usize].other.unwrap();
            let p = events[order[0] as usize].point;
            divide_segment(events, queue, target, p);
        }
        return 2;
    }

    if right_coincide {
        // share the right endpoint: split the earlier-starting segment
        let p = events[order[1] as usize].point;
        divide_segment(events, queue, order[0], p);
        return 3;
    }

    // neither endpoint coincides
    let e0_other = events[order[0] as usize].other.unwrap();
    if e0_other != order[3] {
        // partial overlap, neither contains the other
        let p1 = events[order[1] as usize].point;
        let p2 = events[order[2] as usize].point;
        divide_segment(events, queue, order[0], p1);
        divide_segment(events, queue, order[1], p2);
    } else {
        // one segment fully contains the other
        let p1 = events[order[1] as usize].point;
        let p2 = events[order[2] as usize].point;
        divide_segment(events, queue, order[0], p1);
        let target = events[order[3] as usize].other.unwrap();
        divide_segment(events, queue, target, p2);
    }
    3
}

fn divide_segment(events: &mut Vec<Event>, queue: &mut Vec<u32>, e: u32, p: Coord) {
    let other = events[e as usize].other.unwrap();
    let subject = events[e as usize].subject;

    // right event for the left part
    let r = new_event(events, p, subject);
    // left event for the right part
    let l = new_event(events, p, subject);

    events[r as usize].left = false;
    events[r as usize].other = Some(e);
    events[l as usize].left = true;
    events[l as usize].other = Some(other);

    events[other as usize].other = Some(l);
    events[e as usize].other = Some(r);

    // keep left/right consistent if the split point ordering flipped
    if compare_events(events, l, other) == Ordering::Greater {
        events[other as usize].left = true;
        events[l as usize].left = false;
    }
    // (compare_events(e, r) == Greater is a rare numerical case; leave as-is.)

    heap_push(queue, events, l);
    heap_push(queue, events, r);
}

// ---- connecting result edges into contours ----------------------------------

fn connect_edges(events: &[Event], sorted: &[u32]) -> Vec<Polygon> {
    // events that are part of the result
    let mut result: Vec<u32> = sorted
        .iter()
        .copied()
        .filter(|&e| {
            let ev = &events[e as usize];
            (ev.left && ev.in_result)
                || (!ev.left && events[ev.other.unwrap() as usize].in_result)
        })
        .collect();
    result.sort_by(|&a, &b| compare_events(events, a, b));

    let n = result.len();
    if n == 0 {
        return Vec::new();
    }
    // position of each event index within `result` (arena index → result slot),
    // replacing the per-call HashMap with a flat lookup.
    let mut index_of = vec![u32::MAX; events.len()];
    for (i, &e) in result.iter().enumerate() {
        index_of[e as usize] = i as u32;
    }
    let other_pos: Vec<usize> = result
        .iter()
        .map(|&e| {
            let other = events[e as usize].other.unwrap();
            let p = index_of[other as usize];
            if p == u32::MAX { usize::MAX } else { p as usize }
        })
        .collect();
    if other_pos.contains(&usize::MAX) {
        // result events not paired up; bail (shouldn't happen for valid input)
        return Vec::new();
    }

    let mut visited = vec![false; n];
    let mut rings: Vec<Vec<Coord>> = Vec::new();

    for i in 0..n {
        if visited[i] {
            continue;
        }
        let mut ring: Vec<Coord> = vec![events[result[i] as usize].point];
        let mut pos = i;
        loop {
            visited[pos] = true;
            pos = other_pos[pos]; // far end of the current segment
            visited[pos] = true;
            ring.push(events[result[pos] as usize].point);
            pos = next_pos(events, &result, pos, &visited, i);
            if pos == i || pos >= n {
                break;
            }
        }
        if ring.len() >= 4 {
            ring.pop(); // last point duplicates the first
            rings.push(ring);
        }
    }

    assemble(rings)
}

/// The next unvisited result event sharing the current point (the other segment
/// incident there); returns `start` when the contour closes.
fn next_pos(events: &[Event], result: &[u32], pos: usize, visited: &[bool], start: usize) -> usize {
    let n = result.len();
    let p0 = events[result[pos] as usize].point;
    let mut np = pos + 1;
    while np < n && events[result[np] as usize].point == p0 {
        if !visited[np] {
            return np;
        }
        np += 1;
    }
    np = pos;
    while np > 0 {
        np -= 1;
        if events[result[np] as usize].point != p0 {
            break;
        }
        if !visited[np] {
            return np;
        }
    }
    start
}

/// Turn raw closed rings into polygons. The Martinez output can wind either way,
/// so classify by *containment depth* rather than signed area: a ring contained
/// in an even number of others is a shell, odd is a hole of its immediate
/// container. Orientation is then normalised (shells CCW, holes CW).
fn assemble(rings: Vec<Vec<Coord>>) -> Vec<Polygon> {
    let n = rings.len();
    if n == 0 {
        return Vec::new();
    }
    let reps: Vec<Coord> = rings.iter().map(|r| interior_point(r)).collect();
    let depth: Vec<usize> = (0..n)
        .map(|i| {
            (0..n)
                .filter(|&j| j != i && point_in_ring(reps[i], &rings[j]))
                .count()
        })
        .collect();

    // even-depth rings become shells; collect their holes (odd-depth rings whose
    // immediate container is that shell).
    let mut polys: Vec<Polygon> = Vec::new();
    for i in 0..n {
        if !depth[i].is_multiple_of(2) {
            continue; // hole; handled by its parent shell
        }
        let holes: Vec<LineString> = (0..n)
            .filter(|&h| depth[h] == depth[i] + 1 && immediate_parent(&rings, &reps, &depth, h) == Some(i))
            .map(|h| LineString::new(oriented(rings[h].clone(), false)))
            .collect();
        polys.push(Polygon::new(
            LineString::new(oriented(rings[i].clone(), true)),
            holes,
        ));
    }
    polys
}

/// The innermost ring that contains hole `h` at one level shallower.
fn immediate_parent(
    rings: &[Vec<Coord>],
    reps: &[Coord],
    depth: &[usize],
    h: usize,
) -> Option<usize> {
    let mut best: Option<usize> = None;
    for j in 0..rings.len() {
        if j != h && depth[j] + 1 == depth[h] && point_in_ring(reps[h], &rings[j]) {
            best = match best {
                Some(b) if ring_area2(&rings[b]).abs() <= ring_area2(&rings[j]).abs() => Some(b),
                _ => Some(j),
            };
        }
    }
    best
}

/// Reverse the ring if its winding doesn't match the requested orientation
/// (`ccw = true` for shells, `false` for holes).
fn oriented(mut ring: Vec<Coord>, ccw: bool) -> Vec<Coord> {
    if (ring_area2(&ring) > 0.0) != ccw {
        ring.reverse();
    }
    ring
}

/// A point strictly inside the ring: step inward from the bottom-most (hence
/// convex) vertex along its edge bisector.
fn interior_point(ring: &[Coord]) -> Coord {
    let v = ring_vertices(ring);
    let n = v.len();
    if n < 3 {
        return v.first().copied().unwrap_or(Coord::new(0.0, 0.0));
    }
    let mut k = 0;
    for i in 1..n {
        if v[i].y < v[k].y || (v[i].y == v[k].y && v[i].x < v[k].x) {
            k = i;
        }
    }
    let cur = v[k];
    let e1 = v[(k + n - 1) % n] - cur;
    let e2 = v[(k + 1) % n] - cur;
    let (u1, l1) = unit(e1);
    let (u2, l2) = unit(e2);
    let mut bis = Coord::new(u1.x + u2.x, u1.y + u2.y);
    let blen = (bis.x * bis.x + bis.y * bis.y).sqrt();
    bis = if blen < 1e-12 {
        Coord::new(-u1.y, u1.x) // edges nearly opposite: use a perpendicular
    } else {
        Coord::new(bis.x / blen, bis.y / blen)
    };
    let step = l1.min(l2) * 1e-3;
    let p = Coord::new(cur.x + bis.x * step, cur.y + bis.y * step);
    if point_in_ring(p, ring) {
        p
    } else {
        Coord::new(cur.x - bis.x * step, cur.y - bis.y * step)
    }
}

fn unit(c: Coord) -> (Coord, f64) {
    let l = (c.x * c.x + c.y * c.y).sqrt();
    if l < 1e-12 {
        (Coord::new(0.0, 0.0), 0.0)
    } else {
        (Coord::new(c.x / l, c.y / l), l)
    }
}

fn ring_area2(ring: &[Coord]) -> f64 {
    let v = ring_vertices(ring);
    let n = v.len();
    if n < 3 {
        return 0.0;
    }
    let mut s = 0.0;
    for i in 0..n {
        let a = v[i];
        let b = v[(i + 1) % n];
        s += a.x * b.y - b.x * a.y;
    }
    s
}

fn point_in_ring(p: Coord, ring: &[Coord]) -> bool {
    let v = ring_vertices(ring);
    let n = v.len();
    if n < 3 {
        return false;
    }
    let mut inside = false;
    let mut j = n - 1;
    for i in 0..n {
        let (pi, pj) = (v[i], v[j]);
        if (pi.y > p.y) != (pj.y > p.y)
            && p.x < (pj.x - pi.x) * (p.y - pi.y) / (pj.y - pi.y) + pi.x
        {
            inside = !inside;
        }
        j = i;
    }
    inside
}
