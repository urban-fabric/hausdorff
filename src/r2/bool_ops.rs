//! `r2::bool_ops` — polygon boolean operations (intersection, union,
//! difference, symmetric difference) via the Martinez–Rueda–Feito sweep-line
//! algorithm. Pure, no external deps.
//!
//! Ported from the reference implementation by F. Martínez et al. (2009/2013);
//! the [`BooleanOps`] trait mirrors `geo`'s so call sites migrate by import.
//!
//! **Status: work in progress.** Correct for disjoint and nested polygons; the
//! crossing/overlapping cases (segment subdivision + sweep ordering) are not yet
//! right. Until complete, [`crate::compat`] keeps `geo`'s `BooleanOps` as the
//! blessed path, and this module is intentionally not re-exported at the `r2`
//! root.

use std::cell::RefCell;
use std::cmp::Ordering;
use std::collections::BinaryHeap;
use std::rc::Rc;

use super::{Coord, LineString, MultiPolygon, Polygon};

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

type Link = Rc<RefCell<Event>>;

struct Event {
    point: Coord,
    left: bool,
    other: Option<Link>,
    subject: bool,
    edge_type: EdgeType,
    in_out: bool,
    other_in_out: bool,
    in_result: bool,
    prev_in_result: Option<Link>,
    // connection bookkeeping
    result_in_out: bool,
    pos: usize,
    contour_id: i32,
}

impl Event {
    fn new(point: Coord, subject: bool) -> Link {
        Rc::new(RefCell::new(Event {
            point,
            left: false,
            other: None,
            subject,
            edge_type: EdgeType::Normal,
            in_out: false,
            other_in_out: false,
            in_result: false,
            prev_in_result: None,
            result_in_out: false,
            pos: 0,
            contour_id: 0,
        }))
    }
}

fn other_point(e: &Link) -> Coord {
    e.borrow().other.as_ref().unwrap().borrow().point
}

/// Is `p` below the segment `(a → b)` (a is the left endpoint)?  Uses the signed
/// area of the triangle.
fn signed_area(a: Coord, b: Coord, c: Coord) -> f64 {
    (a.x - c.x) * (b.y - c.y) - (b.x - c.x) * (a.y - c.y)
}

/// Order for the event queue: leftmost-lowest first; at a shared point, right
/// events before left events, then the lower segment first.
fn compare_events(e1: &Link, e2: &Link) -> Ordering {
    let (a, b) = (e1.borrow(), e2.borrow());
    if a.point.x > b.point.x {
        return Ordering::Greater;
    }
    if a.point.x < b.point.x {
        return Ordering::Less;
    }
    if a.point.y != b.point.y {
        return if a.point.y > b.point.y {
            Ordering::Greater
        } else {
            Ordering::Less
        };
    }
    // same point
    if a.left != b.left {
        return if a.left {
            Ordering::Greater
        } else {
            Ordering::Less
        };
    }
    drop((a, b));
    // both left or both right at the same point: order by the other endpoint
    let s = signed_area(e1.borrow().point, other_point(e1), other_point(e2));
    if s != 0.0 {
        // e1 below e2 -> e1 first
        return if s > 0.0 {
            Ordering::Greater
        } else {
            Ordering::Less
        };
    }
    if e1.borrow().subject == e2.borrow().subject {
        Ordering::Equal
    } else if !e1.borrow().subject {
        Ordering::Less
    } else {
        Ordering::Greater
    }
}

/// Min-heap wrapper (BinaryHeap is a max-heap, so invert the order).
struct QueueItem(Link);
impl PartialEq for QueueItem {
    fn eq(&self, other: &Self) -> bool {
        compare_events(&self.0, &other.0) == Ordering::Equal
    }
}
impl Eq for QueueItem {}
impl Ord for QueueItem {
    fn cmp(&self, other: &Self) -> Ordering {
        compare_events(&self.0, &other.0).reverse()
    }
}
impl PartialOrd for QueueItem {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

/// Sweep-line segment order: which left-event's segment is lower at the sweep.
fn compare_segments(le1: &Link, le2: &Link) -> Ordering {
    if Rc::ptr_eq(le1, le2) {
        return Ordering::Equal;
    }
    let (p1, p1o) = (le1.borrow().point, other_point(le1));
    let (p2, p2o) = (le2.borrow().point, other_point(le2));

    if signed_area(p1, p1o, p2) != 0.0 || signed_area(p1, p1o, p2o) != 0.0 {
        // not collinear
        if p1 == p2 {
            // same left endpoint: order by the right endpoints
            return if below(p1, p1o, p2o) {
                Ordering::Less
            } else {
                Ordering::Greater
            };
        }
        if compare_events(le1, le2) == Ordering::Greater {
            // le1 is processed after le2: is le2's segment below le1's left point?
            return if !below(p2, p2o, p1) {
                Ordering::Less
            } else {
                Ordering::Greater
            };
        }
        return if below(p1, p1o, p2) {
            Ordering::Less
        } else {
            Ordering::Greater
        };
    }
    // collinear segments
    if le1.borrow().subject == le2.borrow().subject {
        if p1 == p2 {
            // identical: break ties consistently by pointer
            (Rc::as_ptr(le1) as usize).cmp(&(Rc::as_ptr(le2) as usize))
        } else {
            compare_events(le1, le2)
        }
    } else if le1.borrow().subject {
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
    let mut queue: BinaryHeap<QueueItem> = BinaryHeap::new();
    for p in subject.iter() {
        fill_polygon(p, true, &mut queue);
    }
    for p in clipping.iter() {
        fill_polygon(p, false, &mut queue);
    }

    let sorted = subdivide(queue, op);
    let contours = connect_edges(&sorted);
    MultiPolygon::new(contours)
}

fn fill_polygon(p: &Polygon, subject: bool, queue: &mut BinaryHeap<QueueItem>) {
    fill_ring(&p.exterior().0, subject, queue);
    for hole in p.interiors() {
        fill_ring(&hole.0, subject, queue);
    }
}

fn fill_ring(coords: &[Coord], subject: bool, queue: &mut BinaryHeap<QueueItem>) {
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
        let e1 = Event::new(a, subject);
        let e2 = Event::new(b, subject);
        e1.borrow_mut().other = Some(e2.clone());
        e2.borrow_mut().other = Some(e1.clone());
        // the leftmost endpoint is the "left" event
        if compare_events(&e1, &e2) == Ordering::Less {
            e1.borrow_mut().left = true;
            e2.borrow_mut().left = false;
        } else {
            e1.borrow_mut().left = false;
            e2.borrow_mut().left = true;
        }
        queue.push(QueueItem(e1));
        queue.push(QueueItem(e2));
    }
}

fn ring_vertices(ring: &[Coord]) -> &[Coord] {
    if ring.len() >= 2 && ring[0] == ring[ring.len() - 1] {
        &ring[..ring.len() - 1]
    } else {
        ring
    }
}

fn subdivide(mut queue: BinaryHeap<QueueItem>, op: Operation) -> Vec<Link> {
    let mut sorted: Vec<Link> = Vec::new();
    // sweep-line status: left-events ordered by compare_segments
    let mut sweep: Vec<Link> = Vec::new();

    while let Some(QueueItem(event)) = queue.pop() {
        sorted.push(event.clone());

        if event.borrow().left {
            let pos = sweep_insert(&mut sweep, &event);
            let prev = if pos > 0 { Some(sweep[pos - 1].clone()) } else { None };
            let next = sweep.get(pos + 1).cloned();

            compute_fields(&event, prev.as_ref(), op);

            if let Some(next) = next.as_ref() {
                if possible_intersection(&event, next, &mut queue) == 2 {
                    compute_fields(&event, prev.as_ref(), op);
                    compute_fields(next, Some(&event), op);
                }
            }
            if let Some(prev) = prev.as_ref() {
                if possible_intersection(prev, &event, &mut queue) == 2 {
                    let prev_prev = sweep_prev_of(&sweep, prev);
                    compute_fields(prev, prev_prev.as_ref(), op);
                    compute_fields(&event, Some(prev), op);
                }
            }
        } else {
            // right event: operate on its left counterpart
            let left = event.borrow().other.clone().unwrap();
            if let Some(pos) = sweep.iter().position(|e| Rc::ptr_eq(e, &left)) {
                let prev = if pos > 0 { Some(sweep[pos - 1].clone()) } else { None };
                let next = sweep.get(pos + 1).cloned();
                sweep.remove(pos);
                if let (Some(prev), Some(next)) = (prev, next) {
                    possible_intersection(&prev, &next, &mut queue);
                }
            }
        }
    }
    sorted
}

fn sweep_insert(sweep: &mut Vec<Link>, event: &Link) -> usize {
    let pos = sweep.partition_point(|e| compare_segments(e, event) == Ordering::Less);
    sweep.insert(pos, event.clone());
    pos
}

fn sweep_prev_of(sweep: &[Link], event: &Link) -> Option<Link> {
    let pos = sweep.iter().position(|e| Rc::ptr_eq(e, event))?;
    if pos > 0 {
        Some(sweep[pos - 1].clone())
    } else {
        None
    }
}

fn is_vertical(e: &Link) -> bool {
    e.borrow().point.x == other_point(e).x
}

fn compute_fields(event: &Link, prev: Option<&Link>, op: Operation) {
    {
        let mut ev = event.borrow_mut();
        match prev {
            None => {
                ev.in_out = false;
                ev.other_in_out = true;
            }
            Some(prev) => {
                let pr = prev.borrow();
                if ev.subject == pr.subject {
                    ev.in_out = !pr.in_out;
                    ev.other_in_out = pr.other_in_out;
                } else {
                    ev.in_out = !pr.other_in_out;
                    ev.other_in_out = if is_vertical(prev) { !pr.in_out } else { pr.in_out };
                }
            }
        }
    }
    let res = in_result(event, op);
    event.borrow_mut().in_result = res;
    // prev_in_result chain
    let pir = match prev {
        Some(prev) => {
            if !prev.borrow().in_result || is_vertical(prev) {
                prev.borrow().prev_in_result.clone()
            } else {
                Some(prev.clone())
            }
        }
        None => None,
    };
    event.borrow_mut().prev_in_result = pir;
}

fn in_result(event: &Link, op: Operation) -> bool {
    let ev = event.borrow();
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
fn intersect_segments(
    a1: Coord,
    a2: Coord,
    b1: Coord,
    b2: Coord,
) -> (u8, Coord, Coord) {
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

fn possible_intersection(e1: &Link, e2: &Link, queue: &mut BinaryHeap<QueueItem>) -> u8 {
    let a1 = e1.borrow().point;
    let a2 = other_point(e1);
    let b1 = e2.borrow().point;
    let b2 = other_point(e2);

    let (n, ip1, ip2) = intersect_segments(a1, a2, b1, b2);

    if n == 0 {
        return 0;
    }
    if n == 1 && (a1 == b1 || a2 == b2 || a1 == b2 || a2 == b1) {
        // intersect only at a shared endpoint
        return 0;
    }

    if n == 1 {
        if a1 != ip1 && a2 != ip1 {
            divide_segment(e1, ip1, queue);
        }
        if b1 != ip1 && b2 != ip1 {
            divide_segment(e2, ip1, queue);
        }
        return 1;
    }

    // n == 2: overlapping collinear segments
    let mut events: Vec<(Coord, bool)> = Vec::new(); // (point, is_left_endpoint_of_some_segment)
    let left_coincide = a1 == b1;
    let right_coincide = a2 == b2;

    if !left_coincide {
        if compare_coords(a1, b1) == Ordering::Greater {
            divide_segment(e2, a1, queue);
        } else {
            divide_segment(e1, b1, queue);
        }
    }
    if !right_coincide {
        if compare_coords(a2, b2) == Ordering::Greater {
            divide_segment(e2, b2, queue);
        } else {
            divide_segment(e1, a2, queue);
        }
    }
    let _ = (ip1, ip2, &mut events);

    // mark edge types for the overlapping pair
    if left_coincide {
        // the two segments now share their left endpoint; one is non-contributing
        let same = e1.borrow().in_out == e2.borrow().in_out;
        e2.borrow_mut().edge_type = EdgeType::NonContributing;
        e1.borrow_mut().edge_type = if same {
            EdgeType::SameTransition
        } else {
            EdgeType::DifferentTransition
        };
    }
    2
}

fn compare_coords(a: Coord, b: Coord) -> Ordering {
    if a.x != b.x {
        a.x.partial_cmp(&b.x).unwrap()
    } else {
        a.y.partial_cmp(&b.y).unwrap()
    }
}

fn divide_segment(e: &Link, p: Coord, queue: &mut BinaryHeap<QueueItem>) {
    let other = e.borrow().other.clone().unwrap();
    let subject = e.borrow().subject;

    // right event for the left part
    let r = Event::new(p, subject);
    // left event for the right part
    let l = Event::new(p, subject);

    r.borrow_mut().left = false;
    r.borrow_mut().other = Some(e.clone());
    l.borrow_mut().left = true;
    l.borrow_mut().other = Some(other.clone());

    other.borrow_mut().other = Some(l.clone());
    e.borrow_mut().other = Some(r.clone());

    // keep left/right consistent if the split point ordering flipped
    if compare_events(&l, &other) == Ordering::Greater {
        other.borrow_mut().left = true;
        l.borrow_mut().left = false;
    }
    if compare_events(e, &r) == Ordering::Greater {
        // rare numerical case; leave as-is
    }

    queue.push(QueueItem(l));
    queue.push(QueueItem(r));
}

// ---- connecting result edges into contours ----------------------------------

fn connect_edges(sorted: &[Link]) -> Vec<Polygon> {
    // events that are part of the result
    let mut result: Vec<Link> = sorted
        .iter()
        .filter(|e| {
            let b = e.borrow();
            (b.left && b.in_result)
                || (!b.left && b.other.as_ref().unwrap().borrow().in_result)
        })
        .cloned()
        .collect();
    result.sort_by(compare_events);

    let n = result.len();
    if n == 0 {
        return Vec::new();
    }
    // index of each result event, then the index of each event's other endpoint
    let mut index_of = std::collections::HashMap::with_capacity(n);
    for (i, e) in result.iter().enumerate() {
        index_of.insert(Rc::as_ptr(e) as usize, i);
    }
    let other_pos: Vec<usize> = result
        .iter()
        .map(|e| {
            let other = e.borrow().other.clone().unwrap();
            index_of
                .get(&(Rc::as_ptr(&other) as usize))
                .copied()
                .unwrap_or(usize::MAX)
        })
        .collect();
    if other_pos.iter().any(|&p| p == usize::MAX) {
        // result events not paired up; bail (shouldn't happen for valid input)
        return Vec::new();
    }

    let mut visited = vec![false; n];
    let mut rings: Vec<Vec<Coord>> = Vec::new();

    for i in 0..n {
        if visited[i] {
            continue;
        }
        let mut ring: Vec<Coord> = vec![result[i].borrow().point];
        let mut pos = i;
        loop {
            visited[pos] = true;
            pos = other_pos[pos]; // far end of the current segment
            visited[pos] = true;
            ring.push(result[pos].borrow().point);
            pos = next_pos(&result, pos, &visited, i);
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
fn next_pos(result: &[Link], pos: usize, visited: &[bool], start: usize) -> usize {
    let n = result.len();
    let p0 = result[pos].borrow().point;
    let mut np = pos + 1;
    while np < n && result[np].borrow().point == p0 {
        if !visited[np] {
            return np;
        }
        np += 1;
    }
    np = pos;
    while np > 0 {
        np -= 1;
        if result[np].borrow().point != p0 {
            break;
        }
        if !visited[np] {
            return np;
        }
    }
    start
}

/// Turn raw closed rings into polygons: positive-area rings are shells, the rest
/// become holes of the shell that contains them.
fn assemble(rings: Vec<Vec<Coord>>) -> Vec<Polygon> {
    let mut shells: Vec<Vec<Coord>> = Vec::new();
    let mut holes: Vec<Vec<Coord>> = Vec::new();
    for r in rings {
        if ring_area2(&r) >= 0.0 {
            shells.push(r);
        } else {
            holes.push(r);
        }
    }
    let mut polys: Vec<Polygon> = shells
        .into_iter()
        .map(|s| Polygon::new(LineString::new(s), vec![]))
        .collect();

    for hole in holes {
        // attach to the first shell that contains the hole's first vertex
        let probe = hole[0];
        if let Some(poly) = polys
            .iter_mut()
            .find(|p| point_in_ring(probe, &p.exterior().0))
        {
            let (ext, mut ints) = std::mem::replace(poly, Polygon::new(LineString::new(vec![]), vec![])).into_inner();
            ints.push(LineString::new(hole));
            *poly = Polygon::new(ext, ints);
        } else {
            polys.push(Polygon::new(LineString::new(hole), vec![]));
        }
    }
    polys
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
