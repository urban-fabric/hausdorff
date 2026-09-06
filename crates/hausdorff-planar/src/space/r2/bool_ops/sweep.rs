use super::super::{Coord, CoordsIter, MultiPolygon, Polygon};
use super::event::{EdgeType, Event, compare_events, compare_segments, heap_pop, heap_push};
use super::event::{new_event, other_point};
use super::intersection::possible_intersection;
use super::{Operation, connect::connect_edges, ring_vertices};

pub(super) fn run(subject: &MultiPolygon, clipping: &MultiPolygon, op: Operation) -> MultiPolygon {
    let _span = crate::profile::span("bool_ops::run");
    use super::super::winding::{Direction, Orient};
    let subject = subject.orient(Direction::Default);
    let clipping = clipping.orient(Direction::Default);

    let cap = 2 * (subject.coords_count() + clipping.coords_count()) + 16;
    let mut events: Vec<Event> = Vec::with_capacity(cap);
    let mut queue: Vec<u32> = Vec::with_capacity(cap);
    {
        let _fill = crate::profile::span("fill");
        for p in subject.iter() {
            fill_polygon(&mut events, &mut queue, p, true);
        }
        for p in clipping.iter() {
            fill_polygon(&mut events, &mut queue, p, false);
        }
    }

    let sorted = {
        let _subdivide = crate::profile::span("subdivide");
        subdivide(&mut events, &mut queue, op)
    };
    let contours = {
        let _connect = crate::profile::span("connect_edges");
        connect_edges(&events, &sorted, op)
    };
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
        // The ring runs a -> b, so it leads away from `e1` and into `e2`.
        events[e1 as usize].ring_forward = true;
        events[e2 as usize].ring_forward = false;
        if compare_events(events, e1, e2) == std::cmp::Ordering::Less {
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

fn subdivide(events: &mut Vec<Event>, queue: &mut Vec<u32>, op: Operation) -> Vec<u32> {
    let mut sorted: Vec<u32> = Vec::new();
    let mut sweep: Vec<u32> = Vec::new();

    while let Some(event) = heap_pop(queue, events) {
        sorted.push(event);

        if events[event as usize].left {
            let pos = sweep_insert(&mut sweep, events, event);
            let prev = if pos > 0 { Some(sweep[pos - 1]) } else { None };
            let next = sweep.get(pos + 1).copied();

            compute_fields(events, event, prev, op);

            let mut divided = false;
            if let Some(next) = next {
                match possible_intersection(events, queue, event, next) {
                    0 => {}
                    2 => {
                        compute_fields(events, event, prev, op);
                        compute_fields(events, next, Some(event), op);
                    }
                    _ => divided = true,
                }
            }
            if let Some(prev) = prev {
                match possible_intersection(events, queue, prev, event) {
                    0 => {}
                    2 => {
                        let prev_prev = sweep_prev_of(&sweep, prev);
                        compute_fields(events, prev, prev_prev, op);
                        compute_fields(events, event, Some(prev), op);
                    }
                    _ => divided = true,
                }
            }

            // This segment was placed against a neighbour that has since been
            // divided, and it is the *undivided* neighbour it was compared with.
            // Where the division landed on this segment's own start — a vertex of
            // one polygon sitting on an edge of the other — the comparison was
            // decided by a determinant that is only zero once the split exists,
            // and until then reads as a couple of ULPs to one side. Placing it
            // again now that the neighbour ends where they actually meet is what
            // keeps the sweep's order and its subdivision telling the same story.
            if divided {
                if let Some(at) = sweep.iter().position(|&e| e == event) {
                    sweep.remove(at);
                }
                let pos = sweep_insert(&mut sweep, events, event);
                let prev = if pos > 0 { Some(sweep[pos - 1]) } else { None };
                let next = sweep.get(pos + 1).copied();
                compute_fields(events, event, prev, op);
                // Moving it gives it neighbours it was never tested against, and
                // an untested pair is a crossing the sweep will not split.
                if let Some(next) = next {
                    possible_intersection(events, queue, event, next);
                }
                if let Some(prev) = prev {
                    possible_intersection(events, queue, prev, event);
                }
            }
        } else {
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
    let pos =
        sweep.partition_point(|&e| compare_segments(events, e, event) == std::cmp::Ordering::Less);
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
