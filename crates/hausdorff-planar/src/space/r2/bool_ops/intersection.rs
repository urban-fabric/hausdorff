use std::cmp::Ordering;

use super::super::Coord;
use super::event::{EdgeType, Event, compare_events, heap_push, new_event, other_point};

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

    let kross_e_va = e.x * va.y - e.y * va.x;
    if kross_e_va != 0.0 {
        return (0, a1, a1);
    }
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

pub(super) fn possible_intersection(
    events: &mut Vec<Event>,
    queue: &mut Vec<u32>,
    e1: u32,
    e2: u32,
) -> u8 {
    let a1 = events[e1 as usize].point;
    let a2 = other_point(events, e1);
    let b1 = events[e2 as usize].point;
    let b2 = other_point(events, e2);

    let (n, ip1, _ip2) = intersect_segments(a1, a2, b1, b2);

    if n == 0 {
        return 0;
    }
    if n == 1 && (a1 == b1 || a2 == b2 || a1 == b2 || a2 == b1) {
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

    if events[e1 as usize].subject == events[e2 as usize].subject {
        return 0;
    }

    let e1o = events[e1 as usize].other.unwrap();
    let e2o = events[e2 as usize].other.unwrap();
    let mut order: Vec<u32> = Vec::new();
    let left_coincide = a1 == b1;
    let right_coincide = a2 == b2;

    if left_coincide {
    } else if compare_events(events, e1, e2) == Ordering::Greater {
        order.push(e2);
        order.push(e1);
    } else {
        order.push(e1);
        order.push(e2);
    }

    if right_coincide {
    } else if compare_events(events, e1o, e2o) == Ordering::Greater {
        order.push(e2o);
        order.push(e1o);
    } else {
        order.push(e1o);
        order.push(e2o);
    }

    if left_coincide {
        events[e2 as usize].edge_type = EdgeType::NonContributing;
        let same = events[e2 as usize].in_out == events[e1 as usize].in_out;
        events[e1 as usize].edge_type = if same {
            EdgeType::SameTransition
        } else {
            EdgeType::DifferentTransition
        };
        if !right_coincide {
            let target = events[order[1] as usize].other.unwrap();
            let p = events[order[0] as usize].point;
            divide_segment(events, queue, target, p);
        }
        return 2;
    }

    if right_coincide {
        let p = events[order[1] as usize].point;
        divide_segment(events, queue, order[0], p);
        return 3;
    }

    let e0_other = events[order[0] as usize].other.unwrap();
    if e0_other != order[3] {
        let p1 = events[order[1] as usize].point;
        let p2 = events[order[2] as usize].point;
        divide_segment(events, queue, order[0], p1);
        divide_segment(events, queue, order[1], p2);
    } else {
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

    let r = new_event(events, p, subject);
    let l = new_event(events, p, subject);

    events[r as usize].left = false;
    events[r as usize].other = Some(e);
    events[l as usize].left = true;
    events[l as usize].other = Some(other);

    events[other as usize].other = Some(l);
    events[e as usize].other = Some(r);

    if compare_events(events, l, other) == Ordering::Greater {
        events[other as usize].left = true;
        events[l as usize].left = false;
    }

    heap_push(queue, events, l);
    heap_push(queue, events, r);
}
