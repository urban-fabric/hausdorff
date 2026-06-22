use std::cmp::Ordering;

use super::super::Coord;

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub(super) enum EdgeType {
    Normal,
    NonContributing,
    SameTransition,
    DifferentTransition,
}

/// A sweep event in the arena. `other`/`prev_in_result` are arena indices.
#[derive(Clone, Copy)]
pub(super) struct Event {
    pub(super) point: Coord,
    pub(super) left: bool,
    pub(super) other: Option<u32>,
    pub(super) subject: bool,
    pub(super) edge_type: EdgeType,
    pub(super) in_out: bool,
    pub(super) other_in_out: bool,
    pub(super) in_result: bool,
    pub(super) prev_in_result: Option<u32>,
}

/// Push a fresh event into the arena and return its index.
pub(super) fn new_event(events: &mut Vec<Event>, point: Coord, subject: bool) -> u32 {
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
pub(super) fn other_point(events: &[Event], e: u32) -> Coord {
    events[events[e as usize].other.unwrap() as usize].point
}

/// Is `p` below the segment `(a -> b)` (a is the left endpoint)? Uses the signed
/// area of the triangle.
fn signed_area(a: Coord, b: Coord, c: Coord) -> f64 {
    (a.x - c.x) * (b.y - c.y) - (b.x - c.x) * (a.y - c.y)
}

/// Order for the event queue: leftmost-lowest first; at a shared point, right
/// events before left events, then the lower segment first.
pub(super) fn compare_events(events: &[Event], e1: u32, e2: u32) -> Ordering {
    let a = &events[e1 as usize];
    let b = &events[e2 as usize];
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
    if a.left != b.left {
        return if a.left {
            Ordering::Greater
        } else {
            Ordering::Less
        };
    }
    let s = signed_area(a.point, other_point(events, e1), other_point(events, e2));
    if s != 0.0 {
        return if s > 0.0 {
            Ordering::Less
        } else {
            Ordering::Greater
        };
    }
    if a.subject != b.subject {
        return if a.subject {
            Ordering::Less
        } else {
            Ordering::Greater
        };
    }
    e1.cmp(&e2)
}

/// Push event index `idx` and sift it up to its place.
pub(super) fn heap_push(heap: &mut Vec<u32>, events: &[Event], idx: u32) {
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
pub(super) fn heap_pop(heap: &mut Vec<u32>, events: &[Event]) -> Option<u32> {
    if heap.is_empty() {
        return None;
    }
    let top = heap.swap_remove(0);
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
pub(super) fn compare_segments(events: &[Event], le1: u32, le2: u32) -> Ordering {
    if le1 == le2 {
        return Ordering::Equal;
    }
    let (p1, p1o) = (events[le1 as usize].point, other_point(events, le1));
    let (p2, p2o) = (events[le2 as usize].point, other_point(events, le2));

    if signed_area(p1, p1o, p2) != 0.0 || signed_area(p1, p1o, p2o) != 0.0 {
        if p1 == p2 {
            return if below(p1, p1o, p2o) {
                Ordering::Less
            } else {
                Ordering::Greater
            };
        }
        if compare_events(events, le1, le2) == Ordering::Greater {
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
    let (s1, s2) = (events[le1 as usize].subject, events[le2 as usize].subject);
    if s1 == s2 {
        if p1 == p2 {
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
