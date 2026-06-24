use super::super::{Coord, Polygon};
use super::assemble::assemble;
use super::event::{Event, compare_events};

pub(super) fn connect_edges(events: &[Event], sorted: &[u32]) -> Vec<Polygon> {
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
            pos = other_pos[pos];
            visited[pos] = true;
            ring.push(events[result[pos] as usize].point);
            pos = next_pos(events, &result, pos, &visited, i);
            if pos == i || pos >= n {
                break;
            }
        }
        if ring.len() >= 4 {
            ring.pop();
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
