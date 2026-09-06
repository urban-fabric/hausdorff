//! Splitting a bag of segments into an arrangement: pieces that meet only at
//! shared endpoints and never cross.
//!
//! Every endpoint the arrangement produces is drawn from one shared pool, so
//! two pieces that meet do so at *bit-identical* coordinates. The tracing phase
//! keys vertices on those bits, and nothing downstream has to carry a
//! tolerance.

use std::collections::HashMap;

use super::super::Coord;

/// Vertices closer than this, relative to their own magnitude, are the same
/// vertex. Chosen to sit above intersection round-off and far below any
/// distance a caller could intend — the same reasoning as the boolean sweep's
/// crossing tolerance.
const TOL_ULPS: f64 = 8.0;

fn same(a: Coord, b: Coord) -> bool {
    let close = |u: f64, v: f64| {
        let tol = TOL_ULPS * f64::EPSILON * u.abs().max(v.abs()).max(1.0);
        (u - v).abs() <= tol
    };
    close(a.x, b.x) && close(a.y, b.y)
}

/// The shared vertex pool. Points are looked up in a coarse grid so snapping
/// stays linear rather than comparing against every vertex already seen.
pub(super) struct Nodes {
    points: Vec<Coord>,
    grid: HashMap<(i64, i64), Vec<usize>>,
    cell: f64,
}

impl Nodes {
    pub(super) fn new(extent: f64) -> Self {
        // A cell far larger than the snapping tolerance but small enough that
        // buckets stay short; only the 9 cells around a query are searched.
        let cell = (extent * 1e-6).max(f64::MIN_POSITIVE);
        Nodes {
            points: Vec::new(),
            grid: HashMap::new(),
            cell,
        }
    }

    fn bucket(&self, c: Coord) -> (i64, i64) {
        (
            (c.x / self.cell).floor() as i64,
            (c.y / self.cell).floor() as i64,
        )
    }

    /// The index of the pooled vertex equal to `c`, inserting it if new.
    pub(super) fn intern(&mut self, c: Coord) -> usize {
        let (bx, by) = self.bucket(c);
        for dx in -1..=1 {
            for dy in -1..=1 {
                if let Some(bucket) = self.grid.get(&(bx + dx, by + dy)) {
                    for &i in bucket {
                        if same(self.points[i], c) {
                            return i;
                        }
                    }
                }
            }
        }
        let i = self.points.len();
        self.points.push(c);
        self.grid.entry((bx, by)).or_default().push(i);
        i
    }

    pub(super) fn coord(&self, i: usize) -> Coord {
        self.points[i]
    }
}

fn cross(o: Coord, a: Coord, b: Coord) -> f64 {
    (a.x - o.x) * (b.y - o.y) - (a.y - o.y) * (b.x - o.x)
}

/// Where along `a -> b` the point `p` falls, as a fraction. `p` is assumed
/// collinear with the segment.
fn param(a: Coord, b: Coord, p: Coord) -> f64 {
    let (dx, dy) = (b.x - a.x, b.y - a.y);
    if dx.abs() >= dy.abs() {
        if dx == 0.0 { 0.0 } else { (p.x - a.x) / dx }
    } else if dy == 0.0 {
        0.0
    } else {
        (p.y - a.y) / dy
    }
}

/// Split points on `i`, expressed as parameters along it, from its interaction
/// with `j`.
fn splits_between(i: (Coord, Coord), j: (Coord, Coord)) -> Vec<f64> {
    let (a, b) = i;
    let (c, d) = j;
    let r = Coord::new(b.x - a.x, b.y - a.y);
    let s = Coord::new(d.x - c.x, d.y - c.y);
    let denom = r.x * s.y - r.y * s.x;

    // Collinear: the overlap's endpoints split this segment.
    if denom == 0.0 {
        if cross(a, b, c) != 0.0 {
            return Vec::new();
        }
        let mut out = Vec::new();
        for p in [c, d] {
            let t = param(a, b, p);
            if (0.0..=1.0).contains(&t) {
                out.push(t);
            }
        }
        return out;
    }

    let qp = Coord::new(c.x - a.x, c.y - a.y);
    let t = (qp.x * s.y - qp.y * s.x) / denom;
    let u = (qp.x * r.y - qp.y * r.x) / denom;
    if (0.0..=1.0).contains(&t) && (0.0..=1.0).contains(&u) {
        vec![t]
    } else {
        Vec::new()
    }
}

/// Cut `segments` into non-crossing pieces on a shared vertex pool.
///
/// Returns the pieces as index pairs into `nodes`. Pieces whose ends intern to
/// the same vertex are dropped: they carry no length and would only add
/// degenerate edges to the trace.
pub(super) fn arrange(segments: &[(Coord, Coord)], nodes: &mut Nodes) -> Vec<(usize, usize)> {
    let n = segments.len();
    let mut cuts: Vec<Vec<f64>> = vec![Vec::new(); n];
    for i in 0..n {
        for j in 0..n {
            if i == j {
                continue;
            }
            cuts[i].extend(splits_between(segments[i], segments[j]));
        }
    }

    let mut out = Vec::new();
    for (i, &(a, b)) in segments.iter().enumerate() {
        let mut ts = std::mem::take(&mut cuts[i]);
        ts.push(0.0);
        ts.push(1.0);
        ts.retain(|t| t.is_finite() && (0.0..=1.0).contains(t));
        ts.sort_by(f64::total_cmp);

        let mut prev: Option<usize> = None;
        for t in ts {
            let p = Coord::new(a.x + (b.x - a.x) * t, a.y + (b.y - a.y) * t);
            let v = nodes.intern(p);
            if let Some(u) = prev
                && u != v
            {
                out.push((u, v));
            }
            if prev.is_none_or(|u| u != v) {
                prev = Some(v);
            }
        }
    }
    out
}
