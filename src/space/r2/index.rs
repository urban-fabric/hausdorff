//! `r2::index` — a static **STR-tree** (Sort-Tile-Recursive bulk-loaded R-tree),
//! the spatial index behind fast range queries, spatial joins, and nearest-
//! neighbour search.
//!
//! This is the data structure professional engines (GEOS, PostGIS, `geo`'s
//! `rstar`) use to avoid the `O(n·m)` brute force of pairwise spatial predicates.
//! It indexes `(Rect, T)` items — bounding box plus an arbitrary payload (a
//! geometry, an id, a row) — built once and queried many times:
//!
//! - [`query`](StrTree::query): every item whose box overlaps a query box
//!   (range query / the candidate step of a spatial join).
//! - [`nearest`](StrTree::nearest): the item whose box is closest to a point
//!   (exact 1-NN for point data; a candidate for general geometries).
//!
//! STR packing sorts items into a near-square grid of tiles, giving a balanced,
//! cache-friendly tree built in `O(n log n)`.

use super::{Coord, Rect};

/// A static, bulk-loaded STR-tree over `(Rect, T)` items.
#[derive(Debug, Clone)]
pub struct StrTree<T> {
    root: Option<Node<T>>,
    len: usize,
    node_capacity: usize,
}

#[derive(Debug, Clone)]
enum Node<T> {
    Leaf { bbox: Rect, items: Vec<(Rect, T)> },
    Internal { bbox: Rect, children: Vec<Node<T>> },
}

impl<T> Node<T> {
    fn bbox(&self) -> Rect {
        match self {
            Node::Leaf { bbox, .. } | Node::Internal { bbox, .. } => *bbox,
        }
    }
}

impl<T> StrTree<T> {
    /// Compiled fallback node capacity, matching the bundled `hausdorff.toml`.
    pub const DEFAULT_NODE_CAPACITY: usize = 16;

    /// Bulk-load an index from `(bounding box, payload)` items, using the node
    /// capacity from [`hausdorff.toml`](crate::config) (overridable per machine).
    #[must_use]
    pub fn new(items: Vec<(Rect, T)>) -> Self {
        Self::with_node_capacity(items, crate::config::config().index.node_capacity)
    }

    /// Bulk-load with an explicit node capacity (`>= 1`).
    #[must_use]
    pub fn with_node_capacity(items: Vec<(Rect, T)>, node_capacity: usize) -> Self {
        let cap = node_capacity.max(1);
        let len = items.len();
        if items.is_empty() {
            return Self { root: None, len, node_capacity: cap };
        }
        // Bottom-up: pack items into leaves, then pack levels until one root.
        let mut level = str_pack(items, cap, |(r, _)| *r, |chunk| {
            Node::Leaf { bbox: union(chunk.iter().map(|(r, _)| *r)), items: chunk }
        });
        while level.len() > 1 {
            level = str_pack(level, cap, Node::bbox, |chunk| {
                Node::Internal { bbox: union(chunk.iter().map(Node::bbox)), children: chunk }
            });
        }
        Self { root: level.into_iter().next(), len, node_capacity: cap }
    }

    /// Number of indexed items.
    #[must_use]
    pub fn len(&self) -> usize {
        self.len
    }

    /// Whether the index holds no items.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.len == 0
    }

    /// The per-node item capacity the tree was built with.
    #[must_use]
    pub fn node_capacity(&self) -> usize {
        self.node_capacity
    }

    /// The payloads of every item whose bounding box overlaps `query` (touching
    /// counts) — the range-query / spatial-join candidate set.
    #[must_use]
    pub fn query(&self, query: Rect) -> Vec<&T> {
        let mut out = Vec::new();
        if let Some(root) = &self.root {
            query_node(root, query, &mut out);
        }
        out
    }

    /// Run `f` on every item whose box overlaps `query`, without allocating a
    /// result vector (for large candidate sets).
    pub fn query_each(&self, query: Rect, mut f: impl FnMut(&Rect, &T)) {
        if let Some(root) = &self.root {
            for_each_overlapping(root, query, &mut f);
        }
    }

    /// The payload whose bounding box is closest to `point` (by squared
    /// box-distance), or `None` if the index is empty. Exact 1-nearest-neighbour
    /// for point data; for areal items it is the nearest *by bounding box* — use
    /// it as the candidate, then refine with the exact geometry distance.
    #[must_use]
    pub fn nearest(&self, point: Coord) -> Option<&T> {
        let mut best: Option<(f64, &T)> = None;
        if let Some(root) = &self.root {
            nearest_node(root, point, &mut best);
        }
        best.map(|(_, t)| t)
    }
}

// ---- STR bulk packing ------------------------------------------------------

/// Sort-tile-recursive grouping: arrange `items` into a near-square grid of
/// `cap`-sized tiles (sort by box-centre x into vertical slices, then by y
/// within each slice) and fold each tile into a node via `group`.
fn str_pack<T, E>(
    mut items: Vec<E>,
    cap: usize,
    bbox: impl Fn(&E) -> Rect,
    mut group: impl FnMut(Vec<E>) -> Node<T>,
) -> Vec<Node<T>> {
    let n = items.len();
    let tile_count = n.div_ceil(cap);
    let slice_count = (tile_count as f64).sqrt().ceil() as usize;
    let slice_len = (slice_count * cap).max(cap);

    items.sort_by(|a, b| center_x(&bbox(a)).total_cmp(&center_x(&bbox(b))));
    let mut start = 0;
    while start < n {
        let end = (start + slice_len).min(n);
        items[start..end].sort_by(|a, b| center_y(&bbox(a)).total_cmp(&center_y(&bbox(b))));
        start = end;
    }

    let mut out = Vec::with_capacity(tile_count);
    let mut iter = items.into_iter();
    loop {
        let chunk: Vec<E> = iter.by_ref().take(cap).collect();
        if chunk.is_empty() {
            break;
        }
        out.push(group(chunk));
    }
    out
}

// ---- query / nearest traversal ---------------------------------------------

fn query_node<'a, T>(node: &'a Node<T>, q: Rect, out: &mut Vec<&'a T>) {
    if !rects_overlap(node.bbox(), q) {
        return;
    }
    match node {
        Node::Leaf { items, .. } => {
            for (r, t) in items {
                if rects_overlap(*r, q) {
                    out.push(t);
                }
            }
        }
        Node::Internal { children, .. } => {
            for c in children {
                query_node(c, q, out);
            }
        }
    }
}

fn for_each_overlapping<T>(node: &Node<T>, q: Rect, f: &mut impl FnMut(&Rect, &T)) {
    if !rects_overlap(node.bbox(), q) {
        return;
    }
    match node {
        Node::Leaf { items, .. } => {
            for (r, t) in items {
                if rects_overlap(*r, q) {
                    f(r, t);
                }
            }
        }
        Node::Internal { children, .. } => {
            for c in children {
                for_each_overlapping(c, q, f);
            }
        }
    }
}

fn nearest_node<'a, T>(node: &'a Node<T>, point: Coord, best: &mut Option<(f64, &'a T)>) {
    // Prune whole subtrees that can't beat the current best.
    if let Some((bd, _)) = best
        && point_rect_dist_sq(point, node.bbox()) >= *bd
    {
        return;
    }
    match node {
        Node::Leaf { items, .. } => {
            for (r, t) in items {
                let d = point_rect_dist_sq(point, *r);
                if best.is_none_or(|(bd, _)| d < bd) {
                    *best = Some((d, t));
                }
            }
        }
        Node::Internal { children, .. } => {
            // Descend the closest child first so the bound tightens early.
            let mut order: Vec<&Node<T>> = children.iter().collect();
            order.sort_by(|a, b| {
                point_rect_dist_sq(point, a.bbox()).total_cmp(&point_rect_dist_sq(point, b.bbox()))
            });
            for c in order {
                nearest_node(c, point, best);
            }
        }
    }
}

// ---- rectangle geometry ----------------------------------------------------

fn center_x(r: &Rect) -> f64 {
    (r.min().x + r.max().x) * 0.5
}

fn center_y(r: &Rect) -> f64 {
    (r.min().y + r.max().y) * 0.5
}

fn union(rects: impl Iterator<Item = Rect>) -> Rect {
    let mut lo = Coord::new(f64::INFINITY, f64::INFINITY);
    let mut hi = Coord::new(f64::NEG_INFINITY, f64::NEG_INFINITY);
    for r in rects {
        let (rl, rh) = (r.min(), r.max());
        lo = Coord::new(lo.x.min(rl.x), lo.y.min(rl.y));
        hi = Coord::new(hi.x.max(rh.x), hi.y.max(rh.y));
    }
    Rect::new(lo, hi)
}

fn rects_overlap(a: Rect, b: Rect) -> bool {
    a.min().x <= b.max().x
        && b.min().x <= a.max().x
        && a.min().y <= b.max().y
        && b.min().y <= a.max().y
}

fn point_rect_dist_sq(p: Coord, r: Rect) -> f64 {
    let dx = (r.min().x - p.x).max(p.x - r.max().x).max(0.0);
    let dy = (r.min().y - p.y).max(p.y - r.max().y).max(0.0);
    dx * dx + dy * dy
}
