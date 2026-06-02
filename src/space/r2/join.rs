//! `space::r2::join` — spatial joins over columns of geometries.
//!
//! A build-and-probe join on top of [`StrTree`]: index one side's bounding
//! boxes once, then for each geometry on the other side find the matches whose
//! boxes overlap. [`spatial_join_candidates`] returns those candidate
//! `(left, right)` index pairs — the broad phase a caller refines with an exact
//! predicate — while [`spatial_join`] folds the predicate into the probe and
//! returns only the pairs that satisfy it.
//!
//! `right` is the indexed (build) side: pass the larger or more-frequently-reused
//! side there. Pairs come out grouped by left index. The probe is sequential;
//! split the left column into chunks to parallelize across them downstream.

use super::{BoundingRect, Geometry, Intersects, Rect, StrTree};

/// Index `right`'s bounding boxes (skipping any empty geometry) for probing.
fn build_index(right: &[Geometry]) -> StrTree<usize> {
    let items: Vec<(Rect, usize)> = right
        .iter()
        .enumerate()
        .filter_map(|(j, g)| g.bounding_rect().map(|r| (r, j)))
        .collect();
    StrTree::new(items)
}

/// Broad-phase spatial join: every `(i, j)` such that `left[i]`'s bounding box
/// overlaps `right[j]`'s — the candidate set an exact predicate then refines.
/// Box-disjoint pairs are excluded (they cannot interact), so this is the same
/// set the `O(n·m)` nested-box-test loop would produce. Pairs are grouped by
/// left index.
#[must_use]
pub fn spatial_join_candidates(left: &[Geometry], right: &[Geometry]) -> Vec<(usize, usize)> {
    let tree = build_index(right);
    let mut out = Vec::new();
    for (i, lg) in left.iter().enumerate() {
        if let Some(q) = lg.bounding_rect() {
            tree.query_each(q, |_, &j| out.push((i, j)));
        }
    }
    out
}

/// Spatial join with an exact `predicate`: every `(i, j)` whose bounding boxes
/// overlap **and** for which `predicate(&left[i], &right[j])` holds. The box
/// filter is the broad phase, so the predicate runs only on candidate pairs, not
/// all `n·m`. Pass e.g. `|a, b| a.intersects(b)`, an `ST_DWithin` test, or a
/// `relate`-based predicate. Pairs are grouped by left index.
#[must_use]
pub fn spatial_join<F>(left: &[Geometry], right: &[Geometry], predicate: F) -> Vec<(usize, usize)>
where
    F: Fn(&Geometry, &Geometry) -> bool,
{
    let tree = build_index(right);
    let mut out = Vec::new();
    for (i, lg) in left.iter().enumerate() {
        if let Some(q) = lg.bounding_rect() {
            tree.query_each(q, |_, &j| {
                if predicate(lg, &right[j]) {
                    out.push((i, j));
                }
            });
        }
    }
    out
}

/// One indexed geometry: its bounding box, the geometry, and a caller payload.
struct Indexed<P> {
    rect: Rect,
    geom: Geometry,
    payload: P,
}

/// A reusable spatial index over payload-tagged geometries — **build once, probe
/// many.** Unlike [`spatial_join`] (which rebuilds the tree each call over
/// `&[Geometry]` and returns positional indices), this keeps the bulk-loaded
/// [`StrTree`] and carries an arbitrary payload per geometry (e.g. a source row
/// id), so a broadcast join can index the small side once and stream batches of
/// the large side against it. Empty geometries (no bounding box) are dropped at
/// build, so positions need not line up — hence the payloads.
pub struct GeometryIndex<P> {
    entries: Vec<Indexed<P>>,
    tree: StrTree<u32>,
}

impl<P: Copy> GeometryIndex<P> {
    /// Index `(geometry, payload)` items; geometries with no bounding box are
    /// dropped.
    pub fn build(items: impl IntoIterator<Item = (Geometry, P)>) -> Self {
        let entries: Vec<Indexed<P>> = items
            .into_iter()
            .filter_map(|(geom, payload)| {
                geom.bounding_rect().map(|rect| Indexed { rect, geom, payload })
            })
            .collect();
        let tree =
            StrTree::new(entries.iter().enumerate().map(|(i, e)| (e.rect, i as u32)).collect());
        Self { entries, tree }
    }

    /// Number of indexed geometries.
    #[must_use]
    pub fn len(&self) -> usize {
        self.entries.len()
    }

    /// Whether the index holds no geometries.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    /// `(probe_payload, indexed_payload)` for every probe geometry whose bounding
    /// box overlaps an indexed geometry **and** for which `predicate(probe,
    /// indexed)` holds. The index envelopes are the broad phase; `predicate` is
    /// the exact refine.
    #[must_use]
    pub fn join<Q, F>(&self, probe: &[(Geometry, Q)], predicate: F) -> Vec<(Q, P)>
    where
        Q: Copy,
        F: Fn(&Geometry, &Geometry) -> bool,
    {
        let mut out = Vec::new();
        for (pgeom, pq) in probe {
            if let Some(q) = pgeom.bounding_rect() {
                self.tree.query_each(q, |_, &i| {
                    let e = &self.entries[i as usize];
                    if predicate(pgeom, &e.geom) {
                        out.push((*pq, e.payload));
                    }
                });
            }
        }
        out
    }
}

/// Nested-loop spatial join (no index) over payload-tagged geometries: every
/// `(left_payload, right_payload)` whose bounding boxes overlap **and** for which
/// `predicate(left, right)` holds. Produces the same candidate set as the
/// [`StrTree`]-indexed join — for small or lightly-probed inputs, building the
/// tree does not pay off and this brute-force scan is cheaper. Empty geometries
/// (no bounding box) are skipped.
#[must_use]
pub fn nested_loop_join<P, Q, F>(
    left: &[(Geometry, P)],
    right: &[(Geometry, Q)],
    predicate: F,
) -> Vec<(P, Q)>
where
    P: Copy,
    Q: Copy,
    F: Fn(&Geometry, &Geometry) -> bool,
{
    let right_rects: Vec<Option<Rect>> = right.iter().map(|(g, _)| g.bounding_rect()).collect();
    let mut out = Vec::new();
    for (lg, lp) in left {
        let Some(la) = lg.bounding_rect() else {
            continue;
        };
        for ((rg, rq), rb) in right.iter().zip(&right_rects) {
            let Some(rb) = rb else { continue };
            if la.intersects(rb) && predicate(lg, rg) {
                out.push((*lp, *rq));
            }
        }
    }
    out
}

#[cfg(test)]
mod payload_join_tests {
    use super::*;
    use crate::space::r2::{Coord, LineString, Point, Polygon};

    fn square(x: f64, y: f64) -> Geometry {
        Geometry::Polygon(Polygon::new(
            LineString::from(vec![(x, y), (x + 1.0, y), (x + 1.0, y + 1.0), (x, y + 1.0), (x, y)]),
            vec![],
        ))
    }

    #[test]
    fn index_and_nested_loop_agree_with_payloads() {
        // Build side tagged with row ids 10, 20, 30; probe side with 1, 2.
        let build =
            vec![(square(0.0, 0.0), 10u32), (square(5.0, 5.0), 20), (square(-0.5, -0.5), 30)];
        let probe = vec![
            (Geometry::Point(Point(Coord::new(0.4, 0.4))), 1u32),
            (Geometry::Point(Point(Coord::new(100.0, 100.0))), 2),
        ];
        let pred = |a: &Geometry, b: &Geometry| crate::space::r2::SpatialPredicate::Intersects.eval(a, b);

        let idx = GeometryIndex::build(build.clone());
        assert_eq!(idx.len(), 3);
        let mut indexed: Vec<(u32, u32)> = idx.join(&probe, pred);
        // Nested-loop with sides swapped to match (probe_payload, build_payload).
        let mut nested: Vec<(u32, u32)> = nested_loop_join(&probe, &build, pred);
        indexed.sort_unstable();
        nested.sort_unstable();
        assert_eq!(indexed, nested);
        // Point (0.4,0.4) is inside squares 10 and 30; point 2 hits nothing.
        assert_eq!(indexed, vec![(1, 10), (1, 30)]);
    }
}
