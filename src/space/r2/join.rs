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

use super::{BoundingRect, Geometry, Rect, StrTree};

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
