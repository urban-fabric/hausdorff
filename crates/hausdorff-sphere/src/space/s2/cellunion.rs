//! `s2::CellUnion` — a canonical region represented as a set of S2 cells,
//! ported from golang/geo's `cellunion.go`. The cell ids are kept *normalized*:
//! sorted, non-overlapping (no cell contains another), and with every group of
//! four complete siblings replaced by their parent, repeated up the hierarchy.
//! That canonical form makes equality a plain `Vec` comparison and lets the
//! containment/intersection queries run on sorted leaf-id ranges.

use hausdorff_base::earth::EARTH_RADIUS_M;

use super::measures::point_area;
use super::{S2CellId, S2Point};

/// A normalized set of S2 cells (see the module docs for the invariants).
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct CellUnion(Vec<S2CellId>);

impl From<Vec<S2CellId>> for CellUnion {
    /// Builds a normalized union from arbitrary cell ids: invalid ids are
    /// dropped; the rest are sorted, deduplicated, reduced to maximal cells
    /// (cells contained in another cell of the set are removed), and any four
    /// complete siblings are replaced by their parent, repeating up the
    /// hierarchy.
    ///
    /// Proof-backed by `riemann` `Hausdorff/CellNesting`: the id-range `contains`
    /// test implies the contained cell's whole range nests (`nesting`/`nesting'`,
    /// via the odd-multiple-of-`4^e` divisibility invariant), so dropping a
    /// contained cell loses no leaf coverage (`covers_of_contains`) — the
    /// soundness of the maximal-cell reduction.
    fn from(mut ids: Vec<S2CellId>) -> Self {
        ids.retain(|id| id.is_valid());
        ids.sort_unstable();
        let mut out: Vec<S2CellId> = Vec::with_capacity(ids.len());
        for mut id in ids {
            if out.last().is_some_and(|last| last.contains(id)) {
                continue;
            }
            while out.last().is_some_and(|&last| id.contains(last)) {
                out.pop();
            }
            while out.len() >= 3
                && are_siblings(
                    out[out.len() - 3],
                    out[out.len() - 2],
                    out[out.len() - 1],
                    id,
                )
            {
                out.truncate(out.len() - 3);
                id = id.immediate_parent();
            }
            out.push(id);
        }
        Self(out)
    }
}

/// True iff `a`, `b`, `c`, `d` are the four children of one parent cell, in
/// sorted order (golang/geo's `areSiblings`).
fn are_siblings(a: S2CellId, b: S2CellId, c: S2CellId, d: S2CellId) -> bool {
    if a.as_u64() ^ b.as_u64() ^ c.as_u64() != d.as_u64() {
        return false;
    }
    // Mask out the two child-position bits and everything below them; four
    // siblings agree on all remaining (parent) bits.
    let mut mask = (d.as_u64() & d.as_u64().wrapping_neg()) << 1;
    mask = !(mask + (mask << 1));
    let masked = d.as_u64() & mask;
    a.as_u64() & mask == masked
        && b.as_u64() & mask == masked
        && c.as_u64() & mask == masked
        && !d.is_face()
}

impl CellUnion {
    /// The normalized cell ids, sorted ascending.
    #[must_use]
    pub fn cell_ids(&self) -> &[S2CellId] {
        &self.0
    }

    /// Number of cells in the union.
    #[must_use]
    pub fn len(&self) -> usize {
        self.0.len()
    }

    /// True iff the union contains no cells.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    /// Iterates over the cells in sorted order.
    pub fn iter(&self) -> impl Iterator<Item = S2CellId> + '_ {
        self.0.iter().copied()
    }

    /// True iff the union fully contains `id` (i.e. some cell of the union is
    /// `id` or an ancestor of it). `O(log n)` on the sorted ids.
    #[must_use]
    pub fn contains_cell(&self, id: S2CellId) -> bool {
        let i = self.0.partition_point(|&c| c < id);
        if i < self.0.len() && self.0[i].range_min() <= id {
            return true;
        }
        i > 0 && self.0[i - 1].range_max() >= id
    }

    /// True iff the union and `id` share any cell (either contains the other).
    /// `O(log n)` on the sorted ids.
    #[must_use]
    pub fn intersects_cell(&self, id: S2CellId) -> bool {
        let i = self.0.partition_point(|&c| c < id);
        if i < self.0.len() && self.0[i].range_min() <= id.range_max() {
            return true;
        }
        i > 0 && self.0[i - 1].range_max() >= id.range_min()
    }

    /// True iff the union contains the unit-sphere point `p`, tested via the
    /// leaf cell containing `p`.
    #[must_use]
    pub fn contains_point(&self, p: S2Point) -> bool {
        self.contains_cell(S2CellId::from_point(p))
    }

    /// The normalized union of the two cell sets.
    #[must_use]
    pub fn union(&self, other: &CellUnion) -> CellUnion {
        let mut ids = self.0.clone();
        ids.extend_from_slice(&other.0);
        CellUnion::from(ids)
    }

    /// The normalized intersection of the two cell sets. Because both inputs
    /// are normalized, two intersecting cells always nest, so the result is the
    /// smaller cell of every nesting pair (a sorted merge, `O(n + m)`).
    #[must_use]
    pub fn intersection(&self, other: &CellUnion) -> CellUnion {
        let (a, b) = (&self.0, &other.0);
        let mut out = Vec::new();
        let (mut i, mut j) = (0, 0);
        while i < a.len() && j < b.len() {
            if a[i].range_max() < b[j].range_min() {
                i += 1;
            } else if b[j].range_max() < a[i].range_min() {
                j += 1;
            } else if b[j].contains(a[i]) {
                out.push(a[i]);
                i += 1;
            } else {
                out.push(b[j]);
                j += 1;
            }
        }
        CellUnion::from(out)
    }

    /// The normalized difference `self ∖ other`: every cell of `self` is
    /// recursively subdivided until each piece is either disjoint from `other`
    /// (kept) or contained in it (dropped).
    #[must_use]
    pub fn difference(&self, other: &CellUnion) -> CellUnion {
        let mut out = Vec::new();
        for &id in &self.0 {
            difference_cell(id, other, &mut out);
        }
        CellUnion::from(out)
    }

    /// The exact surface area of the union in square metres on the WGS84 mean
    /// sphere: the sum over the cells of the spherical area of their projected
    /// quadrilaterals, scaled by the Earth radius squared.
    #[must_use]
    pub fn area_m2(&self) -> f64 {
        let steradians: f64 = self.0.iter().map(|&id| cell_area_steradians(id)).sum();
        steradians * EARTH_RADIUS_M * EARTH_RADIUS_M
    }
}

/// Appends the cells of `id ∖ other` to `out` (golang/geo's
/// `cellUnionDifferenceInternal`).
fn difference_cell(id: S2CellId, other: &CellUnion, out: &mut Vec<S2CellId>) {
    if !other.intersects_cell(id) {
        out.push(id);
    } else if !other.contains_cell(id) {
        for child in id.children() {
            difference_cell(child, other, out);
        }
    }
}

/// Exact spherical area of a single cell in steradians: the cell's four
/// vertices split into two triangles.
fn cell_area_steradians(id: S2CellId) -> f64 {
    let [v0, v1, v2, v3] = id.vertices();
    point_area(v0, v1, v2) + point_area(v0, v2, v3)
}

#[cfg(test)]
mod tests {
    use std::f64::consts::PI;

    use super::*;
    use hausdorff_base::earth::GeoPoint;

    /// A small deterministic PRNG (SplitMix-style LCG) so the randomized
    /// invariant tests need no external crate and stay reproducible.
    struct Rng(u64);

    impl Rng {
        fn next_u64(&mut self) -> u64 {
            self.0 = self
                .0
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            self.0 >> 1
        }

        fn f64(&mut self) -> f64 {
            (self.next_u64() >> 10) as f64 / (1u64 << 53) as f64
        }

        fn cell(&mut self, min_level: u32, max_level: u32) -> S2CellId {
            let lat = self.f64() * 170.0 - 85.0;
            let lon = self.f64() * 360.0 - 180.0;
            let level = min_level + (self.next_u64() % u64::from(max_level - min_level + 1)) as u32;
            S2CellId::from_geo_point_at_level(GeoPoint::from_lat_lon(lat, lon), level)
        }

        fn cells(&mut self, n: usize) -> Vec<S2CellId> {
            (0..n).map(|_| self.cell(2, 12)).collect()
        }
    }

    fn assert_normalized(u: &CellUnion) {
        let ids = u.cell_ids();
        for w in ids.windows(2) {
            assert!(w[0] < w[1], "not sorted: {} !< {}", w[0], w[1]);
            assert!(
                w[0].range_max() < w[1].range_min(),
                "overlapping cells {} and {}",
                w[0],
                w[1]
            );
        }
        for w in ids.windows(4) {
            assert!(
                !are_siblings(w[0], w[1], w[2], w[3]),
                "four uncollapsed siblings starting at {}",
                w[0]
            );
        }
    }

    #[test]
    fn normalization_sorts_dedupes_and_drops_contained() {
        let cell = S2CellId::from_geo_point_at_level(GeoPoint::from_lat_lon(10.0, 20.0), 8);
        let child = cell.children()[2];
        let other = S2CellId::from_geo_point_at_level(GeoPoint::from_lat_lon(-30.0, 100.0), 9);
        let u = CellUnion::from(vec![other, cell, child, cell, child.children()[1]]);
        assert_eq!(
            u.cell_ids(),
            &[std::cmp::min(cell, other), std::cmp::max(cell, other)][..]
        );
        assert_normalized(&u);
    }

    #[test]
    fn four_siblings_collapse_to_parent_recursively() {
        let cell = S2CellId::from_geo_point_at_level(GeoPoint::from_lat_lon(48.0, 2.0), 10);
        let grandchildren: Vec<S2CellId> =
            cell.children().iter().flat_map(|c| c.children()).collect();
        assert_eq!(grandchildren.len(), 16);
        let u = CellUnion::from(grandchildren);
        assert_eq!(u.cell_ids(), &[cell][..]);
    }

    #[test]
    fn face_cells_do_not_collapse_and_invalid_ids_are_dropped() {
        let mut ids: Vec<S2CellId> = (0..6).map(S2CellId::from_face).collect();
        ids.push(S2CellId::from_u64(0));
        ids.push(S2CellId::from_u64(u64::MAX));
        let u = CellUnion::from(ids);
        assert_eq!(u.len(), 6);
        assert_normalized(&u);
        let sphere = 4.0 * PI * EARTH_RADIUS_M * EARTH_RADIUS_M;
        assert!(
            (u.area_m2() - sphere).abs() / sphere < 1e-9,
            "six faces should cover the sphere: {} vs {}",
            u.area_m2(),
            sphere
        );
    }

    #[test]
    fn cell_area_children_sum_to_parent() {
        let cell = S2CellId::from_geo_point_at_level(GeoPoint::from_lat_lon(35.0, 139.0), 5);
        let parent = CellUnion::from(vec![cell]).area_m2();
        let children = CellUnion::from(cell.children().to_vec()).area_m2();
        assert_eq!(CellUnion::from(cell.children().to_vec()).len(), 1);
        assert!(
            (parent - children).abs() / parent < 1e-9,
            "parent {parent} vs children {children}"
        );
    }

    #[test]
    fn randomized_union_with_own_children_normalizes_to_parent_set() {
        let mut rng = Rng(42);
        for _ in 0..20 {
            let cells = rng.cells(30);
            let base = CellUnion::from(cells.clone());
            let mut with_children = cells.clone();
            for c in &cells {
                with_children.extend(c.children());
            }
            let round_trip = CellUnion::from(with_children);
            assert_eq!(base, round_trip);
            assert_normalized(&base);
        }
    }

    #[test]
    fn randomized_set_operation_identities() {
        let mut rng = Rng(7);
        for _ in 0..20 {
            let a = CellUnion::from(rng.cells(25));
            let b = CellUnion::from(rng.cells(25));

            let a_or_b = a.union(&b);
            let a_and_b = a.intersection(&b);
            let a_minus_b = a.difference(&b);
            for u in [&a_or_b, &a_and_b, &a_minus_b] {
                assert_normalized(u);
            }

            for &c in a.cell_ids() {
                assert!(a_or_b.contains_cell(c), "A∪B must contain every cell of A");
            }
            for &c in b.cell_ids() {
                assert!(a_or_b.contains_cell(c), "A∪B must contain every cell of B");
            }
            for &c in a_and_b.cell_ids() {
                assert!(a.contains_cell(c), "A∩B must be inside A");
                assert!(b.contains_cell(c), "A∩B must be inside B");
            }
            assert_eq!(a_and_b, b.intersection(&a), "intersection commutes");
            for &c in a_minus_b.cell_ids() {
                assert!(a.contains_cell(c), "A∖B must be inside A");
                assert!(!b.intersects_cell(c), "(A∖B)∩B must be empty");
            }
            assert!(a_minus_b.intersection(&b).is_empty());
            assert_eq!(a_minus_b.union(&a_and_b).union(&a_or_b), a_or_b);
        }
    }

    #[test]
    fn contains_point_via_leaf_cell() {
        let inside = GeoPoint::from_lat_lon(52.52, 13.40);
        let cell = S2CellId::from_geo_point_at_level(inside, 8);
        let u = CellUnion::from(vec![cell]);
        assert!(u.contains_point(S2Point::from_geo_point(inside)));
        let outside = GeoPoint::from_lat_lon(-30.0, -60.0);
        assert!(!u.contains_point(S2Point::from_geo_point(outside)));
        assert!(u.contains_cell(cell.children()[3]));
        assert!(u.intersects_cell(cell.parent(4)));
        assert!(!u.contains_cell(cell.parent(4)));
    }

    #[test]
    fn empty_union_behaves() {
        let empty = CellUnion::from(vec![]);
        assert!(empty.is_empty());
        assert_eq!(empty.len(), 0);
        assert_eq!(empty.area_m2(), 0.0);
        let cell = S2CellId::from_face(2);
        assert!(!empty.contains_cell(cell));
        assert!(!empty.intersects_cell(cell));
        let u = CellUnion::from(vec![cell]);
        assert_eq!(empty.union(&u), u);
        assert!(empty.intersection(&u).is_empty());
        assert_eq!(u.difference(&empty), u);
        assert!(u.difference(&u).is_empty());
    }
}
