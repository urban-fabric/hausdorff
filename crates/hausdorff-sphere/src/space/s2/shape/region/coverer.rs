use super::*;

/// Approximates a [`Region`] by a set of S2 cells whose union contains it.
///
/// The covering is built top-down from the six face cells: a cell that does not
/// intersect the region is dropped; a cell at `max_level` — or at/above
/// `min_level` and fully inside the region — is emitted; otherwise it is
/// subdivided into its four children. `max_cells` is a **soft** budget: once
/// reached, cells at/above `min_level` are emitted coarsely instead of refined,
/// so the result stays a valid (super-)cover rather than honouring the budget
/// exactly (unlike golang/geo's priority-queue coverer).
///
/// [`RegionCoverer::adaptive_covering`] and
/// [`RegionCoverer::interior_covering`] run the golang/geo-style priority-queue
/// refinement instead, producing multi-level normalized [`CellUnion`]s under
/// the same three parameters.
#[derive(Clone, Copy, Debug)]
pub struct RegionCoverer {
    /// Coarsest level a covering cell may have.
    pub min_level: u32,
    /// Finest level a covering cell may have (refinement stops here). Values
    /// above the leaf level (30) are treated as 30.
    pub max_level: u32,
    /// Soft cap on the number of cells.
    pub max_cells: usize,
}

impl Default for RegionCoverer {
    fn default() -> Self {
        Self {
            min_level: 0,
            max_level: MAX_LEVEL,
            max_cells: 8,
        }
    }
}

impl RegionCoverer {
    /// A coverer that emits cells at exactly `level`: every level-`level` cell
    /// that intersects the region, with no coarsening and no cell budget. This is
    /// the fixed-resolution covering a partitioner uses to bucket a geometry.
    /// `level` is clamped to the leaf level (30).
    #[must_use]
    pub fn at_level(level: u32) -> Self {
        let level = level.min(MAX_LEVEL);
        Self {
            min_level: level,
            max_level: level,
            max_cells: usize::MAX,
        }
    }

    /// The covering cell ids for `region`.
    #[must_use]
    pub fn covering<R: Region + ?Sized>(&self, region: &R) -> Vec<S2CellId> {
        let mut out = Vec::new();
        for face in 0..6u8 {
            self.cover_cell(region, S2CellId::from_face(face), &mut out);
        }
        out
    }

    fn cover_cell<R: Region + ?Sized>(&self, region: &R, cell: S2CellId, out: &mut Vec<S2CellId>) {
        if !region.intersects_cell(cell) {
            return;
        }
        let level = cell.level();
        let max_level = self.max_level.min(MAX_LEVEL);
        let at_budget = out.len() >= self.max_cells;
        if level >= max_level
            || (level >= self.min_level && (at_budget || region.contains_cell(cell)))
        {
            out.push(cell);
            return;
        }
        for child in cell.children() {
            self.cover_cell(region, child, out);
        }
    }

    /// An adaptive (multi-level) covering of `region`, golang/geo-style:
    /// starting from the six face cells, the intersecting cell whose children
    /// intersect the region the most is repeatedly subdivided until expanding
    /// it would exceed the `max_cells` budget, at which point cells are
    /// emitted as-is. Cells fully inside the region (at or below `min_level`)
    /// are emitted without further subdivision; refinement never passes
    /// `max_level`. Returns a normalized [`CellUnion`].
    ///
    /// The budget is honoured by the generation step but is still soft at the
    /// margins: covering more than `max_cells` intersecting face cells, or a
    /// `min_level` that forces subdivision, can exceed it. Normalization may
    /// also merge four complete siblings emitted at `min_level` into their
    /// parent, so a result cell can sit one level above `min_level`.
    #[must_use]
    pub fn adaptive_covering<R: Region + ?Sized>(&self, region: &R) -> CellUnion {
        self.adaptive(region, false)
    }

    /// An adaptive *interior* covering of `region`: only cells fully contained
    /// in the region are emitted, so the union of the result is a subset of the
    /// region (the dual of [`RegionCoverer::adaptive_covering`]). Cells that
    /// merely intersect are subdivided down to `max_level` and dropped if they
    /// never become contained; the result can be empty for thin regions or
    /// tight budgets. Returns a normalized [`CellUnion`].
    #[must_use]
    pub fn interior_covering<R: Region + ?Sized>(&self, region: &R) -> CellUnion {
        self.adaptive(region, true)
    }

    fn adaptive<R: Region + ?Sized>(&self, region: &R, interior: bool) -> CellUnion {
        let min_level = self.min_level.min(MAX_LEVEL);
        let max_level = self.max_level.min(MAX_LEVEL).max(min_level);
        let mut result: Vec<S2CellId> = Vec::new();
        let mut heap: BinaryHeap<Candidate> = BinaryHeap::new();

        for face in 0..6u8 {
            add_candidate(
                region,
                S2CellId::from_face(face),
                min_level,
                max_level,
                interior,
                &mut result,
                &mut heap,
            );
        }
        while let Some(candidate) = heap.pop() {
            let expand = candidate.cell.level() < min_level
                || result.len() + heap.len() + candidate.children.len() <= self.max_cells;
            if expand {
                for child in candidate.children {
                    add_candidate(
                        region,
                        child,
                        min_level,
                        max_level,
                        interior,
                        &mut result,
                        &mut heap,
                    );
                }
            } else if !interior {
                result.push(candidate.cell);
            }
        }
        CellUnion::from(result)
    }
}

/// A cell awaiting subdivision in the adaptive coverer, with its
/// region-intersecting children precomputed.
struct Candidate {
    cell: S2CellId,
    children: Vec<S2CellId>,
}

impl Candidate {
    /// Subdivision priority: most intersecting children first, then coarser
    /// cells, then the smaller id for determinism.
    fn key(&self) -> (usize, std::cmp::Reverse<u32>, std::cmp::Reverse<u64>) {
        (
            self.children.len(),
            std::cmp::Reverse(self.cell.level()),
            std::cmp::Reverse(self.cell.as_u64()),
        )
    }
}

impl PartialEq for Candidate {
    fn eq(&self, other: &Self) -> bool {
        self.key() == other.key()
    }
}

impl Eq for Candidate {}

impl PartialOrd for Candidate {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Candidate {
    fn cmp(&self, other: &Self) -> Ordering {
        self.key().cmp(&other.key())
    }
}

/// Routes `cell` for the adaptive coverer: dropped if it misses the region,
/// emitted if terminal (fully contained at/below `min_level`, or — for a
/// non-interior covering — at `max_level`), otherwise queued with its
/// intersecting children.
fn add_candidate<R: Region + ?Sized>(
    region: &R,
    cell: S2CellId,
    min_level: u32,
    max_level: u32,
    interior: bool,
    result: &mut Vec<S2CellId>,
    heap: &mut BinaryHeap<Candidate>,
) {
    if !region.intersects_cell(cell) {
        return;
    }
    let level = cell.level();
    if level >= min_level && region.contains_cell(cell) {
        result.push(cell);
        return;
    }
    if level >= max_level {
        if !interior {
            result.push(cell);
        }
        return;
    }
    let children: Vec<S2CellId> = cell
        .children()
        .iter()
        .copied()
        .filter(|&c| region.intersects_cell(c))
        .collect();
    if !children.is_empty() {
        heap.push(Candidate { cell, children });
    }
}
