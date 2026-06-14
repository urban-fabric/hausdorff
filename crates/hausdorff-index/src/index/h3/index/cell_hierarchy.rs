//! [`H3Cell`] hierarchy operations: parent/child navigation
//! (`parent`, `children_at_resolution`, `children_count`, `center_child`,
//! `child_position`, `child_at`) and set compaction
//! (`compact_cells`, `uncompact_size`, `uncompact_cells`).

use std::collections::BTreeSet;

use super::*;

impl H3Cell {
    /// The parent cell at `resolution` (coarser-or-equal), or `None`.
    pub fn parent(self, resolution: u8) -> Option<Self> {
        (resolution <= self.resolution())
            .then(|| Self(set_unused(set_resolution(self.0, resolution), resolution)))
    }

    /// All descendant cells at `resolution` (finer-or-equal); empty if coarser
    /// or out of range.
    pub fn children_at_resolution(self, resolution: u8) -> Vec<Self> {
        if resolution < self.resolution() || resolution > MAX_RESOLUTION {
            return Vec::new();
        }
        children::Children::new(self.0, self.resolution(), self.is_pentagon(), resolution)
            .map(Self)
            .collect()
    }

    /// Number of descendant cells at `resolution`.
    pub fn children_count(self, resolution: u8) -> u64 {
        if resolution > MAX_RESOLUTION {
            return 0;
        }
        children::children_count(self.resolution(), resolution, self.is_pentagon())
    }

    /// The center descendant at `resolution`.
    pub fn center_child(self, resolution: u8) -> Option<Self> {
        if resolution < self.resolution() || resolution > MAX_RESOLUTION {
            return None;
        }
        children::Children::new(self.0, self.resolution(), self.is_pentagon(), resolution)
            .next()
            .map(Self)
    }

    /// Position of this cell within its `parent_resolution` parent's ordered
    /// child list. Returns `None` when `parent_resolution` is finer than this
    /// cell.
    pub fn child_position(self, parent_resolution: u8) -> Option<u64> {
        if parent_resolution > self.resolution() || parent_resolution > MAX_RESOLUTION {
            return None;
        }

        let parent_is_pentagon = self.parent(parent_resolution)?.is_pentagon();
        let mut position = 0u64;

        if parent_is_pentagon {
            for res in (parent_resolution + 1)..=self.resolution() {
                let parent = self.parent(res - 1).expect("valid parent");
                let parent_is_pentagon = parent.is_pentagon();
                let mut digit = get_direction(self.0, res);

                if parent_is_pentagon && digit > 0 {
                    digit -= 1;
                }
                if digit == 0 {
                    continue;
                }

                let diff = self.resolution() - res;
                let hex_count = children::hexagon_children_count(diff);
                let count0 = if parent_is_pentagon {
                    children::pentagon_children_count(diff)
                } else {
                    hex_count
                };
                position += u64::from(digit - 1) * hex_count + count0;
            }
        } else {
            for res in (parent_resolution + 1)..=self.resolution() {
                let diff = self.resolution() - res;
                let hex_count = children::hexagon_children_count(diff);
                let digit = get_direction(self.0, res);
                position += u64::from(digit) * hex_count;
            }
        }

        Some(position)
    }

    /// Child cell at `position` within this cell's ordered children at
    /// `resolution`.
    pub fn child_at(self, mut position: u64, resolution: u8) -> Option<Self> {
        if resolution < self.resolution()
            || resolution > MAX_RESOLUTION
            || position >= self.children_count(resolution)
        {
            return None;
        }

        let mut child = set_resolution(self.0, resolution);
        let mut current_res = self.resolution();

        if self.is_pentagon() {
            for res in (self.resolution() + 1)..=resolution {
                current_res = res;
                let diff = resolution - res;
                let pentagon_count = children::pentagon_children_count(diff);
                if position < pentagon_count {
                    child = set_direction(child, 0, res);
                } else {
                    let count = children::hexagon_children_count(diff);
                    position -= pentagon_count;
                    child = set_direction(child, (position / count) as u8 + 2, res);
                    position %= count;
                    break;
                }
            }
        }

        for res in (current_res + 1)..=resolution {
            let diff = resolution - res;
            let count = children::hexagon_children_count(diff);
            child = set_direction(child, (position / count) as u8, res);
            position %= count;
        }

        Self::try_new(child)
    }

    /// Compresses a unique, same-resolution cell set by replacing full child
    /// branches with their parent cells.
    pub fn compact_cells(cells: &mut Vec<Self>) -> Result<(), H3CompactionError> {
        let Some(first) = cells.first() else {
            return Ok(());
        };
        let resolution = first.resolution();
        if cells.iter().any(|cell| cell.resolution() != resolution) {
            return Err(H3CompactionError::HeterogeneousResolution);
        }

        cells.sort_unstable();
        if cells.windows(2).any(|window| window[0] == window[1]) {
            return Err(H3CompactionError::DuplicateInput);
        }

        let mut set = cells.iter().copied().collect::<BTreeSet<_>>();
        for res in (1..=resolution).rev() {
            let parents = set
                .iter()
                .filter(|cell| cell.resolution() == res)
                .filter_map(|cell| cell.parent(res - 1))
                .collect::<BTreeSet<_>>();

            for parent in parents {
                let children = parent.children_at_resolution(res);
                if children.iter().all(|child| set.contains(child)) {
                    for child in children {
                        set.remove(&child);
                    }
                    set.insert(parent);
                }
            }
        }

        *cells = set.into_iter().collect();
        Ok(())
    }

    /// Exact output size for [`H3Cell::uncompact_cells`].
    pub fn uncompact_size(
        compacted: impl IntoIterator<Item = Self>,
        resolution: u8,
    ) -> Option<u64> {
        if resolution > MAX_RESOLUTION {
            return None;
        }
        let mut size = 0u64;
        for cell in compacted {
            if cell.resolution() > resolution {
                return None;
            }
            size = size.checked_add(cell.children_count(resolution))?;
        }
        Some(size)
    }

    /// Expands compacted cells to the requested `resolution`.
    pub fn uncompact_cells(
        compacted: impl IntoIterator<Item = Self>,
        resolution: u8,
    ) -> Option<Vec<Self>> {
        let compacted = compacted.into_iter().collect::<Vec<_>>();
        let size = Self::uncompact_size(compacted.iter().copied(), resolution)?;
        let mut out = Vec::with_capacity(usize::try_from(size).ok()?);
        for cell in compacted {
            if cell.resolution() > resolution {
                return None;
            }
            out.extend(cell.children_at_resolution(resolution));
        }
        Some(out)
    }
}
