use super::super::coords::{MAX_LEVEL, MAX_SIZE};
use super::{S2CellId, size_ij};

impl S2CellId {
    /// The four cells adjacent across this cell's edges (down, right, up,
    /// left), all at this cell's level and all distinct.
    #[must_use]
    pub fn edge_neighbors(self) -> [Self; 4] {
        let level = self.level();
        let size = size_ij(level);
        let (f, i, j, _) = self.face_ij_orientation();
        [
            Self::from_face_ij_wrap(f, i, j - size).parent(level),
            Self::from_face_ij_wrap(f, i + size, j).parent(level),
            Self::from_face_ij_wrap(f, i, j + size).parent(level),
            Self::from_face_ij_wrap(f, i - size, j).parent(level),
        ]
    }

    /// The neighbors of the closest vertex to this cell at `level` (normally
    /// four, but three at the 8 cube vertices), including this cell's ancestor.
    #[must_use]
    pub fn vertex_neighbors(self, level: u32) -> Vec<Self> {
        let half_size = size_ij(level + 1);
        let size = half_size << 1;
        let (f, i, j, _) = self.face_ij_orientation();

        let (ioffset, isame) = if i & half_size != 0 {
            (size, i + size < MAX_SIZE)
        } else {
            (-size, i - size >= 0)
        };
        let (joffset, jsame) = if j & half_size != 0 {
            (size, j + size < MAX_SIZE)
        } else {
            (-size, j - size >= 0)
        };

        let mut results = vec![
            self.parent(level),
            Self::from_face_ij_same(f, i + ioffset, j, isame).parent(level),
            Self::from_face_ij_same(f, i, j + joffset, jsame).parent(level),
        ];
        if isame || jsame {
            results.push(
                Self::from_face_ij_same(f, i + ioffset, j + joffset, isame && jsame).parent(level),
            );
        }
        results
    }

    /// All neighbors of this cell at `level` (cells whose boundary intersects
    /// this cell's but whose interior does not). Returns empty if
    /// `level < self.level()` or `level > MAX_LEVEL`.
    #[must_use]
    pub fn all_neighbors(self, level: u32) -> Vec<Self> {
        let mut neighbors = Vec::new();
        if level < self.level() || level > MAX_LEVEL {
            return neighbors;
        }

        let (face, mut i, mut j, _) = self.face_ij_orientation();
        let size = size_ij(self.level());
        i &= -size;
        j &= -size;

        let nbr_size = size_ij(level);
        let mut k = -nbr_size;
        loop {
            let same_face;
            if k < 0 {
                same_face = j + k >= 0;
            } else if k >= size {
                same_face = j + k < MAX_SIZE;
            } else {
                same_face = true;
                neighbors.push(
                    Self::from_face_ij_same(face, i + k, j - nbr_size, j - size >= 0).parent(level),
                );
                neighbors.push(
                    Self::from_face_ij_same(face, i + k, j + size, j + size < MAX_SIZE)
                        .parent(level),
                );
            }
            neighbors.push(
                Self::from_face_ij_same(face, i - nbr_size, j + k, same_face && i - size >= 0)
                    .parent(level),
            );
            neighbors.push(
                Self::from_face_ij_same(face, i + size, j + k, same_face && i + size < MAX_SIZE)
                    .parent(level),
            );

            if k >= size {
                break;
            }
            k += nbr_size;
        }
        neighbors
    }
}
