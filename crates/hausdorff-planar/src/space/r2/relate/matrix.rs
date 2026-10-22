//! The DE-9IM intersection matrix and the OGC predicate reads over it. The
//! algorithm in [`super`] fills the matrix; the predicates here interpret it.

/// Location of a point relative to a geometry.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub(super) enum Loc {
    Interior = 0,
    Boundary = 1,
    Exterior = 2,
}

/// A DE-9IM intersection matrix: `m[a][b]` is the dimension (`-1` = none, else
/// `0`/`1`/`2`) of `a`-part(A) ∩ `b`-part(B), with rows/cols ordered
/// Interior, Boundary, Exterior.
///
/// Proof-backed by `riemann` `Hausdorff/Relate`: the DE-9IM predicate algebra
/// over the 3×3 matrix — `intersects_iff_not_disjoint`, the operand-swap
/// dualities (`contains_iff_transpose_within`, `covers_iff_transpose_coveredby`),
/// symmetry of `disjoint`/`intersects`/`equalTopo`/`touches`/`overlaps`, and the
/// implication lattice (`equal ⟹ within ∧ contains`, `within ⟹ coveredBy ⟹
/// intersects`, `touches ⟹ intersects`). The matrix-filling noding stays
/// solver-side.
#[derive(Clone, Copy, Debug)]
pub struct IntersectionMatrix {
    m: [[i8; 3]; 3],
}

pub(super) const I: usize = 0;
pub(super) const B: usize = 1;
pub(super) const E: usize = 2;

impl IntersectionMatrix {
    pub(super) fn empty() -> Self {
        Self { m: [[-1; 3]; 3] }
    }

    #[inline]
    pub(super) fn set_at_least(&mut self, a: usize, b: usize, dim: i8) {
        if dim > self.m[a][b] {
            self.m[a][b] = dim;
        }
    }

    #[inline]
    fn get(&self, a: usize, b: usize) -> i8 {
        self.m[a][b]
    }

    #[inline]
    fn is_t(&self, a: usize, b: usize) -> bool {
        self.m[a][b] >= 0
    }
    #[inline]
    fn is_f(&self, a: usize, b: usize) -> bool {
        self.m[a][b] < 0
    }

    /// Dimension (`-1` = none, else `0`/`1`/`2`) of Interior(A) ∩ Interior(B).
    #[must_use]
    pub fn interior_interior_dim(&self) -> i8 {
        self.m[I][I]
    }

    /// Dimension (`-1` = none, else `0`/`1`/`2`) of Boundary(A) ∩ Boundary(B).
    #[must_use]
    pub fn boundary_boundary_dim(&self) -> i8 {
        self.m[B][B]
    }

    #[must_use]
    pub fn is_disjoint(&self) -> bool {
        self.is_f(I, I) && self.is_f(I, B) && self.is_f(B, I) && self.is_f(B, B)
    }

    #[must_use]
    pub fn is_intersects(&self) -> bool {
        !self.is_disjoint()
    }

    /// `A` within `B`: `T*F**F***`.
    #[must_use]
    pub fn is_within(&self) -> bool {
        self.is_t(I, I) && self.is_f(I, E) && self.is_f(B, E)
    }

    /// `A` contains `B` — within with the operands swapped.
    #[must_use]
    pub fn is_contains(&self) -> bool {
        self.is_t(I, I) && self.is_f(E, I) && self.is_f(E, B)
    }

    /// `A` covers `B`: `T*****FF*` | `*T****FF*` | `***T**FF*` | `****T*FF*`.
    #[must_use]
    pub fn is_covers(&self) -> bool {
        (self.is_t(I, I) || self.is_t(I, B) || self.is_t(B, I) || self.is_t(B, B))
            && self.is_f(E, I)
            && self.is_f(E, B)
    }

    /// `A` covered by `B` — covers with the operands swapped.
    #[must_use]
    pub fn is_coveredby(&self) -> bool {
        (self.is_t(I, I) || self.is_t(I, B) || self.is_t(B, I) || self.is_t(B, B))
            && self.is_f(I, E)
            && self.is_f(B, E)
    }

    /// `A` equals `B` (same point set): `T*F**FFF*`.
    #[must_use]
    pub fn is_equal_topo(&self) -> bool {
        self.is_t(I, I) && self.is_f(I, E) && self.is_f(E, I) && self.is_f(B, E) && self.is_f(E, B)
    }

    /// `A` touches `B`: interiors disjoint, some boundary contact.
    /// `FT*******` | `F**T*****` | `F***T****`.
    #[must_use]
    pub fn is_touches(&self) -> bool {
        self.is_f(I, I) && (self.is_t(B, B) || self.is_t(I, B) || self.is_t(B, I))
    }

    /// `A` crosses `B` — dimension-dependent (OGC).
    #[allow(clippy::nonminimal_bool)]
    #[must_use]
    pub fn is_crosses(&self, dim_a: i8, dim_b: i8) -> bool {
        if (dim_a == 0 && dim_b == 1) || (dim_a == 0 && dim_b == 2) || (dim_a == 1 && dim_b == 2) {
            self.is_t(I, I) && self.is_t(I, E)
        } else if (dim_a == 1 && dim_b == 0)
            || (dim_a == 2 && dim_b == 0)
            || (dim_a == 2 && dim_b == 1)
        {
            self.is_t(I, I) && self.is_t(E, I)
        } else if dim_a == 1 && dim_b == 1 {
            self.get(I, I) == 0
        } else {
            false
        }
    }

    /// `A` overlaps `B` — dimension-dependent (OGC). Only defined for equal dims.
    #[must_use]
    pub fn is_overlaps(&self, dim_a: i8, dim_b: i8) -> bool {
        if dim_a != dim_b {
            return false;
        }
        if dim_a == 0 || dim_a == 2 {
            self.is_t(I, I) && self.is_t(I, E) && self.is_t(E, I)
        } else if dim_a == 1 {
            self.get(I, I) == 1 && self.is_t(I, E) && self.is_t(E, I)
        } else {
            false
        }
    }

    /// Every cell, rows and columns ordered Interior, Boundary, Exterior.
    ///
    /// The named predicates above each read a fixed handful of cells. A caller
    /// that wants a predicate this type does not name — the OGC `relate`
    /// operation takes an arbitrary pattern — needs the whole matrix, and
    /// cannot recover it from the booleans, because each of those is a
    /// many-to-one function of the nine dimensions.
    #[must_use]
    pub fn cells(&self) -> [[i8; 3]; 3] {
        self.m
    }

    /// The matrix as its nine-character DE-9IM string, row-major from
    /// Interior/Interior, using the OGC spelling: `F` for an empty
    /// intersection and `0`/`1`/`2` for one of that dimension.
    ///
    /// `T` never appears in a rendered matrix. It is a *pattern* character
    /// meaning "any dimension", so it belongs to [`Self::matches`] and not
    /// here; a renderer that emitted it would be discarding the dimension it
    /// had already computed.
    #[must_use]
    pub fn de9im(&self) -> String {
        let mut rendered = String::with_capacity(9);
        for row in self.m {
            for dim in row {
                rendered.push(match dim {
                    0 => '0',
                    1 => '1',
                    2 => '2',
                    // Any negative dimension is "no intersection"; the filler
                    // only ever writes -1, and clamping here keeps a future
                    // sentinel from rendering as a stray digit.
                    _ => 'F',
                });
            }
        }
        rendered
    }

    /// Whether this matrix satisfies a nine-character DE-9IM pattern.
    ///
    /// The pattern alphabet is the OGC one, per cell:
    ///
    /// * `0`, `1`, `2` — the intersection has exactly that dimension
    /// * `T` — the cells intersect, at any dimension
    /// * `F` — the cells do not intersect
    /// * `*` — this cell is not constrained
    ///
    /// Returns `None` for a pattern that is not nine of those characters,
    /// rather than treating a malformed pattern as unsatisfied: a caller that
    /// cannot tell a typo from a false has no way to report the typo, and a
    /// mistyped pattern that silently matches nothing is the failure mode this
    /// avoids. Lowercase is accepted for `t` and `f`.
    #[must_use]
    pub fn matches(&self, pattern: &str) -> Option<bool> {
        let mut cells = self.m.iter().flatten();
        let mut characters = 0usize;
        let mut satisfied = true;
        for expected in pattern.chars() {
            let dim = *cells.next()?;
            characters += 1;
            let cell_ok = match expected {
                '*' => true,
                'T' | 't' => dim >= 0,
                'F' | 'f' => dim < 0,
                '0' => dim == 0,
                '1' => dim == 1,
                '2' => dim == 2,
                _ => return None,
            };
            satisfied &= cell_ok;
        }
        (characters == 9).then_some(satisfied)
    }
}
