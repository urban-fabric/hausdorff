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
}
