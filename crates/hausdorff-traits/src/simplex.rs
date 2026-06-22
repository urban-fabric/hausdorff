//! Point-in-simplex containment via a barycentric orientation partition — after
//! the `riemann` Lean `Barycentric` abstraction (`Geometry/Containment`).
//!
//! Point-in-triangle (R2) and point-in-tetrahedron (R3) are the *same* test: the
//! simplex's signed `whole` orientation, plus the per-vertex `parts` (the query
//! point swapped into each vertex) that **sum to the whole** (the barycentric
//! partition of unity). The containment certificate — `whole` and every part
//! sharing one nonzero sign ⟹ the weights `partᵢ/whole` are all positive and sum
//! to one ⟹ the point is strictly inside — is machine-checked **once** in
//! `riemann` (`Barycentric.weights_pos_sum_one`) and instantiated at
//! `triangleBary` (R2) and `tetraBary` (R3). The `Triangle` and `Tetrahedron`
//! impls (in `hausdorff-planar`) and `SphereTriangle` (in `hausdorff-sphere`)
//! mirror those instances, so the inside-test is written once over [`Simplex`].

use euclid::numeric::{Boolean as BooleanMask, Real};

use crate::orientation::Boolean;

#[inline]
fn mask_and<M: BooleanMask>(a: M, b: M) -> M {
    a.and(&b)
}

#[inline]
fn mask_or<M: BooleanMask>(a: M, b: M) -> M {
    a.or(&b)
}

/// A barycentric (orientation) partition of a query point against a simplex: the
/// simplex's signed `whole` orientation and the per-vertex `parts` (the query
/// point swapped into each vertex), which sum to `whole`.
#[derive(Clone, Debug, PartialEq)]
pub struct Barycentric<T> {
    /// The full signed orientation of the simplex (twice the area in R2, six
    /// times the volume in R3).
    pub whole: T,
    /// The per-vertex sub-orientations — one per vertex, the query point swapped
    /// in. They sum to [`whole`](Self::whole).
    pub parts: Vec<T>,
}

impl<T> Barycentric<T>
where
    T: Real,
{
    /// `true` iff the query point is **strictly inside** — `whole` and every part
    /// share one nonzero sign. Then the weights `partᵢ/whole` are all positive
    /// and sum to one (riemann `Barycentric.weights_pos_sum_one`); a degenerate
    /// simplex (`whole == 0`) contains nothing.
    #[must_use]
    pub fn strictly_inside(&self) -> Boolean<T> {
        let zero = T::zero();
        let mut all_positive = self.whole.gt(&zero);
        let mut all_negative = self.whole.lt(&zero);

        for part in &self.parts {
            all_positive = mask_and(all_positive, part.gt(&zero));
            all_negative = mask_and(all_negative, part.lt(&zero));
        }

        mask_or(all_positive, all_negative)
    }

    /// The barycentric weights `partᵢ / whole` — the convex coordinates of the
    /// query point. Inside the simplex they are all positive and sum to one.
    #[must_use]
    pub fn weights(&self) -> Vec<T> {
        self.parts
            .iter()
            .cloned()
            .map(|p| p / self.whole.clone())
            .collect()
    }
}

/// A simplex whose point-containment is decided by a barycentric orientation
/// partition — point-in-triangle (R2) and point-in-tetrahedron (R3) as one test.
pub trait Simplex {
    /// The point type of the ambient space.
    type Point: Clone;

    /// The scalar field the orientation partition is measured in.
    type Scalar: Real;

    /// The barycentric decomposition of `p` against this simplex (riemann
    /// `triangleBary` / `tetraBary`).
    fn barycentric(&self, p: Self::Point) -> Barycentric<Self::Scalar>;

    /// `true` iff `p` is strictly inside — generic over the simplex dimension.
    #[must_use]
    fn contains(&self, p: Self::Point) -> Boolean<Self::Scalar> {
        self.barycentric(p).strictly_inside()
    }
}
