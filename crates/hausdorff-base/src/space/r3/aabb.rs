use super::*;

#[inline]
fn mask_and<M: BooleanMask>(a: M, b: M) -> M {
    a.and(&b)
}

/// An **axis-aligned 3-D bounding box**, stored by its `min` (lowest x, y, z) and
/// `max` (highest) corners — the ℝ³ member of the Euclidean box ladder alongside
/// [`r1::Interval`](crate::space::r1::Interval) (1-D) and `r2::Rect` (2-D, in
/// `hausdorff-planar`), and the box an R-tree over `Vector` points prunes with.
///
/// Membership and the overlap / containment pruning tests are per-axis interval
/// arithmetic, proof-backed by `riemann` `ExactArith/AABB` (`AABB3`):
/// [`overlaps`](Self::overlaps) by `overlaps_of_mem3` (a shared point forces an
/// overlap, so a non-overlapping pair is safe to prune — `prune_sound3`),
/// [`contains_box`](Self::contains_box) by `contains_of_subset3` (containment is
/// monotone) with `subset_trans3` carrying the leaf⊆…⊆root invariant.
#[derive(Clone, Copy, PartialEq, Debug)]
pub struct Aabb<T = f64> {
    min: Vector<T>,
    max: Vector<T>,
}

impl<T> Aabb<T>
where
    T: Real,
{
    /// The box spanning two opposite corners, given in any order (normalized to
    /// per-axis `min`/`max`).
    #[must_use]
    pub fn new(c1: Vector<T>, c2: Vector<T>) -> Self {
        Self {
            min: Vector::new(
                c1.x.clone().min(c2.x.clone()),
                c1.y.clone().min(c2.y.clone()),
                c1.z.clone().min(c2.z.clone()),
            ),
            max: Vector::new(c1.x.max(c2.x), c1.y.max(c2.y), c1.z.max(c2.z)),
        }
    }

    /// The box from already-ordered corners (`min ≤ max` per axis assumed).
    #[must_use]
    pub fn from_min_max(min: Vector<T>, max: Vector<T>) -> Self {
        Self { min, max }
    }

    /// The degenerate box containing only the point `p`.
    #[must_use]
    pub fn from_point(p: Vector<T>) -> Self {
        Self {
            min: p.clone(),
            max: p,
        }
    }

    /// The lower corner (minimum x, y, z).
    #[inline]
    #[must_use]
    pub fn min(&self) -> Vector<T> {
        self.min.clone()
    }

    /// The upper corner (maximum x, y, z).
    #[inline]
    #[must_use]
    pub fn max(&self) -> Vector<T> {
        self.max.clone()
    }

    /// The centre point.
    #[inline]
    #[must_use]
    pub fn center(&self) -> Vector<T> {
        self.min.clone().add(self.max.clone()).mul(T::half())
    }

    /// True iff `p` lies in the box (edges inclusive) — the per-axis range test
    /// (`riemann` `AABB3.contains`).
    #[inline]
    #[must_use]
    pub fn contains(&self, p: Vector<T>) -> <T as Compare>::Mask {
        mask_and(
            mask_and(
                mask_and(self.min.x.le(&p.x), p.x.le(&self.max.x)),
                mask_and(self.min.y.le(&p.y), p.y.le(&self.max.y)),
            ),
            mask_and(self.min.z.le(&p.z), p.z.le(&self.max.z)),
        )
    }

    /// True iff this box and `other` overlap (share any point, faces inclusive) —
    /// the per-axis interval-overlap test. A non-overlapping pair shares no point,
    /// so a spatial index may prune it (`riemann` `AABB3.overlaps_of_mem3` /
    /// `prune_sound3`).
    #[inline]
    #[must_use]
    pub fn overlaps(&self, other: &Aabb<T>) -> <T as Compare>::Mask {
        mask_and(
            mask_and(
                mask_and(self.min.x.le(&other.max.x), other.min.x.le(&self.max.x)),
                mask_and(self.min.y.le(&other.max.y), other.min.y.le(&self.max.y)),
            ),
            mask_and(self.min.z.le(&other.max.z), other.min.z.le(&self.max.z)),
        )
    }

    /// True iff this box fully contains `other` (`riemann` `AABB3.subset`, with
    /// `contains_of_subset3` ⟹ every point of `other` is in `self`).
    #[inline]
    #[must_use]
    pub fn contains_box(&self, other: &Aabb<T>) -> <T as Compare>::Mask {
        mask_and(
            mask_and(
                mask_and(self.min.x.le(&other.min.x), other.max.x.le(&self.max.x)),
                mask_and(self.min.y.le(&other.min.y), other.max.y.le(&self.max.y)),
            ),
            mask_and(self.min.z.le(&other.min.z), other.max.z.le(&self.max.z)),
        )
    }

    /// The smallest box containing both `self` and `other`.
    #[must_use]
    pub fn union(&self, other: &Aabb<T>) -> Aabb<T> {
        Aabb {
            min: Vector::new(
                self.min.x.clone().min(other.min.x.clone()),
                self.min.y.clone().min(other.min.y.clone()),
                self.min.z.clone().min(other.min.z.clone()),
            ),
            max: Vector::new(
                self.max.x.clone().max(other.max.x.clone()),
                self.max.y.clone().max(other.max.y.clone()),
                self.max.z.clone().max(other.max.z.clone()),
            ),
        }
    }

    /// The box grown by `margin` on every side (a negative `margin` shrinks it).
    #[must_use]
    pub fn expanded(&self, margin: T) -> Aabb<T> {
        let d = Vector::new(margin.clone(), margin.clone(), margin);
        Aabb {
            min: self.min.clone().sub(d.clone()),
            max: self.max.clone().add(d),
        }
    }
}
