//! Planar **isometries** — the rigid-motion subgroup of the affine group.
//!
//! An [`Isometry`] is a distance-preserving planar map: a rotation, reflection,
//! translation, or any composition of these. It is a thin, invariant-carrying
//! wrapper around [`AffineTransform`] whose wrapped matrix is **always rigid**
//! ([`AffineTransform::is_rigid`]). The rigid maps are closed under composition
//! and inversion and contain the identity, so they form a subgroup — that
//! structure is surfaced here as a type so APIs that *require* a metric-preserving
//! transform can ask for one statically instead of trusting an arbitrary
//! `AffineTransform`.

use super::{AffineTransform, Coord};

/// A planar **isometry** (rigid motion) — a distance- and angle-preserving map
/// of the plane onto itself: rotations, reflections, translations, and their
/// compositions. Wraps an [`AffineTransform`] whose linear part has orthonormal
/// columns; the wrapper's constructors and operations preserve that invariant,
/// so an `Isometry` is *always* rigid.
///
/// The isometries are a **subgroup** of the affine group: closed under
/// [`compose`](Isometry::compose) (rigid ∘ rigid is rigid), containing
/// [`identity`](Isometry::identity), and closed under [`inverse`](Isometry::inverse)
/// (the inverse of a rigid map is rigid).
///
/// Proof-backed by `riemann` `Hausdorff/AffineIsometry`: `act_isometry`
/// (a map with orthonormal columns preserves squared distance —
/// `sqDist(M·p, M·q) = sqDist(p, q)`) and `rotation_isometry` (rotations are a
/// special case). The subgroup (group) laws are `Hausdorff/AffineMap`:
/// `id_comp`/`comp_id` (identity), `comp_assoc` (associativity), and
/// `comp_inv`/`inv_comp` (inverses).
#[derive(Copy, Clone, PartialEq, Debug)]
pub struct Isometry(AffineTransform);

impl Default for Isometry {
    fn default() -> Self {
        Self::identity()
    }
}

impl Isometry {
    /// The identity isometry — the map fixing every point.
    ///
    /// The unit of the subgroup (`Hausdorff/AffineMap.id_comp`/`comp_id`).
    #[must_use]
    pub fn identity() -> Self {
        Self(AffineTransform::identity())
    }

    /// A **rotation** about the origin by `radians` (counter-clockwise when
    /// positive). Rigid by `Hausdorff/AffineIsometry.rotation_isometry`.
    #[must_use]
    pub fn rotation(radians: f64) -> Self {
        Self(AffineTransform::rotate(
            radians.to_degrees(),
            Coord::new(0.0, 0.0),
        ))
    }

    /// A **translation** by `(dx, dy)`. Translations are isometries (the linear
    /// part is the identity, trivially orthonormal).
    #[must_use]
    pub fn translation(dx: f64, dy: f64) -> Self {
        Self(AffineTransform::translate(dx, dy))
    }

    /// Wrap an [`AffineTransform`] as an isometry, returning `Some` **iff** it is
    /// [rigid](AffineTransform::is_rigid) (orthonormal-column linear part) and
    /// `None` otherwise — e.g. a non-uniform scale or a shear is rejected.
    ///
    /// This is the constructor that enforces the type's invariant; the precise
    /// distance-preservation it guarantees is `Hausdorff/AffineIsometry.act_isometry`.
    #[must_use]
    pub fn from_affine(t: AffineTransform) -> Option<Self> {
        t.is_rigid().then_some(Self(t))
    }

    /// **Compose** two isometries: the result applies `self` first, then `other`
    /// — the same order as [`AffineTransform::compose`]. Rigid ∘ rigid is rigid,
    /// so the isometries are closed under composition (the subgroup is closed;
    /// `Hausdorff/AffineMap.act_comp`).
    #[must_use]
    pub fn compose(&self, other: &Self) -> Self {
        Self(self.0.compose(&other.0))
    }

    /// The **inverse** isometry, undoing `self`. A rigid map has `|det| = 1 ≠ 0`,
    /// so it is always invertible and the inverse is itself rigid — hence this
    /// returns an `Isometry` (not an `Option`).
    ///
    /// `Hausdorff/AffineMap.comp_inv`/`inv_comp` (composing with the inverse
    /// yields the identity).
    #[must_use]
    pub fn inverse(&self) -> Self {
        // A rigid transform is never singular (|det| = 1), so `inverse` is `Some`.
        Self(
            self.0
                .inverse()
                .expect("a rigid transform is invertible (|det| = 1)"),
        )
    }

    /// Apply the isometry to a coordinate.
    #[must_use]
    pub fn apply(&self, p: Coord) -> Coord {
        self.0.apply(p)
    }

    /// The underlying [`AffineTransform`] (borrowed).
    #[must_use]
    pub fn as_affine(&self) -> &AffineTransform {
        &self.0
    }

    /// Consume the isometry, returning the underlying [`AffineTransform`].
    #[must_use]
    pub fn into_affine(self) -> AffineTransform {
        self.0
    }
}

/// Isometries compose with `*`, applying `self` then `rhs` — the same order as
/// [`Isometry::compose`] and [`AffineTransform`]'s own `Mul`. With the
/// [`identity`](Isometry::identity) unit and [`inverse`](Isometry::inverse) this
/// is the subgroup's group operation (`Hausdorff/AffineMap`).
impl std::ops::Mul for Isometry {
    type Output = Isometry;

    fn mul(self, rhs: Isometry) -> Isometry {
        self.compose(&rhs)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-12;
    const PI: f64 = std::f64::consts::PI;

    /// Squared Euclidean distance between two coordinates.
    fn sq_dist(p: Coord, q: Coord) -> f64 {
        let (dx, dy) = (p.x - q.x, p.y - q.y);
        dx * dx + dy * dy
    }

    fn close(p: Coord, q: Coord) -> bool {
        (p.x - q.x).abs() < EPS && (p.y - q.y).abs() < EPS
    }

    #[test]
    fn isometry_preserves_squared_distance() {
        let p = Coord::new(2.0, -3.0);
        let q = Coord::new(-1.0, 5.0);
        let before = sq_dist(p, q);
        // A rotation, a translation, and a rotation-then-translation composite.
        for iso in [
            Isometry::rotation(0.7),
            Isometry::translation(4.0, -2.0),
            Isometry::rotation(1.3).compose(&Isometry::translation(-5.0, 8.0)),
        ] {
            let after = sq_dist(iso.apply(p), iso.apply(q));
            assert!((after - before).abs() < EPS, "{after} vs {before}");
        }
    }

    #[test]
    fn from_affine_accepts_rigid_rejects_non_uniform_scale() {
        // A rotation is rigid → Some.
        assert!(
            Isometry::from_affine(AffineTransform::rotate(37.0, Coord::new(0.0, 0.0))).is_some()
        );
        // A pure translation is rigid → Some.
        assert!(Isometry::from_affine(AffineTransform::translate(3.0, -1.0)).is_some());
        // A reflection (negative uniform scale, |det| = 1) is rigid → Some.
        assert!(
            Isometry::from_affine(AffineTransform::scale(-1.0, 1.0, Coord::new(0.0, 0.0)))
                .is_some()
        );
        // A non-uniform scale is NOT rigid → None.
        assert!(
            Isometry::from_affine(AffineTransform::scale(2.0, 1.0, Coord::new(0.0, 0.0))).is_none()
        );
        // A shear is NOT rigid → None.
        assert!(
            Isometry::from_affine(AffineTransform::skew(30.0, 0.0, Coord::new(0.0, 0.0))).is_none()
        );
    }

    #[test]
    fn compose_and_inverse_round_trip_to_identity() {
        let p = Coord::new(1.5, -2.25);
        let iso = Isometry::rotation(0.9).compose(&Isometry::translation(7.0, -3.0));
        // i ∘ i⁻¹ is the identity on points.
        let round = iso.compose(&iso.inverse()).apply(p);
        assert!(close(round, p), "{round:?}");
        // i⁻¹ ∘ i too.
        let round2 = iso.inverse().compose(&iso).apply(p);
        assert!(close(round2, p), "{round2:?}");
        // The `*` operator agrees with `compose`.
        assert!(close(
            (iso * iso.inverse()).apply(p),
            iso.compose(&iso.inverse()).apply(p)
        ));
    }

    #[test]
    fn rotation_by_theta_then_negative_theta_is_identity() {
        let p = Coord::new(3.0, 4.0);
        let id = Isometry::rotation(0.6).compose(&Isometry::rotation(-0.6));
        assert!(close(id.apply(p), p));
        // A full turn returns every point to itself.
        let full = Isometry::rotation(2.0 * PI);
        assert!(close(full.apply(p), p));
        // The identity isometry fixes points and is its own inverse.
        assert!(close(Isometry::identity().apply(p), p));
        assert!(close(Isometry::identity().inverse().apply(p), p));
    }

    #[test]
    fn composition_stays_rigid_a_subgroup() {
        // Every constructed isometry's wrapped transform is rigid (closure).
        let a = Isometry::rotation(0.4);
        let b = Isometry::translation(2.0, -6.0);
        assert!(a.compose(&b).as_affine().is_rigid());
        assert!(a.inverse().as_affine().is_rigid());
        assert!(Isometry::identity().as_affine().is_rigid());
        // into_affine returns the same matrix as_affine borrows.
        assert_eq!(a.into_affine(), *a.as_affine());
    }
}
