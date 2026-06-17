//! Point-in-spherical-triangle as the sphere's instance of the shared
//! [`Simplex`](hausdorff_planar::algebra::Simplex) containment trait — the same
//! barycentric orientation partition the planar `Triangle` and the ℝ³
//! `Tetrahedron` use, now over [`S2Point`] vertices and the spherical triple
//! product.
//!
//! A [`SphereTriangle`] is three vertices on the unit sphere. Its
//! [`barycentric`](SphereTriangle::barycentric) decomposition takes
//! `whole = orient(a, b, c)` (the scalar triple product `(a × b) · c`, the
//! [`Orientation`](hausdorff_planar::algebra::Orientation) instance from
//! [`predicates`](super::predicates)) and the three parts swap the query point
//! into each vertex. The strictly-inside verdict is then the generic same-sign
//! test, so `contains` is the spherical mirror of point-in-triangle.
//!
//! Proof-backed by `riemann`: the generic same-sign containment certificate
//! ([`Barycentric`]'s `Geometry/Containment` partition, whose
//! `Barycentric.weights_pos_sum_one` is the *affine* reading) instantiated on the
//! sphere as `sphereInTriangle` (`Hausdorff/SphereOrientation`,
//! `sphereInTriangle_iff` — strictly inside iff
//! `0 < triple a b p ∧ 0 < triple b c p ∧ 0 < triple c a p`), whose left-of-all
//! region is convex (`sphereLeftOfAll_convex`). The same-sign verdict transfers
//! verbatim; the convex-coordinate (sum-to-one) reading of the weights is
//! planar-only — on the sphere the triple products do not partition into one.

use hausdorff_base::space::s1::Angle;
use hausdorff_planar::algebra::{Barycentric, Orientation, Simplex};

use super::{Cap, S2Point};

/// A **spherical triangle** held as its three vertices `(a, b, c)` on the unit
/// sphere — the vertex-carrying companion to [`SphericalTriangle`], whose sides
/// are arc lengths. This form carries the orientation information point-in-triangle
/// containment needs, so it is the sphere's [`Simplex`].
#[derive(Clone, Copy, PartialEq, Debug)]
pub struct SphereTriangle(pub S2Point, pub S2Point, pub S2Point);

impl SphereTriangle {
    /// The triangle through three vertices.
    #[must_use]
    pub fn new(a: S2Point, b: S2Point, c: S2Point) -> Self {
        Self(a, b, c)
    }

    /// The three vertices `(a, b, c)`.
    #[must_use]
    pub fn vertices(&self) -> (S2Point, S2Point, S2Point) {
        (self.0, self.1, self.2)
    }

    /// The **circumcenter** — the pole of the small circle through all three
    /// vertices, equidistant in central angle from each — or `None` when the
    /// vertices are degenerate (a coincident pair, so the chords give no plane).
    ///
    /// Closed form: the cross product of two chords, `(a − b) × (b − c)`, is
    /// orthogonal to each chord and so equidistant (equal dot ⇒ equal central
    /// angle) from all three vertices. It is oriented to the triangle's own side
    /// (`dot ≥ 0`, the near pole) so the circumscribed circle is the *minor* one
    /// the vertices bound. Proof-backed by `riemann`
    /// `Hausdorff/SphereCircle.sphere_circumcenter` — the spherical mirror of the
    /// planar [`Triangle::circumcenter`](hausdorff_planar::space::r2::Triangle::circumcenter)
    /// / `ExactArith/Circle.circumcenter_exists`, with an explicit cross-product
    /// pole instead of a Cramer solve.
    #[must_use]
    pub fn circumcenter(&self) -> Option<S2Point> {
        let SphereTriangle(a, b, c) = *self;
        let (va, vb, vc) = (a.as_vector(), b.as_vector(), c.as_vector());
        let pole = va.sub(vb).cross(vb.sub(vc)); // (a − b) × (b − c)
        if pole.norm2() == 0.0 {
            return None; // degenerate — chords are parallel (a coincident pair)
        }
        // Orient to the triangle's side so the circumscribed cap is the minor cap.
        let pole = if pole.dot(va) >= 0.0 { pole } else { pole.mul(-1.0) };
        Some(S2Point::from_direction(pole))
    }

    /// The **circumscribed cap** — the spherical cap whose boundary passes
    /// through all three vertices: centred at the [`circumcenter`](Self::circumcenter)
    /// with angular radius the central angle to any vertex. `None` if the
    /// vertices are degenerate. The minor circumscribing cap, a ready bounding
    /// volume for spatial indexing (the spherical analog of
    /// [`Tetrahedron::circumsphere`](hausdorff_planar::algebra::Tetrahedron::circumsphere)).
    #[must_use]
    pub fn circumscribed_cap(&self) -> Option<Cap> {
        let SphereTriangle(a, _, _) = *self;
        self.circumcenter()
            .map(|center| Cap::from_center_angle(center, Angle::from_radians(center.distance(a))))
    }
}

impl Simplex for SphereTriangle {
    type Point = S2Point;
    type Scalar = f64;

    /// `sphereInTriangle`: `whole = orient(a, b, c)` (the spherical triple
    /// product); each part swaps `p` into one vertex (riemann
    /// `Hausdorff/SphereOrientation.sphereInTriangle_iff`).
    fn barycentric(&self, p: S2Point) -> Barycentric {
        let SphereTriangle(a, b, c) = *self;
        Barycentric {
            whole: S2Point::orient(a, b, c),
            parts: vec![
                S2Point::orient(p, b, c),
                S2Point::orient(a, p, c),
                S2Point::orient(a, b, p),
            ],
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use hausdorff_base::earth::GeoPoint;

    fn pt(lat: f64, lon: f64) -> S2Point {
        S2Point::from_geo_point(GeoPoint::from_lat_lon(lat, lon))
    }

    #[test]
    fn point_in_spherical_triangle() {
        // Three mid-latitude points 120° apart in longitude; the north-polar cap
        // they bound contains a near-pole point but not a southern one.
        let t = SphereTriangle::new(pt(10.0, 0.0), pt(10.0, 120.0), pt(10.0, 240.0));
        assert!(t.contains(pt(80.0, 0.0))); // deep inside (near the pole)
        assert!(!t.contains(pt(-80.0, 0.0))); // the antipodal cap → outside
        // A vertex lies on the boundary, never strictly inside.
        assert!(!t.contains(pt(10.0, 0.0)));
    }

    #[test]
    fn barycentric_parts_share_the_whole_sign_inside() {
        let t = SphereTriangle::new(pt(10.0, 0.0), pt(10.0, 120.0), pt(10.0, 240.0));
        let inside = pt(80.0, 0.0);
        let bary = t.barycentric(inside);
        // The containment certificate (`sphereInTriangle_iff`): the whole and
        // every part share one nonzero sign — that, not a planar
        // partition-of-unity, is what makes the point strictly inside on the
        // sphere. (Unlike the affine plane the spherical triple products do not
        // sum to the whole, so the convex-coordinate reading is planar-only.)
        assert!(bary.strictly_inside());
        let sign = bary.whole.signum();
        assert!(bary.parts.iter().all(|&p| p.signum() == sign && p != 0.0));
        // A point in the antipodal cap is strictly outside.
        let outside = t.barycentric(pt(-80.0, 0.0));
        assert!(!outside.strictly_inside());
    }

    #[test]
    fn agrees_with_generic_in_triangle() {
        use hausdorff_planar::algebra::in_triangle;
        let (a, b, c) = (pt(10.0, 0.0), pt(10.0, 120.0), pt(10.0, 240.0));
        let t = SphereTriangle::new(a, b, c);
        for q in [pt(80.0, 0.0), pt(-80.0, 0.0), pt(30.0, 60.0)] {
            assert_eq!(t.contains(q), in_triangle(a, b, c, q));
        }
    }

    #[test]
    fn circumcircle_is_equidistant_and_bounds_the_triangle() {
        use std::f64::consts::FRAC_PI_2;
        // Three points at latitude 10° N, 120° apart: their circumcircle is that
        // small circle, so the circumcentre is the north pole.
        let (a, b, c) = (pt(10.0, 0.0), pt(10.0, 120.0), pt(10.0, 240.0));
        let t = SphereTriangle::new(a, b, c);
        let center = t.circumcenter().expect("non-degenerate");
        assert!(center.distance(pt(90.0, 0.0)) < 1e-9); // the north pole
        // Equidistant from all three vertices (sphere_circumcenter).
        let r = center.distance(a);
        assert!((center.distance(b) - r).abs() < 1e-12);
        assert!((center.distance(c) - r).abs() < 1e-12);
        // The circumscribed cap carries that radius and contains its centre.
        let cap = t.circumscribed_cap().unwrap();
        assert!((cap.radius().radians() - r).abs() < 1e-12);
        assert!(cap.contains_point(center));
        // Three points on one great circle (the equator) circumscribe to a
        // radius-π/2 cap whose pole is the great circle's pole.
        let gc = SphereTriangle::new(pt(0.0, 0.0), pt(0.0, 90.0), pt(0.0, 180.0));
        assert!((gc.circumscribed_cap().unwrap().radius().radians() - FRAC_PI_2).abs() < 1e-9);
        // Degenerate: a repeated vertex gives parallel chords → None.
        assert!(SphereTriangle::new(a, a, c).circumcenter().is_none());
        assert!(SphereTriangle::new(a, a, c).circumscribed_cap().is_none());
    }
}
