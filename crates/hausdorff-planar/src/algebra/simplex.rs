//! Point-in-simplex instances — point-in-triangle (R2, [`Triangle`]) and
//! point-in-tetrahedron (R3, [`Tetrahedron`]) — plus the tetrahedron's
//! circumsphere. The [`Simplex`] / [`Barycentric`] trait and the `whole`/`parts`
//! containment certificate live upstream in [`hausdorff_traits`],
//! machine-checked in `riemann` (`Barycentric.weights_pos_sum_one`, instantiated
//! `triangleBary` / `tetraBary`).
//!
//! For R2 the alloc-free [`in_triangle`](hausdorff_traits::in_triangle)
//! remains the fast path; [`Simplex`] additionally exposes the barycentric
//! decomposition (the [`weights`](Barycentric::weights), useful for
//! interpolation/shading) and the 3-D tetrahedron case.

use hausdorff_traits::{Barycentric, Orientation, Simplex};

use crate::space::r2::{Coord, Triangle};
use hausdorff_base::space::r3::Vector;

impl Simplex for Triangle {
    type Point = Coord;
    type Scalar = f64;

    /// `triangleBary`: `whole = orient(a,b,c)`; each part swaps `p` into one vertex.
    fn barycentric(&self, p: Coord) -> Barycentric {
        let Triangle(a, b, c) = *self;
        Barycentric {
            whole: Coord::orient(a, b, c),
            parts: vec![
                Coord::orient(p, b, c),
                Coord::orient(a, p, c),
                Coord::orient(a, b, p),
            ],
        }
    }
}

/// A tetrahedron `(a, b, c, d)` in ℝ³ — the R3 simplex, mirror of [`Triangle`].
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Tetrahedron(pub Vector, pub Vector, pub Vector, pub Vector);

/// The 4-point signed volume `((b−a)×(c−a))·(d−a)` (riemann `orient4`) — the
/// [`Orientation`] triple product on the edge vectors.
#[inline]
fn orient4(a: Vector, b: Vector, c: Vector, d: Vector) -> f64 {
    Vector::orient(b.sub(a), c.sub(a), d.sub(a))
}

impl Simplex for Tetrahedron {
    type Point = Vector;
    type Scalar = f64;

    /// `tetraBary`: `whole = orient4(a,b,c,d)`; each part swaps `q` into one vertex.
    fn barycentric(&self, q: Vector) -> Barycentric {
        let Tetrahedron(a, b, c, d) = *self;
        Barycentric {
            whole: orient4(a, b, c, d),
            parts: vec![
                orient4(q, b, c, d),
                orient4(a, q, c, d),
                orient4(a, b, q, d),
                orient4(a, b, c, q),
            ],
        }
    }
}

impl Tetrahedron {
    /// The **circumcenter** — the point equidistant from all four vertices —
    /// or `None` when the vertices are coplanar (zero signed volume), so no
    /// sphere passes through them.
    ///
    /// Built from the three perpendicular-bisector planes by Cramer's rule, in
    /// the closed vector form `a + (‖u‖²(v×w) + ‖v‖²(w×u) + ‖w‖²(u×v)) / (2·u·(v×w))`
    /// with edges `u = b−a`, `v = c−a`, `w = d−a` (the denominator is `2·orient4`).
    /// Proof-backed by `riemann` `ExactArith/Circle.circumsphere_exists`: when
    /// `orient4 ≠ 0` such an equidistant centre exists (the R3 analog of the
    /// planar `circumcenter_exists`).
    #[must_use]
    pub fn circumcenter(self) -> Option<Vector> {
        let Tetrahedron(a, b, c, d) = self;
        let (u, v, w) = (b.sub(a), c.sub(a), d.sub(a));
        let den = 2.0 * u.dot(v.cross(w)); // 2·orient4(a,b,c,d)
        if den == 0.0 {
            return None; // coplanar — no circumsphere
        }
        let num = v
            .cross(w)
            .mul(u.norm2())
            .add(w.cross(u).mul(v.norm2()))
            .add(u.cross(v).mul(w.norm2()));
        Some(a.add(num.mul(1.0 / den)))
    }

    /// The **circumsphere** as `(centre, radius)`, or `None` if the vertices are
    /// coplanar. The radius is the centre's distance to any vertex (equal to all
    /// four by `circumsphere_exists`).
    #[must_use]
    pub fn circumsphere(self) -> Option<(Vector, f64)> {
        let Tetrahedron(a, _, _, _) = self;
        self.circumcenter().map(|o| (o, o.sub(a).norm()))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn triangle_containment_and_weights() {
        let t = Triangle(
            Coord::new(0.0, 0.0),
            Coord::new(1.0, 0.0),
            Coord::new(0.0, 1.0),
        );
        assert!(t.contains(Coord::new(0.25, 0.25))); // inside
        assert!(!t.contains(Coord::new(1.0, 1.0))); // outside
        assert!(!t.contains(Coord::new(0.5, 0.0))); // on an edge → not strict

        // The centroid's weights are all 1/3 and sum to one.
        let w = t.barycentric(Coord::new(1.0 / 3.0, 1.0 / 3.0)).weights();
        assert!(w.iter().all(|wi| (wi - 1.0 / 3.0).abs() < 1e-12));
        assert!((w.iter().sum::<f64>() - 1.0).abs() < 1e-12);
    }

    #[test]
    fn tetrahedron_containment_and_weights() {
        let t = Tetrahedron(
            Vector::new(0.0, 0.0, 0.0),
            Vector::new(1.0, 0.0, 0.0),
            Vector::new(0.0, 1.0, 0.0),
            Vector::new(0.0, 0.0, 1.0),
        );
        assert!(t.contains(Vector::new(0.25, 0.25, 0.25))); // inside (x+y+z < 1)
        assert!(!t.contains(Vector::new(1.0, 1.0, 1.0))); // outside
        assert!(!t.contains(Vector::new(0.5, 0.5, 0.5))); // beyond the far face

        // An interior point's weights are positive and sum to one.
        let w = t.barycentric(Vector::new(0.1, 0.2, 0.3)).weights();
        assert!(w.iter().all(|&wi| wi > 0.0));
        assert!((w.iter().sum::<f64>() - 1.0).abs() < 1e-12);
    }

    #[test]
    fn tetrahedron_circumsphere_is_equidistant() {
        let verts = [
            Vector::new(0.0, 0.0, 0.0),
            Vector::new(1.0, 0.0, 0.0),
            Vector::new(0.0, 1.0, 0.0),
            Vector::new(0.0, 0.0, 1.0),
        ];
        let t = Tetrahedron(verts[0], verts[1], verts[2], verts[3]);
        let (center, radius) = t.circumsphere().expect("non-degenerate");
        // The unit-corner tetra circumcentre is (0.5, 0.5, 0.5), radius √3/2.
        assert!((center.sub(Vector::new(0.5, 0.5, 0.5))).norm() < 1e-12);
        assert!((radius - 3.0_f64.sqrt() / 2.0).abs() < 1e-12);
        // Equidistant from all four vertices (circumsphere_exists).
        for v in verts {
            assert!((center.sub(v).norm() - radius).abs() < 1e-12);
        }

        // A translated, scaled, irregular tetra is still equidistant.
        let q = Tetrahedron(
            Vector::new(2.0, -1.0, 3.0),
            Vector::new(5.0, 0.0, 3.5),
            Vector::new(2.5, 4.0, 2.0),
            Vector::new(3.0, 1.0, 7.0),
        );
        let (qc, qr) = q.circumsphere().unwrap();
        for v in [q.0, q.1, q.2, q.3] {
            assert!((qc.sub(v).norm() - qr).abs() < 1e-9);
        }

        // Coplanar vertices have no circumsphere.
        let flat = Tetrahedron(
            Vector::new(0.0, 0.0, 0.0),
            Vector::new(1.0, 0.0, 0.0),
            Vector::new(0.0, 1.0, 0.0),
            Vector::new(1.0, 1.0, 0.0), // all z = 0
        );
        assert!(flat.circumcenter().is_none());
    }
}
