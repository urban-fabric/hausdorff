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
//! `triangleBary` (R2) and `tetraBary` (R3). [`Triangle`] and [`Tetrahedron`]
//! mirror those instances, so the inside-test is written once over [`Simplex`].
//!
//! For R2 the alloc-free [`in_triangle`](super::in_triangle) remains the fast
//! path; [`Simplex`] additionally exposes the barycentric decomposition (the
//! [`weights`](Barycentric::weights), useful for interpolation/shading) and the
//! 3-D tetrahedron case.

use super::Orientation;
use crate::space::r2::{Coord, Triangle};
use hausdorff_base::space::r3::Vector;

/// A barycentric (orientation) partition of a query point against a simplex: the
/// simplex's signed `whole` orientation and the per-vertex `parts` (the query
/// point swapped into each vertex), which sum to `whole`.
#[derive(Clone, Debug, PartialEq)]
pub struct Barycentric {
    /// The full signed orientation of the simplex (twice the area in R2, six
    /// times the volume in R3).
    pub whole: f64,
    /// The per-vertex sub-orientations — one per vertex, the query point swapped
    /// in. They sum to [`whole`](Self::whole).
    pub parts: Vec<f64>,
}

impl Barycentric {
    /// `true` iff the query point is **strictly inside** — `whole` and every part
    /// share one nonzero sign. Then the weights `partᵢ/whole` are all positive
    /// and sum to one (riemann `Barycentric.weights_pos_sum_one`); a degenerate
    /// simplex (`whole == 0`) contains nothing.
    #[must_use]
    pub fn strictly_inside(&self) -> bool {
        (self.whole > 0.0 && self.parts.iter().all(|&x| x > 0.0))
            || (self.whole < 0.0 && self.parts.iter().all(|&x| x < 0.0))
    }

    /// The barycentric weights `partᵢ / whole` — the convex coordinates of the
    /// query point. Inside the simplex they are all positive and sum to one.
    #[must_use]
    pub fn weights(&self) -> Vec<f64> {
        self.parts.iter().map(|&p| p / self.whole).collect()
    }
}

/// A simplex whose point-containment is decided by a barycentric orientation
/// partition — point-in-triangle (R2) and point-in-tetrahedron (R3) as one test.
pub trait Simplex {
    /// The point type of the ambient space.
    type Point: Copy;

    /// The barycentric decomposition of `p` against this simplex (riemann
    /// `triangleBary` / `tetraBary`).
    fn barycentric(&self, p: Self::Point) -> Barycentric;

    /// `true` iff `p` is strictly inside — generic over the simplex dimension.
    #[must_use]
    fn contains(&self, p: Self::Point) -> bool {
        self.barycentric(p).strictly_inside()
    }
}

impl Simplex for Triangle {
    type Point = Coord;

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
}
