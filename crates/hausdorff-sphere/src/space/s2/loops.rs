//! `s2::Loop` — a simple closed loop on the **unit sphere**, i.e. a spherical
//! polygon boundary. By convention the interior is the region to the *left* of
//! the directed edges, so a counter-clockwise loop encloses its interior.
//!
//! Everything here is generic unit-sphere geometry: [`area`](Loop::area) is in
//! **steradians** and [`perimeter`](Loop::perimeter) is an [`Angle`], both
//! radius-independent. Convert to metres/m² for a particular planet with
//! [`hausdorff_base::earth::sphere::Sphere`]. Computed from unit vectors — no projection, no
//! pole/antimeridian special cases.

use std::f64::consts::PI;

use hausdorff_base::space::s1::Angle;

use super::S2Point;
use super::edges::{Crossing, crossing_sign};
use super::measures::turn_angle;

/// A closed loop of vertices on the unit sphere.
#[derive(Clone, Debug, PartialEq)]
pub struct Loop {
    vertices: Vec<S2Point>,
}

impl Loop {
    /// A loop over `vertices`. The closing edge from the last vertex back to the
    /// first is implicit; a trailing vertex equal to the first is dropped.
    #[must_use]
    pub fn new(mut vertices: Vec<S2Point>) -> Self {
        if vertices.len() >= 2 && vertices[0] == *vertices.last().unwrap() {
            vertices.pop();
        }
        Self { vertices }
    }

    /// The loop vertices (no repeated closing vertex).
    #[must_use]
    pub fn vertices(&self) -> &[S2Point] {
        &self.vertices
    }

    /// The number of vertices (= the number of edges).
    #[must_use]
    pub fn num_vertices(&self) -> usize {
        self.vertices.len()
    }

    /// The directed edges `(vᵢ, vᵢ₊₁)`, wrapping at the end.
    fn edges(&self) -> impl Iterator<Item = (S2Point, S2Point)> + '_ {
        let n = self.vertices.len();
        (0..n).map(move |i| (self.vertices[i], self.vertices[(i + 1) % n]))
    }

    /// Area of the loop interior in **steradians**, in `[0, 4π]` — the region to
    /// the left of the directed edges. Robust for loops of any size (including
    /// great-circle loops) via the Gauss–Bonnet theorem: `area = 2π − Σ turn`.
    ///
    /// Relative accuracy degrades for sub-square-metre loops (the turning-angle
    /// method has a small absolute error per vertex); for tiny triangles use the
    /// [`super::point_area`] measures directly.
    #[must_use]
    pub fn area(&self) -> f64 {
        let n = self.vertices.len();
        if n < 3 {
            return 0.0;
        }
        let mut turning = 0.0;
        for i in 0..n {
            let prev = self.vertices[(i + n - 1) % n];
            let cur = self.vertices[i];
            let next = self.vertices[(i + 1) % n];
            turning += turn_angle(prev, cur, next).radians();
        }
        (2.0 * PI - turning).clamp(0.0, 4.0 * PI)
    }

    /// True iff `p` lies in the loop interior. Uses the spherical winding number
    /// (the signed solid angle the loop subtends at `p`): for a counter-clockwise
    /// loop an interior point sees ≈ +2π and an exterior point ≈ 0 or −2π.
    /// Robust away from the boundary; results exactly on an edge/vertex are not
    /// guaranteed.
    #[must_use]
    pub fn contains_point(&self, p: S2Point) -> bool {
        if self.vertices.len() < 3 {
            return false;
        }
        let pv = p.vector();
        let mut sum = 0.0;
        for (a, b) in self.edges() {
            let c1 = pv.cross(a.vector());
            let c2 = pv.cross(b.vector());
            sum += c1.cross(c2).dot(pv).atan2(c1.dot(c2));
        }
        sum > PI
    }

    /// Perimeter — the sum of the great-circle edge lengths — as an [`Angle`].
    /// Multiply by a planet radius (see [`hausdorff_base::earth::sphere::Sphere`]) for metres.
    #[must_use]
    pub fn perimeter(&self) -> Angle {
        let total: f64 = self
            .edges()
            .map(|(a, b)| a.vector().angle(b.vector()).radians())
            .sum();
        Angle::from_radians(total)
    }

    /// True iff this loop and `other` overlap — their boundaries cross, or one
    /// loop contains the other (the nested case). `O(n·m)` in the edge counts.
    #[must_use]
    pub fn intersects(&self, other: &Loop) -> bool {
        if self.vertices.len() < 3 || other.vertices.len() < 3 {
            return false;
        }
        for (a, b) in self.edges() {
            for (c, d) in other.edges() {
                if crossing_sign(a, b, c, d) == Crossing::Cross {
                    return true;
                }
            }
        }
        self.contains_point(other.vertices[0]) || other.contains_point(self.vertices[0])
    }
}
