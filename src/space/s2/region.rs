//! `s2::Region` — the common interface over spherical regions ([`Cap`], [`Loop`],
//! [`Polygon`]), mirroring golang/geo's `Region`. A region can test point
//! containment and report a bounding [`Cap`]; the bounding cap is the hook a
//! future cell coverer uses to approximate a region by S2 cells.

use crate::space::r3::Vector;
use crate::space::s1::Angle;

use super::{Cap, Loop, Polygon, S2Point};

/// A region on the unit sphere.
pub trait Region {
    /// True iff the region contains `p`.
    fn contains_point(&self, p: S2Point) -> bool;

    /// A spherical cap that bounds the region (not necessarily the tightest).
    fn bounding_cap(&self) -> Cap;
}

/// The bounding cap of a set of points: centered on their normalized vector
/// mean, with radius equal to the farthest point's angular distance.
fn bounding_cap_of(points: &[S2Point]) -> Cap {
    if points.is_empty() {
        return Cap::empty(S2Point::from_vector(Vector::new(1.0, 0.0, 0.0)));
    }
    let mut sum = Vector::new(0.0, 0.0, 0.0);
    for p in points {
        sum = sum.add(p.vector());
    }
    let mean = sum.normalize();
    // Antipodal points can cancel to the zero vector; fall back to the first.
    let center = if mean.norm2() == 0.0 {
        points[0]
    } else {
        S2Point::from_vector(mean)
    };
    let radius = points
        .iter()
        .map(|p| center.vector().angle(p.vector()).radians())
        .fold(0.0_f64, f64::max);
    Cap::from_center_angle(center, Angle::from_radians(radius))
}

impl Region for Cap {
    fn contains_point(&self, p: S2Point) -> bool {
        Cap::contains_point(*self, p)
    }
    fn bounding_cap(&self) -> Cap {
        *self
    }
}

impl Region for Loop {
    fn contains_point(&self, p: S2Point) -> bool {
        Loop::contains_point(self, p)
    }
    fn bounding_cap(&self) -> Cap {
        bounding_cap_of(self.vertices())
    }
}

impl Region for Polygon {
    fn contains_point(&self, p: S2Point) -> bool {
        Polygon::contains_point(self, p)
    }
    fn bounding_cap(&self) -> Cap {
        bounding_cap_of(self.shell().vertices())
    }
}
