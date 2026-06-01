//! Minimal ℝ³ vector backing the S2 cube-sphere projection — a private
//! analogue of golang/geo's `r3.Vector` carrying only what the cell-indexing
//! core needs.

#[derive(Clone, Copy, PartialEq, Debug)]
pub(super) struct R3Vector {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

/// The three axes of ℝ³.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub(super) enum Axis {
    X,
    Y,
    Z,
}

impl R3Vector {
    pub(super) const fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z }
    }

    pub(super) fn dot(self, o: Self) -> f64 {
        self.x * o.x + self.y * o.y + self.z * o.z
    }

    pub(super) fn norm2(self) -> f64 {
        self.dot(self)
    }

    /// A unit vector in the same direction, or the zero vector if `self` is
    /// zero (matching `r3.Vector.Normalize`).
    pub(super) fn normalize(self) -> Self {
        let n2 = self.norm2();
        if n2 == 0.0 {
            return Self::new(0.0, 0.0, 0.0);
        }
        let inv = 1.0 / n2.sqrt();
        Self::new(self.x * inv, self.y * inv, self.z * inv)
    }

    /// The axis of the largest-magnitude component. Ties resolve toward `Z`
    /// then `Y` (strict `>`), exactly as `r3.Vector.LargestComponent`.
    pub(super) fn largest_component(self) -> Axis {
        let (tx, ty, tz) = (self.x.abs(), self.y.abs(), self.z.abs());
        if tx > ty {
            if tx > tz { Axis::X } else { Axis::Z }
        } else if ty > tz {
            Axis::Y
        } else {
            Axis::Z
        }
    }
}
