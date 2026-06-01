//! `s1` — one-dimensional spherical geometry: the [`Angle`] type (a measure on
//! the circle / great-circle distance), mirroring golang/geo's `s1` package.
//!
//! An [`Angle`] is a thin `f64` radians newtype that keeps angular quantities
//! from being confused with bare scalars, and converts cleanly to/from degrees.
//! Distances on the unit sphere are angles; multiply by an Earth radius (see
//! [`crate::geo::EARTH_RADIUS_M`]) to get metres.

use std::ops::{Add, Neg, Sub};

/// An angle, stored in **radians**.
#[derive(Clone, Copy, PartialEq, PartialOrd, Debug, Default)]
pub struct Angle(f64);

impl Angle {
    /// An angle of `radians`.
    #[must_use]
    pub const fn from_radians(radians: f64) -> Self {
        Self(radians)
    }

    /// An angle of `degrees`.
    #[must_use]
    pub fn from_degrees(degrees: f64) -> Self {
        Self(degrees.to_radians())
    }

    /// The angle in radians.
    #[inline]
    #[must_use]
    pub fn radians(self) -> f64 {
        self.0
    }

    /// The angle in degrees.
    #[inline]
    #[must_use]
    pub fn degrees(self) -> f64 {
        self.0.to_degrees()
    }

    /// The absolute value.
    #[inline]
    #[must_use]
    pub fn abs(self) -> Self {
        Self(self.0.abs())
    }

    /// Normalized into `(-π, π]`.
    #[must_use]
    pub fn normalized(self) -> Self {
        let mut r = self.0 % (2.0 * std::f64::consts::PI);
        if r <= -std::f64::consts::PI {
            r += 2.0 * std::f64::consts::PI;
        } else if r > std::f64::consts::PI {
            r -= 2.0 * std::f64::consts::PI;
        }
        Self(r)
    }
}

impl Add for Angle {
    type Output = Angle;
    fn add(self, rhs: Angle) -> Angle {
        Angle(self.0 + rhs.0)
    }
}

impl Sub for Angle {
    type Output = Angle;
    fn sub(self, rhs: Angle) -> Angle {
        Angle(self.0 - rhs.0)
    }
}

impl Neg for Angle {
    type Output = Angle;
    fn neg(self) -> Angle {
        Angle(-self.0)
    }
}
