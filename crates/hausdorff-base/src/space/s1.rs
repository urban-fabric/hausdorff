//! `s1` — one-dimensional spherical geometry: angular measure. The home for
//! angles in the workspace.
//!
//! Angles live in hausdorff rather than `euclid` on purpose: unlike metres or
//! seconds they carry **no physical unit**, only a choice of representation
//! (radians or degrees) — so they are pure geometry. [`Angle`] is the canonical
//! radians-backed type the geometry computes in (a great-circle distance on the
//! unit sphere *is* an angle); [`AngleDegrees`] is the degrees representation for
//! I/O and human-facing values. Mirrors golang/geo's `s1.Angle` and subsumes
//! `euclid`'s `AngleRadians`/`AngleDegrees`.

#![allow(unsafe_code)]

use std::ops::{Add, Div, Mul, Neg, Sub};

use bytemuck::{Pod, Zeroable};
use serde::{Deserialize, Serialize};

/// An angle in **radians** — the canonical angular measure for the geometry
/// code. POD layout + transparent serde, so it stores and serializes exactly
/// like the bare `f64` it wraps.
#[repr(transparent)]
#[derive(
    Copy, Clone, PartialEq, PartialOrd, Debug, Default, Pod, Zeroable, Serialize, Deserialize,
)]
#[serde(transparent)]
pub struct Angle(f64);

impl Angle {
    /// Zero radians.
    pub const ZERO: Self = Self(0.0);
    /// π radians (a straight angle).
    pub const PI: Self = Self(std::f64::consts::PI);
    /// 2π radians (a full turn).
    pub const TAU: Self = Self(std::f64::consts::TAU);
    /// A positive-infinite angle (a sentinel for "no/unbounded distance").
    pub const INFINITY: Self = Self(f64::INFINITY);

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
    pub const fn radians(self) -> f64 {
        self.0
    }

    /// The angle in degrees as a bare `f64`.
    #[inline]
    #[must_use]
    pub fn degrees(self) -> f64 {
        self.0.to_degrees()
    }

    /// The angle as a typed [`AngleDegrees`].
    #[must_use]
    pub fn to_degrees(self) -> AngleDegrees {
        AngleDegrees(self.0.to_degrees())
    }

    /// The radians value narrowed to `f32`.
    #[must_use]
    pub fn to_f32_radians(self) -> f32 {
        self.0 as f32
    }

    /// The absolute value.
    #[inline]
    #[must_use]
    pub fn abs(self) -> Self {
        Self(self.0.abs())
    }

    /// `sin` / `cos` / `tan` of the angle (radians need no extra conversion).
    #[inline]
    #[must_use]
    pub fn sin(self) -> f64 {
        self.0.sin()
    }
    #[inline]
    #[must_use]
    pub fn cos(self) -> f64 {
        self.0.cos()
    }
    #[inline]
    #[must_use]
    pub fn tan(self) -> f64 {
        self.0.tan()
    }

    /// Normalized into `(-π, π]`.
    #[must_use]
    pub fn normalized(self) -> Self {
        let mut r = self.0 % std::f64::consts::TAU;
        if r <= -std::f64::consts::PI {
            r += std::f64::consts::TAU;
        } else if r > std::f64::consts::PI {
            r -= std::f64::consts::TAU;
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

impl Mul<f64> for Angle {
    type Output = Angle;
    fn mul(self, rhs: f64) -> Angle {
        Angle(self.0 * rhs)
    }
}

impl Div<f64> for Angle {
    type Output = Angle;
    fn div(self, rhs: f64) -> Angle {
        Angle(self.0 / rhs)
    }
}

/// An angle in **degrees** — the representation for I/O and human-facing values.
/// Convert to [`Angle`] (radians) before doing trigonometry.
#[repr(transparent)]
#[derive(
    Copy, Clone, PartialEq, PartialOrd, Debug, Default, Pod, Zeroable, Serialize, Deserialize,
)]
#[serde(transparent)]
pub struct AngleDegrees(f64);

impl AngleDegrees {
    /// Zero degrees.
    pub const ZERO: Self = Self(0.0);

    /// An angle of `degrees`.
    #[must_use]
    pub const fn new(degrees: f64) -> Self {
        Self(degrees)
    }

    /// The value in degrees.
    #[inline]
    #[must_use]
    pub const fn as_degrees(self) -> f64 {
        self.0
    }

    /// As a radians-backed [`Angle`].
    #[must_use]
    pub fn to_radians(self) -> Angle {
        Angle(self.0.to_radians())
    }
}

impl From<Angle> for AngleDegrees {
    fn from(a: Angle) -> Self {
        a.to_degrees()
    }
}

impl From<AngleDegrees> for Angle {
    fn from(d: AngleDegrees) -> Self {
        d.to_radians()
    }
}
