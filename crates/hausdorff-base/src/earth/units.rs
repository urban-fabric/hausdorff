//! Concrete Earth-facing unit carriers.
//!
//! Current `euclid` owns the unit contracts, while consumer crates own concrete
//! carriers. These wrappers preserve Hausdorff's public metre-typed API.

use std::ops::{Add, Div, Mul, Neg, Sub};

use euclid::{Dimension, LinearQuantity, Measure, ScalarQuantity};
use serde::{Deserialize, Serialize};

/// Length in metres.
#[repr(transparent)]
#[derive(Clone, Copy, PartialEq, PartialOrd, Debug, Default, Serialize, Deserialize)]
#[serde(transparent)]
pub struct LengthMeters<T = f64>(T);

impl<T> LengthMeters<T> {
    /// Wrap a raw metre value.
    #[must_use]
    pub const fn new(value: T) -> Self {
        Self(value)
    }

    /// Unwrap the raw metre value.
    #[must_use]
    pub fn value(self) -> T {
        self.0
    }

    /// Borrow the raw metre value.
    #[must_use]
    pub const fn as_value(&self) -> &T {
        &self.0
    }
}

/// Area in square metres.
#[repr(transparent)]
#[derive(Clone, Copy, PartialEq, PartialOrd, Debug, Default, Serialize, Deserialize)]
#[serde(transparent)]
pub struct AreaSquareMeters<T = f64>(T);

impl<T> AreaSquareMeters<T> {
    /// Wrap a raw square-metre value.
    #[must_use]
    pub const fn new(value: T) -> Self {
        Self(value)
    }

    /// Unwrap the raw square-metre value.
    #[must_use]
    pub fn value(self) -> T {
        self.0
    }

    /// Borrow the raw square-metre value.
    #[must_use]
    pub const fn as_value(&self) -> &T {
        &self.0
    }
}

impl<T> Add for LengthMeters<T>
where
    T: Add<Output = T>,
{
    type Output = Self;

    fn add(self, rhs: Self) -> Self {
        Self(self.0 + rhs.0)
    }
}

impl<T> Sub for LengthMeters<T>
where
    T: Sub<Output = T>,
{
    type Output = Self;

    fn sub(self, rhs: Self) -> Self {
        Self(self.0 - rhs.0)
    }
}

impl<T> Neg for LengthMeters<T>
where
    T: Neg<Output = T>,
{
    type Output = Self;

    fn neg(self) -> Self {
        Self(-self.0)
    }
}

impl<T> Mul<T> for LengthMeters<T>
where
    T: Mul<Output = T>,
{
    type Output = Self;

    fn mul(self, rhs: T) -> Self {
        Self(self.0 * rhs)
    }
}

impl<T> Div<T> for LengthMeters<T>
where
    T: Div<Output = T>,
{
    type Output = Self;

    fn div(self, rhs: T) -> Self {
        Self(self.0 / rhs)
    }
}

impl<T> Div for LengthMeters<T>
where
    T: Div<Output = T>,
{
    type Output = T;

    fn div(self, rhs: Self) -> T {
        self.0 / rhs.0
    }
}

impl<T> Mul for LengthMeters<T>
where
    T: Mul<Output = T>,
{
    type Output = AreaSquareMeters<T>;

    fn mul(self, rhs: Self) -> AreaSquareMeters<T> {
        AreaSquareMeters(self.0 * rhs.0)
    }
}

impl<T> Add for AreaSquareMeters<T>
where
    T: Add<Output = T>,
{
    type Output = Self;

    fn add(self, rhs: Self) -> Self {
        Self(self.0 + rhs.0)
    }
}

impl<T> Sub for AreaSquareMeters<T>
where
    T: Sub<Output = T>,
{
    type Output = Self;

    fn sub(self, rhs: Self) -> Self {
        Self(self.0 - rhs.0)
    }
}

impl<T> Neg for AreaSquareMeters<T>
where
    T: Neg<Output = T>,
{
    type Output = Self;

    fn neg(self) -> Self {
        Self(-self.0)
    }
}

impl<T> Mul<T> for AreaSquareMeters<T>
where
    T: Mul<Output = T>,
{
    type Output = Self;

    fn mul(self, rhs: T) -> Self {
        Self(self.0 * rhs)
    }
}

impl<T> Div<T> for AreaSquareMeters<T>
where
    T: Div<Output = T>,
{
    type Output = Self;

    fn div(self, rhs: T) -> Self {
        Self(self.0 / rhs)
    }
}

impl<T> Div<LengthMeters<T>> for AreaSquareMeters<T>
where
    T: Div<Output = T>,
{
    type Output = LengthMeters<T>;

    fn div(self, rhs: LengthMeters<T>) -> LengthMeters<T> {
        LengthMeters(self.0 / rhs.0)
    }
}

impl Measure for LengthMeters {
    type Repr = f64;

    const DIMENSION: Dimension = Dimension::Length;
    const UNIT_NAME: &'static str = "m";

    fn raw(&self) -> f64 {
        self.0
    }

    fn from_raw(v: f64) -> Self {
        Self(v)
    }
}

impl ScalarQuantity for LengthMeters {
    fn min(self, other: Self) -> Self {
        Self(self.0.min(other.0))
    }

    fn max(self, other: Self) -> Self {
        Self(self.0.max(other.0))
    }

    fn clamp(self, lo: Self, hi: Self) -> Self {
        Self(self.0.clamp(lo.0, hi.0))
    }

    fn floor(self) -> Self {
        Self(self.0.floor())
    }

    fn ceil(self) -> Self {
        Self(self.0.ceil())
    }

    fn round(self) -> Self {
        Self(self.0.round())
    }

    fn trunc(self) -> Self {
        Self(self.0.trunc())
    }

    fn is_finite(&self) -> bool {
        self.0.is_finite()
    }

    fn is_nan(&self) -> bool {
        self.0.is_nan()
    }

    fn is_infinite(&self) -> bool {
        self.0.is_infinite()
    }

    fn approx_eq(&self, other: &Self, tolerance: &Self) -> bool {
        (self.0 - other.0).abs() <= tolerance.0.abs()
    }
}

impl LinearQuantity for LengthMeters {
    fn add(self, other: Self) -> Self {
        self + other
    }

    fn sub(self, other: Self) -> Self {
        self - other
    }

    fn neg(self) -> Self {
        -self
    }

    fn abs(self) -> Self {
        Self(self.0.abs())
    }

    fn signum(self) -> f64 {
        self.0.signum()
    }

    fn ratio(self, other: Self) -> f64 {
        self / other
    }

    fn lerp(self, other: Self, t: f64) -> Self {
        Self(self.0 * (1.0 - t) + other.0 * t)
    }

    fn scale(self, factor: f64) -> Self {
        self * factor
    }

    fn scaled(self, k: f64) -> Self {
        self * k
    }
}

impl Measure for AreaSquareMeters {
    type Repr = f64;

    const DIMENSION: Dimension = Dimension::Area;
    const UNIT_NAME: &'static str = "m^2";

    fn raw(&self) -> f64 {
        self.0
    }

    fn from_raw(v: f64) -> Self {
        Self(v)
    }
}

impl ScalarQuantity for AreaSquareMeters {
    fn min(self, other: Self) -> Self {
        Self(self.0.min(other.0))
    }

    fn max(self, other: Self) -> Self {
        Self(self.0.max(other.0))
    }

    fn clamp(self, lo: Self, hi: Self) -> Self {
        Self(self.0.clamp(lo.0, hi.0))
    }

    fn floor(self) -> Self {
        Self(self.0.floor())
    }

    fn ceil(self) -> Self {
        Self(self.0.ceil())
    }

    fn round(self) -> Self {
        Self(self.0.round())
    }

    fn trunc(self) -> Self {
        Self(self.0.trunc())
    }

    fn is_finite(&self) -> bool {
        self.0.is_finite()
    }

    fn is_nan(&self) -> bool {
        self.0.is_nan()
    }

    fn is_infinite(&self) -> bool {
        self.0.is_infinite()
    }

    fn approx_eq(&self, other: &Self, tolerance: &Self) -> bool {
        (self.0 - other.0).abs() <= tolerance.0.abs()
    }
}

impl LinearQuantity for AreaSquareMeters {
    fn add(self, other: Self) -> Self {
        self + other
    }

    fn sub(self, other: Self) -> Self {
        self - other
    }

    fn neg(self) -> Self {
        -self
    }

    fn abs(self) -> Self {
        Self(self.0.abs())
    }

    fn signum(self) -> f64 {
        self.0.signum()
    }

    fn ratio(self, other: Self) -> f64 {
        self.0 / other.0
    }

    fn lerp(self, other: Self, t: f64) -> Self {
        Self(self.0 * (1.0 - t) + other.0 * t)
    }

    fn scale(self, factor: f64) -> Self {
        self * factor
    }

    fn scaled(self, k: f64) -> Self {
        self * k
    }
}
