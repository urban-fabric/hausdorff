//! The atomic planar types: a 2-D [`Coord`]inate and a [`Point`].
//!
//! Both are `f64`-scalar and `Copy`. `Coord` is the bare `(x, y)` pair used
//! throughout the planar geometry types; `Point` is the single-coordinate
//! geometry (a `Coord` newtype) that carries geometry semantics.

use std::ops::{Add, Div, Mul, Neg, Sub};

use serde::{Deserialize, Serialize};

/// A 2-D coordinate `(x, y)` in an unspecified Cartesian plane. The planar
/// counterpart to the spherical [`crate::geo::GeoPoint`]; by GIS convention an
/// interchange maps `x = lon`, `y = lat` (see [`crate::convert`]).
#[derive(Copy, Clone, PartialEq, Debug, Default, Serialize, Deserialize)]
pub struct Coord {
    pub x: f64,
    pub y: f64,
}

impl Coord {
    /// A coordinate at `(x, y)`.
    #[must_use]
    pub const fn new(x: f64, y: f64) -> Self {
        Self { x, y }
    }

    /// The `(x, y)` pair.
    #[inline]
    #[must_use]
    pub fn x_y(self) -> (f64, f64) {
        (self.x, self.y)
    }

    /// The 2-D cross product `self × rhs` about `origin`
    /// (`(self - origin) × (rhs - origin)`): positive iff `origin → self → rhs`
    /// turns counter-clockwise. The orientation kernel the planar algorithms
    /// share.
    #[inline]
    #[must_use]
    pub fn cross(self, origin: Coord, rhs: Coord) -> f64 {
        (self.x - origin.x) * (rhs.y - origin.y) - (self.y - origin.y) * (rhs.x - origin.x)
    }

    /// Dot product of the two coordinates treated as vectors from the origin.
    #[inline]
    #[must_use]
    pub fn dot(self, rhs: Coord) -> f64 {
        self.x * rhs.x + self.y * rhs.y
    }

    /// True iff both components are finite.
    #[inline]
    #[must_use]
    pub fn is_finite(self) -> bool {
        self.x.is_finite() && self.y.is_finite()
    }
}

impl From<(f64, f64)> for Coord {
    fn from((x, y): (f64, f64)) -> Self {
        Self { x, y }
    }
}

impl From<[f64; 2]> for Coord {
    fn from([x, y]: [f64; 2]) -> Self {
        Self { x, y }
    }
}

impl From<Coord> for (f64, f64) {
    fn from(c: Coord) -> Self {
        (c.x, c.y)
    }
}

impl Add for Coord {
    type Output = Coord;
    fn add(self, rhs: Coord) -> Coord {
        Coord::new(self.x + rhs.x, self.y + rhs.y)
    }
}

impl Sub for Coord {
    type Output = Coord;
    fn sub(self, rhs: Coord) -> Coord {
        Coord::new(self.x - rhs.x, self.y - rhs.y)
    }
}

impl Mul<f64> for Coord {
    type Output = Coord;
    fn mul(self, rhs: f64) -> Coord {
        Coord::new(self.x * rhs, self.y * rhs)
    }
}

impl Div<f64> for Coord {
    type Output = Coord;
    fn div(self, rhs: f64) -> Coord {
        Coord::new(self.x / rhs, self.y / rhs)
    }
}

impl Neg for Coord {
    type Output = Coord;
    fn neg(self) -> Coord {
        Coord::new(-self.x, -self.y)
    }
}

/// A single planar point — a [`Coord`] with geometry semantics. Mirrors
/// `geo_types::Point` so consumer migration is a mechanical import swap.
#[derive(Copy, Clone, PartialEq, Debug, Default, Serialize, Deserialize)]
pub struct Point(pub Coord);

impl Point {
    /// The point at `(x, y)`.
    #[must_use]
    pub const fn new(x: f64, y: f64) -> Self {
        Self(Coord::new(x, y))
    }

    /// The x ordinate.
    #[inline]
    #[must_use]
    pub fn x(self) -> f64 {
        self.0.x
    }

    /// The y ordinate.
    #[inline]
    #[must_use]
    pub fn y(self) -> f64 {
        self.0.y
    }

    /// The `(x, y)` pair.
    #[inline]
    #[must_use]
    pub fn x_y(self) -> (f64, f64) {
        self.0.x_y()
    }
}

impl From<Coord> for Point {
    fn from(c: Coord) -> Self {
        Self(c)
    }
}

impl From<Point> for Coord {
    fn from(p: Point) -> Self {
        p.0
    }
}

impl From<(f64, f64)> for Point {
    fn from((x, y): (f64, f64)) -> Self {
        Self(Coord::new(x, y))
    }
}
