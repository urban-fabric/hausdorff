//! Planet-flexible metric layer: turn the generic, unit-less spherical geometry
//! (angles from [`crate::s1`], areas in steradians from [`crate::s2`]) into
//! physical lengths and areas on a sphere of a chosen radius.
//!
//! `s1`/`s2` live on the *unit* circle/sphere and carry no radius; [`crate::geo`]
//! is Earth-specific. [`Sphere`] is the general bridge — pick a radius (Earth,
//! Mars, the Moon, or any custom value) and it converts angles → metres and
//! steradians → m². This is the only place a radius (and hence a physical unit)
//! enters; the geometry itself stays unit-less and planet-agnostic.

use euclid::LengthMeters;

use crate::s1::Angle;

/// A sphere of a given radius — the bridge from unit-less spherical geometry to
/// physical distances and areas.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Sphere {
    radius_m: f64,
}

impl Sphere {
    /// Earth — GRS80 / WGS84 mean radius R₁, matching [`crate::geo::EARTH_RADIUS_M`].
    pub const EARTH: Sphere = Sphere {
        radius_m: crate::geo::EARTH_RADIUS_M,
    };
    /// Mars — IAU mean radius.
    pub const MARS: Sphere = Sphere {
        radius_m: 3_389_500.0,
    };
    /// The Moon — mean radius.
    pub const MOON: Sphere = Sphere {
        radius_m: 1_737_400.0,
    };

    /// A sphere of the given radius.
    #[must_use]
    pub fn new(radius: LengthMeters) -> Self {
        Self {
            radius_m: radius.as_meters(),
        }
    }

    /// A sphere of `radius_m` metres.
    #[must_use]
    pub const fn from_meters(radius_m: f64) -> Self {
        Self { radius_m }
    }

    /// The radius.
    #[must_use]
    pub fn radius(self) -> LengthMeters {
        LengthMeters::new(self.radius_m)
    }

    /// Arc length of a central `angle`: `radius × angle`. The general form of a
    /// great-circle distance — `Sphere::EARTH.length(...)` is the Earth case.
    #[must_use]
    pub fn length(self, angle: Angle) -> LengthMeters {
        LengthMeters::new(self.radius_m * angle.radians())
    }

    /// Surface area subtended by a solid angle in `steradians`: `radius² ×
    /// steradians` (e.g. a full sphere is `4π` steradians).
    #[must_use]
    pub fn area_m2(self, steradians: f64) -> f64 {
        self.radius_m * self.radius_m * steradians
    }
}
