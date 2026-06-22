//! Compass bearing newtype.

use serde::{Deserialize, Serialize};

use euclid::Angle as EuclidAngle;

use super::GeoPoint;

/// A compass bearing in degrees: `0` = North, `90` = East, increasing
/// clockwise. Construction normalizes into `[0, 360)`.
#[derive(Clone, Copy, PartialEq, Debug, Serialize, Deserialize)]
pub struct Bearing(f64);

impl Bearing {
    /// A bearing of `degrees`, normalized into `[0, 360)`.
    #[must_use]
    pub fn from_degrees(degrees: f64) -> Self {
        Self(degrees.rem_euclid(360.0))
    }

    /// A bearing from a generic Euclid angle, normalized into `[0, 360)`.
    #[must_use]
    pub fn from_angle(angle: EuclidAngle) -> Self {
        Self::from_degrees(angle.degrees())
    }

    /// The initial great-circle bearing when travelling from `from` toward `to`.
    ///
    /// Proof-backed by `riemann` `Hausdorff/BearingAzimuth.initialBearing_is_azimuth`:
    /// the reported `atan2(y, x)` is the *literal* compass azimuth — its cosine and
    /// sine are the geodesic tangent's normalised North (`x`) and East (`y`)
    /// components in the local ENU frame (building on `Hausdorff/Bearing`'s proof
    /// that `(y, x)` are those tangent components).
    #[must_use]
    pub fn initial(from: GeoPoint, to: GeoPoint) -> Self {
        let lat1 = from.lat.to_radians();
        let lat2 = to.lat.to_radians();
        let dlon = (to.lon - from.lon).to_radians();
        let y = dlon.sin() * lat2.cos();
        let x = lat1.cos() * lat2.sin() - lat1.sin() * lat2.cos() * dlon.cos();
        Self::from_degrees(y.atan2(x).to_degrees())
    }

    /// The bearing in degrees, in `[0, 360)`.
    #[inline]
    #[must_use]
    pub fn degrees(self) -> f64 {
        self.0
    }

    /// The bearing as a generic Euclid angle.
    #[inline]
    #[must_use]
    pub fn angle(self) -> EuclidAngle {
        EuclidAngle::from_degrees(self.0)
    }

    /// The bearing in radians, in `[0, 2π)`.
    #[inline]
    #[must_use]
    pub fn radians(self) -> f64 {
        self.angle().radians()
    }
}
