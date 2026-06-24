use super::*;

/// A coordinate value tagged by which kind of CRS it belongs to — geographic
/// (lat/lon [`GeoPoint`] + ellipsoidal height), projected planar metres
/// ([`Coord`]), or geocentric ECEF ([`EcefPoint`]). [`transform`] consumes and
/// produces these.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum CrsCoord {
    /// Lat/lon on a geographic CRS, with **ellipsoidal height** in metres. A
    /// 2-D geographic coordinate uses `height = 0`; a 3-D one (e.g. EPSG:4979)
    /// carries the height through the [`transform`] pivot. Build with
    /// [`geographic`](Self::geographic) / [`geographic_with_height`](Self::geographic_with_height).
    Geographic { point: GeoPoint, height: f64 },
    /// Easting/northing on a projected CRS (planar metres; no height).
    Projected(Coord),
    /// X/Y/Z on a geocentric CRS.
    Geocentric(EcefPoint),
}

impl CrsCoord {
    /// A 2-D geographic coordinate — lat/lon at ellipsoidal height 0.
    #[must_use]
    pub fn geographic(point: GeoPoint) -> CrsCoord {
        CrsCoord::Geographic { point, height: 0.0 }
    }

    /// A 3-D geographic coordinate — lat/lon at `height` metres above the ellipsoid.
    #[must_use]
    pub fn geographic_with_height(point: GeoPoint, height: f64) -> CrsCoord {
        CrsCoord::Geographic { point, height }
    }

    /// The geographic lat/lon (panics if this is not a geographic coordinate).
    #[must_use]
    pub fn point(self) -> GeoPoint {
        self.as_geographic().0
    }

    /// The ellipsoidal height in metres (panics if this is not geographic).
    #[must_use]
    pub fn height(self) -> f64 {
        self.as_geographic().1
    }

    /// The projected easting/northing (panics if this is not projected).
    #[must_use]
    pub fn projected(self) -> Coord {
        match self {
            CrsCoord::Projected(c) => c,
            _ => panic!("expected a Projected coordinate for a projected CRS"),
        }
    }

    /// The geocentric ECEF point (panics if this is not geocentric).
    #[must_use]
    pub fn geocentric(self) -> EcefPoint {
        match self {
            CrsCoord::Geocentric(p) => p,
            _ => panic!("expected a Geocentric coordinate for a geocentric CRS"),
        }
    }

    /// (lat/lon, ellipsoidal height) — the pivot's geographic accessor.
    pub(super) fn as_geographic(self) -> (GeoPoint, f64) {
        match self {
            CrsCoord::Geographic { point, height } => (point, height),
            _ => panic!("expected a Geographic coordinate for a geographic CRS"),
        }
    }
}

/// The dimensionality of a geographic CRS: **2-D** (lat/lon, e.g. EPSG:4326 — no
/// height) or **3-D** (lat/lon + ellipsoidal height, e.g. EPSG:4979). The pivot
/// carries height through ECEF regardless; `dim` decides whether the *output* of
/// a geographic CRS retains it (`D3`) or drops it to 0 (`D2`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum GeoDim {
    /// Geographic 2-D — lat/lon only.
    D2,
    /// Geographic 3-D — lat/lon + ellipsoidal height (metres).
    D3,
}
