use super::*;

/// A concrete map projection — the projected-CRS member of the enum, so a [`Crs`]
/// stays `Copy` and match-able (no `dyn Projection`). Phase 2 adds the conic /
/// stereographic / equal-area members.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum ProjectionKind {
    /// Spherical Web Mercator (EPSG:3857).
    WebMercator(WebMercator),
    /// Transverse Mercator (UTM and national grids).
    TransverseMercator(TransverseMercator),
    /// Lambert Conformal Conic, 2 standard parallels (national/regional grids).
    LambertConformalConic(LambertConformalConic),
    /// Polar Stereographic (the polar regions).
    PolarStereographic(PolarStereographic),
    /// Albers Equal-Area Conic (statistical/choropleth maps).
    AlbersEqualArea(AlbersEqualArea),
    /// Local tangent-plane (topocentric ENU) — a site-local East/North system.
    LocalTangentPlane(LocalTangentPlane),
}

impl ProjectionKind {
    /// Geographic (lat/lon) → projected planar metres.
    #[must_use]
    pub fn project(self, p: GeoPoint) -> Coord {
        match self {
            ProjectionKind::WebMercator(m) => m.project(p),
            ProjectionKind::TransverseMercator(m) => m.project(p),
            ProjectionKind::LambertConformalConic(m) => m.project(p),
            ProjectionKind::PolarStereographic(m) => m.project(p),
            ProjectionKind::AlbersEqualArea(m) => m.project(p),
            ProjectionKind::LocalTangentPlane(m) => m.project(p),
        }
    }

    /// Projected planar metres → geographic (lat/lon).
    #[must_use]
    pub fn inverse(self, c: Coord) -> GeoPoint {
        match self {
            ProjectionKind::WebMercator(m) => m.inverse(c),
            ProjectionKind::TransverseMercator(m) => m.inverse(c),
            ProjectionKind::LambertConformalConic(m) => m.inverse(c),
            ProjectionKind::PolarStereographic(m) => m.inverse(c),
            ProjectionKind::AlbersEqualArea(m) => m.inverse(c),
            ProjectionKind::LocalTangentPlane(m) => m.inverse(c),
        }
    }
}
