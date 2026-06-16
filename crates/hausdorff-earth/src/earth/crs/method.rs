//! Shared **projection-method registry** for the WKT2 and PROJJSON codecs. Both
//! use the EPSG method names ("Transverse Mercator", …) and parameter names
//! ("Latitude of natural origin", …), so the method ↔ [`ProjectionKind`] mapping —
//! parse and emit — lives here once rather than in each codec.

use super::ProjectionKind;
use crate::earth::projection::{
    AlbersEqualArea, LambertConformalConic, PolarStereographic, TransverseMercator, WebMercator,
};
use hausdorff_base::earth::Ellipsoid;

/// The unit class of a projection parameter (for serialization).
pub(super) enum ParamUnit {
    /// An angle (degrees).
    Angle,
    /// A length (metres).
    Length,
    /// A dimensionless scale factor.
    Scale,
}

/// Build a [`ProjectionKind`] from an EPSG `method` name and a `param`-by-name
/// lookup — the parse side shared by the WKT2 and PROJJSON readers.
pub(super) fn projection_from_method(
    method: &str,
    ellipsoid: Ellipsoid,
    param: impl Fn(&str) -> Option<f64>,
) -> Option<ProjectionKind> {
    Some(match method.to_ascii_lowercase().as_str() {
        "transverse mercator" => ProjectionKind::TransverseMercator(TransverseMercator {
            ellipsoid,
            central_meridian: param("Longitude of natural origin")?,
            lat_origin: param("Latitude of natural origin").unwrap_or(0.0),
            k0: param("Scale factor at natural origin").unwrap_or(1.0),
            false_easting: param("False easting").unwrap_or(0.0),
            false_northing: param("False northing").unwrap_or(0.0),
        }),
        "lambert conic conformal (2sp)" => {
            ProjectionKind::LambertConformalConic(LambertConformalConic {
                ellipsoid,
                lat1: param("Latitude of 1st standard parallel")?,
                lat2: param("Latitude of 2nd standard parallel")?,
                lat0: param("Latitude of false origin").unwrap_or(0.0),
                lon0: param("Longitude of false origin")?,
                false_easting: param("Easting at false origin").unwrap_or(0.0),
                false_northing: param("Northing at false origin").unwrap_or(0.0),
            })
        }
        "lambert conic conformal (1sp)" => {
            let lat0 = param("Latitude of natural origin")?;
            ProjectionKind::LambertConformalConic(LambertConformalConic {
                ellipsoid,
                lat1: lat0,
                lat2: lat0,
                lat0,
                lon0: param("Longitude of natural origin")?,
                false_easting: param("False easting").unwrap_or(0.0),
                false_northing: param("False northing").unwrap_or(0.0),
            })
        }
        "albers equal area" => ProjectionKind::AlbersEqualArea(AlbersEqualArea {
            ellipsoid,
            lat1: param("Latitude of 1st standard parallel")?,
            lat2: param("Latitude of 2nd standard parallel")?,
            lat0: param("Latitude of false origin").unwrap_or(0.0),
            lon0: param("Longitude of false origin")?,
            false_easting: param("Easting at false origin").unwrap_or(0.0),
            false_northing: param("Northing at false origin").unwrap_or(0.0),
        }),
        "polar stereographic (variant b)" => {
            let lat_ts = param("Latitude of standard parallel")?;
            ProjectionKind::PolarStereographic(PolarStereographic {
                ellipsoid,
                lat_ts,
                lon0: param("Longitude of origin").unwrap_or(0.0),
                north: lat_ts > 0.0,
                false_easting: param("False easting").unwrap_or(0.0),
                false_northing: param("False northing").unwrap_or(0.0),
            })
        }
        "popular visualisation pseudo mercator" | "mercator (variant a)"
        | "mercator (variant b)" => {
            ProjectionKind::WebMercator(WebMercator { radius: ellipsoid.a })
        }
        _ => return None,
    })
}

/// The EPSG method name, method code, and `(name, value, unit)` parameter list for
/// a projection — the emit side shared by the WKT2 and PROJJSON serializers.
/// `None` for a projection with no standard method (the local tangent plane).
pub(super) fn conversion_params(
    p: ProjectionKind,
) -> Option<(&'static str, u32, Vec<(&'static str, f64, ParamUnit)>)> {
    use ParamUnit::{Angle, Length, Scale};
    Some(match p {
        ProjectionKind::TransverseMercator(m) => (
            "Transverse Mercator",
            9807,
            vec![
                ("Latitude of natural origin", m.lat_origin, Angle),
                ("Longitude of natural origin", m.central_meridian, Angle),
                ("Scale factor at natural origin", m.k0, Scale),
                ("False easting", m.false_easting, Length),
                ("False northing", m.false_northing, Length),
            ],
        ),
        ProjectionKind::LambertConformalConic(m) => (
            "Lambert Conic Conformal (2SP)",
            9802,
            vec![
                ("Latitude of false origin", m.lat0, Angle),
                ("Longitude of false origin", m.lon0, Angle),
                ("Latitude of 1st standard parallel", m.lat1, Angle),
                ("Latitude of 2nd standard parallel", m.lat2, Angle),
                ("Easting at false origin", m.false_easting, Length),
                ("Northing at false origin", m.false_northing, Length),
            ],
        ),
        ProjectionKind::AlbersEqualArea(m) => (
            "Albers Equal Area",
            9822,
            vec![
                ("Latitude of false origin", m.lat0, Angle),
                ("Longitude of false origin", m.lon0, Angle),
                ("Latitude of 1st standard parallel", m.lat1, Angle),
                ("Latitude of 2nd standard parallel", m.lat2, Angle),
                ("Easting at false origin", m.false_easting, Length),
                ("Northing at false origin", m.false_northing, Length),
            ],
        ),
        ProjectionKind::PolarStereographic(m) => (
            "Polar Stereographic (variant B)",
            9829,
            vec![
                ("Latitude of standard parallel", m.lat_ts, Angle),
                ("Longitude of origin", m.lon0, Angle),
                ("False easting", m.false_easting, Length),
                ("False northing", m.false_northing, Length),
            ],
        ),
        ProjectionKind::WebMercator(_) => (
            "Popular Visualisation Pseudo Mercator",
            1024,
            vec![
                ("Latitude of natural origin", 0.0, Angle),
                ("Longitude of natural origin", 0.0, Angle),
                ("False easting", 0.0, Length),
                ("False northing", 0.0, Length),
            ],
        ),
        ProjectionKind::LocalTangentPlane(_) => return None,
    })
}
