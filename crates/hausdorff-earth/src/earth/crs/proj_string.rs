//! Parsing **PROJ-strings** (`+proj=… +datum=…`) into a [`Crs`].
//!
//! [`Crs::from_proj_string`] reads the PROJ pipeline-free single-CRS form — the
//! `+proj`/`+datum`/`+ellps`/`+towgs84` parameters PROJ emits for a CRS — and
//! builds the matching [`Crs`]. Every PROJ projection parameter maps onto a field
//! this crate already models (`+lat_0`→`lat_origin`, `+lon_0`→`central_meridian`,
//! `+k`→`k0`, `+x_0`/`+y_0`→false easting/northing, `+lat_1`/`+lat_2`→standard
//! parallels), so parsing is a direct translation, not a re-derivation.
//!
//! Supported `+proj`: `longlat` (geographic), `geocent` (geocentric), `merc`
//! (Web Mercator), `utm`, `tmerc` (Transverse Mercator), `lcc` (Lambert Conformal
//! Conic), `aea` (Albers), `stere` (**polar** stereographic only). The datum comes
//! from `+datum`, or `+ellps` (+ optional `+towgs84` Helmert), or an explicit
//! `+a`/`+b`/`+f`/`+rf` ellipsoid; `+init=epsg:NNNN` defers to [`Crs::from_epsg`].
//! Units are assumed metres. The PROJ→PROJ-string conversion is lossy (it drops
//! the 2-D/3-D distinction, axis order, etc.), so this is best-effort by design.

use super::{Crs, Datum, DatumShift, GeoDim, Helmert, ProjectionKind, datum};
use crate::earth::projection::{
    AlbersEqualArea, LambertConformalConic, PolarStereographic, TransverseMercator, WebMercator,
};
use hausdorff_base::earth::Ellipsoid;

impl Crs {
    /// Parse a **PROJ-string** into a [`Crs`], or `None` if it is unrecognised or
    /// uses a projection/datum this crate does not model.
    ///
    /// The resulting CRS carries no EPSG code (`epsg() == None`) — it is defined by
    /// its parameters. `+init=epsg:NNNN` is the exception: it defers to
    /// [`from_epsg`](Self::from_epsg), keeping the code.
    ///
    /// ```
    /// # use hausdorff_earth::earth::crs::Crs;
    /// let bng = Crs::from_proj_string(
    ///     "+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 \
    ///      +x_0=400000 +y_0=-100000 +ellps=airy +units=m +no_defs",
    /// );
    /// assert!(matches!(bng, Some(Crs::Projected { .. })));
    /// ```
    #[must_use]
    pub fn from_proj_string(s: &str) -> Option<Crs> {
        let p = ProjParams::parse(s);
        // `+init=epsg:NNNN` — defer to the registry (and keep the code).
        if let Some(init) = p.get("init") {
            let code = init.strip_prefix("epsg:").or_else(|| init.strip_prefix("EPSG:"))?;
            return Crs::from_epsg(code.parse().ok()?);
        }
        let datum = p.datum()?;
        let ell = datum.ellipsoid;
        let projected = |projection| {
            Some(Crs::Projected {
                datum,
                projection,
                epsg: None,
            })
        };
        match p.get("proj")? {
            "longlat" | "latlong" | "latlon" | "lonlat" => Some(Crs::Geographic {
                datum,
                epsg: None,
                dim: GeoDim::D2,
            }),
            "geocent" => Some(Crs::Geocentric { datum, epsg: None }),
            "merc" => projected(ProjectionKind::WebMercator(WebMercator {
                radius: p.f("a").unwrap_or(ell.a),
            })),
            "utm" => {
                let zone = p.f("zone")? as u8;
                if !(1..=TransverseMercator::UTM_ZONES).contains(&zone) {
                    return None;
                }
                let mut tm = TransverseMercator::utm(zone, !p.flag("south"));
                tm.ellipsoid = ell;
                projected(ProjectionKind::TransverseMercator(tm))
            }
            "tmerc" => projected(ProjectionKind::TransverseMercator(TransverseMercator {
                ellipsoid: ell,
                central_meridian: p.f("lon_0").unwrap_or(0.0),
                lat_origin: p.f("lat_0").unwrap_or(0.0),
                k0: p.k0(),
                false_easting: p.f("x_0").unwrap_or(0.0),
                false_northing: p.f("y_0").unwrap_or(0.0),
            })),
            "lcc" => {
                let lat1 = p.f("lat_1")?;
                projected(ProjectionKind::LambertConformalConic(LambertConformalConic {
                    ellipsoid: ell,
                    lat1,
                    lat2: p.f("lat_2").unwrap_or(lat1),
                    lat0: p.f("lat_0").unwrap_or(lat1),
                    lon0: p.f("lon_0").unwrap_or(0.0),
                    false_easting: p.f("x_0").unwrap_or(0.0),
                    false_northing: p.f("y_0").unwrap_or(0.0),
                }))
            }
            "aea" => projected(ProjectionKind::AlbersEqualArea(AlbersEqualArea {
                ellipsoid: ell,
                lat1: p.f("lat_1")?,
                lat2: p.f("lat_2")?,
                lat0: p.f("lat_0").unwrap_or(0.0),
                lon0: p.f("lon_0").unwrap_or(0.0),
                false_easting: p.f("x_0").unwrap_or(0.0),
                false_northing: p.f("y_0").unwrap_or(0.0),
            })),
            "stere" => {
                // Polar stereographic only — the latitude of origin is a pole.
                let lat0 = p.f("lat_0")?;
                if lat0.abs() < 89.9 {
                    return None;
                }
                projected(ProjectionKind::PolarStereographic(PolarStereographic {
                    ellipsoid: ell,
                    lat_ts: p.f("lat_ts").unwrap_or(lat0),
                    lon0: p.f("lon_0").unwrap_or(0.0),
                    north: lat0 > 0.0,
                    false_easting: p.f("x_0").unwrap_or(0.0),
                    false_northing: p.f("y_0").unwrap_or(0.0),
                }))
            }
            _ => None,
        }
    }
}

/// A parsed bag of PROJ `+key[=value]` parameters (order preserved; first match
/// wins). Borrows from the input string.
struct ProjParams<'a> {
    items: Vec<(&'a str, Option<&'a str>)>,
}

impl<'a> ProjParams<'a> {
    fn parse(s: &'a str) -> Self {
        let items = s
            .split_whitespace()
            .filter_map(|tok| {
                let tok = tok.strip_prefix('+').unwrap_or(tok);
                if tok.is_empty() {
                    return None;
                }
                Some(match tok.split_once('=') {
                    Some((k, v)) => (k, Some(v)),
                    None => (tok, None),
                })
            })
            .collect();
        ProjParams { items }
    }

    /// The value of `+key=value`, or `None` if absent or a bare flag.
    fn get(&self, key: &str) -> Option<&'a str> {
        self.items.iter().find(|(k, _)| *k == key).and_then(|(_, v)| *v)
    }

    /// Whether a bare `+key` flag (e.g. `+south`) is present.
    fn flag(&self, key: &str) -> bool {
        self.items.iter().any(|(k, _)| *k == key)
    }

    /// A `+key=value` parsed as `f64`.
    fn f(&self, key: &str) -> Option<f64> {
        self.get(key)?.parse().ok()
    }

    /// The central-meridian scale factor `k0` (`+k_0`, else `+k`, else 1).
    fn k0(&self) -> f64 {
        self.f("k_0").or_else(|| self.f("k")).unwrap_or(1.0)
    }

    /// Resolve the [`Datum`] from `+datum`, or `+ellps` (+ optional `+towgs84`),
    /// or an explicit `+a`/`+b`/`+f`/`+rf` ellipsoid (assumed WGS84-aligned).
    fn datum(&self) -> Option<Datum> {
        if let Some(name) = self.get("datum") {
            return Some(match name {
                "WGS84" => Datum::WGS84,
                // NAD83 / GRS80 are treated as a zero shift from WGS84 (as in the registry).
                "NAD83" | "GRS80" | "ETRS89" => Datum::GRS80,
                "NAD27" => Datum::nad27(),
                "OSGB36" => Datum::osgb36(),
                _ => return None,
            });
        }
        let ellipsoid = self.ellipsoid()?;
        let to_wgs84 = match self.towgs84() {
            Some(h) => DatumShift::Helmert(h),
            None => DatumShift::Helmert(Helmert::IDENTITY),
        };
        Some(Datum { ellipsoid, to_wgs84 })
    }

    /// Resolve the reference ellipsoid from `+ellps`, or an explicit
    /// `+a` with `+b` / `+f` / `+rf`; defaults to WGS84 if nothing is given.
    fn ellipsoid(&self) -> Option<Ellipsoid> {
        if let Some(name) = self.get("ellps") {
            return datum::ellipsoid_by_proj_name(name);
        }
        let Some(a) = self.f("a") else {
            return Some(Ellipsoid::WGS84);
        };
        let f = if let Some(b) = self.f("b") {
            (a - b) / a
        } else if let Some(rf) = self.f("rf") {
            1.0 / rf
        } else {
            self.f("f")?
        };
        Some(Ellipsoid { a, f })
    }

    /// The `+towgs84=` Helmert transform: 3 params (translation only) or 7
    /// (translation, rotation in arc-seconds, scale in ppm). `None` if absent.
    fn towgs84(&self) -> Option<Helmert> {
        let v: Vec<f64> = self
            .get("towgs84")?
            .split(',')
            .filter_map(|x| x.parse().ok())
            .collect();
        match v.as_slice() {
            [tx, ty, tz] => Some(Helmert::bursa_wolf(*tx, *ty, *tz, 0.0, 0.0, 0.0, 0.0)),
            [tx, ty, tz, rx, ry, rz, s] => {
                Some(Helmert::bursa_wolf(*tx, *ty, *tz, *rx, *ry, *rz, *s))
            }
            _ => None,
        }
    }
}
