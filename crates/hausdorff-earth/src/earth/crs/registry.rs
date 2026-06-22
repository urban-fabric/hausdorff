use super::*;

impl Crs {
    /// The CRS for a standard EPSG `code`, or `None` if unrecognised.
    ///
    /// Recognised:
    /// - **Geographic 2D**: `4326` (WGS84), `4258` (ETRS89), `4269` (NAD83); the
    ///   Helmert-shifted `4277` (OSGB36), `4230` (ED50), `4267` (NAD27); the
    ///   **NTv2 grid-shifted** `4314` (DHDN), `4275` (NTF), `4272` (NZGD49).
    /// - **Geographic 3D** (lat/lon + ellipsoidal height): `4979` (WGS84),
    ///   `4937` (ETRS89).
    /// - **Geocentric**: `4978` (WGS84 ECEF).
    /// - **Web Mercator**: `3857`/`900913`.
    /// - **Transverse Mercator / UTM**: the WGS84 UTM zones `32601..=32660`
    ///   (north) / `32701..=32760` (south); the ETRS89 zones `25828..=25838`,
    ///   NAD83 zones `26901..=26923`, GDA94/MGA zones `28348..=28358`; `2193`
    ///   (NZGD2000 / NZTM2000), `2180` (ETRS89 / Poland CS92), `27700` (British
    ///   National Grid).
    /// - **Lambert Conformal Conic**: `2154` (RGF93 Lambert-93), `3978` (NAD83 /
    ///   Canada Atlas Lambert).
    /// - **Albers Equal-Area**: `5070` (NAD83 Conus Albers), `3577` (GDA94 /
    ///   Australian Albers).
    /// - **Polar Stereographic**: `3995` (Arctic) / `3031` (Antarctic); `3413` /
    ///   `3976` (WGS84 NSIDC Sea-Ice North / South).
    ///
    /// (NAD83/ETRS89/GDA94/RGF93 share the GRS80 ellipsoid and are treated as a
    /// zero shift from WGS84.)
    #[must_use]
    pub fn from_epsg(code: u32) -> Option<Crs> {
        match code {
            4326 => Some(Crs::wgs84()),
            4978 => Some(Crs::Geocentric {
                datum: Datum::WGS84,
                epsg: Some(4978),
            }),
            3857 | 900_913 => Some(Crs::Projected {
                datum: Datum::WGS84,
                projection: ProjectionKind::WebMercator(WebMercator::EPSG_3857),
                epsg: Some(3857),
            }),
            32_601..=32_660 => Some(Self::utm((code - 32_600) as u8, true, code)),
            32_701..=32_760 => Some(Self::utm((code - 32_700) as u8, false, code)),
            5070 => Some(Crs::Projected {
                datum: Datum::GRS80,
                projection: ProjectionKind::AlbersEqualArea(AlbersEqualArea {
                    ellipsoid: Ellipsoid::GRS80,
                    lat1: 29.5,
                    lat2: 45.5,
                    lat0: 23.0,
                    lon0: -96.0,
                    false_easting: 0.0,
                    false_northing: 0.0,
                }),
                epsg: Some(5070),
            }),
            2154 => Some(Crs::Projected {
                datum: Datum::GRS80,
                projection: ProjectionKind::LambertConformalConic(LambertConformalConic {
                    ellipsoid: Ellipsoid::GRS80,
                    lat1: 44.0,
                    lat2: 49.0,
                    lat0: 46.5,
                    lon0: 3.0,
                    false_easting: 700_000.0,
                    false_northing: 6_600_000.0,
                }),
                epsg: Some(2154),
            }),
            3995 => Some(Self::polar_stereographic(71.0, 0.0, true, 3995)),
            3031 => Some(Self::polar_stereographic(-71.0, 0.0, false, 3031)),
            // Geographic 3-D — lat/lon + ellipsoidal height (height retained).
            4979 => Some(Crs::wgs84_3d()), // WGS84 3D
            4937 => Some(Crs::Geographic {
                datum: Datum::GRS80, // ETRS89 3D
                epsg: Some(4937),
                dim: GeoDim::D3,
            }),
            // Datum-shifted geographic CRSes (non-WGS84 datums).
            4277 => Some(Crs::Geographic {
                datum: Datum::osgb36(),
                epsg: Some(4277),
                dim: GeoDim::D2,
            }),
            4230 => Some(Crs::Geographic {
                datum: Datum::ed50(),
                epsg: Some(4230),
                dim: GeoDim::D2,
            }),
            4267 => Some(Crs::Geographic {
                datum: Datum::nad27(),
                epsg: Some(4267),
                dim: GeoDim::D2,
            }),
            // NTv2 grid-shifted geographic datums (~cm accuracy, embedded grids).
            4314 => Some(Crs::Geographic {
                datum: Datum::dhdn(), // German DHDN — BETA2007
                epsg: Some(4314),
                dim: GeoDim::D2,
            }),
            4275 => Some(Crs::Geographic {
                datum: Datum::ntf(), // French NTF — ntf_r93
                epsg: Some(4275),
                dim: GeoDim::D2,
            }),
            4272 => Some(Crs::Geographic {
                datum: Datum::nzgd49(), // New Zealand NZGD49 — nzgd2kgrid0005
                epsg: Some(4272),
                dim: GeoDim::D2,
            }),
            // OSGB36 / British National Grid (Transverse Mercator on Airy 1830).
            27_700 => Some(Crs::Projected {
                datum: Datum::osgb36(),
                projection: ProjectionKind::TransverseMercator(TransverseMercator {
                    ellipsoid: Datum::osgb36().ellipsoid, // Airy 1830
                    central_meridian: -2.0,
                    lat_origin: 49.0,
                    k0: 0.999_601_271_7,
                    false_easting: 400_000.0,
                    false_northing: -100_000.0,
                }),
                epsg: Some(27_700),
            }),
            // GRS80 geographic bases (NAD83 / ETRS89 — zero shift from WGS84).
            4269 => Some(Crs::Geographic {
                datum: Datum::GRS80,
                epsg: Some(4269),
                dim: GeoDim::D2,
            }),
            4258 => Some(Crs::Geographic {
                datum: Datum::GRS80,
                epsg: Some(4258),
                dim: GeoDim::D2,
            }),
            // ETRS89 / UTM zones 28N..=38N (GRS80) — Europe.
            25_828..=25_838 => Some(Self::utm_on(
                Datum::GRS80,
                (code - 25_800) as u8,
                true,
                code,
            )),
            // NAD83 / UTM zones 1N..=23N (GRS80) — North America.
            26_901..=26_923 => Some(Self::utm_on(
                Datum::GRS80,
                (code - 26_900) as u8,
                true,
                code,
            )),
            // GDA94 / MGA zones 48..=58 (south, GRS80) — Australia.
            28_348..=28_358 => Some(Self::utm_on(
                Datum::GRS80,
                (code - 28_300) as u8,
                false,
                code,
            )),
            // NZGD2000 / New Zealand Transverse Mercator 2000.
            2193 => Some(Crs::Projected {
                datum: Datum::GRS80,
                projection: ProjectionKind::TransverseMercator(TransverseMercator {
                    ellipsoid: Ellipsoid::GRS80,
                    central_meridian: 173.0,
                    lat_origin: 0.0,
                    k0: 0.9996,
                    false_easting: 1_600_000.0,
                    false_northing: 10_000_000.0,
                }),
                epsg: Some(2193),
            }),
            // ETRS89 / Poland CS92 (note k0 = 0.9993, a single-zone TM).
            2180 => Some(Crs::Projected {
                datum: Datum::GRS80,
                projection: ProjectionKind::TransverseMercator(TransverseMercator {
                    ellipsoid: Ellipsoid::GRS80,
                    central_meridian: 19.0,
                    lat_origin: 0.0,
                    k0: 0.9993,
                    false_easting: 500_000.0,
                    false_northing: -5_300_000.0,
                }),
                epsg: Some(2180),
            }),
            // NAD83 / Canada Atlas Lambert.
            3978 => Some(Crs::Projected {
                datum: Datum::GRS80,
                projection: ProjectionKind::LambertConformalConic(LambertConformalConic {
                    ellipsoid: Ellipsoid::GRS80,
                    lat1: 49.0,
                    lat2: 77.0,
                    lat0: 49.0,
                    lon0: -95.0,
                    false_easting: 0.0,
                    false_northing: 0.0,
                }),
                epsg: Some(3978),
            }),
            // GDA94 / Australian Albers.
            3577 => Some(Crs::Projected {
                datum: Datum::GRS80,
                projection: ProjectionKind::AlbersEqualArea(AlbersEqualArea {
                    ellipsoid: Ellipsoid::GRS80,
                    lat1: -18.0,
                    lat2: -36.0,
                    lat0: 0.0,
                    lon0: 132.0,
                    false_easting: 0.0,
                    false_northing: 0.0,
                }),
                epsg: Some(3577),
            }),
            // WGS84 / NSIDC Sea Ice Polar Stereographic (North off −45°, South off 0°).
            3413 => Some(Self::polar_stereographic(70.0, -45.0, true, 3413)),
            3976 => Some(Self::polar_stereographic(-70.0, 0.0, false, 3976)),
            _ => None,
        }
    }

    fn utm(zone: u8, north: bool, epsg: u32) -> Crs {
        Self::utm_on(Datum::WGS84, zone, north, epsg)
    }

    /// A UTM `zone` projected CRS on an arbitrary `datum` — the standard UTM
    /// Transverse Mercator parameters (6° zones, k0 = 0.9996), but on the datum's
    /// own ellipsoid (e.g. GRS80 for the ETRS89 / NAD83 / GDA94 zones).
    fn utm_on(datum: Datum, zone: u8, north: bool, epsg: u32) -> Crs {
        let mut tm = TransverseMercator::utm(zone, north);
        tm.ellipsoid = datum.ellipsoid;
        Crs::Projected {
            datum,
            projection: ProjectionKind::TransverseMercator(tm),
            epsg: Some(epsg),
        }
    }

    fn polar_stereographic(lat_ts: f64, lon0: f64, north: bool, epsg: u32) -> Crs {
        Crs::Projected {
            datum: Datum::WGS84,
            projection: ProjectionKind::PolarStereographic(PolarStereographic {
                ellipsoid: Ellipsoid::WGS84,
                lat_ts,
                lon0,
                north,
                false_easting: 0.0,
                false_northing: 0.0,
            }),
            epsg: Some(epsg),
        }
    }

    /// The EPSG code, if this CRS carries one.
    #[must_use]
    pub fn epsg(self) -> Option<u32> {
        match self {
            Crs::Geographic { epsg, .. }
            | Crs::Geocentric { epsg, .. }
            | Crs::Projected { epsg, .. } => epsg,
        }
    }

    /// The datum this CRS is referenced to.
    #[must_use]
    pub fn datum(self) -> Datum {
        match self {
            Crs::Geographic { datum, .. }
            | Crs::Geocentric { datum, .. }
            | Crs::Projected { datum, .. } => datum,
        }
    }

    /// The CRS's **area of use** as a geographic (WGS84 lon/lat degrees) bounding
    /// box — the extent over which it is intended to be used — or `None` for a
    /// CRS with no standard code (e.g. a [`local_tangent_plane`](Self::local_tangent_plane)).
    ///
    /// For the UTM-family CRSes this is the nominal 6° zone band (the projection's
    /// valid extent, the right domain for a distortion check); for the others it
    /// is the EPSG-published extent. A box spanning the antimeridian is returned
    /// with `min_lon > max_lon` (a wrapping box) — use [`contains`](Self::contains),
    /// which handles the wrap, rather than [`BBox::contains`] for membership.
    #[must_use]
    pub fn area_of_use(self) -> Option<BBox> {
        let code = self.epsg()?;
        Some(match code {
            // Global frames.
            4326 | 4978 | 4979 => BBox::new(-90.0, -180.0, 90.0, 180.0),
            3857 | 900_913 => BBox::new(-85.06, -180.0, 85.06, 180.0),
            // UTM-family: the 6° zone band (north 0..84, south −80..0).
            32_601..=32_660 => utm_band((code - 32_600) as u8, true),
            32_701..=32_760 => utm_band((code - 32_700) as u8, false),
            25_828..=25_838 => utm_band((code - 25_800) as u8, true),
            26_901..=26_923 => utm_band((code - 26_900) as u8, true),
            28_348..=28_358 => utm_band((code - 28_300) as u8, false),
            // EPSG-published extents (S, W, N, E from the EPSG registry).
            5070 => BBox::new(24.41, -124.79, 49.38, -66.91),
            2154 => BBox::new(41.15, -9.86, 51.56, 10.38),
            3978 => BBox::new(38.21, -141.01, 86.46, -40.73),
            3577 => BBox::new(-43.7, 112.85, -9.86, 153.69),
            2193 => BBox::new(-47.33, 166.37, -34.1, 178.63),
            2180 => BBox::new(49.0, 14.14, 55.93, 24.15),
            3995 | 3413 => BBox::new(60.0, -180.0, 90.0, 180.0),
            3031 | 3976 => BBox::new(-90.0, -180.0, -60.0, 180.0),
            4258 | 4937 => BBox::new(32.88, -16.1, 84.73, 40.18),
            4230 => BBox::new(25.71, -16.1, 84.73, 48.61),
            4314 => BBox::new(47.27, 5.86, 55.09, 13.84),
            4275 => BBox::new(41.31, -4.87, 51.14, 9.63),
            4272 => BBox::new(-47.65, 165.87, -33.89, 179.27),
            4277 | 27_700 => BBox::new(49.75, -9.01, 61.01, 2.01),
            // Antimeridian-wrapping extents (min_lon > max_lon).
            4267 => BBox::new(7.15, 167.65, 83.17, -47.74),
            4269 => BBox::new(14.92, 167.65, 86.45, -40.73),
            _ => return None,
        })
    }

    /// True iff the WGS84 geographic point `p` lies within this CRS's
    /// [`area_of_use`](Self::area_of_use) (inclusive), handling an
    /// antimeridian-wrapping extent correctly. A CRS with no declared area of use
    /// contains every finite point; non-finite points are never contained.
    ///
    /// The check is a coarse domain guard — `p` is taken as WGS84 lon/lat and the
    /// datum-shift difference (≤ ~100 m) is negligible against the degree-scale
    /// bounds.
    #[must_use]
    pub fn contains(self, p: GeoPoint) -> bool {
        if !p.lat.is_finite() || !p.lon.is_finite() {
            return false;
        }
        let Some(b) = self.area_of_use() else {
            return true;
        };
        let lat_ok = p.lat >= b.min_lat && p.lat <= b.max_lat;
        let lon_ok = if b.min_lon <= b.max_lon {
            p.lon >= b.min_lon && p.lon <= b.max_lon
        } else {
            // Box crosses the antimeridian: valid east of min_lon OR west of max_lon.
            p.lon >= b.min_lon || p.lon <= b.max_lon
        };
        lat_ok && lon_ok
    }
}

/// The area-of-use bounding box for UTM `zone` (1..=60): the 6° longitude band
/// `[-180 + 6(zone−1), -180 + 6·zone]`, latitude `0..84` (north) or `−80..0`
/// (south).
fn utm_band(zone: u8, north: bool) -> BBox {
    let west = -180.0 + 6.0 * (f64::from(zone) - 1.0);
    let (south, north_lat) = if north { (0.0, 84.0) } else { (-80.0, 0.0) };
    BBox::new(south, west, north_lat, west + 6.0)
}
