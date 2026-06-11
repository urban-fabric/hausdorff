//! [`Ellipsoid`] — a reference ellipsoid of revolution (or sphere): the
//! datum/planet model shared by the map projections
//! ([`crate::earth::projection`]) and the ellipsoidal (Karney) geodesics
//! ([`Ellipsoid::inverse`] / [`Ellipsoid::direct`]). The spherical counterpart
//! on the metric side is [`crate::earth::sphere::Sphere`].

use euclid::LengthMeters;

use crate::earth::karney::{self, Karney, LineKernel};
use crate::earth::{Bearing, GeoPoint};

/// A reference ellipsoid (or sphere) — the datum/planet a projection or
/// geodesic computation is built on, defined by its semi-major axis and
/// flattening.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Ellipsoid {
    /// Semi-major (equatorial) axis, metres.
    pub a: f64,
    /// Flattening (`0` for a perfect sphere).
    pub f: f64,
}

impl Ellipsoid {
    /// WGS84 (the GPS datum).
    pub const WGS84: Ellipsoid = Ellipsoid {
        a: 6_378_137.0,
        f: 1.0 / 298.257_223_563,
    };
    /// GRS80 (almost identical to WGS84).
    pub const GRS80: Ellipsoid = Ellipsoid {
        a: 6_378_137.0,
        f: 1.0 / 298.257_222_101,
    };
    /// A perfect sphere of the WGS84 mean radius.
    pub const SPHERE: Ellipsoid = Ellipsoid {
        a: 6_371_008.8,
        f: 0.0,
    };
    /// Mars (IAU 2009).
    pub const MARS: Ellipsoid = Ellipsoid {
        a: 3_396_190.0,
        f: 1.0 / 169.894_447_223_612,
    };

    /// Semi-minor (polar) axis `b = a(1 − f)`, metres.
    #[must_use]
    pub fn b(self) -> f64 {
        self.a * (1.0 - self.f)
    }

    /// First eccentricity squared `e² = f(2 − f)`.
    #[must_use]
    pub fn e2(self) -> f64 {
        self.f * (2.0 - self.f)
    }

    /// Second eccentricity squared `e′² = e² / (1 − e²)`.
    #[must_use]
    pub fn ep2(self) -> f64 {
        self.e2() / (1.0 - self.e2())
    }

    /// Third flattening `n = f / (2 − f)`.
    #[must_use]
    pub fn n(self) -> f64 {
        self.f / (2.0 - self.f)
    }

    /// Run `f` with this ellipsoid's geodesic solver — the cached one for
    /// WGS84 (the common case), a fresh one otherwise.
    fn with_solver<R>(self, f: impl FnOnce(&Karney) -> R) -> R {
        if self == Self::WGS84 {
            f(karney::wgs84())
        } else {
            f(&Karney::new(self))
        }
    }

    /// The inverse geodesic problem (Karney): shortest-path distance and the
    /// forward bearings at both endpoints. Accurate to round-off (~15 nm on
    /// Earth) for any pair, including near-antipodal points. Results are NaN
    /// for a non-finite or out-of-range latitude.
    #[must_use]
    pub fn inverse(self, from: GeoPoint, to: GeoPoint) -> GeodesicInverse {
        self.with_solver(|k| {
            let (distance_m, azi1, azi2) = k.inverse(from.lat, from.lon, to.lat, to.lon);
            GeodesicInverse {
                distance_m,
                initial_bearing: Bearing::from_degrees(azi1),
                final_bearing: Bearing::from_degrees(azi2),
            }
        })
    }

    /// Geodesic distance between two points on this ellipsoid, in metres —
    /// the distance component of [`inverse`](Self::inverse). The ellipsoidal
    /// counterpart of [`crate::earth::haversine_m`] (which can differ from it
    /// by up to ~0.5% on Earth).
    #[must_use]
    pub fn distance_m(self, from: GeoPoint, to: GeoPoint) -> f64 {
        self.with_solver(|k| k.inverse(from.lat, from.lon, to.lat, to.lon).0)
    }

    /// [`distance_m`](Self::distance_m) as a typed length.
    #[must_use]
    pub fn distance(self, from: GeoPoint, to: GeoPoint) -> LengthMeters {
        LengthMeters::new(self.distance_m(from, to))
    }

    /// The direct geodesic problem (Karney): travel `distance_m` metres from
    /// `origin` along the initial compass `bearing`; returns the destination
    /// and the forward bearing there. Inverse-pair to
    /// [`inverse`](Self::inverse).
    #[must_use]
    pub fn direct(self, origin: GeoPoint, bearing: Bearing, distance_m: f64) -> GeodesicDirect {
        self.with_solver(|k| {
            let (lat2, lon2, azi2) =
                k.direct(origin.lat, origin.lon, bearing.degrees(), distance_m);
            GeodesicDirect {
                point: GeoPoint::from_lat_lon(lat2, lon2),
                final_bearing: Bearing::from_degrees(azi2),
            }
        })
    }

    /// The destination point of [`direct`](Self::direct) — the ellipsoidal
    /// counterpart of [`crate::earth::destination_m`].
    #[must_use]
    pub fn destination_m(self, origin: GeoPoint, bearing: Bearing, distance_m: f64) -> GeoPoint {
        self.direct(origin, bearing, distance_m).point
    }

    /// [`destination_m`](Self::destination_m) with a typed length.
    #[must_use]
    pub fn destination(
        self,
        origin: GeoPoint,
        bearing: Bearing,
        distance: LengthMeters,
    ) -> GeoPoint {
        self.destination_m(origin, bearing, distance.as_meters())
    }

    /// The geodesic line from `origin` with initial compass `bearing` — the
    /// origin-dependent setup of the direct problem done once, so each
    /// [`GeodesicLine::position_m`] along it costs roughly a third of a full
    /// [`direct`](Self::direct). Use this to place many waypoints on one
    /// geodesic (sampling, densification, interpolation).
    #[must_use]
    pub fn line(self, origin: GeoPoint, bearing: Bearing) -> GeodesicLine {
        let kernel = self.with_solver(|k| k.line(origin.lat, origin.lon, bearing.degrees()));
        GeodesicLine { kernel }
    }

    /// Geodesic perimeter (metres) and **signed** area (m²,
    /// counter-clockwise positive) of the polygon ring through `ring`'s
    /// vertices — Karney's ellipsoidal polygon area, exact to round-off where
    /// the spherical
    /// `ChamberlainDuquetteArea` (in `hausdorff-planar`'s `space::r2`)
    /// is off by up to ~0.5%. The ring closes itself; an already-closed ring
    /// is fine. A clockwise ring gives a negative area; for rings that may
    /// enclose more than half the planet use
    /// [`ring_perimeter_area_unsigned_m`](Self::ring_perimeter_area_unsigned_m).
    #[must_use]
    pub fn ring_perimeter_area_m(self, ring: &[GeoPoint]) -> (f64, f64) {
        self.ring_perimeter_area_inner(ring, true)
    }

    /// Geodesic perimeter (metres) and **unsigned** enclosed area (m²) of the
    /// ring through `ring`'s vertices, read as counter-clockwise — correct
    /// even for rings enclosing more than half the planet (where the signed
    /// form wraps negative).
    #[must_use]
    pub fn ring_perimeter_area_unsigned_m(self, ring: &[GeoPoint]) -> (f64, f64) {
        self.ring_perimeter_area_inner(ring, false)
    }

    fn ring_perimeter_area_inner(self, ring: &[GeoPoint], signed: bool) -> (f64, f64) {
        self.ring_perimeter_area_winding_m(ring, true, signed)
    }

    /// The winding-flexible form of [`ring_perimeter_area_m`](Self::ring_perimeter_area_m):
    /// `ccw_positive` picks which winding counts as positive (holes are
    /// conventionally clockwise). Workspace-internal (hence `#[doc(hidden)]`) —
    /// `hausdorff-planar`'s `space::r2` `GeodesicArea` trait builds polygon
    /// semantics on it.
    #[doc(hidden)]
    pub fn ring_perimeter_area_winding_m(
        self,
        ring: &[GeoPoint],
        ccw_positive: bool,
        signed: bool,
    ) -> (f64, f64) {
        let latlon: Vec<(f64, f64)> = ring.iter().map(|p| (p.lat, p.lon)).collect();
        self.with_solver(|k| k.ring_perimeter_area(&latlon, ccw_positive, signed))
    }

    /// The ellipsoid's total surface area, m² (`4π` times the authalic
    /// radius squared).
    #[must_use]
    pub fn total_area_m2(self) -> f64 {
        self.with_solver(Karney::total_area)
    }
}

/// Solution of the inverse geodesic problem ([`Ellipsoid::inverse`]): the
/// shortest-path distance and the forward bearing at each endpoint.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GeodesicInverse {
    /// Shortest-path (geodesic) distance, metres.
    pub distance_m: f64,
    /// Forward bearing at the start point.
    pub initial_bearing: Bearing,
    /// Forward bearing at the end point (the direction of travel on arrival,
    /// not the back-azimuth).
    pub final_bearing: Bearing,
}

/// Solution of the direct geodesic problem ([`Ellipsoid::direct`]): the
/// destination and the forward bearing on arrival.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GeodesicDirect {
    /// The destination point.
    pub point: GeoPoint,
    /// Forward bearing at the destination.
    pub final_bearing: Bearing,
}

/// A geodesic line ([`Ellipsoid::line`]): an origin and initial bearing with
/// the direct-problem series precomputed, so points can be placed at any
/// distance along it cheaply. Distances may be negative (backwards along the
/// line).
pub struct GeodesicLine {
    kernel: LineKernel,
}

impl GeodesicLine {
    /// The point `distance_m` metres along the line.
    #[must_use]
    pub fn point_m(&self, distance_m: f64) -> GeoPoint {
        let (lat, lon, _) = self.kernel.position(distance_m);
        GeoPoint::from_lat_lon(lat, lon)
    }

    /// The point `distance_m` metres along the line, with the forward bearing
    /// there.
    #[must_use]
    pub fn position_m(&self, distance_m: f64) -> GeodesicDirect {
        let (lat, lon, azi2) = self.kernel.position(distance_m);
        GeodesicDirect {
            point: GeoPoint::from_lat_lon(lat, lon),
            final_bearing: Bearing::from_degrees(azi2),
        }
    }

    /// [`position_m`](Self::position_m) with a typed length.
    #[must_use]
    pub fn position(&self, distance: LengthMeters) -> GeodesicDirect {
        self.position_m(distance.as_meters())
    }
}
