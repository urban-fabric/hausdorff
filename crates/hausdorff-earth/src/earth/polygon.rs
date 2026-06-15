//! `polygon` — WGS84 polygons measured **natively on the curved Earth**, with no
//! map projection: a [`GeoPolygon`] (exterior ring + holes) and a
//! [`GeoMultiPolygon`], the "OSM without a projection" geometry types.
//!
//! Each measure is computed in the model that makes it *exact*:
//! - **area** — through the **authalic (equal-area) sphere**
//!   ([`AuthalicSphere`](crate::earth::sphere::AuthalicSphere)): the φ→ξ warp is
//!   area-preserving, so the proven spherical-excess kernel (`riemann`
//!   Hausdorff/SphereArea) yields the true datum area;
//! - **perimeter** — the **geodesic** boundary length (Karney, exact to round-off);
//! - **containment** — the spherical winding-number point-in-polygon
//!   ([`s2::Loop::contains_point`](hausdorff_sphere::space::s2::Loop::contains_point),
//!   great-circle edges on the unit sphere), robust away from the boundary.
//!
//! Rings may be given open or closed, in either winding order — areas use
//! unsigned magnitudes and containment normalizes each ring to counter-clockwise.

use crate::earth::sphere::AuthalicSphere;
use crate::earth::{EcefPoint, Ellipsoid, GeoPoint};
use hausdorff_sphere::space::s2::edges::{Crossing, crossing_sign};
use hausdorff_sphere::space::s2::{Loop, Polygon as S2Polygon, S2Point, polygon_signed_area_vos};

/// A polygon on the WGS84 ellipsoid (or any [`Ellipsoid`]): an **exterior** ring
/// of lat/lon vertices with zero or more **interior** (hole) rings. Measured
/// natively on the curved Earth — see the [module docs](self).
#[derive(Clone, Debug, PartialEq)]
pub struct GeoPolygon {
    ellipsoid: Ellipsoid,
    exterior: Vec<GeoPoint>,
    interiors: Vec<Vec<GeoPoint>>,
}

impl GeoPolygon {
    /// A **WGS84** polygon with the given `exterior` ring and `interiors` (holes).
    #[must_use]
    pub fn new(exterior: Vec<GeoPoint>, interiors: Vec<Vec<GeoPoint>>) -> Self {
        Self::on(Ellipsoid::WGS84, exterior, interiors)
    }

    /// A WGS84 polygon with no holes.
    #[must_use]
    pub fn from_exterior(exterior: Vec<GeoPoint>) -> Self {
        Self::new(exterior, Vec::new())
    }

    /// A polygon on an arbitrary `ellipsoid` (e.g. [`Ellipsoid::MARS`]).
    #[must_use]
    pub fn on(
        ellipsoid: Ellipsoid,
        exterior: Vec<GeoPoint>,
        interiors: Vec<Vec<GeoPoint>>,
    ) -> Self {
        Self {
            ellipsoid,
            exterior,
            interiors,
        }
    }

    /// The exterior ring vertices.
    #[must_use]
    pub fn exterior(&self) -> &[GeoPoint] {
        &self.exterior
    }

    /// The interior (hole) rings.
    #[must_use]
    pub fn interiors(&self) -> &[Vec<GeoPoint>] {
        &self.interiors
    }

    /// The reference ellipsoid this polygon is measured on.
    #[must_use]
    pub fn ellipsoid(&self) -> Ellipsoid {
        self.ellipsoid
    }

    /// **Exact ellipsoidal area (m²)** — the exterior area minus the hole areas,
    /// each computed through the authalic sphere (the proven spherical-excess
    /// path). Clamped at `0`.
    #[must_use]
    pub fn area_m2(&self) -> f64 {
        let au = AuthalicSphere::of(self.ellipsoid);
        let holes: f64 = self.interiors.iter().map(|h| au.polygon_area_m2(h)).sum();
        (au.polygon_area_m2(&self.exterior) - holes).max(0.0)
    }

    /// **Geodesic boundary length (m)** — the exterior perimeter plus every hole
    /// perimeter, each summed as Karney geodesic edge lengths.
    #[must_use]
    pub fn perimeter_m(&self) -> f64 {
        ring_geodesic_len(self.ellipsoid, &self.exterior)
            + self
                .interiors
                .iter()
                .map(|h| ring_geodesic_len(self.ellipsoid, h))
                .sum::<f64>()
    }

    /// True iff `p` lies inside the exterior and outside every hole — the
    /// spherical winding-number test (great-circle edges; the ellipsoid's
    /// flattening does not affect which side of an edge a point is on, away from
    /// the boundary).
    #[must_use]
    pub fn contains(&self, p: GeoPoint) -> bool {
        self.as_s2().contains_point(S2Point::from_geo_point(p))
    }

    /// True iff this polygon and `other` overlap — their exterior boundaries cross
    /// or one exterior contains the other (holes are not subtracted from the
    /// overlap test). The great-circle edge-crossing test is proof-backed by
    /// `riemann` Hausdorff/SphereSegmentCross (`crossing_sign` soundness).
    #[must_use]
    pub fn intersects(&self, other: &GeoPolygon) -> bool {
        ccw_loop(&self.exterior).intersects(&ccw_loop(&other.exterior))
    }

    /// True iff `other` lies **entirely within** this polygon: every vertex of
    /// `other`'s exterior is inside, and no boundary edges cross (so all of
    /// `other` is inside, not just its vertices). With holes it also requires
    /// `other` not to cross any hole boundary. Sound for simple (non-self-crossing)
    /// polygons; the edge-crossing test is proof-backed by `riemann`
    /// Hausdorff/SphereSegmentCross.
    #[must_use]
    pub fn contains_polygon(&self, other: &GeoPolygon) -> bool {
        let target = ccw_loop(&other.exterior);
        other.exterior.iter().all(|&v| self.contains(v))
            && !loops_cross(&ccw_loop(&self.exterior), &target)
            && self
                .interiors
                .iter()
                .all(|h| !loops_cross(&ccw_loop(h), &target))
    }

    /// The spherical [`Polygon`](S2Polygon) of this shape — exterior + holes as
    /// orientation-normalized (counter-clockwise) unit-sphere loops.
    fn as_s2(&self) -> S2Polygon {
        S2Polygon::new(
            ccw_loop(&self.exterior),
            self.interiors.iter().map(|h| ccw_loop(h)).collect(),
        )
    }
}

/// A collection of [`GeoPolygon`]s read as one shape — `Σ` of the parts (a
/// point is inside iff it is inside any member; area / perimeter sum).
#[derive(Clone, Debug, PartialEq)]
pub struct GeoMultiPolygon {
    polygons: Vec<GeoPolygon>,
}

impl GeoMultiPolygon {
    /// A multipolygon of `polygons`.
    #[must_use]
    pub fn new(polygons: Vec<GeoPolygon>) -> Self {
        Self { polygons }
    }

    /// The member polygons.
    #[must_use]
    pub fn polygons(&self) -> &[GeoPolygon] {
        &self.polygons
    }

    /// Total exact ellipsoidal area (m²) — the sum over members.
    #[must_use]
    pub fn area_m2(&self) -> f64 {
        self.polygons.iter().map(GeoPolygon::area_m2).sum()
    }

    /// Total geodesic boundary length (m) — the sum over members.
    #[must_use]
    pub fn perimeter_m(&self) -> f64 {
        self.polygons.iter().map(GeoPolygon::perimeter_m).sum()
    }

    /// True iff `p` lies in any member polygon.
    #[must_use]
    pub fn contains(&self, p: GeoPoint) -> bool {
        self.polygons.iter().any(|poly| poly.contains(p))
    }
}

/// Convert a WGS84 lat/lon `ring` to **unit-sphere** points — the spherical
/// coordinate form of an OSM ring ([`S2Point::from_geo_point`]).
#[must_use]
pub fn geo_ring_to_s2(ring: &[GeoPoint]) -> Vec<S2Point> {
    ring.iter().map(|&p| S2Point::from_geo_point(p)).collect()
}

/// Convert a WGS84 lat/lon `ring` to **ECEF** points on `ellipsoid`'s surface
/// (ellipsoidal height 0) — the elliptical coordinate form of an OSM ring
/// ([`Ellipsoid::geodetic_to_ecef`]).
#[must_use]
pub fn geo_ring_to_ecef(ring: &[GeoPoint], ellipsoid: Ellipsoid) -> Vec<EcefPoint> {
    ring.iter().map(|&p| ellipsoid.geodetic_to_ecef(p, 0.0)).collect()
}

/// A ring as a counter-clockwise unit-sphere [`Loop`] — reversed if its signed
/// spherical area is negative (clockwise), so containment reads the intended
/// (smaller) interior regardless of input winding.
fn ccw_loop(ring: &[GeoPoint]) -> Loop {
    let mut pts = geo_ring_to_s2(ring);
    if polygon_signed_area_vos(&pts) < 0.0 {
        pts.reverse();
    }
    Loop::new(pts)
}

/// True iff any edge of loop `a` properly crosses any edge of loop `b` — an
/// interior great-circle crossing (`s2::edges::crossing_sign == Cross`, the test
/// proven sound in `riemann` Hausdorff/SphereSegmentCross). `O(n·m)`.
fn loops_cross(a: &Loop, b: &Loop) -> bool {
    let (va, vb) = (a.vertices(), b.vertices());
    let (na, nb) = (va.len(), vb.len());
    if na < 2 || nb < 2 {
        return false;
    }
    (0..na).any(|i| {
        let (pa, qa) = (va[i], va[(i + 1) % na]);
        (0..nb).any(|j| crossing_sign(pa, qa, vb[j], vb[(j + 1) % nb]) == Crossing::Cross)
    })
}

/// Closed geodesic length of `ring` (metres) — Karney edge lengths around the
/// ring, including the closing edge. `0` for fewer than three distinct vertices.
fn ring_geodesic_len(ellipsoid: Ellipsoid, ring: &[GeoPoint]) -> f64 {
    // Drop a trailing vertex equal to the first (an already-closed ring).
    let v: &[GeoPoint] = match ring {
        [first, .., last] if first == last => &ring[..ring.len() - 1],
        _ => ring,
    };
    if v.len() < 3 {
        return 0.0;
    }
    (0..v.len())
        .map(|i| ellipsoid.distance_m(v[i], v[(i + 1) % v.len()]))
        .sum()
}

#[cfg(test)]
mod tests;
