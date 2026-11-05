//! Map-projection kernels over any [`Real`] scalar — geographic
//! (lat/lon degrees) ↔ projected planar metres, written once so the same
//! source runs on `f64`, `f32`, and tensor/expression backends.
//!
//! These are the backend-neutral peers of the concrete `Projection` impls in
//! `hausdorff-earth`: that crate keeps the `Ellipsoid`/`GeoPoint`/`Coord`-typed
//! API and delegates the arithmetic here, exactly as the earth geodesic
//! wrappers delegate to [`crate::earth`]. The per-point coordinate is the lazy
//! scalar `T`; an ellipsoid's constants (the eccentricity, the Krüger series,
//! the false origin) are `f64` — they depend only on the datum, not on the
//! point being transformed. The f64 convenience entry points lift those
//! constants internally; typed entry points accept backend-typed constants so
//! expression/tensor callers can lift once and reuse them across many points.
//!
//! Provided: [`WebMercator`](web_mercator_forward) (spherical Mercator) and the
//! Karney–Krüger [`TransverseMercator`](transverse_mercator_forward) (the basis
//! of UTM). Following Karney 2011, *Transverse Mercator with an accuracy of a
//! few nanometers*, and Snyder, *Map Projections — A Working Manual* (USGS
//! PP1395, 1987).

use euclid_numeric::Real;

use crate::earth::LatLon;

mod web_mercator;
pub use web_mercator::{
    WebMercatorConstants, web_mercator_forward, web_mercator_forward_typed, web_mercator_inverse,
    web_mercator_inverse_typed, web_mercator_project_at, web_mercator_project_at_typed,
};

mod transverse_mercator;
pub use transverse_mercator::{
    KrugerSeries, KrugerSeriesT, transverse_mercator_forward, transverse_mercator_forward_typed,
    transverse_mercator_inverse, transverse_mercator_inverse_typed, transverse_mercator_project_at,
    transverse_mercator_project_at_typed,
};

mod auxiliary;
pub use auxiliary::{
    authalic_q, central_lon_delta, conformal_t, lat_from_authalic_q, lat_from_conformal_t,
    m_factor, tauf, taupf,
};

mod families;
pub use families::{
    AlbersConstants, AlbersConstantsT, LccConstants, LccConstantsT, PolarStereoConstants,
    PolarStereoConstantsT, albers_forward, albers_forward_typed, albers_inverse,
    albers_inverse_typed, lcc_forward, lcc_forward_typed, lcc_inverse, lcc_inverse_typed,
    polar_stereographic_forward, polar_stereographic_forward_typed, polar_stereographic_inverse,
    polar_stereographic_inverse_typed,
};

#[cfg(test)]
mod tests;
