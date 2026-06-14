//! Geospatial foundations shared by every hausdorff crate: the WGS84 lat/lon
//! primitives and geodesy kernels ([`earth`]) and the unit-less
//! circle/vector substrates ([`space`] — `s1` angles, `r3` vectors).
//!
//! This is the thin slice most consumers need (`GeoPoint`, `haversine_m`,
//! `BBox`); the planar/spherical geometry, the indexes, and the codecs live in
//! the sibling `hausdorff-*` crates, all re-assembled by the `hausdorff`
//! façade crate under the original module paths.

pub mod earth;
pub mod space;
pub mod web;
