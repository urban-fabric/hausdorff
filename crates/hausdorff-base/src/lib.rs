//! Geospatial foundations shared by every hausdorff crate: the WGS84 lat/lon
//! primitives and geodesy kernels ([`earth`]) and the unit-less
//! circle/vector substrates ([`space`] — `s1` angles, `r3` vectors).
//!
//! This is the thin slice most consumers need (`GeoPoint`, `haversine_m`,
//! `BBox`); the planar/spherical geometry, the indexes, and the codecs live in
//! the sibling `hausdorff-*` crates, all re-assembled by the `hausdorff`
//! façade crate under the original module paths.
//!
//! [`algebra`] implements the shared `hausdorff-geometry-traits` abstractions
//! (orientation / vector space / interpolation / region pruning) for these
//! substrate types, so the generic geometry algorithms run on `r1` / `s1` / `r3`
//! the same as on the planar and spherical spaces.

pub mod algebra;
pub mod earth;
pub mod space;
pub mod web;
