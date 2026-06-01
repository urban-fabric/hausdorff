// `deny` rather than `forbid` so the `bytemuck` POD derives in `geo` can locally
// `#![allow(unsafe_code)]`; the H3 implementation stays unsafe-free under deny.
#![deny(unsafe_code)]

//! Internal geospatial primitives and index implementations for UrbanFabric.

pub mod geo;
pub mod h3;
pub mod s2;

pub use geo::{FixedGeoPoint, GeoPoint};
