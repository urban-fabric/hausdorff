//! `r2_format` — serialization of planar [`crate::r2::Geometry`] to and from the
//! standard interchange formats (the "r2-format" layer).
//!
//! WKT/WKB/GeoJSON are 2-D Cartesian geometry formats: their coordinates map
//! directly to [`crate::r2::Coord`] regardless of CRS (a lon/lat GeoJSON is an
//! `r2::Geometry` you then interpret geographically via [`crate::convert`] /
//! [`crate::projection`]). Native, self-contained codecs that let the workspace
//! deprecate the external `wkt` / `wkb` / `geojson` crates.

pub mod geojson;
pub mod wkb;
pub mod wkt;
