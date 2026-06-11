//! `format` — serialization of planar [`hausdorff_planar::space::r2::Geometry`] to and from
//! the standard interchange formats (WKB, WKT, GeoJSON).
//!
//! WKT/WKB/GeoJSON are 2-D Cartesian geometry formats: their coordinates map
//! directly to [`hausdorff_planar::space::r2::Coord`] regardless of CRS (a lon/lat GeoJSON is an
//! `r2::Geometry` you then interpret geographically via `hausdorff-earth`'s `convert` /
//! `earth::projection`). Native, self-contained codecs that let the workspace
//! deprecate the external `wkt` / `wkb` / `geojson` crates.

pub mod geojson;
pub mod shapefile;
pub mod wkb;
pub mod wkt;
