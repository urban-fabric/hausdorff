//! Codec round-trips: WKB, WKT, GeoJSON, shapefile.
//!
//! The oracle is **the round trip itself** — encode, decode, compare — which needs
//! no fixtures and covers every geometry variant the model has. Where an external
//! byte format is involved the suites also pin known-good encodings, because a
//! round trip alone cannot catch two inverse bugs that cancel.
//!
//! Hostile and malformed input is *not* here: it is in
//! [`robustness`](super::robustness), because what it asserts is a different kind
//! of claim — that the process survives — and because a failure there can take the
//! whole test binary down rather than report itself.

mod geojson;
mod geojson_features;
mod shapefile;
mod wkb;
mod wkt;
