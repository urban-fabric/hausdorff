//! Planar (Cartesian) geometry — the in-house replacement for `geo`/`geo-types`.
//!
//! `f64`-scalar geometry types ([`Coord`], [`Point`], [`Line`], [`LineString`],
//! [`Polygon`], [`MultiPolygon`], …) and the algorithms over them, mirroring
//! `geo`'s shapes and trait names so consumers migrate by swapping the import
//! path. This is the planar sibling of the spherical WGS84 types in
//! [`crate::geo`]; move between the two with [`crate::convert`].
//!
//! Ported from georust `geo` / `geo-types` (dual MIT/Apache-2.0); see the crate
//! `NOTICE`.

#[macro_use]
mod macros;

mod coord;
mod geometry;
mod line_string;
mod polygon;
mod primitives;

pub use coord::{Coord, Point};
pub use geometry::Geometry;
pub use line_string::{LineString, MultiLineString};
pub use polygon::{MultiPoint, MultiPolygon, Polygon};
pub use primitives::{Line, Rect, Triangle};
