//! Public-API tests for `hausdorff::space::r2`, relocated out of `src/`.

#[path = "r2/coord.rs"]
mod coord;

#[path = "r2/primitives.rs"]
mod primitives;

#[path = "r2/line_string.rs"]
mod line_string;

#[path = "r2/polygon.rs"]
mod polygon;

#[path = "r2/geometry.rs"]
mod geometry;

#[path = "r2/winding.rs"]
mod winding;

#[path = "r2/area.rs"]
mod area;

#[path = "r2/bounding_rect.rs"]
mod bounding_rect;

#[path = "r2/dimensions.rs"]
mod dimensions;

#[path = "r2/centroid.rs"]
mod centroid;

#[path = "r2/contains.rs"]
mod contains;

#[path = "r2/convex_hull.rs"]
mod convex_hull;

#[path = "r2/macros.rs"]
mod macros;

#[path = "r2/geometry_dispatch.rs"]
mod geometry_dispatch;
