//! Small fixed-vertex planar geometries: [`Line`] (a segment), [`Triangle`],
//! and the axis-aligned [`Rect`]. Each can expand to a [`Polygon`] for the
//! algorithms that only speak in rings.

use serde::{Deserialize, Serialize};

use super::{Coord, LineString, Point, Polygon};

mod fixed;
pub use fixed::{Line, Rect, Triangle};

mod ellipse;
pub use ellipse::{ConfidenceEllipse, Ellipse};

#[cfg(test)]
mod tests;
