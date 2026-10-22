//! Planar buffering: the set of points within a distance of a geometry.
//!
//! # Why this is a union of parts rather than an offset curve
//!
//! The textbook construction walks the boundary emitting a parallel curve,
//! mitring at convex turns and arcing at reflex ones, then resolves the
//! self-intersections the curve makes wherever the distance exceeds the local
//! feature size. That last step is the whole difficulty: the curve folds back
//! on itself at narrow necks and across short segments, and deciding which
//! folds are holes and which are spurious is where offset-curve implementations
//! accumulate their special cases.
//!
//! A buffer is also exactly a Minkowski sum of the geometry with a disc, and a
//! Minkowski sum distributes over union: the buffer of a shape is the union of
//! the buffers of its parts. So this emits one simple, always-valid polygon per
//! part — a disc per vertex, a rectangle per segment, the polygon itself for
//! areal input — and dissolves them with [`unary_union`]. The sweep that
//! implements the union already resolves overlaps correctly, so the hard case
//! is handled by code that exists and is tested rather than by a second
//! topology routine written for this operation.
//!
//! The cost is more work for the union than a traced curve would need. The
//! gain is that every intermediate is a convex polygon with no self-contact,
//! which is the shape boolean operations are most reliable on.
//!
//! # Joins and caps
//!
//! Round only, and that is a deliberate consequence of the construction above:
//! placing a disc at every vertex *is* a round join and a round cap, so the
//! style falls out of the method rather than being a parameter it has to
//! branch on. Mitred and flat variants need the traced curve, because they are
//! statements about how to leave a vertex rather than about what to place on
//! it. Round is the default in PostGIS, GEOS and PostGIS's callers.
//!
//! The disc is a regular polygon, so a buffer is an approximation whose error
//! is bounded by the sagitta of one of its arcs. [`BufferStyle::quadrant_segments`]
//! sets how many segments approximate a quarter turn; the default of 8 gives a
//! radial error below 0.2% of the distance, which is the same default GEOS
//! uses.
//!
//! # Negative distances
//!
//! A negative distance erodes rather than dilates, and it is only meaningful
//! for areal input: a curve or a point has no interior to retreat into, so the
//! result is empty. Erosion is the region whose distance to the exterior
//! exceeds the magnitude, which is the input minus a dilation of its own
//! boundary — one buffer and one difference, both of which exist here.

use std::f64::consts::TAU;

use super::bool_ops::unary_union;
use super::{BooleanOps as _, Geometry, LineString, MultiPolygon, Polygon};

/// How finely a buffer approximates its arcs.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct BufferStyle {
    quadrant_segments: u32,
}

/// Arcs per quarter turn when the caller does not choose. Matches GEOS, whose
/// choice is the one most published coordinates were produced under.
const DEFAULT_QUADRANT_SEGMENTS: u32 = 8;

/// Below this many segments a quarter turn is not an arc in any useful sense,
/// and a 1-segment quadrant turns a disc into a square — which would silently
/// make a round join a mitred one.
const MIN_QUADRANT_SEGMENTS: u32 = 2;

/// A ceiling on arc refinement. Each vertex contributes `4 * n` coordinates to
/// the union's input, so an unbounded `n` is an unbounded amount of work driven
/// by one argument.
const MAX_QUADRANT_SEGMENTS: u32 = 256;

impl Default for BufferStyle {
    fn default() -> Self {
        Self {
            quadrant_segments: DEFAULT_QUADRANT_SEGMENTS,
        }
    }
}

impl BufferStyle {
    /// A style with an explicit arc refinement, clamped to a sane range.
    ///
    /// Clamped rather than rejected: the argument controls accuracy, not
    /// meaning, so every value has a reasonable interpretation and refusing
    /// one would turn a cosmetic choice into an error path for the caller.
    #[must_use]
    pub fn with_quadrant_segments(segments: u32) -> Self {
        Self {
            quadrant_segments: segments.clamp(MIN_QUADRANT_SEGMENTS, MAX_QUADRANT_SEGMENTS),
        }
    }

    #[must_use]
    pub fn quadrant_segments(self) -> u32 {
        self.quadrant_segments
    }

    /// Coordinates around a full turn.
    fn turn_segments(self) -> u32 {
        self.quadrant_segments.saturating_mul(4)
    }
}

/// The set of points within `distance` of `geometry`.
///
/// A positive distance dilates, a negative distance erodes, and zero returns
/// the geometry's own areal part — which is the limit of both and is empty for
/// anything without area.
///
/// Returns an empty `MultiPolygon` rather than an error for input a buffer is
/// not defined on, so a column of mixed geometries buffers row by row without
/// one row failing the batch. A non-finite distance is likewise empty: there is
/// no region it describes.
#[must_use]
pub fn buffer(geometry: &Geometry, distance: f64, style: BufferStyle) -> MultiPolygon {
    if !distance.is_finite() {
        return MultiPolygon::new(Vec::new());
    }
    if distance == 0.0 {
        return areal_part(geometry);
    }
    if distance < 0.0 {
        return erode(geometry, -distance, style);
    }
    let mut parts: Vec<Polygon> = Vec::new();
    // Areal input contributes its own interior: dilation only adds to a
    // region, and the per-boundary parts below cover the collar, not the
    // middle.
    parts.extend(areal_part(geometry).0);
    collect_dilation_parts(geometry, distance, style, &mut parts);
    unary_union(parts.iter())
}

/// Erosion: the region at least `distance` inside the geometry.
fn erode(geometry: &Geometry, distance: f64, style: BufferStyle) -> MultiPolygon {
    let interior = areal_part(geometry);
    if interior.0.is_empty() {
        // A point or a curve has no interior to retreat into. Empty is the
        // honest answer rather than the input unchanged.
        return MultiPolygon::new(Vec::new());
    }
    // Every ring, dilated outward, is the collar to remove. Removing it from
    // the interior leaves exactly the points further than `distance` from the
    // boundary.
    let mut collar: Vec<Polygon> = Vec::new();
    for polygon in &interior.0 {
        for ring in std::iter::once(polygon.exterior()).chain(polygon.interiors()) {
            ring_dilation_parts(ring, distance, style, &mut collar);
        }
    }
    if collar.is_empty() {
        return interior;
    }
    interior.difference(&unary_union(collar.iter()))
}

/// The geometry's polygonal components, which a dilation grows from and an
/// erosion shrinks within.
fn areal_part(geometry: &Geometry) -> MultiPolygon {
    match geometry {
        Geometry::Polygon(polygon) => MultiPolygon::new(vec![polygon.clone()]),
        Geometry::MultiPolygon(polygons) => polygons.clone(),
        Geometry::Rect(rect) => MultiPolygon::new(vec![rect.to_polygon()]),
        Geometry::Triangle(triangle) => MultiPolygon::new(vec![triangle.to_polygon()]),
        Geometry::GeometryCollection(collection) => {
            let mut held: Vec<Polygon> = Vec::new();
            for member in collection {
                held.extend(areal_part(member).0);
            }
            MultiPolygon::new(held)
        }
        _ => MultiPolygon::new(Vec::new()),
    }
}

/// One convex part per vertex and per segment of every component of
/// `geometry`.
fn collect_dilation_parts(
    geometry: &Geometry,
    distance: f64,
    style: BufferStyle,
    out: &mut Vec<Polygon>,
) {
    match geometry {
        Geometry::Point(p) => out.push(disc(p.0.x, p.0.y, distance, style)),
        Geometry::MultiPoint(points) => {
            for p in points {
                out.push(disc(p.0.x, p.0.y, distance, style));
            }
        }
        Geometry::Line(segment) => {
            chain_dilation_parts(&[segment.start, segment.end], distance, style, out);
        }
        Geometry::LineString(line) => {
            ring_dilation_parts(line, distance, style, out);
        }
        Geometry::MultiLineString(lines) => {
            for line in &lines.0 {
                ring_dilation_parts(line, distance, style, out);
            }
        }
        Geometry::Polygon(polygon) => {
            for ring in std::iter::once(polygon.exterior()).chain(polygon.interiors()) {
                ring_dilation_parts(ring, distance, style, out);
            }
        }
        Geometry::MultiPolygon(polygons) => {
            for polygon in polygons {
                for ring in std::iter::once(polygon.exterior()).chain(polygon.interiors()) {
                    ring_dilation_parts(ring, distance, style, out);
                }
            }
        }
        Geometry::Rect(rect) => {
            let polygon = rect.to_polygon();
            ring_dilation_parts(polygon.exterior(), distance, style, out);
        }
        Geometry::Triangle(triangle) => {
            let polygon = triangle.to_polygon();
            ring_dilation_parts(polygon.exterior(), distance, style, out);
        }
        Geometry::GeometryCollection(collection) => {
            for member in collection {
                collect_dilation_parts(member, distance, style, out);
            }
        }
    }
}

fn ring_dilation_parts(
    ring: &LineString,
    distance: f64,
    style: BufferStyle,
    out: &mut Vec<Polygon>,
) {
    let coords: Vec<super::Coord> = ring.coords().copied().collect();
    chain_dilation_parts(&coords, distance, style, out);
}

/// A disc at every vertex and a rectangle along every segment.
///
/// The discs are what make the joins and caps round, and they also cover the
/// degenerate cases for free: a chain of one vertex is a point buffer, and a
/// zero-length segment contributes a rectangle of no extent that the union
/// discards while its endpoint discs remain.
fn chain_dilation_parts(
    coords: &[super::Coord],
    distance: f64,
    style: BufferStyle,
    out: &mut Vec<Polygon>,
) {
    for coord in coords {
        if coord.x.is_finite() && coord.y.is_finite() {
            out.push(disc(coord.x, coord.y, distance, style));
        }
    }
    for pair in coords.windows(2) {
        let (from, to) = (pair[0], pair[1]);
        if !(from.x.is_finite() && from.y.is_finite() && to.x.is_finite() && to.y.is_finite()) {
            continue;
        }
        let (dx, dy) = (to.x - from.x, to.y - from.y);
        let length = dx.hypot(dy);
        if length == 0.0 {
            continue;
        }
        // The unit normal, scaled to the distance: offsetting both endpoints
        // by ±this gives the segment's rectangle.
        let (nx, ny) = (-dy / length * distance, dx / length * distance);
        out.push(Polygon::new(
            LineString::from(vec![
                (from.x + nx, from.y + ny),
                (to.x + nx, to.y + ny),
                (to.x - nx, to.y - ny),
                (from.x - nx, from.y - ny),
                (from.x + nx, from.y + ny),
            ]),
            Vec::new(),
        ));
    }
}

/// A regular polygon approximating the disc of `radius` about a point.
fn disc(x: f64, y: f64, radius: f64, style: BufferStyle) -> Polygon {
    let segments = style.turn_segments();
    let mut coords = Vec::with_capacity(segments as usize + 1);
    // Half a step of rotation, so no vertex lands at exactly 0, 90, 180 or
    // 270 degrees. Those are the directions a segment's offset rectangle has
    // its edges along, and a disc vertex landing exactly on a rectangle edge
    // is a degenerate tangency for the union sweep rather than a crossing.
    // Rotating the disc keeps it the same shape and the same inscribed area
    // while making every boundary meeting transversal.
    let half_step = TAU / (2.0 * f64::from(segments));
    for step in 0..segments {
        let angle = half_step + TAU * f64::from(step) / f64::from(segments);
        coords.push((x + radius * angle.cos(), y + radius * angle.sin()));
    }
    // Close the ring on the exact first coordinate rather than on a recomputed
    // one, so the closure is bitwise rather than approximate.
    if let Some(first) = coords.first().copied() {
        coords.push(first);
    }
    Polygon::new(LineString::from(coords), Vec::new())
}

#[cfg(test)]
mod tests;
