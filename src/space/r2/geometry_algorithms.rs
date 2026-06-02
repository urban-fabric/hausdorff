//! Trait dispatch over the [`Geometry`] enum, so the planar algorithms run on a
//! heterogeneous geometry (not only the concrete types) — the native equivalent
//! of calling `geo`'s traits on a `geo::Geometry`.

use super::convex_hull::convex_hull_coords;
use super::{
    Area, BoundingRect, Centroid, ChamberlainDuquetteArea, ConvexHull, Coord, CoordsIter, Euclidean,
    Geometry, Length, LineString, MultiPoint, MultiPolygon, Point, Polygon, Rect, Simplify, Winding,
};

impl ChamberlainDuquetteArea for Geometry {
    fn chamberlain_duquette_signed_area(&self) -> f64 {
        match self {
            Geometry::Polygon(g) => g.chamberlain_duquette_signed_area(),
            Geometry::MultiPolygon(g) => g.chamberlain_duquette_signed_area(),
            Geometry::Rect(g) => g.to_polygon().chamberlain_duquette_signed_area(),
            Geometry::Triangle(g) => g.to_polygon().chamberlain_duquette_signed_area(),
            Geometry::GeometryCollection(gs) => {
                gs.iter().map(ChamberlainDuquetteArea::chamberlain_duquette_signed_area).sum()
            }
            _ => 0.0,
        }
    }
}

impl Area for Geometry {
    fn signed_area(&self) -> f64 {
        match self {
            Geometry::Polygon(g) => g.signed_area(),
            Geometry::MultiPolygon(g) => g.signed_area(),
            Geometry::Rect(g) => g.signed_area(),
            Geometry::Triangle(g) => g.signed_area(),
            Geometry::GeometryCollection(gs) => gs.iter().map(Area::signed_area).sum(),
            // points and curves have no area
            Geometry::Point(_)
            | Geometry::Line(_)
            | Geometry::LineString(_)
            | Geometry::MultiPoint(_)
            | Geometry::MultiLineString(_) => 0.0,
        }
    }

    fn unsigned_area(&self) -> f64 {
        // A collection's unsigned area is the sum of its parts' *unsigned* areas,
        // not `|Σ signed|` — members of opposite winding must not cancel.
        match self {
            Geometry::GeometryCollection(gs) => gs.iter().map(Area::unsigned_area).sum(),
            other => other.signed_area().abs(),
        }
    }
}

impl BoundingRect for Geometry {
    type Output = Option<Rect>;

    fn bounding_rect(&self) -> Option<Rect> {
        let (mut min_x, mut min_y) = (f64::INFINITY, f64::INFINITY);
        let (mut max_x, mut max_y) = (f64::NEG_INFINITY, f64::NEG_INFINITY);
        let mut any = false;
        for c in self.coords_iter() {
            if !c.is_finite() {
                continue;
            }
            any = true;
            min_x = min_x.min(c.x);
            min_y = min_y.min(c.y);
            max_x = max_x.max(c.x);
            max_y = max_y.max(c.y);
        }
        any.then(|| Rect::new(Coord::new(min_x, min_y), Coord::new(max_x, max_y)))
    }
}

impl ConvexHull for Geometry {
    fn convex_hull(&self) -> Polygon {
        let coords: Vec<Coord> = self.coords_iter().collect();
        Polygon::new(LineString::new(convex_hull_coords(&coords)), vec![])
    }
}

impl Centroid for Geometry {
    fn centroid(&self) -> Option<Point> {
        match self {
            Geometry::Point(g) => g.centroid(),
            Geometry::LineString(g) => g.centroid(),
            Geometry::Polygon(g) => g.centroid(),
            Geometry::MultiPolygon(g) => g.centroid(),
            Geometry::MultiPoint(g) => g.centroid(),
            Geometry::Rect(g) => Some(g.center()),
            Geometry::Triangle(g) => g.to_polygon().centroid(),
            Geometry::Line(l) => Some(Point((l.start + l.end) / 2.0)),
            Geometry::MultiLineString(mls) => {
                // Length-weighted mean of the component centroids (a long line
                // pulls the centroid more than a short one).
                let mut acc = Coord::new(0.0, 0.0);
                let mut total = 0.0;
                for ls in mls.iter() {
                    if let Some(c) = ls.centroid() {
                        let len = Euclidean.length(ls);
                        acc = acc + c.0 * len;
                        total += len;
                    }
                }
                if total > 0.0 {
                    Some(Point(acc / total))
                } else {
                    // every component is zero-length: fall back to their mean
                    let pts: Vec<Point> = mls.iter().filter_map(Centroid::centroid).collect();
                    MultiPoint::new(pts).centroid()
                }
            }
            Geometry::GeometryCollection(_) => {
                // The mean of every coordinate in the collection.
                let mut sum = Coord::new(0.0, 0.0);
                let mut n = 0.0;
                for c in self.coords_iter() {
                    sum = sum + c;
                    n += 1.0;
                }
                (n > 0.0).then(|| Point(sum / n))
            }
        }
    }
}

impl Simplify for Geometry {
    fn simplify(&self, epsilon: f64) -> Self {
        match self {
            Geometry::LineString(g) => Geometry::LineString(g.simplify(epsilon)),
            Geometry::Polygon(g) => Geometry::Polygon(g.simplify(epsilon)),
            Geometry::MultiLineString(g) => Geometry::MultiLineString(g.simplify(epsilon)),
            Geometry::MultiPolygon(g) => Geometry::MultiPolygon(g.simplify(epsilon)),
            // recurse into a collection, simplifying each member
            Geometry::GeometryCollection(gs) => {
                Geometry::GeometryCollection(gs.iter().map(|g| g.simplify(epsilon)).collect())
            }
            // nothing to simplify on points / lines / rects / triangles
            other => other.clone(),
        }
    }
}

/// OGC-style accessors and classifiers over the [`Geometry`] enum — the heterogeneous
/// equivalents of the per-type methods, so a column of mixed geometry can be queried
/// without first matching the variant. Indices are **0-based** (idiomatic Rust); the
/// 1-based OGC `ST_*` convention is applied by the caller (e.g. galois's `wkb_*`).
/// "Wrong geometry kind" returns `None`, mirroring SQL null.
impl Geometry {
    /// The exterior ring of a polygon as a [`LineString`]; `None` otherwise
    /// (`ST_ExteriorRing`).
    pub fn exterior_ring(&self) -> Option<LineString> {
        match self {
            Geometry::Polygon(p) => Some(p.exterior().clone()),
            _ => None,
        }
    }

    /// The `n`-th interior ring (hole) of a polygon, 0-indexed; `None` for
    /// non-polygons or out of range (`ST_InteriorRingN`).
    pub fn interior_ring(&self, n: usize) -> Option<LineString> {
        match self {
            Geometry::Polygon(p) => p.interiors().get(n).cloned(),
            _ => None,
        }
    }

    /// Number of interior rings (holes) of a polygon; `None` otherwise
    /// (`ST_NumInteriorRings`).
    pub fn num_interior_rings(&self) -> Option<usize> {
        match self {
            Geometry::Polygon(p) => Some(p.interiors().len()),
            _ => None,
        }
    }

    /// The `n`-th sub-geometry, 0-indexed: the nth member of a multi-geometry or
    /// collection, or the geometry itself when `n == 0` for a single geometry;
    /// `None` out of range (`ST_GeometryN`).
    pub fn geometry_n(&self, n: usize) -> Option<Geometry> {
        match self {
            Geometry::MultiPoint(mp) => mp.iter().nth(n).map(|p| Geometry::Point(*p)),
            Geometry::MultiLineString(mls) => mls.iter().nth(n).cloned().map(Geometry::LineString),
            Geometry::MultiPolygon(mp) => mp.iter().nth(n).cloned().map(Geometry::Polygon),
            Geometry::GeometryCollection(c) => c.get(n).cloned(),
            single => (n == 0).then(|| single.clone()),
        }
    }

    /// Total number of sub-geometries: the member count of a multi-geometry or
    /// collection, else `1` (`ST_NumGeometries`).
    pub fn num_geometries(&self) -> usize {
        match self {
            Geometry::MultiPoint(mp) => mp.0.len(),
            Geometry::MultiLineString(mls) => mls.0.len(),
            Geometry::MultiPolygon(mp) => mp.0.len(),
            Geometry::GeometryCollection(gc) => gc.len(),
            _ => 1,
        }
    }

    /// The `n`-th vertex of a `LineString`, 0-indexed, as a [`Point`]; `None` for
    /// non-lines or out of range (`ST_PointN`).
    pub fn point_n(&self, n: usize) -> Option<Point> {
        match self {
            Geometry::LineString(ls) => ls.0.get(n).map(|c| Point(*c)),
            _ => None,
        }
    }

    /// The first vertex of a `LineString`; `None` for non-lines or an empty line
    /// (`ST_StartPoint`).
    pub fn start_point(&self) -> Option<Point> {
        match self {
            Geometry::LineString(ls) => ls.0.first().map(|c| Point(*c)),
            _ => None,
        }
    }

    /// The last vertex of a `LineString`; `None` for non-lines or an empty line
    /// (`ST_EndPoint`).
    pub fn end_point(&self) -> Option<Point> {
        match self {
            Geometry::LineString(ls) => ls.0.last().map(|c| Point(*c)),
            _ => None,
        }
    }

    /// Whether a curve is closed (first vertex equals last): a `LineString`
    /// directly, a `MultiLineString` iff every part is; `None` for non-curves
    /// (`ST_IsClosed`).
    pub fn is_closed(&self) -> Option<bool> {
        match self {
            Geometry::LineString(ls) => Some(ls.is_closed()),
            Geometry::MultiLineString(mls) => Some(mls.iter().all(LineString::is_closed)),
            _ => None,
        }
    }

    /// Whether a polygon's exterior ring is wound clockwise; `None` for
    /// non-polygons (`ST_IsPolygonCW`).
    pub fn is_polygon_cw(&self) -> Option<bool> {
        match self {
            Geometry::Polygon(p) => Some(p.exterior().is_cw()),
            _ => None,
        }
    }

    /// Whether a polygon's exterior ring is wound counter-clockwise; `None` for
    /// non-polygons (`ST_IsPolygonCCW`).
    pub fn is_polygon_ccw(&self) -> Option<bool> {
        match self {
            Geometry::Polygon(p) => Some(p.exterior().is_ccw()),
            _ => None,
        }
    }

    /// OGC geometry-type name (`"Point"`, `"LineString"`, …) (`ST_GeometryType`).
    pub fn geometry_type(&self) -> &'static str {
        match self {
            Geometry::Point(_) => "Point",
            Geometry::Line(_) => "Line",
            Geometry::LineString(_) => "LineString",
            Geometry::Polygon(_) => "Polygon",
            Geometry::MultiPoint(_) => "MultiPoint",
            Geometry::MultiLineString(_) => "MultiLineString",
            Geometry::MultiPolygon(_) => "MultiPolygon",
            Geometry::GeometryCollection(_) => "GeometryCollection",
            Geometry::Rect(_) => "Rect",
            Geometry::Triangle(_) => "Triangle",
        }
    }

    /// View the areal parts of a geometry as a [`MultiPolygon`] (a single polygon
    /// wrapped; `Rect`/`Triangle` converted); `None` for non-areal geometries. The
    /// areal coercion the boolean-overlay ops consume.
    pub fn as_multi_polygon(&self) -> Option<MultiPolygon> {
        match self {
            Geometry::Polygon(p) => Some(MultiPolygon::new(vec![p.clone()])),
            Geometry::MultiPolygon(mp) => Some(mp.clone()),
            Geometry::Rect(rect) => Some(MultiPolygon::new(vec![rect.to_polygon()])),
            Geometry::Triangle(triangle) => Some(MultiPolygon::new(vec![triangle.to_polygon()])),
            _ => None,
        }
    }
}
