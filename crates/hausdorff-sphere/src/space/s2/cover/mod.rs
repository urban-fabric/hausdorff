//! Cover an r2 [`Geometry`] with S2 cells — the geometry-enum dispatch over
//! [`RegionCoverer`]. The `Series`/WKB plumbing stays in the caller; this is
//! the pure geometry → cells core. [`cover_r2_geometry_at_level`] emits a
//! fixed-level cover; [`cover_r2_geometry`] emits an adaptive multi-level
//! [`CellUnion`].

use std::collections::BTreeSet;

use hausdorff_base::earth::GeoPoint;
use hausdorff_planar::space::r2::{Coord, Geometry, LineString, Polygon};

use super::coords::MAX_LEVEL;
use super::{
    CellUnion, Loop, Polygon as S2Polygon, Polyline, Rect, RegionCoverer, S2CellId, S2Point,
};

/// S2 cells at `level` covering `geom`:
/// - **points** → the containing cell;
/// - **areal** geometries → the [`RegionCoverer`] cover (every level-`level` cell
///   that intersects the region — boundary-complete, not a centroid sample);
/// - **lines** → the [`Polyline`] boundary cover (every level-`level` cell the
///   chain of great-circle edges passes through);
/// - **collections** recurse.
///
/// Returns `None` when `level` exceeds the leaf level (30), or when a polygon
/// ring or line has a non-finite or out-of-range coordinate. An out-of-range
/// *point* coordinate is skipped (a point is a single sample, not a malformed
/// ring); a degenerate hole (fewer than three vertices, hence no area) is
/// dropped.
#[must_use]
pub fn cover_r2_geometry_at_level(geom: &Geometry, level: u32) -> Option<Vec<S2CellId>> {
    if level > MAX_LEVEL {
        return None;
    }
    let mut cells = BTreeSet::new();
    walk_parts(
        geom,
        &mut FixedCover {
            level,
            out: &mut cells,
        },
    )?;
    Some(cells.into_iter().collect())
}

/// An adaptive multi-level cover of `geom` driven by `coverer`
/// (see [`RegionCoverer::adaptive_covering`]): points map to their cell at
/// `coverer.max_level`, areal geometries to the adaptive region cover, lines to
/// the adaptive [`Polyline`] boundary cover, and collections to the normalized
/// union of their parts. Coordinate validation matches
/// [`cover_r2_geometry_at_level`].
#[must_use]
pub fn cover_r2_geometry(geom: &Geometry, coverer: &RegionCoverer) -> Option<CellUnion> {
    let mut cells = Vec::new();
    walk_parts(
        geom,
        &mut AdaptiveCover {
            coverer,
            point_level: coverer.max_level.min(MAX_LEVEL),
            out: &mut cells,
        },
    )?;
    Some(CellUnion::from(cells))
}

/// The antimeridian-correct lat/lng bounding [`Rect`] of an r2 (lon/lat)
/// [`Geometry`] — the wrap-aware peer of the geometry's planar `bounding_rect`.
/// Every finite, in-range vertex is folded in via [`Rect::add_point`], so a ring
/// straddling ±180° yields a tight *inverted* longitude span instead of the
/// globe-spanning `[-180°, 180°]` a min/max box produces. `None` when the geometry
/// has no valid vertex (empty, or every coordinate non-finite / out of range).
///
/// This is a **vertex** bound (see [`Rect::from_geo_points`]): the latitude span is
/// the vertex min/max, matching the planar `bounding_rect` away from the poles,
/// while the longitude span is wrap-aware. Points recurse through collections.
#[must_use]
pub fn bound_r2_geometry(geom: &Geometry) -> Option<Rect> {
    let mut pts = Vec::new();
    walk_parts(geom, &mut VertexCollector { out: &mut pts });
    (!pts.is_empty()).then(|| Rect::from_geo_points(&pts))
}

/// The great-circle-**arc**-aware lat/lng bound of an r2 geometry — like
/// [`bound_r2_geometry`], but each edge's latitude span covers its great-circle
/// arc (which can ride higher than either endpoint; see [`super::arc::rect_bound`]),
/// so the bound conservatively contains the geometry *as the spherical predicates
/// interpret it*. This is the bound a spherical-predicate broad phase must use: a
/// vertex-only box can miss an intersection where two edges bow toward each other
/// between their endpoints. `None` when the geometry has no valid vertex.
#[must_use]
pub fn bound_r2_geometry_arc(geom: &Geometry) -> Option<Rect> {
    let mut arcs = ArcBounder {
        bound: Rect::empty(),
    };
    walk_parts(geom, &mut arcs);
    (!arcs.bound.is_empty()).then_some(arcs.bound)
}

/// A sink for the primitive parts of an r2 [`Geometry`], visited in document
/// order by [`walk_parts`].
///
/// Every geometry→S2 pass here (fixed-level cover, adaptive cover, vertex bound,
/// arc bound) shares the exact same [`Geometry`]-enum walk: `Multi*` iterate,
/// [`Geometry::Rect`]/[`Geometry::Triangle`] fold in as their `to_polygon()`, and
/// [`Geometry::GeometryCollection`] recurses. Hoisting that skeleton into one
/// walker keeps the four passes from each re-deriving it (and re-deriving the
/// exact Rect/Triangle canonicalization), so a new geometry variant is handled in
/// one place instead of four. Each pass only says what to do with the three
/// leaf kinds — a point, a coordinate chain, or a polygon.
///
/// Callbacks return [`Option`] so a validation failure short-circuits the whole
/// walk (mirroring the `?` in the cover passes); an infallible pass returns
/// `Some(())`.
trait GeometryParts {
    /// A single lon/lat coordinate (a `Point`, or one member of a `MultiPoint`).
    fn point(&mut self, lon: f64, lat: f64) -> Option<()>;
    /// A chain of lon/lat coordinates that forms a line (a `Line`'s two
    /// endpoints, or one `LineString`).
    fn line(&mut self, coords: &[Coord]) -> Option<()>;
    /// One polygon, exterior plus holes (also how `Rect`/`Triangle` arrive).
    fn polygon(&mut self, poly: &Polygon) -> Option<()>;
}

/// Walk `geom`'s primitive parts, dispatching each to `parts` — the single copy
/// of the [`Geometry`]-enum plumbing every pass in this module shares. Returns
/// `None` as soon as a callback does.
fn walk_parts<V: GeometryParts>(geom: &Geometry, parts: &mut V) -> Option<()> {
    match geom {
        Geometry::Point(p) => parts.point(p.x(), p.y())?,
        Geometry::MultiPoint(mp) => {
            for p in mp.iter() {
                parts.point(p.x(), p.y())?;
            }
        }
        Geometry::Line(l) => parts.line(&[l.start, l.end])?,
        Geometry::LineString(ls) => parts.line(&ls.0)?,
        Geometry::MultiLineString(mls) => {
            for ls in mls.iter() {
                parts.line(&ls.0)?;
            }
        }
        Geometry::Polygon(poly) => parts.polygon(poly)?,
        Geometry::MultiPolygon(mp) => {
            for poly in mp.iter() {
                parts.polygon(poly)?;
            }
        }
        Geometry::Rect(r) => parts.polygon(&r.to_polygon())?,
        Geometry::Triangle(t) => parts.polygon(&t.to_polygon())?,
        Geometry::GeometryCollection(gs) => {
            for g in gs {
                walk_parts(g, parts)?;
            }
        }
    }
    Some(())
}

/// The valid lon/lat vertices of one ring or line, in order.
fn ring_geo_points(ls: &LineString) -> Vec<GeoPoint> {
    let mut out = Vec::new();
    for c in ls.coords() {
        push_coord(c.x, c.y, &mut out);
    }
    out
}

/// Fold one vertex sequence into `bound`: each consecutive pair as a great-circle
/// arc bound, or a lone point directly.
fn fold_arc_seq(pts: &[GeoPoint], bound: Rect) -> Rect {
    match pts {
        [] => bound,
        [p] => bound.add_point(*p),
        _ => pts.windows(2).fold(bound, |acc, w| {
            acc.union(&super::arc::rect_bound(
                S2Point::from_geo_point(w[0]),
                S2Point::from_geo_point(w[1]),
            ))
        }),
    }
}

fn fold_polygon_arcs(poly: &Polygon, bound: Rect) -> Rect {
    let mut b = std::iter::once(poly.exterior())
        .chain(poly.interiors())
        .fold(bound, |acc, ring| fold_arc_seq(&ring_geo_points(ring), acc));
    // A polygon whose interior contains a pole must bound that pole: the pole is
    // interior (and at every longitude), not on any edge, so the arc fold over
    // the boundary misses it — which would let a broad phase drop a point/region
    // at the pole that the polygon actually contains.
    if let Some(shell) = ring_loop(poly.exterior()) {
        for lat in [90.0_f64, -90.0] {
            let pole = S2Point::from_geo_point(GeoPoint::from_lat_lon(lat, 0.0));
            if shell.contains_point(pole) {
                b = b.union(&Rect::from_degrees(lat, -180.0, lat, 180.0));
            }
        }
    }
    b
}

/// Accumulates the great-circle-arc lat/lng bound of a geometry (see
/// [`bound_r2_geometry_arc`]): points and line chains fold in as arc sequences,
/// polygons via [`fold_polygon_arcs`] (which also captures an enclosed pole).
struct ArcBounder {
    bound: Rect,
}

impl GeometryParts for ArcBounder {
    fn point(&mut self, lon: f64, lat: f64) -> Option<()> {
        self.bound = fold_arc_seq(&ring_geo_points_of_coords([(lon, lat)]), self.bound);
        Some(())
    }
    fn line(&mut self, coords: &[Coord]) -> Option<()> {
        let pts = ring_geo_points_of_coords(coords.iter().map(|c| (c.x, c.y)));
        self.bound = fold_arc_seq(&pts, self.bound);
        Some(())
    }
    fn polygon(&mut self, poly: &Polygon) -> Option<()> {
        self.bound = fold_polygon_arcs(poly, self.bound);
        Some(())
    }
}

/// Valid lon/lat points from raw `(x, y)` coordinate pairs, in order.
fn ring_geo_points_of_coords<I: IntoIterator<Item = (f64, f64)>>(coords: I) -> Vec<GeoPoint> {
    let mut out = Vec::new();
    for (x, y) in coords {
        push_coord(x, y, &mut out);
    }
    out
}

/// Appends every finite, in-range lon/lat vertex of a geometry — the vertex set
/// behind [`bound_r2_geometry`]. Points, line vertices, and polygon-ring
/// vertices all contribute alike.
struct VertexCollector<'a> {
    out: &'a mut Vec<GeoPoint>,
}

impl GeometryParts for VertexCollector<'_> {
    fn point(&mut self, lon: f64, lat: f64) -> Option<()> {
        push_coord(lon, lat, self.out);
        Some(())
    }
    fn line(&mut self, coords: &[Coord]) -> Option<()> {
        for c in coords {
            push_coord(c.x, c.y, self.out);
        }
        Some(())
    }
    fn polygon(&mut self, poly: &Polygon) -> Option<()> {
        push_polygon(poly, self.out);
        Some(())
    }
}

fn push_polygon(poly: &Polygon, out: &mut Vec<GeoPoint>) {
    for ring in std::iter::once(poly.exterior()).chain(poly.interiors()) {
        for c in ring.coords() {
            push_coord(c.x, c.y, out);
        }
    }
}

fn push_coord(lon: f64, lat: f64, out: &mut Vec<GeoPoint>) {
    if valid_lonlat(lon, lat) {
        out.push(GeoPoint::from_lon_lat(lon, lat));
    }
}

/// Collects a fixed-level cover into a deduplicating [`BTreeSet`]: points map to
/// their containing cell, lines to their [`Polyline`] cover, polygons to their
/// [`RegionCoverer`] cover — all at `level`. Backs [`cover_r2_geometry_at_level`].
struct FixedCover<'a> {
    level: u32,
    out: &'a mut BTreeSet<S2CellId>,
}

impl GeometryParts for FixedCover<'_> {
    fn point(&mut self, lon: f64, lat: f64) -> Option<()> {
        insert_point(lon, lat, self.level, self.out);
        Some(())
    }
    fn line(&mut self, coords: &[Coord]) -> Option<()> {
        self.out.extend(cover_line_coords(
            coords,
            &RegionCoverer::at_level(self.level),
        )?);
        Some(())
    }
    fn polygon(&mut self, poly: &Polygon) -> Option<()> {
        cover_polygon(poly, self.level, self.out)
    }
}

/// Collects an adaptive multi-level cover into a `Vec` (later normalized into a
/// [`CellUnion`]): points map to their cell at `point_level`, lines to the
/// adaptive [`Polyline`] cover, polygons to the adaptive region cover. Backs
/// [`cover_r2_geometry`].
struct AdaptiveCover<'a> {
    coverer: &'a RegionCoverer,
    /// `coverer.max_level` clamped to the leaf level, cached so each point does
    /// not recompute it.
    point_level: u32,
    out: &'a mut Vec<S2CellId>,
}

impl GeometryParts for AdaptiveCover<'_> {
    fn point(&mut self, lon: f64, lat: f64) -> Option<()> {
        if valid_lonlat(lon, lat) {
            let gp = GeoPoint::from_lon_lat(lon, lat);
            self.out
                .push(S2CellId::from_geo_point_at_level(gp, self.point_level));
        }
        Some(())
    }
    fn line(&mut self, coords: &[Coord]) -> Option<()> {
        self.out.extend(cover_line_coords(coords, self.coverer)?);
        Some(())
    }
    fn polygon(&mut self, poly: &Polygon) -> Option<()> {
        self.out
            .extend(adaptive_polygon(poly, self.coverer)?.iter());
        Some(())
    }
}

/// True iff `(lon, lat)` is a finite, in-range WGS84 coordinate.
pub(super) fn valid_lonlat(lon: f64, lat: f64) -> bool {
    lon.is_finite()
        && lat.is_finite()
        && (-90.0..=90.0).contains(&lat)
        && (-180.0..=180.0).contains(&lon)
}

fn insert_point(lon: f64, lat: f64, level: u32, out: &mut BTreeSet<S2CellId>) {
    if valid_lonlat(lon, lat) {
        let gp = GeoPoint::from_lon_lat(lon, lat);
        out.insert(S2CellId::from_geo_point_at_level(gp, level));
    }
}

fn cover_polygon(poly: &Polygon, level: u32, out: &mut BTreeSet<S2CellId>) -> Option<()> {
    let s2 = polygon_region(poly)?;
    out.extend(RegionCoverer::at_level(level).covering(&s2));
    Some(())
}

fn adaptive_polygon(poly: &Polygon, coverer: &RegionCoverer) -> Option<CellUnion> {
    Some(coverer.adaptive_covering(&polygon_region(poly)?))
}

/// The lon/lat polygon as a spherical [`S2Polygon`]; `None` when a ring has an
/// invalid coordinate.
pub(super) fn polygon_region(poly: &Polygon) -> Option<S2Polygon> {
    for ring in std::iter::once(poly.exterior()).chain(poly.interiors()) {
        if !ring.coords().all(|c| valid_lonlat(c.x, c.y)) {
            return None;
        }
    }
    let shell = ring_loop(poly.exterior())?;
    let holes: Vec<Loop> = poly.interiors().iter().filter_map(ring_loop).collect();
    Some(S2Polygon::new(shell, holes))
}

/// The cells covering a chain of lon/lat coordinates: empty for no
/// coordinates, the containing cell for one, and the [`Polyline`] cover
/// otherwise; `None` when any coordinate is invalid.
fn cover_line_coords(coords: &[Coord], coverer: &RegionCoverer) -> Option<Vec<S2CellId>> {
    if !coords.iter().all(|c| valid_lonlat(c.x, c.y)) {
        return None;
    }
    let level = coverer.max_level.min(MAX_LEVEL);
    match coords {
        [] => Some(Vec::new()),
        [c] => {
            let gp = GeoPoint::from_lon_lat(c.x, c.y);
            Some(vec![S2CellId::from_geo_point_at_level(gp, level)])
        }
        _ => {
            let pts: Vec<S2Point> = coords
                .iter()
                .map(|c| S2Point::from_geo_point(GeoPoint::from_lon_lat(c.x, c.y)))
                .collect();
            let line = Polyline::new(pts);
            if coverer.min_level == coverer.max_level && coverer.max_cells == usize::MAX {
                Some(coverer.covering(&line))
            } else {
                Some(coverer.adaptive_covering(&line).cell_ids().to_vec())
            }
        }
    }
}

/// A lon/lat ring as a counter-clockwise spherical [`Loop`] (interior to the left
/// — the small side); `None` for a ring with fewer than three vertices.
fn ring_loop(ring: &LineString) -> Option<Loop> {
    let pts: Vec<S2Point> = ring
        .coords()
        .map(|c| S2Point::from_geo_point(GeoPoint::from_lon_lat(c.x, c.y)))
        .collect();
    let l = Loop::new(pts);
    if l.num_vertices() < 3 {
        return None;
    }
    if l.area() > std::f64::consts::PI * 2.0 {
        let mut v = l.vertices().to_vec();
        v.reverse();
        Some(Loop::new(v))
    } else {
        Some(l)
    }
}

#[cfg(test)]
mod tests;
