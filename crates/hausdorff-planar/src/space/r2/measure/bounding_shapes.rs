//! Bounding-shape descriptors: the smallest enclosing circle, the smallest
//! enclosing rectangle at *any* angle, and the largest circle that fits inside.
//!
//! The axis-aligned [`BoundingRect`](crate::space::r2::BoundingRect) answers "roughly where
//! is this". These answer questions about *shape*, and they are the standard
//! vocabulary for it — `ST_MinimumBoundingCircle`, `ST_OrientedEnvelope`, and
//! `ST_MaximumInscribedCircle` in PostGIS, the same three in JTS and Shapely.
//!
//! The reason they are worth having rather than approximating from the
//! axis-aligned box: the box is an artifact of the coordinate axes, not of the
//! geometry. A long thin parcel at 45° has a bounding box nearly twice its own
//! area, so "how elongated is this" and "which way does it point" cannot be
//! answered from it at all. The oriented envelope can answer both.
//!
//! ## All three reduce to the convex hull
//!
//! Each of the first two depends only on the hull, which is what makes them cheap
//! here — the hull already exists and is already tested.
//!
//! - **Minimum bounding circle** is determined by two or three hull points, so
//!   Welzl's algorithm needs only those.
//! - **Oriented envelope**: the minimum-area enclosing rectangle always has a side
//!   flush with a hull edge. That is the whole content of rotating calipers — it
//!   turns a search over all angles into a loop over hull edges.
//! - **Maximum inscribed circle** is the one that does *not* reduce to the hull,
//!   because the largest circle inside a concave polygon is generally not inside
//!   its hull's incircle. It needs the polygon itself.

use std::cmp::Ordering;
use std::collections::BinaryHeap;

use crate::space::r2::convex_hull_coords;
use crate::space::r2::{
    Area, Contains, Coord, CoordsIter, Distance, Euclidean, Geometry, LineString, Point, Polygon,
};

/// A circle as a centre and a radius — the result of the two circle descriptors.
///
/// Deliberately not a geometry: a circle is not representable in the linear
/// geometry model, and returning an approximating polygon would bake a segment
/// count into the answer. Use [`Circle::to_polygon`] when a geometry is wanted,
/// where the caller chooses the fidelity.
#[derive(Clone, Copy, PartialEq, Debug)]
pub struct Circle {
    pub centre: Coord,
    pub radius: f64,
}

impl Circle {
    /// Approximate the circle as a regular polygon with `segments` sides.
    ///
    /// The polygon is *inscribed*, so it lies inside the circle and its area is
    /// smaller. Worth knowing when the result feeds an area comparison.
    #[must_use]
    pub fn to_polygon(self, segments: usize) -> Polygon {
        let n = segments.max(3);
        let ring: Vec<Coord> = (0..n)
            .map(|i| {
                let t = std::f64::consts::TAU * (i as f64) / (n as f64);
                Coord::new(
                    self.centre.x + self.radius * t.cos(),
                    self.centre.y + self.radius * t.sin(),
                )
            })
            .collect();
        Polygon::new(LineString::new(ring), vec![])
    }

    /// Whether `p` is inside or on the circle, with a relative tolerance so a
    /// point placed exactly on the rim is not excluded by rounding.
    #[must_use]
    pub fn contains(self, p: Coord) -> bool {
        let dx = p.x - self.centre.x;
        let dy = p.y - self.centre.y;
        dx * dx + dy * dy <= self.radius * self.radius * (1.0 + 1e-12) + 1e-12
    }
}

/// The smallest enclosing rectangle at any orientation.
#[derive(Clone, Copy, PartialEq, Debug)]
pub struct OrientedRect {
    /// The four corners, counter-clockwise.
    pub corners: [Coord; 4],
    /// Length of the longer side.
    pub length: f64,
    /// Length of the shorter side.
    pub width: f64,
    /// Direction of the longer side, radians counter-clockwise from the x axis,
    /// in `[0, π)` — an axis has no head or tail, so the two opposite directions
    /// are the same orientation.
    pub angle: f64,
}

impl OrientedRect {
    /// As a closed polygon ring.
    #[must_use]
    pub fn to_polygon(self) -> Polygon {
        Polygon::new(LineString::new(self.corners.to_vec()), vec![])
    }

    /// `width / length` — 1 for a square, approaching 0 for a sliver. The scale-free
    /// measure of elongation the axis-aligned box cannot provide.
    #[must_use]
    pub fn aspect_ratio(self) -> f64 {
        if self.length <= 0.0 {
            1.0
        } else {
            self.width / self.length
        }
    }
}

/// Smallest enclosing circle, oriented envelope, and largest inscribed circle.
pub trait BoundingShapes {
    /// The smallest circle containing every coordinate, or `None` if there are
    /// none.
    fn minimum_bounding_circle(&self) -> Option<Circle>;

    /// The smallest-area enclosing rectangle at any orientation, or `None` for
    /// fewer than two distinct coordinates.
    fn oriented_envelope(&self) -> Option<OrientedRect>;

    /// The largest circle whose centre lies inside the geometry, found to
    /// `tolerance`. `None` for a geometry with no interior.
    fn maximum_inscribed_circle(&self, tolerance: f64) -> Option<Circle>;
}

/// Welzl's algorithm, iterative rather than recursive.
///
/// The textbook form recurses once per point, which on a 100 000-point input is
/// 100 000 stack frames — the same unbounded-recursion hazard the codecs had. The
/// move-to-front iterative variant below is the standard remedy and is what most
/// production implementations use.
fn welzl(points: &[Coord]) -> Option<Circle> {
    let mut pts: Vec<Coord> = points.iter().copied().filter(|c| c.is_finite()).collect();
    if pts.is_empty() {
        return None;
    }
    // A deterministic shuffle. Welzl's expected-linear bound needs the input in
    // random order, and real geometry arrives sorted — along a boundary, in scan
    // order — which is the worst case. A fixed LCG keeps the result reproducible
    // while breaking that adversarial ordering, so the same input always gives the
    // same circle.
    let mut seed = 0x2545_F491_4F6C_DD1Du64;
    for i in (1..pts.len()).rev() {
        seed = seed.wrapping_mul(6_364_136_223_846_793_005).wrapping_add(1);
        let j = (seed >> 33) as usize % (i + 1);
        pts.swap(i, j);
    }

    let mut circle = Circle {
        centre: pts[0],
        radius: 0.0,
    };
    for i in 0..pts.len() {
        if circle.contains(pts[i]) {
            continue;
        }
        circle = Circle {
            centre: pts[i],
            radius: 0.0,
        };
        for j in 0..i {
            if circle.contains(pts[j]) {
                continue;
            }
            circle = from_two(pts[i], pts[j]);
            for k in 0..j {
                if circle.contains(pts[k]) {
                    continue;
                }
                circle = from_three(pts[i], pts[j], pts[k]).unwrap_or(circle);
            }
        }
    }
    Some(circle)
}

fn from_two(a: Coord, b: Coord) -> Circle {
    let centre = Coord::new((a.x + b.x) / 2.0, (a.y + b.y) / 2.0);
    Circle {
        centre,
        radius: Euclidean.distance(Point(centre), Point(a)),
    }
}

/// The circumcircle of three points, or `None` when they are collinear and so
/// define no circle.
fn from_three(a: Coord, b: Coord, c: Coord) -> Option<Circle> {
    let d = 2.0 * (a.x * (b.y - c.y) + b.x * (c.y - a.y) + c.x * (a.y - b.y));
    if d.abs() < f64::MIN_POSITIVE {
        return None;
    }
    let (a2, b2, c2) = (
        a.x * a.x + a.y * a.y,
        b.x * b.x + b.y * b.y,
        c.x * c.x + c.y * c.y,
    );
    let centre = Coord::new(
        (a2 * (b.y - c.y) + b2 * (c.y - a.y) + c2 * (a.y - b.y)) / d,
        (a2 * (c.x - b.x) + b2 * (a.x - c.x) + c2 * (b.x - a.x)) / d,
    );
    if !centre.is_finite() {
        return None;
    }
    Some(Circle {
        centre,
        radius: Euclidean.distance(Point(centre), Point(a)),
    })
}

/// Rotating calipers over the convex hull.
///
/// Only hull edges are tried, because the minimum-area enclosing rectangle always
/// has one side flush with a hull edge. That is what reduces a continuous search
/// over all angles to a loop over a handful of them, and it is why this is cheap.
fn calipers(points: &[Coord]) -> Option<OrientedRect> {
    let hull = convex_hull_coords(points);
    // `convex_hull_coords` returns a closed ring, so the last point repeats.
    let ring: Vec<Coord> = if hull.len() >= 2 && hull[0] == hull[hull.len() - 1] {
        hull[..hull.len() - 1].to_vec()
    } else {
        hull
    };
    if ring.is_empty() {
        return None;
    }
    if ring.len() == 1 {
        let p = ring[0];
        return Some(OrientedRect {
            corners: [p, p, p, p],
            length: 0.0,
            width: 0.0,
            angle: 0.0,
        });
    }

    let mut best: Option<(f64, OrientedRect)> = None;
    for i in 0..ring.len() {
        let a = ring[i];
        let b = ring[(i + 1) % ring.len()];
        let (dx, dy) = (b.x - a.x, b.y - a.y);
        let len = (dx * dx + dy * dy).sqrt();
        if len <= 0.0 {
            continue;
        }
        // Unit vectors along and across this edge.
        let (ux, uy) = (dx / len, dy / len);
        let (vx, vy) = (-uy, ux);

        // Project every hull vertex onto the two axes; the extents are the
        // rectangle's sides.
        let (mut lo_u, mut hi_u, mut lo_v, mut hi_v) = (
            f64::INFINITY,
            f64::NEG_INFINITY,
            f64::INFINITY,
            f64::NEG_INFINITY,
        );
        for p in &ring {
            let pu = p.x * ux + p.y * uy;
            let pv = p.x * vx + p.y * vy;
            lo_u = lo_u.min(pu);
            hi_u = hi_u.max(pu);
            lo_v = lo_v.min(pv);
            hi_v = hi_v.max(pv);
        }
        let (su, sv) = (hi_u - lo_u, hi_v - lo_v);
        let area = su * sv;
        if best.as_ref().is_some_and(|(b, _)| area >= *b) {
            continue;
        }

        // Back from axis coordinates to the plane.
        let at = |u: f64, v: f64| Coord::new(u * ux + v * vx, u * uy + v * vy);
        let corners = [
            at(lo_u, lo_v),
            at(hi_u, lo_v),
            at(hi_u, hi_v),
            at(lo_u, hi_v),
        ];
        // The reported angle is that of the *longer* side, folded into `[0, π)`
        // because an axis has no direction — 10° and 190° describe the same
        // orientation, and a caller comparing orientations needs one of them.
        let along = if su >= sv { uy.atan2(ux) } else { vy.atan2(vx) };
        let angle = along.rem_euclid(std::f64::consts::PI);
        best = Some((
            area,
            OrientedRect {
                corners,
                length: su.max(sv),
                width: su.min(sv),
                angle,
            },
        ));
    }
    best.map(|(_, r)| r)
}

/// A square cell of the search, carrying the signed clearance at its centre.
///
/// Ordered by [`Cell::bound`] so a `BinaryHeap` yields the most promising cell
/// first — best-first search, not breadth-first, which is what keeps the explored
/// set small.
#[derive(Clone, Copy, Debug)]
struct Cell {
    centre: Coord,
    /// Half the side length, so the cell is `centre ± half` on each axis.
    half: f64,
    /// Clearance at `centre`, negated outside the geometry.
    value: f64,
}

impl Cell {
    /// The largest clearance any point in this cell could have.
    ///
    /// Clearance is distance to the boundary, so it is 1-Lipschitz: moving by `d`
    /// changes it by at most `d`. No point of the cell is further from the centre
    /// than the half-diagonal, so nothing inside can beat `value + half·√2`. That
    /// inequality is the whole basis of the pruning below — it is an upper bound,
    /// not an estimate, so discarding a cell that fails it cannot lose the optimum.
    fn bound(self) -> f64 {
        self.value + self.half * std::f64::consts::SQRT_2
    }
}

impl PartialEq for Cell {
    fn eq(&self, other: &Self) -> bool {
        self.bound() == other.bound()
    }
}
impl Eq for Cell {}
impl Ord for Cell {
    fn cmp(&self, other: &Self) -> Ordering {
        self.bound().total_cmp(&other.bound())
    }
}
impl PartialOrd for Cell {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

/// The largest inscribed circle, by Lipschitz branch and bound over the extent.
///
/// A deliberately simple search rather than the medial-axis construction JTS uses
/// for the exact answer. The honest trade is stated in the public docs: this
/// converges to `tolerance`, where the medial axis is exact immediately but needs a
/// Voronoi diagram of the boundary. Since the answer is used for label placement
/// and clearance checks, where a tolerance is acceptable and a dependency on a
/// whole extra algorithm is not, refinement is the better trade here.
///
/// **Why branch and bound rather than a coarse-to-fine sweep.** The obvious
/// structure — grid the extent, keep the best point, re-grid around *it* at half
/// the spacing — is wrong, and was what this did. Clearance is not unimodal: a
/// polygon with two wide lobes has a local maximum in each, and re-gridding around
/// the winner discards every basin but one. Once the coarse pass favours the wrong
/// lobe the search can never reach the right one, so the error is set by the
/// *initial* spacing and does not shrink when `tolerance` does. A dumbbell with
/// lobes of radius 1.0 and 1.1 returned 1.0 for every tolerance from `1e-3` to
/// `1e-9` — a shortfall of 0.098 against a requested 1e-9.
///
/// Branch and bound keeps the guarantee instead of assuming it. Every cell stays
/// live until its own [`Cell::bound`] proves it cannot hold anything better than
/// the incumbent, so no basin is ever dropped on the strength of one sample. The
/// search stops when the most promising remaining cell cannot beat the incumbent by
/// more than `tolerance`, which makes the returned radius within `tolerance` of the
/// true maximum rather than within one grid step of a local one.
fn inscribed(geom: &Geometry, tolerance: f64) -> Option<Circle> {
    let boundary = boundary_segments(geom);
    if boundary.is_empty() {
        return None;
    }
    let coords = geom.to_coords();
    let (mut lo_x, mut hi_x, mut lo_y, mut hi_y) = (
        f64::INFINITY,
        f64::NEG_INFINITY,
        f64::INFINITY,
        f64::NEG_INFINITY,
    );
    for c in &coords {
        lo_x = lo_x.min(c.x);
        hi_x = hi_x.max(c.x);
        lo_y = lo_y.min(c.y);
        hi_y = hi_y.max(c.y);
    }
    if !(lo_x.is_finite() && hi_x.is_finite() && lo_y.is_finite() && hi_y.is_finite()) {
        return None;
    }

    let clearance = |p: Coord| {
        boundary
            .iter()
            .map(|&(a, b)| segment_distance(p, a, b))
            .fold(f64::INFINITY, f64::min)
    };

    // Negated outside, which keeps the function 1-Lipschitz across the boundary so
    // one bound covers cells that straddle it. A cell outside the geometry is still
    // worth subdividing — it may overlap an interior region the centre missed —
    // but its own centre can never become the incumbent.
    let signed = |p: Coord| {
        let d = clearance(p);
        if geom.contains(&Point(p)) { d } else { -d }
    };

    let tol = tolerance.max(1e-12);
    let span = (hi_x - lo_x).max(hi_y - lo_y);
    if !span.is_finite() || span <= 0.0 {
        // Every coordinate at one spot: no interior to inscribe anything in. The
        // finiteness arm is unreachable given the bounds check above, and is kept
        // so the root cell below can never be built with a non-finite half-width.
        return None;
    }

    let root = Cell {
        centre: Coord::new(f64::midpoint(lo_x, hi_x), f64::midpoint(lo_y, hi_y)),
        half: span / 2.0,
        value: 0.0,
    };
    let root = Cell {
        value: signed(root.centre),
        ..root
    };
    let (mut best_centre, mut best_radius) = (root.centre, root.value);

    let mut heap = BinaryHeap::new();
    heap.push(root);

    // A node budget as well as the tolerance test, because the number of cells
    // within `tol` of optimal is a property of the geometry, not of the algorithm:
    // a shape with a long ridge of near-equal clearance has many. Exhausting the
    // budget degrades the guarantee to "best found", which is why it is generous
    // rather than tight — it is a backstop against a pathological input, not the
    // normal stopping condition.
    let mut budget = 200_000_u32;
    while let Some(cell) = heap.pop() {
        // The heap is ordered by `bound`, so if the most promising cell left cannot
        // beat the incumbent by more than `tol`, none of the others can either.
        if cell.bound() <= best_radius + tol {
            break;
        }
        if budget == 0 {
            break;
        }
        budget -= 1;

        let half = cell.half / 2.0;
        for (sx, sy) in [(-1.0, -1.0), (-1.0, 1.0), (1.0, -1.0), (1.0, 1.0)] {
            let centre = Coord::new(cell.centre.x + sx * half, cell.centre.y + sy * half);
            let value = signed(centre);
            if value > best_radius {
                best_radius = value;
                best_centre = centre;
            }
            heap.push(Cell {
                centre,
                half,
                value,
            });
        }
    }

    (best_radius > 0.0).then_some(Circle {
        centre: best_centre,
        radius: best_radius,
    })
}

/// Every boundary segment of an areal geometry, shells and holes alike.
///
/// Holes count: the largest circle in an annulus is limited by the inner ring as
/// much as the outer one, and omitting holes would place a label in a lake.
fn boundary_segments(geom: &Geometry) -> Vec<(Coord, Coord)> {
    let mut out = Vec::new();
    let mut push_ring = |ring: &LineString| {
        let v = &ring.0;
        for i in 0..v.len().saturating_sub(1) {
            out.push((v[i], v[i + 1]));
        }
        if v.len() >= 3 && v[0] != v[v.len() - 1] {
            out.push((v[v.len() - 1], v[0]));
        }
    };
    let mut push_poly = |p: &Polygon| {
        push_ring(p.exterior());
        for h in p.interiors() {
            push_ring(h);
        }
    };
    match geom {
        Geometry::Polygon(p) => push_poly(p),
        Geometry::MultiPolygon(mp) => mp.iter().for_each(&mut push_poly),
        Geometry::Rect(r) => push_poly(&r.to_polygon()),
        Geometry::Triangle(t) => push_poly(&t.to_polygon()),
        Geometry::GeometryCollection(gs) => {
            for g in gs {
                out.extend(boundary_segments(g));
            }
        }
        _ => {}
    }
    out
}

fn segment_distance(p: Coord, a: Coord, b: Coord) -> f64 {
    let (dx, dy) = (b.x - a.x, b.y - a.y);
    let len2 = dx * dx + dy * dy;
    if len2 <= 0.0 {
        return ((p.x - a.x).powi(2) + (p.y - a.y).powi(2)).sqrt();
    }
    let t = (((p.x - a.x) * dx + (p.y - a.y) * dy) / len2).clamp(0.0, 1.0);
    let (qx, qy) = (a.x + t * dx, a.y + t * dy);
    ((p.x - qx).powi(2) + (p.y - qy).powi(2)).sqrt()
}

impl<T: CoordsIter + ?Sized> BoundingShapes for T
where
    Self: AsGeometry,
{
    fn minimum_bounding_circle(&self) -> Option<Circle> {
        welzl(&self.to_coords())
    }

    fn oriented_envelope(&self) -> Option<OrientedRect> {
        calipers(&self.to_coords())
    }

    fn maximum_inscribed_circle(&self, tolerance: f64) -> Option<Circle> {
        inscribed(&self.as_geometry(), tolerance)
    }
}

/// Bridge so the blanket impl above can reach the areal structure a maximum
/// inscribed circle needs, which `CoordsIter` alone does not expose.
pub trait AsGeometry {
    /// This value as a `Geometry`.
    fn as_geometry(&self) -> Geometry;
}

impl AsGeometry for Geometry {
    fn as_geometry(&self) -> Geometry {
        self.clone()
    }
}

impl AsGeometry for Polygon {
    fn as_geometry(&self) -> Geometry {
        Geometry::Polygon(self.clone())
    }
}

impl AsGeometry for crate::space::r2::MultiPolygon {
    fn as_geometry(&self) -> Geometry {
        Geometry::MultiPolygon(self.clone())
    }
}

impl AsGeometry for LineString {
    fn as_geometry(&self) -> Geometry {
        Geometry::LineString(self.clone())
    }
}

impl AsGeometry for crate::space::r2::MultiPoint {
    fn as_geometry(&self) -> Geometry {
        Geometry::MultiPoint(self.clone())
    }
}

/// Area of an oriented rectangle, for the tests and for callers comparing it to
/// the geometry's own area.
#[must_use]
pub fn oriented_envelope_area(r: OrientedRect) -> f64 {
    r.length * r.width
}

/// The convex hull's area, exposed because "how much of its hull does this shape
/// fill" (solidity) is the other standard shape descriptor and is one division
/// away once the hull is available.
#[must_use]
pub fn hull_area(points: &[Coord]) -> f64 {
    let ring = convex_hull_coords(points);
    if ring.len() < 4 {
        return 0.0;
    }
    Polygon::new(LineString::new(ring), vec![]).unsigned_area()
}
