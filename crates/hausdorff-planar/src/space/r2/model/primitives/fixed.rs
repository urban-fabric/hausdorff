use super::*;

/// A line segment from `start` to `end`.
#[derive(Copy, Clone, PartialEq, Debug, Serialize, Deserialize)]
pub struct Line {
    pub start: Coord,
    pub end: Coord,
}

impl Line {
    /// The segment from `start` to `end`.
    #[must_use]
    pub fn new(start: impl Into<Coord>, end: impl Into<Coord>) -> Self {
        Self {
            start: start.into(),
            end: end.into(),
        }
    }

    /// `end - start`.
    #[inline]
    #[must_use]
    pub fn delta(&self) -> Coord {
        self.end - self.start
    }

    /// Squared Euclidean length — cheaper than [`length`](Self::length) when you
    /// only need to compare.
    #[inline]
    #[must_use]
    pub fn length_squared(&self) -> f64 {
        let d = self.delta();
        d.x * d.x + d.y * d.y
    }

    /// Euclidean length.
    #[inline]
    #[must_use]
    pub fn length(&self) -> f64 {
        self.length_squared().sqrt()
    }

    /// The **perpendicular bisector** of the segment, as `[A, B, C]` for the line
    /// `A·x + B·y = C` — exactly the points equidistant from `start` and `end`.
    ///
    /// Equidistance is *linear*: `(end − start)·p = (|end|² − |start|²)/2` (the
    /// `|p|²` terms cancel), so the bisector's normal `(A, B)` is the chord
    /// `end − start` and it passes through the midpoint. Proof-backed by `riemann`
    /// `ExactArith/Circle.equidistant_iff` (and the inner-product-generic
    /// `Geometry/Metric.equidist_iff`); it is the linear condition the planar
    /// [`Triangle::circumcenter`] solves.
    #[must_use]
    pub fn perpendicular_bisector(&self) -> [f64; 3] {
        let (a, b) = (self.start, self.end);
        let na = b.x - a.x;
        let nb = b.y - a.y;
        let cc = 0.5 * ((b.x * b.x + b.y * b.y) - (a.x * a.x + a.y * a.y));
        [na, nb, cc]
    }
}

/// A triangle defined by its three vertices.
#[derive(Copy, Clone, PartialEq, Debug, Serialize, Deserialize)]
pub struct Triangle(pub Coord, pub Coord, pub Coord);

impl Triangle {
    /// The triangle with the given vertices.
    #[must_use]
    pub fn new(a: Coord, b: Coord, c: Coord) -> Self {
        Self(a, b, c)
    }

    /// The three vertices as an array.
    #[inline]
    #[must_use]
    pub fn to_array(self) -> [Coord; 3] {
        [self.0, self.1, self.2]
    }

    /// This triangle as a closed three-sided [`Polygon`].
    #[must_use]
    pub fn to_polygon(self) -> Polygon {
        Polygon::new(LineString::new(vec![self.0, self.1, self.2]), vec![])
    }

    /// The **circumcenter** — the point equidistant from all three vertices (the
    /// centre of the circumscribed circle) — or `None` when the vertices are
    /// collinear (zero signed area), so no circle passes through them.
    ///
    /// The three sides' perpendicular bisectors meet here; because equidistance
    /// is the *linear* condition `(b−a)·p = (|b|²−|a|²)/2` (the `|p|²` term
    /// cancels), the centre solves a 2×2 system, here in closed Cramer form with
    /// denominator `2·(signed area)`. Proof-backed by `riemann`
    /// `ExactArith/Circle.circumcenter_exists` (a non-degenerate triangle has an
    /// equidistant centre) via `equidistant_iff` — the R2 companion of
    /// `Tetrahedron::circumcenter`'s `circumsphere_exists`, one dimension down.
    #[must_use]
    pub fn circumcenter(self) -> Option<Coord> {
        let Triangle(a, b, c) = self;
        let den = 2.0 * ((b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x));
        if den == 0.0 {
            return None; // collinear — no circumcircle
        }
        let (da, db, dc) = (
            a.x * a.x + a.y * a.y,
            b.x * b.x + b.y * b.y,
            c.x * c.x + c.y * c.y,
        );
        let ux = ((db - da) * (c.y - a.y) - (dc - da) * (b.y - a.y)) / den;
        let uy = ((b.x - a.x) * (dc - da) - (c.x - a.x) * (db - da)) / den;
        Some(Coord::new(ux, uy))
    }

    /// The **circumcircle** as `(centre, radius)`, or `None` if the vertices are
    /// collinear. The radius is the centre's distance to any vertex (equal to all
    /// three by `circumcenter_exists`).
    #[must_use]
    pub fn circumcircle(self) -> Option<(Coord, f64)> {
        let Triangle(a, _, _) = self;
        self.circumcenter()
            .map(|o| (o, (o.x - a.x).hypot(o.y - a.y)))
    }
}

/// An axis-aligned rectangle, stored by its `min` (lowest x, y) and `max`
/// (highest x, y) corners. [`new`](Self::new) accepts the corners in any order
/// and normalizes them.
#[derive(Copy, Clone, PartialEq, Debug, Serialize, Deserialize)]
pub struct Rect {
    min: Coord,
    max: Coord,
}

impl Rect {
    /// The rectangle spanning two opposite corners, given in any order.
    #[must_use]
    pub fn new(c1: impl Into<Coord>, c2: impl Into<Coord>) -> Self {
        let (a, b) = (c1.into(), c2.into());
        Self {
            min: Coord::new(a.x.min(b.x), a.y.min(b.y)),
            max: Coord::new(a.x.max(b.x), a.y.max(b.y)),
        }
    }

    /// The lower corner (minimum x and y).
    #[inline]
    #[must_use]
    pub fn min(self) -> Coord {
        self.min
    }

    /// The upper corner (maximum x and y).
    #[inline]
    #[must_use]
    pub fn max(self) -> Coord {
        self.max
    }

    /// Width (extent along x).
    #[inline]
    #[must_use]
    pub fn width(self) -> f64 {
        self.max.x - self.min.x
    }

    /// Height (extent along y).
    #[inline]
    #[must_use]
    pub fn height(self) -> f64 {
        self.max.y - self.min.y
    }

    /// The center point.
    #[inline]
    #[must_use]
    pub fn center(self) -> Point {
        Point::new(
            (self.min.x + self.max.x) * 0.5,
            (self.min.y + self.max.y) * 0.5,
        )
    }

    /// True iff `coord` lies within the rectangle (inclusive).
    #[inline]
    #[must_use]
    pub fn contains_coord(self, coord: Coord) -> bool {
        coord.x >= self.min.x
            && coord.x <= self.max.x
            && coord.y >= self.min.y
            && coord.y <= self.max.y
    }

    /// True iff this rectangle and `other` overlap — share any point, edges
    /// inclusive — by the per-axis interval-overlap test. A non-overlapping pair
    /// shares no point, so a spatial index may prune it: proof-backed by `riemann`
    /// `ExactArith/AABB.overlaps_of_mem` (a shared point forces an overlap) and
    /// `prune_sound` (the R-tree pruning soundness it underwrites).
    #[inline]
    #[must_use]
    pub fn overlaps(self, other: Rect) -> bool {
        self.min.x <= other.max.x
            && other.min.x <= self.max.x
            && self.min.y <= other.max.y
            && other.min.y <= self.max.y
    }

    /// This rectangle as a closed, counter-clockwise [`Polygon`].
    #[must_use]
    pub fn to_polygon(self) -> Polygon {
        Polygon::new(
            LineString::new(vec![
                Coord::new(self.min.x, self.min.y),
                Coord::new(self.max.x, self.min.y),
                Coord::new(self.max.x, self.max.y),
                Coord::new(self.min.x, self.max.y),
            ]),
            vec![],
        )
    }
}
