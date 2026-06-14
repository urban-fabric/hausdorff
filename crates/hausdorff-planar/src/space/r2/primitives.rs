//! Small fixed-vertex planar geometries: [`Line`] (a segment), [`Triangle`],
//! and the axis-aligned [`Rect`]. Each can expand to a [`Polygon`] for the
//! algorithms that only speak in rings.

use serde::{Deserialize, Serialize};

use super::{Coord, LineString, Point, Polygon};

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

/// An axis-aligned **ellipse** region: every point with
/// `(x − cx)²/a² + (y − cy)²/b² ≤ 1`, centred at `(cx, cy)` with semi-axes `a`
/// (along x) and `b` (along y). The planar bounding-region mirror of the
/// spherical `S2Cap`, and the curved generalization of [`Rect`] for radial
/// containment / spatial pre-filtering.
///
/// Membership ([`contains_coord`](Self::contains_coord)) and convexity are
/// proof-backed in `riemann` `Hausdorff/Ellipse` (`inEllipse`,
/// `inEllipse_convex`): the test is the division-free form
/// `b²(x−cx)² + a²(y−cy)² ≤ a²b²`, and the region is convex (so a segment
/// between two contained points stays contained).
#[derive(Copy, Clone, PartialEq, Debug, Serialize, Deserialize)]
pub struct Ellipse {
    center: Coord,
    a: f64,
    b: f64,
}

impl Ellipse {
    /// An ellipse centred at `center` with semi-axes `a` (along x) and `b`
    /// (along y). A circle of radius `r` is `Ellipse::new(c, r, r)`.
    #[must_use]
    pub fn new(center: impl Into<Coord>, a: f64, b: f64) -> Self {
        Self {
            center: center.into(),
            a,
            b,
        }
    }

    /// The centre point.
    #[inline]
    #[must_use]
    pub fn center(self) -> Coord {
        self.center
    }

    /// Semi-axis along x.
    #[inline]
    #[must_use]
    pub fn semi_x(self) -> f64 {
        self.a
    }

    /// Semi-axis along y.
    #[inline]
    #[must_use]
    pub fn semi_y(self) -> f64 {
        self.b
    }

    /// Area `π·a·b`.
    #[inline]
    #[must_use]
    pub fn area(self) -> f64 {
        std::f64::consts::PI * self.a * self.b
    }

    /// True iff `coord` lies within the ellipse (boundary inclusive).
    ///
    /// The division-free membership test `b²(x−cx)² + a²(y−cy)² ≤ a²b²`,
    /// proof-backed by `riemann` `Hausdorff/Ellipse.inEllipse`; the region it
    /// defines is convex (`inEllipse_convex`).
    #[inline]
    #[must_use]
    pub fn contains_coord(self, coord: Coord) -> bool {
        let dx = coord.x - self.center.x;
        let dy = coord.y - self.center.y;
        let (a2, b2) = (self.a * self.a, self.b * self.b);
        b2 * dx * dx + a2 * dy * dy <= a2 * b2
    }

    /// The tight axis-aligned [`Rect`] bounding this ellipse.
    #[must_use]
    pub fn bounding_rect(self) -> Rect {
        Rect::new(
            Coord::new(self.center.x - self.a, self.center.y - self.b),
            Coord::new(self.center.x + self.a, self.center.y + self.b),
        )
    }

    /// The linear eccentricity `c = √(a_major² − a_minor²)` — the focus offset
    /// from the centre (`0` for a circle).
    #[must_use]
    pub fn focal_distance(self) -> f64 {
        let (major, minor) = (self.a.max(self.b), self.a.min(self.b));
        (major * major - minor * minor).sqrt()
    }

    /// Eccentricity `e = c / a_major ∈ [0, 1)` (`0` for a circle).
    #[must_use]
    pub fn eccentricity(self) -> f64 {
        self.focal_distance() / self.a.max(self.b)
    }

    /// The two foci, on the major axis at `centre ± c` (coincident at the centre
    /// for a circle).
    #[must_use]
    pub fn foci(self) -> (Coord, Coord) {
        let c = self.focal_distance();
        if self.a >= self.b {
            (
                Coord::new(self.center.x - c, self.center.y),
                Coord::new(self.center.x + c, self.center.y),
            )
        } else {
            (
                Coord::new(self.center.x, self.center.y - c),
                Coord::new(self.center.x, self.center.y + c),
            )
        }
    }

    /// The boundary point at eccentric angle `t`: `centre + (a cos t, b sin t)`.
    #[must_use]
    pub fn point_at_angle(self, t: f64) -> Coord {
        let (s, c) = t.sin_cos();
        Coord::new(self.center.x + self.a * c, self.center.y + self.b * s)
    }

    /// The outward unit **normal** at eccentric angle `t` — direction
    /// `(b cos t, a sin t)` normalized. It makes equal angles with the two focal
    /// radii (the reflection / whispering-gallery property), proof-backed by
    /// `riemann` `Hausdorff/Ellipse.ellipse_reflection`.
    #[must_use]
    pub fn normal_at(self, t: f64) -> [f64; 2] {
        let (s, c) = t.sin_cos();
        let (nx, ny) = (self.b * c, self.a * s);
        let len = nx.hypot(ny);
        [nx / len, ny / len]
    }

    /// The unit **tangent** at eccentric angle `t` — direction `(−a sin t,
    /// b cos t)` normalized (the parametric velocity, ⟂ the normal).
    #[must_use]
    pub fn tangent_at(self, t: f64) -> [f64; 2] {
        let (s, c) = t.sin_cos();
        let (tx, ty) = (-self.a * s, self.b * c);
        let len = tx.hypot(ty);
        [tx / len, ty / len]
    }

    /// This ellipse as a closed, counter-clockwise inscribed **`n`-gon** — `n`
    /// vertices on the boundary at equal eccentric-angle steps (`n` clamped to
    /// ≥ 3). Every vertex lies on the ellipse and the polygon ⊆ the (convex)
    /// ellipse, so it is a faithful polygonal stand-in that flows into the
    /// crate's polygon algorithms (area, clipping, predicates, indexing).
    ///
    /// Proof-backed by `riemann` `Hausdorff/Ellipse.ellipse_boundary` (vertices
    /// on the curve) and `ellipse_chord_inside` (each edge stays inside).
    #[must_use]
    pub fn to_polygon(self, n: usize) -> Polygon {
        let n = n.max(3);
        let coords: Vec<Coord> = (0..n)
            .map(|i| self.point_at_angle(std::f64::consts::TAU * (i as f64) / (n as f64)))
            .collect();
        Polygon::new(LineString::new(coords), vec![])
    }
}

/// A **confidence / error ellipse** built from a 2×2 covariance: the Mahalanobis
/// region `(p − μ)ᵀ Σ⁻¹ (p − μ) ≤ k²`, the standard way to draw positional
/// uncertainty (GPS error ellipses, Kalman covariance, scatter summaries).
///
/// `k` is the scale (e.g. `k = 1` for the 1-σ ellipse, `k = √5.991 ≈ 2.448` for
/// the 95% region in 2-D). The principal axes are the eigenvectors of `Σ`, with
/// semi-axes `k√λ`; this is a *rotated* ellipse (unlike axis-aligned [`Ellipse`]).
///
/// Convexity of the region is proof-backed by `riemann`
/// `Hausdorff/Ellipse.inQuad_convex` (a positive-semidefinite quadratic-form
/// sublevel set is convex — true for any covariance).
#[derive(Copy, Clone, PartialEq, Debug, Serialize, Deserialize)]
pub struct ConfidenceEllipse {
    center: Coord,
    sxx: f64,
    sxy: f64,
    syy: f64,
    k: f64,
}

impl ConfidenceEllipse {
    /// The `k`-scaled confidence ellipse of the covariance `[[sxx, sxy], [sxy,
    /// syy]]` centred at `center`.
    #[must_use]
    pub fn from_covariance(center: impl Into<Coord>, sxx: f64, sxy: f64, syy: f64, k: f64) -> Self {
        Self {
            center: center.into(),
            sxx,
            sxy,
            syy,
            k,
        }
    }

    /// The covariance eigenvalues `(λ_max, λ_min)` (variances along the principal
    /// axes), in closed form for the symmetric 2×2.
    #[must_use]
    pub fn eigenvalues(self) -> (f64, f64) {
        let mean = 0.5 * (self.sxx + self.syy);
        let half_diff = 0.5 * (self.sxx - self.syy);
        let r = half_diff.hypot(self.sxy);
        (mean + r, mean - r)
    }

    /// Semi-axis along the major direction, `k√λ_max`.
    #[must_use]
    pub fn semi_major(self) -> f64 {
        self.k * self.eigenvalues().0.max(0.0).sqrt()
    }

    /// Semi-axis along the minor direction, `k√λ_min`.
    #[must_use]
    pub fn semi_minor(self) -> f64 {
        self.k * self.eigenvalues().1.max(0.0).sqrt()
    }

    /// Orientation of the major axis (radians), `½·atan2(2 sxy, sxx − syy)`.
    #[must_use]
    pub fn rotation_rad(self) -> f64 {
        0.5 * (2.0 * self.sxy).atan2(self.sxx - self.syy)
    }

    /// Area `π · semi_major · semi_minor = π k² √(det Σ)`.
    #[must_use]
    pub fn area(self) -> f64 {
        std::f64::consts::PI * self.semi_major() * self.semi_minor()
    }

    /// True iff `coord` is inside the ellipse — the Mahalanobis test cleared by
    /// `det Σ > 0`: `syy·dx² − 2 sxy·dx·dy + sxx·dy² ≤ k²·det Σ` (the convex
    /// quadratic-form region of `riemann` `Hausdorff/Ellipse.inQuad_convex`).
    #[inline]
    #[must_use]
    pub fn contains_coord(self, coord: Coord) -> bool {
        let dx = coord.x - self.center.x;
        let dy = coord.y - self.center.y;
        let det = self.sxx * self.syy - self.sxy * self.sxy;
        self.syy * dx * dx - 2.0 * self.sxy * dx * dy + self.sxx * dy * dy <= self.k * self.k * det
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ellipse_membership_boundary_and_convexity() {
        let e = Ellipse::new(Coord::new(1.0, 2.0), 3.0, 2.0);
        // Centre inside; the two semi-axis tips lie on the (inclusive) boundary.
        assert!(e.contains_coord(Coord::new(1.0, 2.0)));
        assert!(e.contains_coord(Coord::new(4.0, 2.0)));
        assert!(e.contains_coord(Coord::new(1.0, 4.0)));
        // Just outside the boundary is excluded.
        assert!(!e.contains_coord(Coord::new(4.001, 2.0)));
        assert!(!e.contains_coord(Coord::new(1.0, 4.001)));
        // Convexity (inEllipse_convex): the midpoint of two contained points is
        // contained.
        let mid = Coord::new((4.0 + 1.0) * 0.5, (2.0 + 4.0) * 0.5);
        assert!(e.contains_coord(mid));
    }

    #[test]
    fn circle_is_a_degenerate_ellipse() {
        let c = Ellipse::new(Coord::new(0.0, 0.0), 5.0, 5.0);
        assert!(c.contains_coord(Coord::new(3.0, 4.0))); // on the radius-5 circle
        assert!(!c.contains_coord(Coord::new(3.0, 4.1)));
    }

    #[test]
    fn ellipse_area_and_bounds() {
        let e = Ellipse::new(Coord::new(1.0, 2.0), 3.0, 2.0);
        assert!((e.area() - std::f64::consts::PI * 6.0).abs() < 1e-12);
        let r = e.bounding_rect();
        assert_eq!(r.min(), Coord::new(-2.0, 0.0));
        assert_eq!(r.max(), Coord::new(4.0, 4.0));
    }

    #[test]
    fn confidence_ellipse_axis_aligned() {
        // Diagonal covariance → axis-aligned; eigenvalues are the variances.
        let e = ConfidenceEllipse::from_covariance(Coord::new(0.0, 0.0), 4.0, 0.0, 1.0, 1.0);
        let (lmax, lmin) = e.eigenvalues();
        assert!((lmax - 4.0).abs() < 1e-12 && (lmin - 1.0).abs() < 1e-12);
        assert!((e.semi_major() - 2.0).abs() < 1e-12);
        assert!((e.semi_minor() - 1.0).abs() < 1e-12);
        assert!(e.rotation_rad().abs() < 1e-12);
        assert!((e.area() - std::f64::consts::PI * 2.0).abs() < 1e-12);
        // Boundary at the major-axis tip (k·σ_x = 2 along x); centre inside.
        assert!(e.contains_coord(Coord::new(2.0, 0.0)));
        assert!(!e.contains_coord(Coord::new(2.001, 0.0)));
        assert!(e.contains_coord(Coord::new(0.0, 0.0)));
    }

    #[test]
    fn confidence_ellipse_correlated_is_rotated() {
        // Equicorrelated covariance → 45° major axis, eigenvalues 3 and 1.
        let e = ConfidenceEllipse::from_covariance(Coord::new(0.0, 0.0), 2.0, 1.0, 2.0, 1.0);
        let (lmax, lmin) = e.eigenvalues();
        assert!((lmax - 3.0).abs() < 1e-12 && (lmin - 1.0).abs() < 1e-12);
        assert!((e.rotation_rad() - std::f64::consts::FRAC_PI_4).abs() < 1e-12);
        // Area = π·k²·√(det Σ), det = 3.
        assert!((e.area() - std::f64::consts::PI * 3.0_f64.sqrt()).abs() < 1e-12);
    }

    #[test]
    fn ellipse_foci_and_eccentricity() {
        let e = Ellipse::new(Coord::new(0.0, 0.0), 3.0, 2.0); // major along x, c = √5
        let c = 5.0_f64.sqrt();
        let (f1, f2) = e.foci();
        assert!((f1.x + c).abs() < 1e-12 && f1.y.abs() < 1e-12);
        assert!((f2.x - c).abs() < 1e-12 && f2.y.abs() < 1e-12);
        assert!((e.eccentricity() - c / 3.0).abs() < 1e-12);
        // A circle has zero eccentricity and coincident foci at the centre.
        let circle = Ellipse::new(Coord::new(1.0, 1.0), 5.0, 5.0);
        assert!(circle.eccentricity().abs() < 1e-12);
        let (g1, g2) = circle.foci();
        assert_eq!(g1, Coord::new(1.0, 1.0));
        assert_eq!(g2, Coord::new(1.0, 1.0));
    }

    #[test]
    fn ellipse_reflection_property() {
        // The normal makes equal angles with the two focal radii: a ray from one
        // focus reflects to the other (riemann Hausdorff/Ellipse.ellipse_reflection).
        let e = Ellipse::new(Coord::new(0.0, 0.0), 3.0, 2.0);
        let (f1, f2) = e.foci();
        for &t in &[0.3, 1.0, 2.5, 4.0, 5.5] {
            let p = e.point_at_angle(t);
            let n = e.normal_at(t);
            let u1 = [p.x - f1.x, p.y - f1.y];
            let u2 = [p.x - f2.x, p.y - f2.y];
            let cos1 = (n[0] * u1[0] + n[1] * u1[1]) / u1[0].hypot(u1[1]);
            let cos2 = (n[0] * u2[0] + n[1] * u2[1]) / u2[0].hypot(u2[1]);
            assert!((cos1 - cos2).abs() < 1e-12, "t={t}: {cos1} vs {cos2}");
            // Tangent ⟂ normal.
            let tan = e.tangent_at(t);
            assert!((n[0] * tan[0] + n[1] * tan[1]).abs() < 1e-12);
        }
    }

    #[test]
    fn ellipse_to_polygon_inscribes_and_approximates_area() {
        use crate::space::r2::Area;
        let e = Ellipse::new(Coord::new(1.0, -2.0), 3.0, 2.0);
        // The inscribed n-gon's area converges to π·a·b = 6π as n grows.
        let exact = std::f64::consts::PI * 6.0;
        let area = e.to_polygon(4000).unsigned_area();
        assert!((area - exact).abs() / exact < 1e-5, "area {area} vs {exact}");
        // `n` is clamped to ≥ 3 — no panic, still a valid positive-area polygon.
        assert!(e.to_polygon(0).unsigned_area() > 0.0);
    }
}
