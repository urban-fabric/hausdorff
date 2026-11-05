use super::*;

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

    /// The **semi-latus rectum** `ℓ = a_minor² / a_major` — half the focal chord
    /// perpendicular to the major axis.
    #[must_use]
    pub fn semi_latus_rectum(self) -> f64 {
        let (major, minor) = (self.a.max(self.b), self.a.min(self.b));
        minor * minor / major
    }

    /// Distance from the centre to each **directrix** along the major axis,
    /// `a_major / e = a_major² / c` (`+∞` for a circle). A boundary point's
    /// focus-distance is `e ×` its directrix-distance — the focus–directrix
    /// definition, `riemann` `Hausdorff/Ellipse.ellipse_directrix`.
    #[must_use]
    pub fn directrix_offset(self) -> f64 {
        let major = self.a.max(self.b);
        major * major / self.focal_distance()
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

    /// The **tangent line** at eccentric angle `t`, as `[A, B, C]` for the line
    /// `A·x + B·y = C`. Its normal `(A, B) = (b cos t, a sin t)` is the surface
    /// normal direction, and the boundary point lies on it. Proof-backed by
    /// `riemann` `Hausdorff/Ellipse.{ellipse_on_tangent_line, ellipse_tangent_perp}`.
    #[must_use]
    pub fn tangent_line_at(self, t: f64) -> [f64; 3] {
        let (s, c) = t.sin_cos();
        let (na, nb) = (self.b * c, self.a * s); // line normal = ellipse normal dir
        let cc = self.a * self.b + na * self.center.x + nb * self.center.y;
        [na, nb, cc]
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

    #[cfg(test)]
    pub(super) fn covariance(self) -> (f64, f64, f64) {
        (self.sxx, self.sxy, self.syy)
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
        self.syy * dx * dx - 2.0 * self.sxy * dx * dy + self.sxx * dy * dy
            <= self.k * self.k * self.det()
    }

    /// The covariance determinant `sxx·syy − sxy²` (= `λ_max·λ_min`); the area is
    /// `π k² √det`. Proof-backed by `riemann`
    /// `Hausdorff/Ellipse.sym2_eig_prod_eq_det`.
    #[must_use]
    pub fn det(self) -> f64 {
        self.sxx * self.syy - self.sxy * self.sxy
    }

    /// The two drawable **principal semi-axis vectors** — the major axis
    /// `(semi_major·cosθ, semi_major·sinθ)` (`θ = rotation_rad`) and the minor
    /// axis perpendicular to it: the eigenvectors of `Σ` scaled by the semi-axes
    /// `k√λ`. Proof-backed by `riemann` `Hausdorff/Ellipse.sym2_eigenvalue`.
    #[must_use]
    pub fn axes(self) -> ([f64; 2], [f64; 2]) {
        let (s, c) = self.rotation_rad().sin_cos();
        let (maj, min) = (self.semi_major(), self.semi_minor());
        ([maj * c, maj * s], [-min * s, min * c])
    }
}
