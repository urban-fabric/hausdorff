//! `s2::Polygon` — a spherical polygon on the unit sphere: a counter-clockwise
//! shell [`Loop`] with zero or more hole loops removed. The spherical sibling of
//! the planar [`hausdorff_planar::space::r2::Polygon`]; area is in **steradians** and there is no
//! projection.
//!
//! Convention: the shell and every hole are wound counter-clockwise (interior to
//! the left), so each hole loop encloses exactly the region it removes.

use hausdorff_base::space::s1::Angle;

use super::{Loop, Rect, S2Point};

/// A spherical polygon: a shell with holes.
#[derive(Clone, Debug, PartialEq)]
pub struct Polygon {
    shell: Loop,
    holes: Vec<Loop>,
}

impl Polygon {
    /// A polygon bounded by `shell` with `holes` removed (all CCW).
    #[must_use]
    pub fn new(shell: Loop, holes: Vec<Loop>) -> Self {
        Self { shell, holes }
    }

    /// A polygon with no holes.
    #[must_use]
    pub fn from_shell(shell: Loop) -> Self {
        Self {
            shell,
            holes: Vec::new(),
        }
    }

    /// The outer boundary.
    #[must_use]
    pub fn shell(&self) -> &Loop {
        &self.shell
    }

    /// The hole loops.
    #[must_use]
    pub fn holes(&self) -> &[Loop] {
        &self.holes
    }

    /// Area of the polygon interior in **steradians**: the shell area minus the
    /// hole areas (clamped at 0).
    #[must_use]
    pub fn area(&self) -> f64 {
        let holes: f64 = self.holes.iter().map(Loop::area).sum();
        (self.shell.area() - holes).max(0.0)
    }

    /// True iff `p` is inside the shell and outside every hole.
    #[must_use]
    pub fn contains_point(&self, p: S2Point) -> bool {
        self.shell.contains_point(p) && !self.holes.iter().any(|h| h.contains_point(p))
    }

    /// Total boundary length (shell + holes) as an [`Angle`]; multiply by a
    /// planet radius (see [`hausdorff_base::earth::sphere::Sphere`]) for metres.
    #[must_use]
    pub fn perimeter(&self) -> Angle {
        let total = self.shell.perimeter().radians()
            + self
                .holes
                .iter()
                .map(|h| h.perimeter().radians())
                .sum::<f64>();
        Angle::from_radians(total)
    }

    /// A conservative lat/lng rectangle containing the polygon.
    ///
    /// Holes remove area, so the shell's rectangle already bounds the polygon.
    #[must_use]
    pub fn rect_bound(&self) -> Rect {
        self.shell.rect_bound()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use hausdorff_base::earth::GeoPoint;

    fn sp(lat: f64, lon: f64) -> S2Point {
        S2Point::from_geo_point(GeoPoint::from_lat_lon(lat, lon))
    }

    #[test]
    fn rect_bound_is_shell_bound() {
        let shell = Loop::new(vec![
            sp(0.0, 0.0),
            sp(0.0, 10.0),
            sp(10.0, 10.0),
            sp(10.0, 0.0),
        ]);
        let polygon = Polygon::from_shell(shell.clone());
        assert_eq!(polygon.rect_bound(), shell.rect_bound());
    }
}
