//! `impl Region for Rect` — the cell-coverage / containment surface that lets a
//! [`crate::space::s2::RegionCoverer`] cover a rectangle with cells, ported from
//! golang/geo's `Rect.CapBound` and the cell-intersection logic.

use std::f64::consts::{FRAC_PI_2, PI};

use hausdorff_base::earth::GeoPoint;
use hausdorff_base::space::r3::Vector;
use hausdorff_base::space::s1::Angle;

use super::{Rect, Region};
use crate::space::s2::{Cap, S2CellId, S2Point};

impl Region for Rect {
    fn contains_point(&self, p: S2Point) -> bool {
        self.contains_lat_lng(p.to_geo_point())
    }

    /// Port of golang/geo's `Rect.CapBound`: the smaller of a polar cap (around
    /// the nearer pole, reaching the far latitude edge) and — when the
    /// longitude span is at most a half turn — a cap centered on the rectangle
    /// midpoint that reaches all four corners.
    fn bounding_cap(&self) -> Cap {
        if self.is_empty() {
            return Cap::empty(S2Point::from_vector(Vector::new(1.0, 0.0, 0.0)));
        }
        // Padding absorbs the cos() round trip in Cap::contains_point so a
        // boundary point computed at exactly the cap radius stays inside.
        let pad = 1e-14;
        let (pole_z, pole_angle) = if self.lat.lo + self.lat.hi < 0.0 {
            (-1.0, FRAC_PI_2 + self.lat.hi)
        } else {
            (1.0, FRAC_PI_2 - self.lat.lo)
        };
        let pole_cap = Cap::from_center_angle(
            S2Point::from_vector(Vector::new(0.0, 0.0, pole_z)),
            Angle::from_radians(pole_angle + pad),
        );
        if self.lng.length() <= PI {
            let mid_lat = 0.5 * (self.lat.lo + self.lat.hi);
            let center = S2Point::from_geo_point(GeoPoint::from_lat_lon(
                mid_lat.to_degrees(),
                self.lng.center().to_degrees(),
            ));
            let radius = self
                .corner_points()
                .iter()
                .map(|c| center.vector().angle(c.vector()).radians())
                .fold(0.0_f64, f64::max);
            let mid_cap = Cap::from_center_angle(center, Angle::from_radians(radius + pad));
            if mid_cap.height() < pole_cap.height() {
                return mid_cap;
            }
        }
        pole_cap
    }

    fn bounding_rect(&self) -> Rect {
        *self
    }

    fn intersects_cell(&self, cell: S2CellId) -> bool {
        if self.is_empty() {
            return false;
        }
        !self.intersection(&Rect::cell_bound(cell)).is_empty()
    }

    fn contains_cell(&self, cell: S2CellId) -> bool {
        self.contains_rect(&Rect::cell_bound(cell))
    }
}
