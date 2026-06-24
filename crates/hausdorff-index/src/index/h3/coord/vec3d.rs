use super::super::face;
use super::{
    AP7_ROT_RADS, EPSILON, INV_SQRT7_POWERS, ONE_THIRD, RES0_U_GNOMONIC, Vec2d, is_class3,
    to_positive_angle,
};

#[derive(Debug, Clone, Copy, PartialEq)]
pub(in crate::index::h3) struct Vec3d {
    pub(super) x: f64,
    pub(super) y: f64,
    pub(super) z: f64,
}

const NORTH_POLE: Vec3d = Vec3d {
    x: 0.0,
    y: 0.0,
    z: 1.0,
};

impl Vec3d {
    pub(in crate::index::h3) const fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z }
    }

    pub(super) fn norm(&self) -> f64 {
        self.norm_squared().sqrt()
    }

    /// n-vector of a cell center, given its 2D hex coords on `face`.
    /// `res` is the (extended) resolution index; `is_substrate` selects the
    /// substrate grid used for boundary vertices.
    pub(in crate::index::h3) fn from_vec2d(
        value: Vec2d,
        face: usize,
        res: usize,
        is_substrate: bool,
    ) -> Self {
        let r = {
            let mut r = value.magnitude();
            if r < EPSILON {
                return face::CENTER_POINT[face];
            }
            r *= INV_SQRT7_POWERS[res];
            if is_substrate {
                r *= ONE_THIRD;
            }
            (r * RES0_U_GNOMONIC).atan()
        };
        let theta = {
            let mut theta = value.y.atan2(value.x);
            if !is_substrate && is_class3(res) {
                theta = to_positive_angle(theta + AP7_ROT_RADS);
            }
            to_positive_angle(face::AXES_AZ_RADS_CII[face] - theta)
        };
        let center = face::CENTER_POINT[face];
        if r < EPSILON {
            return center;
        }
        let (north, east) = center.tangent_basis();
        let dir = linear_combination(theta.cos(), &north, theta.sin(), &east);
        let mut res_v = linear_combination(r.cos(), &center, r.sin(), &dir);
        res_v.normalize();
        res_v
    }

    /// Closest icosahedron face and the squared euclidean distance to it.
    pub(in crate::index::h3) fn closest_face(self) -> (usize, f64) {
        const MAX_DIST: f64 = 5.0;
        let mut best = (0usize, MAX_DIST);
        for (i, center) in face::CENTER_POINT.iter().enumerate() {
            let d = self.distance_squared(center);
            if d < best.1 {
                best = (i, d);
            }
        }
        best
    }

    /// Azimuth (radians) from `self` to `other`.
    pub(super) fn azimuth(&self, other: &Self) -> f64 {
        let (north, east) = self.tangent_basis();
        let mut proj = linear_combination(1.0, other, -other.dot(self), self);
        proj.normalize();
        proj.dot(&east).atan2(proj.dot(&north))
    }

    fn tangent_basis(&self) -> (Self, Self) {
        let mut north = linear_combination(1.0, &NORTH_POLE, -NORTH_POLE.dot(self), self);
        north.normalize();
        let east = north.cross(self);
        (north, east)
    }

    fn distance_squared(&self, other: &Self) -> f64 {
        linear_combination(1.0, self, -1.0, other).norm_squared()
    }

    fn normalize(&mut self) {
        let norm = self.norm();
        if norm > 0.0 {
            let s = 1.0 / norm;
            self.x *= s;
            self.y *= s;
            self.z *= s;
        } else {
            self.x = 0.0;
            self.y = 0.0;
            self.z = 0.0;
        }
    }

    fn cross(&self, o: &Self) -> Self {
        Self::new(
            self.y.mul_add(o.z, -(self.z * o.y)),
            self.z.mul_add(o.x, -(self.x * o.z)),
            self.x.mul_add(o.y, -(self.y * o.x)),
        )
    }

    fn dot(&self, o: &Self) -> f64 {
        self.x.mul_add(o.x, self.y.mul_add(o.y, self.z * o.z))
    }

    fn norm_squared(&self) -> f64 {
        self.dot(self)
    }
}

#[inline]
fn linear_combination(a: f64, v1: &Vec3d, b: f64, v2: &Vec3d) -> Vec3d {
    Vec3d::new(
        a.mul_add(v1.x, b * v2.x),
        a.mul_add(v1.y, b * v2.y),
        a.mul_add(v1.z, b * v2.z),
    )
}

impl Vec3d {
    /// n-vector from geographic coordinates in radians.
    pub(in crate::index::h3) fn from_lat_lng_rad(lat: f64, lng: f64) -> Self {
        let r = lat.cos();
        Self::new(lng.cos() * r, lng.sin() * r, lat.sin())
    }

    /// Geographic coordinates (radians) of this n-vector: `(lat, lng)`.
    /// Matches h3o's `LatLng::from(Vec3d)` exactly: `asin(z)` / `atan2(y, x)`.
    pub(in crate::index::h3) fn to_lat_lng_rad(self) -> (f64, f64) {
        (self.z.asin(), self.y.atan2(self.x))
    }
}
