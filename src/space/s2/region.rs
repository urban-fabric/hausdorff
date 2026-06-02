//! `s2::Region` — the common interface over spherical regions ([`Cap`], [`Loop`],
//! [`Polygon`]), mirroring golang/geo's `Region`. A region tests point
//! containment, reports a bounding [`Cap`], and answers the cell-level
//! intersect/contains queries that [`RegionCoverer`] uses to approximate the
//! region by a set of S2 cells.

use crate::space::r3::Vector;
use crate::space::s1::Angle;

use super::{Cap, Crossing, Loop, Polygon, S2CellId, S2Point, crossing_sign};

/// The deepest S2 subdivision level (a leaf cell).
const MAX_LEVEL: u32 = 30;

/// A region on the unit sphere.
pub trait Region {
    /// True iff the region contains `p`.
    fn contains_point(&self, p: S2Point) -> bool;

    /// A spherical cap that bounds the region (not necessarily the tightest).
    fn bounding_cap(&self) -> Cap;

    /// True iff the region may intersect `cell`. The coverer biases toward
    /// *over*-covering — it would rather keep a cell than drop one — but this is
    /// a heuristic, not a strict guarantee: it rests on point-in-loop and
    /// edge-crossing predicates whose results on exact boundary coincidences are
    /// not guaranteed, so a cell meeting the region only along a shared boundary
    /// can occasionally be dropped. Robust for typical float lon/lat data away
    /// from such coincidences.
    fn intersects_cell(&self, cell: S2CellId) -> bool;

    /// True iff the region fully contains `cell`. **Strict:** a `false` result
    /// only means "not known to be fully inside", so a coverer keeps subdividing
    /// rather than wrongly emitting a coarse cell that pokes outside the region.
    fn contains_cell(&self, cell: S2CellId) -> bool;
}

// ---- region ↔ cell geometry helpers ----------------------------------------

/// The cell as a counter-clockwise spherical [`Loop`] (interior to the left). The
/// cube→sphere projection can wind a cell either way depending on the face, so we
/// orient by the cell center: a CCW loop contains its own center.
fn cell_loop(cell: S2CellId) -> Loop {
    let verts = cell.vertices();
    let l = Loop::new(verts.to_vec());
    if l.contains_point(cell.center_point()) {
        l
    } else {
        let mut rev = verts.to_vec();
        rev.reverse();
        Loop::new(rev)
    }
}

/// The directed edges of a loop `(vᵢ, vᵢ₊₁)`, wrapping at the end.
fn loop_edges(l: &Loop) -> impl Iterator<Item = (S2Point, S2Point)> + '_ {
    let v = l.vertices();
    let n = v.len();
    (0..n).map(move |i| (v[i], v[(i + 1) % n]))
}

/// True iff any edge of `a` crosses any edge of `b` at a point interior to both
/// (not merely a shared/touching vertex).
fn boundary_crosses(a: &Loop, b: &Loop) -> bool {
    if a.num_vertices() < 2 || b.num_vertices() < 2 {
        return false;
    }
    for (a0, a1) in loop_edges(a) {
        for (b0, b1) in loop_edges(b) {
            if crossing_sign(a0, a1, b0, b1) == Crossing::Cross {
                return true;
            }
        }
    }
    false
}

/// True iff `inner` lies entirely within `outer`: every `inner` vertex is inside
/// `outer` and their boundaries do not cross.
fn loop_contains_loop(outer: &Loop, inner: &Loop) -> bool {
    inner.vertices().iter().all(|&v| outer.contains_point(v)) && !boundary_crosses(outer, inner)
}

/// The bounding cap of a set of points: centered on their normalized vector mean,
/// with radius equal to the farthest point's angular distance.
fn bounding_cap_of(points: &[S2Point]) -> Cap {
    if points.is_empty() {
        return Cap::empty(S2Point::from_vector(Vector::new(1.0, 0.0, 0.0)));
    }
    let mut sum = Vector::new(0.0, 0.0, 0.0);
    for p in points {
        sum = sum.add(p.vector());
    }
    let mean = sum.normalize();
    // Antipodal points can cancel to the zero vector; fall back to the first.
    let center = if mean.norm2() == 0.0 {
        points[0]
    } else {
        S2Point::from_vector(mean)
    };
    let radius = points
        .iter()
        .map(|p| center.vector().angle(p.vector()).radians())
        .fold(0.0_f64, f64::max);
    Cap::from_center_angle(center, Angle::from_radians(radius))
}

impl Region for Cap {
    fn contains_point(&self, p: S2Point) -> bool {
        Cap::contains_point(*self, p)
    }
    fn bounding_cap(&self) -> Cap {
        *self
    }
    fn intersects_cell(&self, cell: S2CellId) -> bool {
        let verts = cell.vertices();
        verts.iter().any(|&v| self.contains_point(v)) || self.intersects(bounding_cap_of(&verts))
    }
    fn contains_cell(&self, cell: S2CellId) -> bool {
        // A cap is convex; if every corner is inside, the cell is (near enough)
        // inside. Approximate — `Cap` is not the coverer's primary region.
        cell.vertices().iter().all(|&v| self.contains_point(v))
    }
}

impl Region for Loop {
    fn contains_point(&self, p: S2Point) -> bool {
        Loop::contains_point(self, p)
    }
    fn bounding_cap(&self) -> Cap {
        bounding_cap_of(self.vertices())
    }
    fn intersects_cell(&self, cell: S2CellId) -> bool {
        // `Loop::intersects` is true on a boundary crossing or either-way nesting.
        self.intersects(&cell_loop(cell))
    }
    fn contains_cell(&self, cell: S2CellId) -> bool {
        loop_contains_loop(self, &cell_loop(cell))
    }
}

impl Region for Polygon {
    fn contains_point(&self, p: S2Point) -> bool {
        Polygon::contains_point(self, p)
    }
    fn bounding_cap(&self) -> Cap {
        bounding_cap_of(self.shell().vertices())
    }
    fn intersects_cell(&self, cell: S2CellId) -> bool {
        let cl = cell_loop(cell);
        // Touches the shell and is not swallowed whole by a hole.
        self.shell().intersects(&cl) && !self.holes().iter().any(|h| loop_contains_loop(h, &cl))
    }
    fn contains_cell(&self, cell: S2CellId) -> bool {
        let cl = cell_loop(cell);
        loop_contains_loop(self.shell(), &cl) && !self.holes().iter().any(|h| h.intersects(&cl))
    }
}

// ---- the coverer ------------------------------------------------------------

/// Approximates a [`Region`] by a set of S2 cells whose union contains it.
///
/// The covering is built top-down from the six face cells: a cell that does not
/// intersect the region is dropped; a cell at `max_level` — or at/above
/// `min_level` and fully inside the region — is emitted; otherwise it is
/// subdivided into its four children. `max_cells` is a **soft** budget: once
/// reached, cells at/above `min_level` are emitted coarsely instead of refined,
/// so the result stays a valid (super-)cover rather than honouring the budget
/// exactly (unlike golang/geo's priority-queue coverer).
#[derive(Clone, Copy, Debug)]
pub struct RegionCoverer {
    /// Coarsest level a covering cell may have.
    pub min_level: u32,
    /// Finest level a covering cell may have (refinement stops here). Values
    /// above the leaf level (30) are treated as 30.
    pub max_level: u32,
    /// Soft cap on the number of cells.
    pub max_cells: usize,
}

impl Default for RegionCoverer {
    fn default() -> Self {
        Self {
            min_level: 0,
            max_level: MAX_LEVEL,
            max_cells: 8,
        }
    }
}

impl RegionCoverer {
    /// A coverer that emits cells at exactly `level`: every level-`level` cell
    /// that intersects the region, with no coarsening and no cell budget. This is
    /// the fixed-resolution covering a partitioner uses to bucket a geometry.
    /// `level` is clamped to the leaf level (30).
    #[must_use]
    pub fn at_level(level: u32) -> Self {
        let level = level.min(MAX_LEVEL);
        Self {
            min_level: level,
            max_level: level,
            max_cells: usize::MAX,
        }
    }

    /// The covering cell ids for `region`.
    #[must_use]
    pub fn covering<R: Region + ?Sized>(&self, region: &R) -> Vec<S2CellId> {
        let mut out = Vec::new();
        for face in 0..6u8 {
            self.cover_cell(region, S2CellId::from_face(face), &mut out);
        }
        out
    }

    fn cover_cell<R: Region + ?Sized>(&self, region: &R, cell: S2CellId, out: &mut Vec<S2CellId>) {
        if !region.intersects_cell(cell) {
            return;
        }
        let level = cell.level();
        // Never refine past a leaf: `S2CellId::children()` is undefined at the
        // leaf level, so clamp the finest level regardless of `self.max_level`.
        let max_level = self.max_level.min(MAX_LEVEL);
        let at_budget = out.len() >= self.max_cells;
        if level >= max_level
            || (level >= self.min_level && (at_budget || region.contains_cell(cell)))
        {
            out.push(cell);
            return;
        }
        for child in cell.children() {
            self.cover_cell(region, child, out);
        }
    }
}

#[cfg(test)]
mod tests {
    use std::f64::consts::PI;

    use super::*;
    use crate::earth::GeoPoint;

    fn pt(lon: f64, lat: f64) -> S2Point {
        S2Point::from_geo_point(GeoPoint::from_lat_lon(lat, lon))
    }

    fn leaf(lon: f64, lat: f64) -> S2CellId {
        S2CellId::from_geo_point(GeoPoint::from_lat_lon(lat, lon))
    }

    /// A small CCW lon/lat box `[lon0,lon1] × [lat0,lat1]` as a spherical loop
    /// (interior to the left — the small side).
    fn square(lon0: f64, lat0: f64, lon1: f64, lat1: f64) -> Loop {
        let l = Loop::new(vec![
            pt(lon0, lat0),
            pt(lon1, lat0),
            pt(lon1, lat1),
            pt(lon0, lat1),
        ]);
        if l.area() > 2.0 * PI {
            let mut v = l.vertices().to_vec();
            v.reverse();
            Loop::new(v)
        } else {
            l
        }
    }

    #[test]
    fn cell_loop_is_oriented_and_contains_its_center() {
        // A handful of cells across different faces: each cell loop must contain
        // the cell center (i.e. it is wound CCW after orientation).
        for &(lon, lat) in &[(0.0, 0.0), (90.0, 10.0), (-120.0, -45.0), (179.0, 80.0)] {
            let cell = S2CellId::from_geo_point_at_level(GeoPoint::from_lat_lon(lat, lon), 6);
            let cl = cell_loop(cell);
            assert!(cl.contains_point(cell.center_point()), "({lon},{lat})");
        }
    }

    #[test]
    fn covering_contains_interior_excludes_exterior() {
        let region = square(-1.0, -1.0, 1.0, 1.0);
        assert!(region.contains_point(pt(0.0, 0.0)), "origin inside the loop");
        assert!(!region.contains_point(pt(20.0, 20.0)), "far point outside");

        let cover = RegionCoverer::at_level(8).covering(&region);
        assert!(!cover.is_empty());
        for &c in &cover {
            assert_eq!(c.level(), 8);
            assert!(region.intersects_cell(c));
        }
        assert!(
            cover.iter().any(|c| c.contains(leaf(0.0, 0.0))),
            "interior point not covered"
        );
        assert!(
            !cover.iter().any(|c| c.contains(leaf(45.0, 45.0))),
            "exterior point wrongly covered"
        );
    }

    #[test]
    fn multi_level_cover_bounds_levels_and_covers_interior() {
        let region = square(-2.0, -2.0, 2.0, 2.0);
        let cover = RegionCoverer {
            min_level: 2,
            max_level: 12,
            max_cells: 64,
        }
        .covering(&region);
        assert!(!cover.is_empty());
        for &c in &cover {
            assert!((2..=12).contains(&c.level()));
            assert!(region.intersects_cell(c));
        }
        assert!(cover.iter().any(|c| c.contains(leaf(0.0, 0.0))));
    }

    #[test]
    fn polygon_hole_is_excluded_from_the_cover() {
        let poly = Polygon::new(
            square(-4.0, -4.0, 4.0, 4.0),
            vec![square(-1.0, -1.0, 1.0, 1.0)],
        );
        let cover = RegionCoverer::at_level(10).covering(&poly);
        assert!(
            cover.iter().any(|c| c.contains(leaf(0.0, 3.0))),
            "ring point should be covered"
        );
        assert!(
            !cover.iter().any(|c| c.contains(leaf(0.0, 0.0))),
            "deep-hole point should not be covered"
        );
    }

    #[test]
    fn at_level_clamps_to_the_leaf() {
        assert_eq!(RegionCoverer::at_level(7).max_level, 7);
        let clamped = RegionCoverer::at_level(40);
        assert_eq!(clamped.min_level, MAX_LEVEL);
        assert_eq!(clamped.max_level, MAX_LEVEL);
    }
}
