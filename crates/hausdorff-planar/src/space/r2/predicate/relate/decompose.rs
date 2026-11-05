//! The [`Decomposed`] view of a geometry — its 0-D points, mod-2 boundary
//! points, 1-D [`Seg`]ments, and 2-D polygons (with a prebuilt
//! [`ArealIndex`]) — and the point [`locate`](Decomposed::locate)
//! that classifies a coordinate against it.

use super::Loc;
use super::geom::dimension;
use super::noding::ArealIndex;
use super::predicates::{locate_areal, point_on_segment_interior};
use crate::space::r2::util::ring_vertices;
use crate::space::r2::{Coord, Geometry, Polygon};

/// A 1-D edge tagged with its location in its *own* geometry: a line/curve edge
/// is part of the Interior; a polygon ring edge is part of the Boundary. Tagging
/// avoids re-deriving that location with floating-point geometry.
#[derive(Clone, Copy)]
pub(super) struct Seg {
    pub(super) a: Coord,
    pub(super) b: Coord,
    pub(super) owner: Loc,
}

/// Decomposed view of a geometry: its 0-D points, boundary points (mod-2 for
/// lines), 1-D segments, and 2-D polygons.
pub(super) struct Decomposed {
    pub(super) points: Vec<Coord>,
    boundary: Vec<Coord>,
    pub(super) segments: Vec<Seg>,
    pub(super) polygons: Vec<Polygon>,
    pub(super) dim: i8,
    /// Ring-edge spatial index for O(log n) areal point location, built once.
    areal_index: Option<ArealIndex>,
    /// True iff any segment is a line/curve edge (owner `Interior`) rather than a
    /// polygon ring edge. When false, `locate`'s linear branch can't match, so it
    /// is skipped instead of scanning every ring edge only to `continue`.
    has_linear: bool,
}

impl Decomposed {
    pub(super) fn new(g: &Geometry) -> Self {
        let mut d = Decomposed {
            points: Vec::new(),
            boundary: Vec::new(),
            segments: Vec::new(),
            polygons: Vec::new(),
            dim: dimension(g),
            areal_index: None,
            has_linear: false,
        };
        d.add(g);
        d.areal_index = ArealIndex::build(&d.polygons);
        d.has_linear = d.segments.iter().any(|s| s.owner == Loc::Interior);
        d
    }

    fn add(&mut self, g: &Geometry) {
        match g {
            Geometry::Point(p) => self.points.push(p.0),
            Geometry::MultiPoint(mp) => self.points.extend(mp.0.iter().map(|p| p.0)),
            Geometry::Line(l) => {
                if l.start != l.end {
                    self.segments.push(Seg {
                        a: l.start,
                        b: l.end,
                        owner: Loc::Interior,
                    });
                    self.add_line_boundary(&[l.start, l.end]);
                } else {
                    self.points.push(l.start);
                }
            }
            Geometry::LineString(ls) => self.add_linestring(&ls.0),
            Geometry::MultiLineString(mls) => {
                for ls in mls.iter() {
                    self.add_linestring(&ls.0);
                }
            }
            Geometry::Polygon(p) => self.add_polygon(p),
            Geometry::MultiPolygon(mp) => {
                for p in mp.iter() {
                    self.add_polygon(p);
                }
            }
            Geometry::Rect(r) => self.add_polygon(&r.to_polygon()),
            Geometry::Triangle(t) => self.add_polygon(&t.to_polygon()),
            Geometry::GeometryCollection(gs) => {
                for sub in gs {
                    self.add(sub);
                }
            }
        }
    }

    fn add_linestring(&mut self, coords: &[Coord]) {
        for w in coords.windows(2) {
            if w[0] != w[1] {
                self.segments.push(Seg {
                    a: w[0],
                    b: w[1],
                    owner: Loc::Interior,
                });
            }
        }
        if coords.len() >= 2 {
            self.add_line_boundary(&[coords[0], coords[coords.len() - 1]]);
        }
    }

    /// Boundary endpoints follow the mod-2 rule: a point that is an endpoint of
    /// an odd number of curve components is on the boundary.
    fn add_line_boundary(&mut self, endpoints: &[Coord]) {
        for &e in endpoints {
            if e.is_finite() {
                self.boundary.push(e);
            }
        }
    }

    fn add_polygon(&mut self, p: &Polygon) {
        ring_segments(&p.exterior().0, &mut self.segments);
        for h in p.interiors() {
            ring_segments(&h.0, &mut self.segments);
        }
        self.polygons.push(p.clone());
    }

    pub(super) fn all_vertices(&self) -> Vec<Coord> {
        let mut v: Vec<Coord> = self.points.clone();
        for s in &self.segments {
            v.push(s.a);
            v.push(s.b);
        }
        v
    }

    /// Effective boundary points after the mod-2 rule.
    pub(super) fn boundary_points(&self) -> Vec<Coord> {
        let mut counts: Vec<(Coord, usize)> = Vec::new();
        for &b in &self.boundary {
            if let Some(c) = counts.iter_mut().find(|(p, _)| *p == b) {
                c.1 += 1;
            } else {
                counts.push((b, 1));
            }
        }
        counts
            .into_iter()
            .filter(|(_, n)| n % 2 == 1)
            .map(|(p, _)| p)
            .collect()
    }

    pub(super) fn locate(&self, p: Coord) -> Loc {
        if !self.polygons.is_empty() {
            let areal = match &self.areal_index {
                Some(idx) => idx.locate(p),
                None => locate_areal(p, &self.polygons),
            };
            match areal {
                Loc::Interior => return Loc::Interior,
                Loc::Boundary => return Loc::Boundary,
                Loc::Exterior => {}
            }
        }
        if self.has_linear {
            let bpts = self.boundary_points();
            if bpts.contains(&p) {
                return Loc::Boundary;
            }
            let areal = !self.polygons.is_empty();
            for s in &self.segments {
                if areal && s.owner == Loc::Boundary {
                    continue;
                }
                if p == s.a || p == s.b || point_on_segment_interior(p, s.a, s.b) {
                    return s.owner;
                }
            }
        }
        if self.points.contains(&p) {
            return Loc::Interior;
        }
        Loc::Exterior
    }

    /// Indexed areal-only location of `p` (ignoring any linear/puntal parts) —
    /// the accelerated equivalent of `locate_areal(p, &self.polygons)`.
    pub(super) fn locate_areal_only(&self, p: Coord) -> Loc {
        match &self.areal_index {
            Some(idx) => idx.locate(p),
            None => locate_areal(p, &self.polygons),
        }
    }
}

fn ring_segments(coords: &[Coord], out: &mut Vec<Seg>) {
    let v = ring_vertices(coords);
    let n = v.len();
    if n < 2 {
        return;
    }
    for i in 0..n {
        let (a, b) = (v[i], v[(i + 1) % n]);
        if a != b {
            out.push(Seg {
                a,
                b,
                owner: Loc::Boundary,
            });
        }
    }
}
