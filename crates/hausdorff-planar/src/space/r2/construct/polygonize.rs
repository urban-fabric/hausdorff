//! `polygonize` — build polygons from a network of lines.
//!
//! The inverse of taking a polygon's boundary, and the operation behind "turn
//! this road centreline network into a block layer" or "turn these surveyed
//! boundary lines into parcels". Cadastral and street data arrives as lines far
//! more often than as polygons, because a line is what a surveyor records and a
//! shared boundary only has to be drawn once.
//!
//! ## What makes it more than ring-finding
//!
//! The input is a *soup*: lines that cross without sharing a vertex, lines that
//! overlap, lines that dangle, duplicates. So the real work is in the
//! arrangement, and this module reuses the one [`make_valid`]
//! already built — cut every segment against every other on a shared vertex pool,
//! then walk. That reuse is the point. An independent noding implementation here
//! would be a second set of tolerance decisions to keep in sync with the first,
//! and the two would drift.
//!
//! What is new is the **classification**, because polygonize answers a different
//! question than repair does. `make_valid` keeps edges under an odd cover, which
//! is a fill rule. Polygonize keeps *every* edge that bounds a face, and then has
//! to say what happened to the edges that bound nothing.
//!
//! The *nesting* is shared too — `make_valid::assemble` sorts rings into shells
//! and holes by containment depth, and that is the same problem here. Worth noting
//! why reuse matters rather than reimplementing it: the subtlety is in choosing a
//! point to test containment *with*. A ring's centroid is inside the ring but can
//! also be inside a different ring that the first merely encloses, which inflates
//! its depth and turns a shell into a hole. `assemble` steps in from the lowest
//! vertex instead, staying next to its own boundary. The first version of this
//! module used a centroid and produced zero polygons for a ring inside a ring.
//!
//! ## Edges that bound nothing are reported, not discarded
//!
//! JTS splits the leftovers into two kinds, and the distinction is diagnostic
//! rather than cosmetic:
//!
//! - a **dangle** is an edge with a loose end — it reaches a vertex nothing else
//!   touches, so it cannot be part of any ring. Usually a line that stops short of
//!   the boundary it was meant to meet.
//! - a **cut edge** (a bridge) has both ends connected but still bounds no face —
//!   removing it would disconnect the network. Usually a spur or a connector drawn
//!   between two otherwise separate loops.
//!
//! Both mean the input did not quite close, and which one it is tells you what to
//! fix: a dangle wants its endpoint extended, a cut edge wants a missing parallel
//! line. Returning them is what makes polygonize usable as a data-quality tool
//! rather than only as a converter — a silent "here are 97 of your 100 parcels" is
//! the worst possible outcome.

use std::collections::{HashMap, HashSet};

use super::make_valid::{self, node, trace};
use crate::space::r2::{Coord, Geometry, LineString, MultiPolygon};

/// The result of polygonizing a line network.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Polygonized {
    /// The faces that closed, with holes nested into their shells.
    pub polygons: MultiPolygon,
    /// Edges with a loose end, which can bound no face.
    pub dangles: Vec<LineString>,
    /// Edges connected at both ends that still bound no face — bridges whose
    /// removal would disconnect the network.
    pub cut_edges: Vec<LineString>,
}

impl Polygonized {
    /// Whether every input edge ended up bounding a face.
    ///
    /// The question a caller converting a parcel layer actually wants answered:
    /// not "did I get polygons" but "did I lose anything".
    #[must_use]
    pub fn is_fully_closed(&self) -> bool {
        self.dangles.is_empty() && self.cut_edges.is_empty()
    }
}

/// Build polygons from a line network.
pub trait Polygonize {
    /// Node the input, trace its faces, and report the edges that bound none.
    fn polygonize(&self) -> Polygonized;
}

/// Collect every segment of any geometry, so a caller can pass a mixed
/// collection without flattening it first.
///
/// Areal inputs contribute their ring edges. That is deliberate rather than an
/// oversight: polygonizing a polygon's own boundary should give the polygon back,
/// and it is the natural way to re-derive faces from a layer that has already been
/// split into overlapping pieces.
fn segments_of(geom: &Geometry, out: &mut Vec<(Coord, Coord)>) {
    let mut ring = |cs: &[Coord]| {
        for w in cs.windows(2) {
            if w[0] != w[1] {
                out.push((w[0], w[1]));
            }
        }
    };
    match geom {
        Geometry::Line(l) => {
            if l.start != l.end {
                out.push((l.start, l.end));
            }
        }
        Geometry::LineString(ls) => ring(&ls.0),
        Geometry::MultiLineString(mls) => mls.iter().for_each(|ls| ring(&ls.0)),
        Geometry::Polygon(p) => {
            ring(&p.exterior().0);
            for h in p.interiors() {
                ring(&h.0);
            }
        }
        Geometry::MultiPolygon(mp) => {
            for p in mp.iter() {
                ring(&p.exterior().0);
                for h in p.interiors() {
                    ring(&h.0);
                }
            }
        }
        Geometry::Rect(r) => ring(&r.to_polygon().exterior().0),
        Geometry::Triangle(t) => ring(&t.to_polygon().exterior().0),
        Geometry::GeometryCollection(gs) => {
            for g in gs {
                segments_of(g, out);
            }
        }
        // A point bounds nothing and contributes no edge.
        Geometry::Point(_) | Geometry::MultiPoint(_) => {}
    }
}

fn polygonize_segments(segments: &[(Coord, Coord)]) -> Polygonized {
    let mut extent = 1.0_f64;
    for &(a, b) in segments {
        for c in [a, b] {
            if !c.x.is_finite() || !c.y.is_finite() {
                // A non-finite coordinate cannot be noded or traced. Returning an
                // empty result rather than erroring keeps the signature simple and
                // matches `make_valid`'s treatment of the same input as unusable.
                return Polygonized::default();
            }
            extent = extent.max(c.x.abs()).max(c.y.abs());
        }
    }
    if segments.is_empty() {
        return Polygonized::default();
    }

    let mut nodes = node::Nodes::new(extent);
    let pieces = node::arrange(segments, &mut nodes);
    if pieces.is_empty() {
        return Polygonized::default();
    }

    // Deduplicate: the same edge drawn twice (a shared parcel boundary digitized
    // from both sides) is one edge of the arrangement, not two. Without this a
    // duplicated boundary would be traced as a zero-area sliver face.
    let mut unique: Vec<(usize, usize)> = pieces
        .iter()
        .map(|&(a, b)| (a.min(b), a.max(b)))
        .collect::<HashSet<_>>()
        .into_iter()
        .collect();
    // A hash set yields a different order each process and the trace must not
    // depend on that.
    unique.sort_unstable();

    // Strip dangles before tracing. An edge reaching a degree-1 vertex cannot
    // bound a face, and removing it can expose another — a chain of spurs peels
    // off one layer per pass — so this iterates to a fixed point.
    let mut live: Vec<bool> = vec![true; unique.len()];
    let mut dangle_edges: Vec<usize> = Vec::new();
    loop {
        let mut degree: HashMap<usize, usize> = HashMap::new();
        for (i, &(a, b)) in unique.iter().enumerate() {
            if !live[i] {
                continue;
            }
            *degree.entry(a).or_insert(0) += 1;
            *degree.entry(b).or_insert(0) += 1;
        }
        let mut peeled = false;
        for (i, &(a, b)) in unique.iter().enumerate() {
            if !live[i] {
                continue;
            }
            if degree.get(&a).copied().unwrap_or(0) < 2 || degree.get(&b).copied().unwrap_or(0) < 2
            {
                live[i] = false;
                dangle_edges.push(i);
                peeled = true;
            }
        }
        if !peeled {
            break;
        }
    }

    let remaining: Vec<(usize, usize)> = unique
        .iter()
        .enumerate()
        .filter(|(i, _)| live[*i])
        .map(|(_, &e)| e)
        .collect();

    let rings = if remaining.is_empty() {
        Vec::new()
    } else {
        trace::rings(&remaining, &nodes).unwrap_or_default()
    };

    // Any surviving edge that no traced ring uses is a cut edge: both ends are
    // connected, so it is not a dangle, yet it bounds no face.
    let mut used: HashSet<(usize, usize)> = HashSet::new();
    for ring in &rings {
        let closed: Vec<Coord> = ring.clone();
        for w in closed.windows(2) {
            if let (Some(u), Some(v)) = (
                find(&nodes, &remaining, w[0]),
                find(&nodes, &remaining, w[1]),
            ) {
                used.insert((u.min(v), u.max(v)));
            }
        }
        if closed.len() >= 2
            && let (Some(u), Some(v)) = (
                find(&nodes, &remaining, closed[closed.len() - 1]),
                find(&nodes, &remaining, closed[0]),
            )
        {
            used.insert((u.min(v), u.max(v)));
        }
    }

    let to_line = |&(a, b): &(usize, usize)| LineString::new(vec![nodes.coord(a), nodes.coord(b)]);
    let cut_edges: Vec<LineString> = remaining
        .iter()
        .filter(|&&(a, b)| !used.contains(&(a.min(b), a.max(b))))
        .map(to_line)
        .collect();
    let dangles: Vec<LineString> = dangle_edges.iter().map(|&i| to_line(&unique[i])).collect();

    Polygonized {
        polygons: make_valid::assemble(rings),
        dangles,
        cut_edges,
    }
}

/// The pool index of a coordinate, by exact match against the edge endpoints.
///
/// The arrangement guarantees traced ring vertices are pooled coordinates, so an
/// exact comparison is correct here and no tolerance is wanted — a tolerance would
/// risk matching a *different* nearby vertex.
fn find(nodes: &node::Nodes, edges: &[(usize, usize)], c: Coord) -> Option<usize> {
    edges
        .iter()
        .flat_map(|&(a, b)| [a, b])
        .find(|&i| nodes.coord(i) == c)
}

impl Polygonize for Geometry {
    fn polygonize(&self) -> Polygonized {
        let mut segments = Vec::new();
        segments_of(self, &mut segments);
        polygonize_segments(&segments)
    }
}

impl Polygonize for [Geometry] {
    fn polygonize(&self) -> Polygonized {
        let mut segments = Vec::new();
        for g in self {
            segments_of(g, &mut segments);
        }
        polygonize_segments(&segments)
    }
}

impl Polygonize for crate::space::r2::MultiLineString {
    fn polygonize(&self) -> Polygonized {
        let mut segments = Vec::new();
        segments_of(&Geometry::MultiLineString(self.clone()), &mut segments);
        polygonize_segments(&segments)
    }
}

impl Polygonize for [LineString] {
    fn polygonize(&self) -> Polygonized {
        let mut segments = Vec::new();
        for ls in self {
            segments_of(&Geometry::LineString(ls.clone()), &mut segments);
        }
        polygonize_segments(&segments)
    }
}
