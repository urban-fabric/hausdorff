//! Delaunay triangulation, and the Voronoi diagram that is its dual.
//!
//! The Delaunay triangulation of a point set is the one where no point lies inside
//! any triangle's circumcircle. That single property is what makes it the default
//! triangulation everywhere: among all triangulations of the same points it
//! maximises the smallest angle, so it avoids the slivers that make interpolation
//! on a TIN surface unstable.
//!
//! Its dual is the **Voronoi diagram** — one cell per input point, containing
//! everything closer to that point than to any other. Nearest-facility questions
//! ("which depot serves this address") are a point-in-polygon lookup once you have
//! it, rather than a scan over every facility.
//!
//! ## Why the dual comes for free, and why that matters
//!
//! A Voronoi cell's vertices are the circumcentres of the Delaunay triangles
//! around its point, in order. So building Voronoi after Delaunay is bookkeeping
//! rather than a second algorithm — which is the reason to do them together
//! instead of implementing Fortune's sweep separately. One geometric predicate to
//! trust, not two.
//!
//! ## The incircle test
//!
//! The whole algorithm is driven by "is this point inside that circumcircle", and
//! that is a *sign* question — so it has the same fragility as orientation, and
//! for the same reason. The determinant is degree 4 in the coordinates, so
//! cancellation is worse than for orientation's degree 2.
//!
//! This uses the `f64` determinant with the coordinates **translated to their own
//! centroid** first. That translation is not cosmetic: the determinant's
//! conditioning depends on coordinate magnitude, and geographic data arrives at
//! magnitudes like `(−122.4, 37.8)` or worse, projected, at `(551000, 4182000)`.
//! Subtracting the centroid brings the working values near zero, which is where
//! `f64` has its digits. The same input untranslated loses roughly
//! `log2(magnitude)` bits of the answer.
//!
//! There is an exact `i128` incircle in [`in_circle_position`](crate::space::r2::in_circle_position), and it is *not*
//! used here. It requires coordinates already placed on an integer grid, and it
//! overflows above roughly `2^42` because of that degree-4 determinant — so
//! feeding it arbitrary `f64` data means choosing a scale factor, and the choice
//! would silently change which triangulation you get. Translating and staying in
//! `f64` is the honest trade. For input genuinely on an integer grid, the exact
//! predicate is the better tool and is available directly.
//!
//! **What that costs**: for cocircular or near-cocircular points — four corners of
//! a square, a regular polygon's vertices — which diagonal gets chosen is
//! arbitrary. Both choices are valid Delaunay triangulations of such a set, since
//! the Delaunay triangulation is only unique when no four points are cocircular.
//! The result is always *a* valid triangulation; it is the tie-break that is
//! unspecified.

use std::collections::{HashMap, HashSet};

use crate::space::r2::{Coord, Geometry, LineString, MultiPolygon, Point, Polygon, Triangle};

/// A Delaunay triangulation of a point set.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Triangulation {
    /// The input points, deduplicated, in the order the triangles index them.
    pub points: Vec<Coord>,
    /// Each triangle as three indices into [`points`](Self::points),
    /// counter-clockwise.
    pub triangles: Vec<[usize; 3]>,
}

impl Triangulation {
    /// The triangles as geometry.
    #[must_use]
    pub fn to_triangles(&self) -> Vec<Triangle> {
        self.triangles
            .iter()
            .map(|&[a, b, c]| Triangle(self.points[a], self.points[b], self.points[c]))
            .collect()
    }

    /// The triangles as a `MultiPolygon`, for callers that want one geometry.
    #[must_use]
    pub fn to_multi_polygon(&self) -> MultiPolygon {
        MultiPolygon::new(
            self.to_triangles()
                .into_iter()
                .map(|t| t.to_polygon())
                .collect(),
        )
    }

    /// The unique undirected edges of the triangulation.
    ///
    /// This is the **Delaunay graph**, and it is the useful output in its own
    /// right: it contains the Euclidean minimum spanning tree and the
    /// nearest-neighbour graph, so proximity questions reduce to walking it.
    #[must_use]
    pub fn edges(&self) -> Vec<(usize, usize)> {
        let mut set: HashSet<(usize, usize)> = HashSet::new();
        for &[a, b, c] in &self.triangles {
            for (u, v) in [(a, b), (b, c), (c, a)] {
                set.insert((u.min(v), u.max(v)));
            }
        }
        let mut out: Vec<(usize, usize)> = set.into_iter().collect();
        // A hash set orders differently each process; callers comparing output
        // across runs need this stable.
        out.sort_unstable();
        out
    }

    /// Whether the Delaunay property holds: no point strictly inside any
    /// triangle's circumcircle.
    ///
    /// Exposed rather than kept in the tests because it is the definition, and a
    /// caller feeding degenerate input may want to check rather than assume. The
    /// tolerance is relative to the circumradius, so cocircular points — which are
    /// *on* the circle, not inside it — do not count as violations.
    #[must_use]
    pub fn is_delaunay(&self) -> bool {
        if self.points.is_empty() {
            return true;
        }
        // Translated to the centroid before any circumcircle is formed, for exactly
        // the reason the builder translates: the circumcentre of three points at
        // magnitude 1e9 is a difference of squares near 1e18, and the check would
        // report violations that are artifacts of its own arithmetic rather than
        // defects in the triangulation it is checking. A checker less accurate than
        // the thing it checks is worse than no checker.
        let n = self.points.len() as f64;
        let shift = Coord::new(
            self.points.iter().map(|p| p.x).sum::<f64>() / n,
            self.points.iter().map(|p| p.y).sum::<f64>() / n,
        );
        let at = |i: usize| Coord::new(self.points[i].x - shift.x, self.points[i].y - shift.y);
        self.triangles.iter().all(|&[a, b, c]| {
            let Some((centre, radius)) = circumcircle(at(a), at(b), at(c)) else {
                // A degenerate triangle has no circumcircle to violate.
                return true;
            };
            let slack = 1e-9 * radius.max(1.0);
            (0..self.points.len()).all(|i| {
                if i == a || i == b || i == c {
                    return true;
                }
                let p = at(i);
                let d = ((p.x - centre.x).powi(2) + (p.y - centre.y).powi(2)).sqrt();
                d >= radius - slack
            })
        })
    }
}

/// A Voronoi diagram: one cell per input point.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Voronoi {
    /// The input points, in the same order as [`cells`](Self::cells).
    pub sites: Vec<Coord>,
    /// The cell for each site, or `None` where it could not be closed.
    ///
    /// `None` is not a failure to compute — it is the honest answer for a site on
    /// the convex hull, whose true cell is **unbounded**. An unbounded region has
    /// no polygon, so the alternatives are to clip it to some arbitrary rectangle
    /// (which invents geometry the caller did not ask for) or to say so. See
    /// [`clipped_to`](Self::clipped_to) for the clip when you want one.
    pub cells: Vec<Option<Polygon>>,
}

impl Voronoi {
    /// Every cell that closed, as a `MultiPolygon`.
    #[must_use]
    pub fn to_multi_polygon(&self) -> MultiPolygon {
        MultiPolygon::new(self.cells.iter().flatten().cloned().collect())
    }

    /// The cells clipped to `bounds`, so hull sites get a cell too.
    ///
    /// The clip is what makes the diagram usable as a lookup table — an unbounded
    /// cell cannot be put in a spatial index. The rectangle is the caller's
    /// choice precisely because it is arbitrary: it should be the extent of the
    /// area the answer is wanted for.
    #[must_use]
    pub fn clipped_to(&self, bounds: crate::space::r2::Rect) -> Vec<Option<Polygon>> {
        use crate::space::r2::BooleanOps;
        let frame = MultiPolygon::new(vec![bounds.to_polygon()]);
        self.cells
            .iter()
            .map(|cell| {
                let cell = cell.as_ref()?;
                let clipped = MultiPolygon::new(vec![cell.clone()]).intersection(&frame);
                clipped.0.into_iter().next()
            })
            .collect()
    }
}

/// The circumcircle of three points, or `None` when they are collinear.
fn circumcircle(a: Coord, b: Coord, c: Coord) -> Option<(Coord, f64)> {
    let d = 2.0 * (a.x * (b.y - c.y) + b.x * (c.y - a.y) + c.x * (a.y - b.y));
    if d == 0.0 || !d.is_finite() {
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
    let r = ((a.x - centre.x).powi(2) + (a.y - centre.y).powi(2)).sqrt();
    r.is_finite().then_some((centre, r))
}

/// Whether `p` is strictly inside the circumcircle of the counter-clockwise
/// triangle `a, b, c`.
///
/// The in-circle determinant, evaluated on coordinates already translated near
/// the origin by the caller — see the module header for why that translation is
/// load-bearing rather than tidy.
fn in_circumcircle(a: Coord, b: Coord, c: Coord, p: Coord) -> bool {
    let (ax, ay) = (a.x - p.x, a.y - p.y);
    let (bx, by) = (b.x - p.x, b.y - p.y);
    let (cx, cy) = (c.x - p.x, c.y - p.y);
    let det = (ax * ax + ay * ay) * (bx * cy - cx * by) - (bx * bx + by * by) * (ax * cy - cx * ay)
        + (cx * cx + cy * cy) * (ax * by - bx * ay);
    det > 0.0
}

fn orient(a: Coord, b: Coord, c: Coord) -> f64 {
    (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x)
}

/// Delaunay-triangulate a point set by Bowyer–Watson incremental insertion.
///
/// The structure: start from one huge triangle containing everything, insert each
/// point by deleting the triangles whose circumcircle contains it and re-filling
/// the resulting cavity, then drop every triangle still touching the super-triangle.
///
/// The super-triangle is sized from the data extent rather than from a fixed
/// constant, because a fixed one is either too small for large coordinates (and
/// silently clips the hull) or so large that its vertices dominate every
/// determinant and destroy the precision of the real ones.
#[must_use]
pub fn delaunay(points: &[Coord]) -> Triangulation {
    // Deduplicate and drop non-finite input. Duplicate sites have no defined
    // Voronoi cell between them and would make the cavity walk ambiguous.
    //
    // Merged within a *relative* tolerance rather than on exact equality. Two
    // points a few ulp apart are distinct by `==` but indistinguishable to the
    // degree-4 in-circle determinant, so keeping both asks the predicate a question
    // it cannot answer and the result is not Delaunay with respect to either. The
    // tolerance scales with the data so it means the same thing at 1e-9 and 1e9.
    let finite: Vec<Coord> = points.iter().copied().filter(|c| c.is_finite()).collect();
    let span = finite
        .iter()
        .map(|p| p.x.abs().max(p.y.abs()))
        .fold(1.0_f64, f64::max);
    let merge_tol = span * 1e-12;
    let mut pts: Vec<Coord> = Vec::new();
    for p in finite {
        if !pts
            .iter()
            .any(|q| (q.x - p.x).abs() <= merge_tol && (q.y - p.y).abs() <= merge_tol)
        {
            pts.push(p);
        }
    }
    if pts.len() < 3 {
        return Triangulation {
            points: pts,
            triangles: Vec::new(),
        };
    }

    // Translate to the centroid before any determinant is evaluated. This is the
    // precision decision described in the module header.
    let n = pts.len() as f64;
    let shift = Coord::new(
        pts.iter().map(|p| p.x).sum::<f64>() / n,
        pts.iter().map(|p| p.y).sum::<f64>() / n,
    );
    let work: Vec<Coord> = pts
        .iter()
        .map(|p| Coord::new(p.x - shift.x, p.y - shift.y))
        .collect();

    let extent = work
        .iter()
        .map(|p| p.x.abs().max(p.y.abs()))
        .fold(1.0_f64, f64::max);
    let big = extent * 1000.0;
    // Three super-triangle vertices, appended after the real points so their
    // indices are known and easy to strip.
    let mut all = work.clone();
    let s0 = all.len();
    all.push(Coord::new(-big * 3.0, -big));
    all.push(Coord::new(big * 3.0, -big));
    all.push(Coord::new(0.0, big * 3.0));

    let mut tris: Vec<[usize; 3]> = vec![[s0, s0 + 1, s0 + 2]];

    for i in 0..work.len() {
        let p = all[i];

        // The cavity is grown by **breadth-first search across triangle
        // adjacency**, not by scanning every triangle independently. That
        // distinction is the whole robustness of this step.
        //
        // Bowyer–Watson assumes the cavity — the triangles whose circumcircle
        // contains `p` — is a connected, star-shaped region, because it then has a
        // single simple boundary cycle that can be re-fanned to `p`. A global scan
        // does not guarantee that. When many points are cocircular (a lattice, a
        // regular polygon, anything on a common circle) the in-circle determinant
        // is near zero for all of them and rounds inconsistently, so the scan picks
        // a set that is disconnected or has a pinched boundary. Re-fanning that
        // produces overlapping or missing triangles: 16 points on a circle gave 22
        // triangles where only 14 can exist.
        //
        // Growing outward from one seed cannot select a disconnected cavity, which
        // removes that failure mode by construction rather than by tightening a
        // tolerance.
        let mut edge_tris: HashMap<(usize, usize), Vec<usize>> = HashMap::new();
        for (t, &[a, b, c]) in tris.iter().enumerate() {
            for (u, v) in [(a, b), (b, c), (c, a)] {
                edge_tris.entry((u.min(v), u.max(v))).or_default().push(t);
            }
        }
        let contains_p = |t: usize, tris: &[[usize; 3]]| {
            let [a, b, c] = tris[t];
            in_circumcircle(all[a], all[b], all[c], p)
        };

        // Seed on the triangle actually containing `p` where there is one; it is
        // always in the cavity, and starting there keeps the growth centred.
        let seed = (0..tris.len())
            .find(|&t| {
                let [a, b, c] = tris[t];
                orient(all[a], all[b], p) >= 0.0
                    && orient(all[b], all[c], p) >= 0.0
                    && orient(all[c], all[a], p) >= 0.0
            })
            .or_else(|| (0..tris.len()).find(|&t| contains_p(t, &tris)));
        let Some(seed) = seed else { continue };

        let mut in_cavity = vec![false; tris.len()];
        in_cavity[seed] = true;
        let mut queue = vec![seed];
        while let Some(t) = queue.pop() {
            let [a, b, c] = tris[t];
            for (u, v) in [(a, b), (b, c), (c, a)] {
                let key = (u.min(v), u.max(v));
                for &nt in edge_tris.get(&key).map_or(&[][..], Vec::as_slice) {
                    if !in_cavity[nt] && contains_p(nt, &tris) {
                        in_cavity[nt] = true;
                        queue.push(nt);
                    }
                }
            }
        }

        // The boundary, as *directed* edges taken from the cavity triangles' own
        // winding. Direction matters: it is what makes the replacement fan inherit
        // the cavity's orientation instead of guessing it per triangle.
        //
        // Then shrink until every boundary edge is visible from `p`. A boundary
        // edge `p` is not strictly left of would produce an inverted or degenerate
        // triangle, which means the cavity is not star-shaped there; dropping the
        // triangle that owns the edge restores it. Each pass removes at least one
        // triangle, so this terminates.
        let boundary = loop {
            let mut directed: HashMap<(usize, usize), ()> = HashMap::new();
            for (t, &[a, b, c]) in tris.iter().enumerate() {
                if !in_cavity[t] {
                    continue;
                }
                for (u, v) in [(a, b), (b, c), (c, a)] {
                    directed.insert((u, v), ());
                }
            }
            // An edge is interior to the cavity when its reverse is also present.
            let mut edges: Vec<(usize, usize)> = directed
                .keys()
                .copied()
                .filter(|&(u, v)| !directed.contains_key(&(v, u)))
                .collect();
            edges.sort_unstable();

            let offender = edges
                .iter()
                .copied()
                .find(|&(u, v)| orient(all[u], all[v], p) <= 0.0);
            match offender {
                None => break edges,
                Some((u, v)) => {
                    // Drop the cavity triangle owning this directed edge.
                    let owner = (0..tris.len()).find(|&t| {
                        if !in_cavity[t] {
                            return false;
                        }
                        let [a, b, c] = tris[t];
                        [(a, b), (b, c), (c, a)].contains(&(u, v))
                    });
                    match owner {
                        Some(t) if t != seed => in_cavity[t] = false,
                        // The seed itself cannot be dropped without emptying the
                        // cavity, so give up on this point rather than corrupt the
                        // mesh. It stays out of the triangulation, which the Euler
                        // count below would reveal.
                        _ => break Vec::new(),
                    }
                }
            }
        };

        if boundary.is_empty() {
            continue;
        }

        // Remove the cavity, highest index first so earlier indices stay valid.
        let mut doomed: Vec<usize> = (0..tris.len()).filter(|&t| in_cavity[t]).collect();
        doomed.sort_unstable_by(|x, y| y.cmp(x));
        for t in doomed {
            tris.swap_remove(t);
        }
        // Re-fan: one triangle per boundary edge, inheriting its direction, so
        // every new triangle is counter-clockwise without a per-triangle test.
        for (u, v) in boundary {
            tris.push([u, v, i]);
        }
    }

    // Drop everything touching the super-triangle, and un-translate.
    let triangles: Vec<[usize; 3]> = tris
        .into_iter()
        .filter(|t| t.iter().all(|&v| v < s0))
        .collect();

    Triangulation {
        points: pts,
        triangles,
    }
}

/// The Voronoi diagram of a point set, built from its Delaunay dual.
///
/// Each cell is the circumcentres of the triangles incident to its site, ordered
/// by angle around the site. A site on the convex hull gets `None`, because its
/// true cell is unbounded — see [`Voronoi::cells`].
#[must_use]
pub fn voronoi(points: &[Coord]) -> Voronoi {
    let tri = delaunay(points);
    let sites = tri.points.clone();
    if tri.triangles.is_empty() {
        return Voronoi {
            cells: vec![None; sites.len()],
            sites,
        };
    }

    // Circumcentre per triangle, and the triangles around each site.
    let centres: Vec<Option<Coord>> = tri
        .triangles
        .iter()
        .map(|&[a, b, c]| {
            circumcircle(tri.points[a], tri.points[b], tri.points[c]).map(|(ctr, _)| ctr)
        })
        .collect();
    let mut incident: Vec<Vec<usize>> = vec![Vec::new(); sites.len()];
    for (t, &[a, b, c]) in tri.triangles.iter().enumerate() {
        for v in [a, b, c] {
            incident[v].push(t);
        }
    }

    // A site is on the hull when one of its incident edges belongs to only one
    // triangle. That edge's cell boundary runs off to infinity.
    let mut edge_use: HashMap<(usize, usize), usize> = HashMap::new();
    for &[a, b, c] in &tri.triangles {
        for (u, v) in [(a, b), (b, c), (c, a)] {
            *edge_use.entry((u.min(v), u.max(v))).or_insert(0) += 1;
        }
    }
    let on_hull: Vec<bool> = (0..sites.len())
        .map(|v| {
            edge_use
                .iter()
                .any(|(&(u, w), &n)| n == 1 && (u == v || w == v))
        })
        .collect();

    let cells = (0..sites.len())
        .map(|v| {
            if on_hull[v] || incident[v].len() < 3 {
                return None;
            }
            let site = sites[v];
            let mut ring: Vec<Coord> = incident[v].iter().filter_map(|&t| centres[t]).collect();
            if ring.len() < 3 {
                return None;
            }
            // Order by angle about the site. The incident triangles arrive in
            // insertion order, which is not the order round the cell.
            ring.sort_by(|p, q| {
                (p.y - site.y)
                    .atan2(p.x - site.x)
                    .total_cmp(&(q.y - site.y).atan2(q.x - site.x))
            });
            ring.dedup();
            if ring.len() < 3 {
                return None;
            }
            let mut closed = ring;
            closed.push(closed[0]);
            Some(Polygon::new(LineString::new(closed), vec![]))
        })
        .collect();

    Voronoi { sites, cells }
}

/// Delaunay-triangulate and build the Voronoi diagram of a geometry's coordinates.
pub trait Triangulate {
    /// The Delaunay triangulation of this geometry's vertices.
    fn delaunay_triangulation(&self) -> Triangulation;
    /// The Voronoi diagram of this geometry's vertices.
    fn voronoi_diagram(&self) -> Voronoi;
}

impl<T: crate::space::r2::CoordsIter + ?Sized> Triangulate for T {
    fn delaunay_triangulation(&self) -> Triangulation {
        delaunay(&self.to_coords())
    }
    fn voronoi_diagram(&self) -> Voronoi {
        voronoi(&self.to_coords())
    }
}

/// The Delaunay triangulation of a geometry, as geometry — the convenience form
/// of `ST_DelaunayTriangles`.
#[must_use]
pub fn delaunay_triangles(geom: &Geometry) -> MultiPolygon {
    delaunay(&crate::space::r2::CoordsIter::to_coords(geom)).to_multi_polygon()
}

/// The Voronoi cells of a geometry's vertices that closed, as geometry.
///
/// Hull sites are absent, for the reason given on [`Voronoi::cells`]. Use
/// [`Voronoi::clipped_to`] when every site needs a cell.
#[must_use]
pub fn voronoi_polygons(geom: &Geometry) -> MultiPolygon {
    voronoi(&crate::space::r2::CoordsIter::to_coords(geom)).to_multi_polygon()
}

/// The site nearest to `p`, by walking the Voronoi cells.
///
/// Here because it is the question a Voronoi diagram is usually built to answer,
/// and doing it from the cells rather than by scanning the sites is the whole
/// point — though with cells already in hand a caller will normally put them in
/// an STR-tree instead.
#[must_use]
pub fn nearest_site(v: &Voronoi, p: Coord) -> Option<usize> {
    use crate::space::r2::Contains;
    for (i, cell) in v.cells.iter().enumerate() {
        if let Some(poly) = cell
            && poly.contains(&Point(p))
        {
            return Some(i);
        }
    }
    // Outside every bounded cell — the honest fallback is the nearest site by
    // distance, which is also the correct answer there.
    v.sites
        .iter()
        .enumerate()
        .min_by(|(_, a), (_, b)| {
            let da = (a.x - p.x).powi(2) + (a.y - p.y).powi(2);
            let db = (b.x - p.x).powi(2) + (b.y - p.y).powi(2);
            da.total_cmp(&db)
        })
        .map(|(i, _)| i)
}
