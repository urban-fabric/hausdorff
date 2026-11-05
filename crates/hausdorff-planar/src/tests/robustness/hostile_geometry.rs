//! Every public planar operation, against geometry that is degenerate,
//! non-finite, or extreme — asserting only that the process survives.
//!
//! The companion to [`codec_hostile`](super::codec_hostile), which does this for
//! the parsers. A caller can hand these algorithms a `NaN` coordinate or a
//! two-point "polygon" without reading any documentation, so the contract that
//! matters first is *no panic*: return an empty result, a `None`, or a garbage
//! number, but do not unwind through the caller. Correctness on such input is a
//! separate question this suite deliberately does not ask.

use crate::space::r2::validation::Validation;
use crate::space::r2::*;
use std::panic::{AssertUnwindSafe, catch_unwind};

/// Run `f`, reporting whether it panicked, with the panic message suppressed so
/// the test output stays readable.
fn survives(f: impl FnOnce()) -> bool {
    let prev = std::panic::take_hook();
    std::panic::set_hook(Box::new(|_| {}));
    let ok = catch_unwind(AssertUnwindSafe(f)).is_ok();
    std::panic::set_hook(prev);
    ok
}

/// A named operation: a label and a thunk that exercises it once.
type Op = (&'static str, Box<dyn FnOnce()>);

fn ring(pts: &[(f64, f64)]) -> LineString {
    LineString::from(pts.to_vec())
}

/// The inputs: each is something a caller can construct without error.
fn hostile() -> Vec<(&'static str, Polygon, LineString, Vec<Coord>)> {
    let n = f64::NAN;
    let i = f64::INFINITY;
    let cases: Vec<(&str, Vec<(f64, f64)>)> = vec![
        ("empty", vec![]),
        ("one point", vec![(1.0, 1.0)]),
        ("two points", vec![(0.0, 0.0), (1.0, 1.0)]),
        ("all identical", vec![(2.0, 2.0); 5]),
        (
            "collinear",
            vec![(0.0, 0.0), (1.0, 0.0), (2.0, 0.0), (3.0, 0.0)],
        ),
        (
            "one NaN",
            vec![(0.0, 0.0), (1.0, 0.0), (n, 1.0), (0.0, 1.0)],
        ),
        ("all NaN", vec![(n, n), (n, n), (n, n)]),
        ("one infinity", vec![(0.0, 0.0), (i, 0.0), (1.0, 1.0)]),
        ("mixed inf", vec![(-i, -i), (i, i), (0.0, 0.0)]),
        ("huge", vec![(1e300, 1e300), (-1e300, 1e300), (0.0, -1e300)]),
        ("subnormal", vec![(0.0, 0.0), (5e-324, 0.0), (0.0, 5e-324)]),
        (
            "zero-area sliver",
            vec![(0.0, 0.0), (1.0, 1e-300), (2.0, 0.0)],
        ),
    ];
    cases
        .into_iter()
        .map(|(name, pts)| {
            let coords: Vec<Coord> = pts.iter().map(|&(x, y)| Coord::new(x, y)).collect();
            (name, Polygon::new(ring(&pts), vec![]), ring(&pts), coords)
        })
        .collect()
}

#[test]
fn no_planar_operation_panics_on_hostile_input() {
    let mut failures: Vec<String> = Vec::new();

    for (name, poly, line, coords) in hostile() {
        let mp = MultiPolygon::new(vec![poly.clone()]);
        let geom = Geometry::Polygon(poly.clone());
        let pts = MultiPoint::new(coords.iter().map(|c| Point(*c)).collect());

        // (label, operation) — every public entry point that takes geometry.
        let ops: Vec<Op> = vec![
            (
                "area",
                Box::new({
                    let p = poly.clone();
                    move || {
                        p.unsigned_area();
                    }
                }),
            ),
            (
                "signed_area",
                Box::new({
                    let p = poly.clone();
                    move || {
                        p.signed_area();
                    }
                }),
            ),
            (
                "centroid",
                Box::new({
                    let p = poly.clone();
                    move || {
                        p.centroid();
                    }
                }),
            ),
            (
                "bounding_rect",
                Box::new({
                    let p = poly.clone();
                    move || {
                        p.bounding_rect();
                    }
                }),
            ),
            (
                "interior_point",
                Box::new({
                    let p = poly.clone();
                    move || {
                        p.interior_point();
                    }
                }),
            ),
            (
                "convex_hull",
                Box::new({
                    let p = poly.clone();
                    move || {
                        p.convex_hull();
                    }
                }),
            ),
            (
                "concave_hull",
                Box::new({
                    let p = pts.clone();
                    move || {
                        p.concave_hull(2.0);
                    }
                }),
            ),
            (
                "min_bounding_circle",
                Box::new({
                    let p = pts.clone();
                    move || {
                        p.minimum_bounding_circle();
                    }
                }),
            ),
            (
                "oriented_envelope",
                Box::new({
                    let p = pts.clone();
                    move || {
                        p.oriented_envelope();
                    }
                }),
            ),
            (
                "inscribed_circle",
                Box::new({
                    let p = poly.clone();
                    move || {
                        p.maximum_inscribed_circle(1e-3);
                    }
                }),
            ),
            (
                "simplify",
                Box::new({
                    let p = poly.clone();
                    move || {
                        p.simplify(0.1);
                    }
                }),
            ),
            (
                "simplify_topology",
                Box::new({
                    let m = mp.clone();
                    move || {
                        m.simplify_topology_preserving(0.1);
                    }
                }),
            ),
            (
                "densify",
                Box::new({
                    let l = line.clone();
                    move || {
                        let _ = Euclidean.densify(&l, 0.5);
                    }
                }),
            ),
            (
                "remove_repeated",
                Box::new({
                    let p = poly.clone();
                    move || {
                        p.remove_repeated_points();
                    }
                }),
            ),
            (
                "make_valid",
                Box::new({
                    let p = poly.clone();
                    move || {
                        let _ = p.make_valid();
                    }
                }),
            ),
            (
                "polygonize",
                Box::new({
                    let l = line.clone();
                    move || {
                        [l].polygonize();
                    }
                }),
            ),
            (
                "delaunay",
                Box::new({
                    let c = coords.clone();
                    move || {
                        let _ = delaunay(&c);
                    }
                }),
            ),
            (
                "voronoi",
                Box::new({
                    let c = coords.clone();
                    move || {
                        let _ = voronoi(&c);
                    }
                }),
            ),
            (
                "nearest_site",
                Box::new({
                    let c = coords.clone();
                    move || {
                        let _ = nearest_site(&voronoi(&c), Coord::new(0.5, 0.5));
                    }
                }),
            ),
            (
                "buffer",
                Box::new({
                    let g = geom.clone();
                    move || {
                        let _ = buffer(&g, 1.0, BufferStyle::default());
                    }
                }),
            ),
            (
                "validation",
                Box::new({
                    let p = poly.clone();
                    move || {
                        let _ = p.validation_errors();
                    }
                }),
            ),
            (
                "is_ccw",
                Box::new({
                    let l = line.clone();
                    move || {
                        l.is_ccw();
                    }
                }),
            ),
            (
                "extremes",
                Box::new({
                    let p = poly.clone();
                    move || {
                        p.extremes();
                    }
                }),
            ),
            (
                "closest_point",
                Box::new({
                    let l = line.clone();
                    move || {
                        l.closest_point(&Point::new(0.5, 0.5));
                    }
                }),
            ),
            (
                "line_locate",
                Box::new({
                    let l = line.clone();
                    move || {
                        l.line_locate_point(&Point::new(0.5, 0.5));
                    }
                }),
            ),
            (
                "line_interpolate",
                Box::new({
                    let l = line.clone();
                    move || {
                        l.line_interpolate_point(0.5);
                    }
                }),
            ),
            (
                "contains",
                Box::new({
                    let p = poly.clone();
                    move || {
                        p.contains(&Point::new(0.5, 0.5));
                    }
                }),
            ),
            (
                "intersects",
                Box::new({
                    let g = geom.clone();
                    move || {
                        g.intersects(&g);
                    }
                }),
            ),
            (
                "relate",
                Box::new({
                    let g = geom.clone();
                    move || {
                        let _ = relate(&g, &g);
                    }
                }),
            ),
            (
                "prepared",
                Box::new({
                    let g = geom.clone();
                    move || {
                        let pg = PreparedGeometry::new(&g);
                        let _ = pg.intersects(&g);
                    }
                }),
            ),
            (
                "union",
                Box::new({
                    let m = mp.clone();
                    move || {
                        m.union(&m);
                    }
                }),
            ),
            (
                "intersection",
                Box::new({
                    let m = mp.clone();
                    move || {
                        m.intersection(&m);
                    }
                }),
            ),
            (
                "difference",
                Box::new({
                    let m = mp.clone();
                    move || {
                        m.difference(&m);
                    }
                }),
            ),
            (
                "unary_union",
                Box::new({
                    let m = mp.clone();
                    move || {
                        unary_union(m.iter());
                    }
                }),
            ),
            (
                "hausdorff",
                Box::new({
                    let g = geom.clone();
                    move || {
                        let _ = hausdorff_distance(&Euclidean, &g, &g);
                    }
                }),
            ),
            (
                "frechet",
                Box::new({
                    let l = line.clone();
                    move || {
                        let g = Geometry::LineString(l);
                        let _ = frechet_distance(&Euclidean, &g, &g);
                    }
                }),
            ),
            (
                "euclidean_dist",
                Box::new({
                    let g = geom.clone();
                    move || {
                        let _ = euclidean_distance(&g, &g);
                    }
                }),
            ),
            (
                "str_tree",
                Box::new({
                    let g = geom.clone();
                    move || {
                        let r = g
                            .bounding_rect()
                            .unwrap_or(Rect::new((0.0, 0.0), (0.0, 0.0)));
                        let t: StrTree<u32> = StrTree::new(vec![(r, 0u32)]);
                        let _ = t.query(Rect::new((0.0, 0.0), (1.0, 1.0))).len();
                    }
                }),
            ),
            (
                "spatial_join",
                Box::new({
                    let g = geom.clone();
                    move || {
                        let l = [g.clone()];
                        let _ = spatial_join_candidates(&l, &l).len();
                        let _ = spatial_join(&l, &l, |a, b| a.intersects(b)).len();
                    }
                }),
            ),
            (
                "affine",
                Box::new({
                    let p = poly.clone();
                    move || {
                        let _ = p
                            .affine_transform(&AffineTransform::new(1.0, 0.5, 0.0, 0.25, 1.0, 0.0));
                    }
                }),
            ),
            (
                "map_coords",
                Box::new({
                    let p = poly.clone();
                    move || {
                        let _ = p.map_coords(|c| c);
                    }
                }),
            ),
            (
                "coords_count",
                Box::new({
                    let p = poly.clone();
                    move || {
                        p.coords_count();
                    }
                }),
            ),
            (
                "dimensions",
                Box::new({
                    let g = geom.clone();
                    move || {
                        g.dimensions();
                    }
                }),
            ),
        ];

        for (op, f) in ops {
            if !survives(f) {
                failures.push(format!("{op} panicked on \"{name}\""));
            }
        }
    }

    assert!(
        failures.is_empty(),
        "{} public operations panicked on input a caller can construct:\n  {}",
        failures.len(),
        failures.join("\n  ")
    );
}
