//! Conformance tests for the planar query operations added alongside the core
//! algorithm surface — `ClosestPoint`, linear referencing
//! (`LineInterpolatePoint` / `LineLocatePoint`), the affine transform stack
//! (`AffineTransform` + `AffineOps` / `Translate` / `Rotate` / `Scale` /
//! `Skew`), and `Extremes` — cross-checked against georust `geo` exactly as
//! `tests/geo_conformance.rs` does for the rest of the surface.

use geo::line_locate_point::LineLocatePoint as _;
use geo::{
    AffineOps as _, ClosestPoint as _, Extremes as _, Rotate as _, Scale as _, Skew as _,
    Translate as _,
};

use hausdorff::space::r2::{
    AffineOps, AffineTransform, Closest, ClosestPoint, CoordsIter, Extremes, LineInterpolatePoint,
    LineLocatePoint, Rotate, Scale, Skew, Translate,
};
use hausdorff::{Coord, Geometry, Line, LineString, MultiPoint, Point, Polygon, Rect, Triangle};

const EPS: f64 = 1e-9;

fn h_line_string(coords: &[(f64, f64)]) -> LineString {
    LineString::new(coords.iter().map(|&(x, y)| Coord::new(x, y)).collect())
}

fn g_line_string(coords: &[(f64, f64)]) -> geo::LineString<f64> {
    geo::LineString(coords.iter().map(|&(x, y)| geo::Coord { x, y }).collect())
}

fn h_polygon(exterior: &[(f64, f64)], holes: &[Vec<(f64, f64)>]) -> Polygon {
    Polygon::new(
        h_line_string(exterior),
        holes.iter().map(|h| h_line_string(h)).collect(),
    )
}

fn g_polygon(exterior: &[(f64, f64)], holes: &[Vec<(f64, f64)>]) -> geo::Polygon<f64> {
    geo::Polygon::new(
        g_line_string(exterior),
        holes.iter().map(|h| g_line_string(h)).collect(),
    )
}

fn assert_points_eq(h: Point, g: geo::Point<f64>, ctx: &str) {
    assert!(
        (h.x() - g.x()).abs() < EPS && (h.y() - g.y()).abs() < EPS,
        "{ctx}: ({}, {}) vs geo ({}, {})",
        h.x(),
        h.y(),
        g.x(),
        g.y()
    );
}

fn assert_closest_eq(h: Closest, g: geo::Closest<f64>, ctx: &str) {
    match (h, g) {
        (Closest::Intersection(hp), geo::Closest::Intersection(gp)) => {
            assert_points_eq(hp, gp, ctx)
        }
        (Closest::SinglePoint(hp), geo::Closest::SinglePoint(gp)) => assert_points_eq(hp, gp, ctx),
        (Closest::Indeterminate, geo::Closest::Indeterminate) => {}
        (h, g) => panic!("{ctx}: variant mismatch {h:?} vs geo {g:?}"),
    }
}

fn query_points() -> Vec<(f64, f64)> {
    vec![
        (0.0, 0.0),
        (4.0, 0.0),
        (2.0, 0.0),
        (2.0, 2.0),
        (5.0, 5.0),
        (-3.0, 1.0),
        (10.0, -10.0),
        (0.5, 3.75),
        (4.0, 4.0),
        (100.0, 0.0),
    ]
}

#[test]
fn closest_point_on_line_matches_geo() {
    let segments = [
        ((0.0, 0.0), (4.0, 0.0)),
        ((0.0, 0.0), (4.0, 4.0)),
        ((-2.0, 1.0), (3.0, -1.0)),
        ((1.0, 1.0), (1.0, 1.0)),
    ];
    for (start, end) in segments {
        let h = Line::new(Coord::new(start.0, start.1), Coord::new(end.0, end.1));
        let g = geo::Line::new(
            geo::Coord {
                x: start.0,
                y: start.1,
            },
            geo::Coord { x: end.0, y: end.1 },
        );
        for (x, y) in query_points() {
            assert_closest_eq(
                h.closest_point(&Point::new(x, y)),
                g.closest_point(&geo::Point::new(x, y)),
                &format!("line {start:?}-{end:?} vs ({x}, {y})"),
            );
        }
    }
}

#[test]
fn closest_point_on_point_and_multipoint_matches_geo() {
    let h_pt = Point::new(2.0, 2.0);
    let g_pt = geo::Point::new(2.0, 2.0);
    let pts = vec![(0.0, 0.0), (4.0, 0.0), (4.0, 4.0), (0.0, 4.0)];
    let h_mp = MultiPoint::from(pts.clone());
    let g_mp: geo::MultiPoint<f64> = pts.iter().map(|&(x, y)| geo::Point::new(x, y)).collect();
    for (x, y) in query_points() {
        let p = Point::new(x, y);
        let gp = geo::Point::new(x, y);
        assert_closest_eq(
            h_pt.closest_point(&p),
            g_pt.closest_point(&gp),
            &format!("point vs ({x}, {y})"),
        );
        assert_closest_eq(
            h_mp.closest_point(&p),
            g_mp.closest_point(&gp),
            &format!("multipoint vs ({x}, {y})"),
        );
    }
    assert_closest_eq(
        MultiPoint::new(vec![]).closest_point(&h_pt),
        geo::MultiPoint::<f64>::new(vec![]).closest_point(&g_pt),
        "empty multipoint",
    );
}

#[test]
fn closest_point_on_line_string_matches_geo() {
    let paths = [
        vec![(0.0, 0.0), (4.0, 0.0), (4.0, 4.0), (0.0, 4.0)],
        vec![(-3.0, -2.0), (1.0, 3.0), (5.0, -1.0), (8.0, 6.0)],
        vec![],
    ];
    for path in &paths {
        let h = h_line_string(path);
        let g = g_line_string(path);
        for (x, y) in query_points() {
            assert_closest_eq(
                h.closest_point(&Point::new(x, y)),
                g.closest_point(&geo::Point::new(x, y)),
                &format!("linestring {path:?} vs ({x}, {y})"),
            );
        }
    }
}

#[test]
fn closest_point_on_polygon_matches_geo() {
    let exterior = vec![(0.0, 0.0), (8.0, 0.0), (8.0, 8.0), (0.0, 8.0)];
    let hole = vec![(3.0, 3.0), (5.0, 3.0), (5.0, 5.0), (3.0, 5.0)];
    let cases: Vec<(Polygon, geo::Polygon<f64>)> = vec![
        (h_polygon(&exterior, &[]), g_polygon(&exterior, &[])),
        (
            h_polygon(&exterior, &[hole.clone()]),
            g_polygon(&exterior, &[hole.clone()]),
        ),
        (h_polygon(&[], &[]), g_polygon(&[], &[])),
    ];
    let queries = [
        (1.0, 1.0),
        (0.0, 0.0),
        (4.0, 0.0),
        (8.0, 8.0),
        (4.0, 3.5),
        (3.0, 4.0),
        (-2.0, 4.0),
        (12.0, 12.0),
        (4.0, -5.0),
    ];
    for (hp, gp) in &cases {
        for (x, y) in queries {
            assert_closest_eq(
                hp.closest_point(&Point::new(x, y)),
                gp.closest_point(&geo::Point::new(x, y)),
                &format!("polygon vs ({x}, {y})"),
            );
        }
    }
}

#[test]
fn closest_point_on_multi_geometries_matches_geo() {
    let part_a = vec![(0.0, 0.0), (2.0, 0.0), (2.0, 2.0), (0.0, 2.0)];
    let part_b = vec![(6.0, 6.0), (9.0, 6.0), (9.0, 9.0), (6.0, 9.0)];
    let h_mp = hausdorff::MultiPolygon::new(vec![h_polygon(&part_a, &[]), h_polygon(&part_b, &[])]);
    let g_mp = geo::MultiPolygon::new(vec![g_polygon(&part_a, &[]), g_polygon(&part_b, &[])]);

    let line_a = vec![(0.0, 0.0), (4.0, 0.0)];
    let line_b = vec![(0.0, 6.0), (4.0, 6.0)];
    let h_mls =
        hausdorff::MultiLineString::new(vec![h_line_string(&line_a), h_line_string(&line_b)]);
    let g_mls = geo::MultiLineString::new(vec![g_line_string(&line_a), g_line_string(&line_b)]);

    for (x, y) in query_points() {
        let p = Point::new(x, y);
        let gp = geo::Point::new(x, y);
        assert_closest_eq(
            h_mp.closest_point(&p),
            g_mp.closest_point(&gp),
            &format!("multipolygon vs ({x}, {y})"),
        );
        assert_closest_eq(
            h_mls.closest_point(&p),
            g_mls.closest_point(&gp),
            &format!("multilinestring vs ({x}, {y})"),
        );
    }
    assert_closest_eq(
        hausdorff::MultiPolygon::new(vec![]).closest_point(&Point::new(1.0, 1.0)),
        geo::MultiPolygon::<f64>::new(vec![]).closest_point(&geo::Point::new(1.0, 1.0)),
        "empty multipolygon",
    );
}

#[test]
fn closest_point_on_rect_triangle_and_geometry_matches_geo() {
    let h_rect = Rect::new(Coord::new(1.0, 1.0), Coord::new(5.0, 3.0));
    let g_rect = geo::Rect::new(geo::Coord { x: 1.0, y: 1.0 }, geo::Coord { x: 5.0, y: 3.0 });
    let h_tri = Triangle::new(
        Coord::new(0.0, 0.0),
        Coord::new(6.0, 0.0),
        Coord::new(0.0, 6.0),
    );
    let g_tri = geo::Triangle::new(
        geo::Coord { x: 0.0, y: 0.0 },
        geo::Coord { x: 6.0, y: 0.0 },
        geo::Coord { x: 0.0, y: 6.0 },
    );
    for (x, y) in query_points() {
        let p = Point::new(x, y);
        let gp = geo::Point::new(x, y);
        assert_closest_eq(
            h_rect.closest_point(&p),
            g_rect.closest_point(&gp),
            &format!("rect vs ({x}, {y})"),
        );
        assert_closest_eq(
            h_tri.closest_point(&p),
            g_tri.closest_point(&gp),
            &format!("triangle vs ({x}, {y})"),
        );
        assert_closest_eq(
            Geometry::Triangle(h_tri).closest_point(&p),
            geo::Geometry::Triangle(g_tri).closest_point(&gp),
            &format!("geometry vs ({x}, {y})"),
        );
    }
}

#[allow(deprecated)]
#[test]
fn line_interpolate_point_matches_geo() {
    use geo::LineInterpolatePoint as _;

    let fractions = [
        -1.0,
        0.0,
        0.1,
        0.25,
        0.5,
        0.75,
        0.9,
        1.0,
        2.0,
        f64::NEG_INFINITY,
        f64::INFINITY,
        f64::NAN,
    ];
    let h_seg = Line::new(Coord::new(-1.0, 0.0), Coord::new(3.0, 2.0));
    let g_seg = geo::Line::new(
        geo::Coord { x: -1.0, y: 0.0 },
        geo::Coord { x: 3.0, y: 2.0 },
    );
    let paths = [
        vec![(-1.0, 0.0), (0.0, 0.0), (0.0, 1.0)],
        vec![(0.0, 0.0), (3.0, 4.0), (3.0, 10.0), (-2.0, 10.0)],
        vec![(1.0, 1.0), (1.0, 1.0)],
        vec![],
    ];
    for &f in &fractions {
        let (h, g) = (
            h_seg.line_interpolate_point(f),
            g_seg.line_interpolate_point(f),
        );
        match (h, g) {
            (Some(hp), Some(gp)) => assert_points_eq(hp, gp, &format!("line interpolate {f}")),
            (None, None) => {}
            (h, g) => panic!("line interpolate {f}: {h:?} vs geo {g:?}"),
        }
        for path in &paths {
            let (h, g) = (
                h_line_string(path).line_interpolate_point(f),
                g_line_string(path).line_interpolate_point(f),
            );
            match (h, g) {
                (Some(hp), Some(gp)) => {
                    assert_points_eq(hp, gp, &format!("ls interpolate {path:?} at {f}"))
                }
                (None, None) => {}
                (h, g) => panic!("ls interpolate {path:?} at {f}: {h:?} vs geo {g:?}"),
            }
        }
    }
}

#[test]
fn line_locate_point_matches_geo() {
    let h_seg = Line::new(Coord::new(-1.0, 0.0), Coord::new(3.0, 2.0));
    let g_seg = geo::Line::new(
        geo::Coord { x: -1.0, y: 0.0 },
        geo::Coord { x: 3.0, y: 2.0 },
    );
    let h_zero = Line::new(Coord::new(1.0, 1.0), Coord::new(1.0, 1.0));
    let g_zero = geo::Line::new(geo::Coord { x: 1.0, y: 1.0 }, geo::Coord { x: 1.0, y: 1.0 });
    let path = vec![(0.0, 0.0), (3.0, 4.0), (3.0, 10.0), (-2.0, 10.0)];
    let h_ls = h_line_string(&path);
    let g_ls = g_line_string(&path);

    let mut queries = query_points();
    queries.push((f64::NAN, 1.0));
    for (x, y) in queries {
        let p = Point::new(x, y);
        let gp = geo::Point::new(x, y);
        for (name, h, g) in [
            (
                "segment",
                h_seg.line_locate_point(&p),
                g_seg.line_locate_point(&gp),
            ),
            (
                "zero",
                h_zero.line_locate_point(&p),
                g_zero.line_locate_point(&gp),
            ),
            (
                "linestring",
                h_ls.line_locate_point(&p),
                g_ls.line_locate_point(&gp),
            ),
        ] {
            match (h, g) {
                (Some(hf), Some(gf)) => assert!(
                    (hf - gf).abs() < EPS,
                    "{name} locate ({x}, {y}): {hf} vs geo {gf}"
                ),
                (None, None) => {}
                (h, g) => panic!("{name} locate ({x}, {y}): {h:?} vs geo {g:?}"),
            }
        }
    }

    assert_eq!(
        h_line_string(&[]).line_locate_point(&Point::new(1.0, 1.0)),
        g_line_string(&[]).line_locate_point(&geo::Point::new(1.0, 1.0)),
    );
}

#[allow(deprecated)]
#[test]
fn locate_interpolate_round_trip_matches_geo() {
    use geo::LineInterpolatePoint as _;

    let path = vec![(0.0, 0.0), (3.0, 4.0), (3.0, 10.0), (-2.0, 10.0)];
    let h_ls = h_line_string(&path);
    let g_ls = g_line_string(&path);
    for i in 0..=20 {
        let f = i as f64 / 20.0;
        let hp = h_ls.line_interpolate_point(f).unwrap();
        let gp = g_ls.line_interpolate_point(f).unwrap();
        assert_points_eq(hp, gp, &format!("round-trip interpolate {f}"));
        let hf = h_ls.line_locate_point(&hp).unwrap();
        let gf = g_ls.line_locate_point(&gp).unwrap();
        assert!(
            (hf - gf).abs() < EPS && (hf - f).abs() < 1e-12,
            "round-trip locate {f}: {hf} vs geo {gf}"
        );
    }
}

fn assert_line_strings_eq(h: &LineString, g: &geo::LineString<f64>, ctx: &str) {
    assert_eq!(h.0.len(), g.0.len(), "{ctx}: coord count");
    for (hc, gc) in h.0.iter().zip(g.0.iter()) {
        assert!(
            (hc.x - gc.x).abs() < EPS && (hc.y - gc.y).abs() < EPS,
            "{ctx}: ({}, {}) vs geo ({}, {})",
            hc.x,
            hc.y,
            gc.x,
            gc.y
        );
    }
}

fn assert_polygons_eq(h: &Polygon, g: &geo::Polygon<f64>, ctx: &str) {
    assert_line_strings_eq(h.exterior(), g.exterior(), ctx);
    assert_eq!(h.interiors().len(), g.interiors().len(), "{ctx}: holes");
    for (hr, gr) in h.interiors().iter().zip(g.interiors()) {
        assert_line_strings_eq(hr, gr, ctx);
    }
}

fn test_polygon_pair() -> (Polygon, geo::Polygon<f64>) {
    let exterior = vec![(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0)];
    let hole = vec![(3.0, 3.0), (5.0, 3.0), (5.0, 5.0), (3.0, 5.0)];
    (
        h_polygon(&exterior, &[hole.clone()]),
        g_polygon(&exterior, &[hole]),
    )
}

fn test_line_string_pair() -> (LineString, geo::LineString<f64>) {
    let path = vec![(0.0, 0.0), (5.0, 1.0), (7.0, 6.0), (-2.0, 4.0)];
    (h_line_string(&path), g_line_string(&path))
}

#[test]
fn affine_transform_application_matches_geo() {
    let transforms: Vec<(AffineTransform, geo::AffineTransform<f64>)> = vec![
        (
            AffineTransform::identity(),
            geo::AffineTransform::identity(),
        ),
        (
            AffineTransform::translate(3.5, -2.25),
            geo::AffineTransform::translate(3.5, -2.25),
        ),
        (
            AffineTransform::scale(2.0, -0.5, Coord::new(1.0, 2.0)),
            geo::AffineTransform::scale(2.0, -0.5, geo::Coord { x: 1.0, y: 2.0 }),
        ),
        (
            AffineTransform::rotate(37.0, Coord::new(-1.0, 4.0)),
            geo::AffineTransform::rotate(37.0, geo::Coord { x: -1.0, y: 4.0 }),
        ),
        (
            AffineTransform::skew(15.0, -30.0, Coord::new(2.0, 2.0)),
            geo::AffineTransform::skew(15.0, -30.0, geo::Coord { x: 2.0, y: 2.0 }),
        ),
        (
            AffineTransform::translate(1.0, 2.0)
                .rotated(45.0, Coord::new(0.0, 0.0))
                .scaled(3.0, 0.25, Coord::new(5.0, 5.0))
                .skewed(10.0, 20.0, Coord::new(-2.0, 3.0)),
            geo::AffineTransform::translate(1.0, 2.0)
                .rotated(45.0, geo::Coord { x: 0.0, y: 0.0 })
                .scaled(3.0, 0.25, geo::Coord { x: 5.0, y: 5.0 })
                .skewed(10.0, 20.0, geo::Coord { x: -2.0, y: 3.0 }),
        ),
    ];
    let (h_poly, g_poly) = test_polygon_pair();
    let (h_ls, g_ls) = test_line_string_pair();
    for (i, (ht, gt)) in transforms.iter().enumerate() {
        assert_eq!(
            [ht.a(), ht.b(), ht.xoff(), ht.d(), ht.e(), ht.yoff()],
            [gt.a(), gt.b(), gt.xoff(), gt.d(), gt.e(), gt.yoff()],
            "transform {i}: matrix"
        );
        assert_polygons_eq(
            &h_poly.affine_transform(ht),
            &g_poly.affine_transform(gt),
            &format!("transform {i} on polygon"),
        );
        assert_line_strings_eq(
            &h_ls.affine_transform(ht),
            &g_ls.affine_transform(gt),
            &format!("transform {i} on linestring"),
        );
        let mut h_mut = h_ls.clone();
        h_mut.affine_transform_mut(ht);
        assert_eq!(h_mut, h_ls.affine_transform(ht), "transform {i} mut");
    }
}

#[test]
fn translate_rotate_scale_skew_match_geo() {
    let (h_poly, g_poly) = test_polygon_pair();
    let (h_ls, g_ls) = test_line_string_pair();

    assert_polygons_eq(
        &h_poly.translate(17.0, -3.0),
        &g_poly.translate(17.0, -3.0),
        "translate polygon",
    );
    assert_line_strings_eq(
        &h_ls.translate(-0.5, 9.0),
        &g_ls.translate(-0.5, 9.0),
        "translate linestring",
    );

    assert_polygons_eq(
        &h_poly.rotate_around_point(-45.0, Point::new(10.0, 0.0)),
        &g_poly.rotate_around_point(-45.0, geo::Point::new(10.0, 0.0)),
        "rotate around point",
    );
    assert_polygons_eq(
        &h_poly.rotate_around_centroid(33.0),
        &g_poly.rotate_around_centroid(33.0),
        "rotate around centroid",
    );
    assert_polygons_eq(
        &h_poly.rotate_around_center(-120.0),
        &g_poly.rotate_around_center(-120.0),
        "rotate around center",
    );
    assert_line_strings_eq(
        &h_ls.rotate_around_centroid(75.0),
        &g_ls.rotate_around_centroid(75.0),
        "rotate linestring around centroid",
    );

    assert_polygons_eq(&h_poly.scale(2.5), &g_poly.scale(2.5), "scale");
    assert_polygons_eq(
        &h_poly.scale_xy(2.0, -4.0),
        &g_poly.scale_xy(2.0, -4.0),
        "scale_xy",
    );
    assert_line_strings_eq(
        &h_ls.scale_around_point(2.0, 4.0, Coord::new(100.0, 100.0)),
        &g_ls.scale_around_point(2.0, 4.0, geo::Coord { x: 100.0, y: 100.0 }),
        "scale_around_point",
    );

    assert_polygons_eq(&h_poly.skew(30.0), &g_poly.skew(30.0), "skew");
    assert_polygons_eq(
        &h_poly.skew_xy(30.0, 12.0),
        &g_poly.skew_xy(30.0, 12.0),
        "skew_xy",
    );
    assert_line_strings_eq(
        &h_ls.skew_around_point(45.0, 10.0, Coord::new(2.0, 2.0)),
        &g_ls.skew_around_point(45.0, 10.0, geo::Coord { x: 2.0, y: 2.0 }),
        "skew_around_point",
    );

    let mut h_mut = h_poly.clone();
    h_mut.rotate_around_centroid_mut(33.0);
    assert_eq!(h_mut, h_poly.rotate_around_centroid(33.0), "rotate mut");
    let mut h_mut = h_ls.clone();
    h_mut.scale_mut(2.5);
    assert_eq!(h_mut, h_ls.scale(2.5), "scale mut");
    let mut h_mut = h_ls.clone();
    h_mut.skew_xy_mut(30.0, 12.0);
    assert_eq!(h_mut, h_ls.skew_xy(30.0, 12.0), "skew mut");
    let mut h_mut = h_ls.clone();
    h_mut.translate_mut(1.0, 2.0);
    assert_eq!(h_mut, h_ls.translate(1.0, 2.0), "translate mut");
}

#[test]
fn extremes_match_geo() {
    let diamond = vec![(1.0, 0.0), (2.0, 1.0), (1.0, 2.0), (0.0, 1.0), (1.0, 0.0)];
    let irregular = vec![
        (1.0, -3.1),
        (3.7, 2.7),
        (0.9, 7.6),
        (-4.8, 6.7),
        (-7.5, 0.9),
    ];

    let h_cases: Vec<Geometry> = vec![
        Geometry::Polygon(h_polygon(&diamond, &[])),
        Geometry::Polygon(h_polygon(&irregular, &[])),
        Geometry::LineString(h_line_string(&irregular)),
        Geometry::MultiPoint(MultiPoint::from(irregular.clone())),
    ];
    let g_cases: Vec<geo::Geometry<f64>> = vec![
        geo::Geometry::Polygon(g_polygon(&diamond, &[])),
        geo::Geometry::Polygon(g_polygon(&irregular, &[])),
        geo::Geometry::LineString(g_line_string(&irregular)),
        geo::Geometry::MultiPoint(
            irregular
                .iter()
                .map(|&(x, y)| geo::Point::new(x, y))
                .collect(),
        ),
    ];

    for (h, g) in h_cases.iter().zip(&g_cases) {
        let ho = h.extremes().unwrap();
        let go = g.extremes().unwrap();
        for (name, he, ge) in [
            ("x_min", ho.x_min, &go.x_min),
            ("y_min", ho.y_min, &go.y_min),
            ("x_max", ho.x_max, &go.x_max),
            ("y_max", ho.y_max, &go.y_max),
        ] {
            assert_eq!(he.index, ge.index, "{name} index on {h:?}");
            assert!(
                (he.coord.x - ge.coord.x).abs() < EPS && (he.coord.y - ge.coord.y).abs() < EPS,
                "{name} coord on {h:?}: {:?} vs geo {:?}",
                he.coord,
                ge.coord
            );
        }
    }

    assert!(MultiPoint::new(vec![]).extremes().is_none());
    assert!(geo::MultiPoint::<f64>::new(vec![]).extremes().is_none());
}

#[test]
fn extremes_ignore_interior_rings_like_geo() {
    let exterior = vec![(0.0, 0.0), (4.0, 0.0), (4.0, 4.0), (0.0, 4.0)];
    let hole = vec![(1.0, 1.0), (3.0, 1.0), (3.0, 3.0), (1.0, 3.0)];
    let h = h_polygon(&exterior, &[hole.clone()]);
    let g = g_polygon(&exterior, &[hole]);
    let ho = h.extremes().unwrap();
    let go = g.extremes().unwrap();
    assert_eq!(ho.x_max.index, go.x_max.index);
    assert_eq!(ho.x_max.coord.x, go.x_max.coord.x);
    assert!(ho.x_max.index < h.exterior().coords_count());
}
