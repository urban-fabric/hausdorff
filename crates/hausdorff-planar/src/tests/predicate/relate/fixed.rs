use super::*;

#[test]
fn polygon_polygon() {
    let a = sq(0.0, 0.0, 4.0, 4.0);
    check2(&a, &sq(2.0, 2.0, 6.0, 6.0));
    check2(&a, &sq(1.0, 1.0, 2.0, 2.0));
    check2(&a, &sq(10.0, 10.0, 12.0, 12.0));
    check2(&a, &sq(4.0, 0.0, 8.0, 4.0));
    check2(&a, &sq(4.0, 4.0, 8.0, 8.0));
    check2(&a, &sq(0.0, 0.0, 4.0, 4.0));
}

#[test]
fn point_polygon() {
    let a = sq(0.0, 0.0, 4.0, 4.0);
    check2(&a, &pt(2.0, 2.0));
    check2(&a, &pt(4.0, 2.0));
    check2(&a, &pt(0.0, 0.0));
    check2(&a, &pt(9.0, 9.0));
}

#[test]
fn point_point() {
    check2(&pt(1.0, 1.0), &pt(1.0, 1.0));
    check2(&pt(1.0, 1.0), &pt(2.0, 2.0));
}

#[test]
fn line_line() {
    check2(&line(0.0, 0.0, 4.0, 4.0), &line(0.0, 4.0, 4.0, 0.0));
    check2(&line(0.0, 0.0, 4.0, 0.0), &line(4.0, 0.0, 8.0, 0.0));
    check2(&line(0.0, 0.0, 4.0, 0.0), &line(2.0, 0.0, 6.0, 0.0));
    check2(&line(0.0, 0.0, 4.0, 0.0), &line(0.0, 1.0, 4.0, 1.0));
    check2(&line(0.0, 0.0, 4.0, 4.0), &line(0.0, 0.0, 4.0, 4.0));
    check2(&line(0.0, 0.0, 4.0, 0.0), &line(2.0, 0.0, 2.0, 4.0));
}

#[test]
fn line_polygon() {
    let a = sq(0.0, 0.0, 4.0, 4.0);
    check2(&a, &line(2.0, 2.0, 3.0, 3.0));
    check2(&a, &line(-1.0, 2.0, 5.0, 2.0));
    check2(&a, &line(4.0, 0.0, 4.0, 4.0));
    check2(&a, &line(5.0, 0.0, 9.0, 4.0));
    check2(&a, &line(2.0, 2.0, 6.0, 2.0));
}

#[test]
fn point_line() {
    let l = line(0.0, 0.0, 4.0, 0.0);
    check2(&l, &pt(2.0, 0.0));
    check2(&l, &pt(0.0, 0.0));
    check2(&l, &pt(2.0, 1.0));
}

#[test]
fn polygon_with_hole() {
    let outer = [
        (0.0, 0.0),
        (10.0, 0.0),
        (10.0, 10.0),
        (0.0, 10.0),
        (0.0, 0.0),
    ];
    let hole = [(3.0, 3.0), (7.0, 3.0), (7.0, 7.0), (3.0, 7.0), (3.0, 3.0)];
    let donut = Geometry::Polygon(Polygon::new(
        LineString::from(outer.to_vec()),
        vec![LineString::from(hole.to_vec())],
    ));
    check2(&donut, &pt(5.0, 5.0));
    check2(&donut, &pt(1.0, 1.0));
    check2(&donut, &sq(4.0, 4.0, 6.0, 6.0));
}

#[test]
fn adversarial_areal_relate_matches_geo() {
    let holed = |shell: &[(f64, f64)], hole: &[(f64, f64)]| {
        Geometry::Polygon(Polygon::new(
            LineString::from(shell.to_vec()),
            vec![LineString::from(hole.to_vec())],
        ))
    };
    let shell = [
        (0.0, 0.0),
        (10.0, 0.0),
        (10.0, 10.0),
        (0.0, 10.0),
        (0.0, 0.0),
    ];
    let hole = [(3.0, 3.0), (7.0, 3.0), (7.0, 7.0), (3.0, 7.0), (3.0, 3.0)];
    let donut = holed(&shell, &hole);

    check2(&sq(0.0, 0.0, 10.0, 10.0), &sq(3.0, 3.0, 7.0, 7.0));
    check2(&sq(0.0, 0.0, 10.0, 1.0), &sq(4.0, -5.0, 5.0, 5.0));
    check2(&donut, &sq(3.0, 3.0, 7.0, 7.0));
    check2(&donut, &sq(5.0, 5.0, 12.0, 12.0));
    check2(&donut, &sq(0.5, 0.5, 2.5, 2.5));
    check2(&donut, &holed(&shell, &hole));
    let donut2 = holed(
        &[
            (5.0, 5.0),
            (15.0, 5.0),
            (15.0, 15.0),
            (5.0, 15.0),
            (5.0, 5.0),
        ],
        &[
            (8.0, 8.0),
            (12.0, 8.0),
            (12.0, 12.0),
            (8.0, 12.0),
            (8.0, 8.0),
        ],
    );
    check2(&donut, &donut2);
}

#[test]
fn multipolygon_cases() {
    let mp = Geometry::MultiPolygon(MultiPolygon::new(vec![
        Polygon::new(
            LineString::from(vec![
                (0.0, 0.0),
                (2.0, 0.0),
                (2.0, 2.0),
                (0.0, 2.0),
                (0.0, 0.0),
            ]),
            vec![],
        ),
        Polygon::new(
            LineString::from(vec![
                (5.0, 5.0),
                (7.0, 5.0),
                (7.0, 7.0),
                (5.0, 7.0),
                (5.0, 5.0),
            ]),
            vec![],
        ),
    ]));
    check2(&mp, &pt(1.0, 1.0));
    check2(&mp, &pt(6.0, 6.0));
    check2(&mp, &pt(3.0, 3.0));
    check2(&mp, &sq(1.0, 1.0, 6.0, 6.0));
}
