use super::*;
use crate::fixed::Lcg;

fn rrect(r: &mut Lcg, g: i64) -> Geometry {
    let (x, y) = (r.n(g), r.n(g));
    let (w, h) = (1.0 + r.n(g / 2), 1.0 + r.n(g / 2));
    sq(x, y, x + w, y + h)
}
fn rline(r: &mut Lcg, g: i64) -> Geometry {
    let (x0, y0) = (r.n(g), r.n(g));
    let (mut x1, mut y1) = (r.n(g), r.n(g));
    if x0 == x1 && y0 == y1 {
        x1 += 1.0;
        y1 += 1.0;
    }
    line(x0, y0, x1, y1)
}
fn rpt(r: &mut Lcg, g: i64) -> Geometry {
    pt(r.n(g), r.n(g))
}

#[test]
fn stress_polygon_polygon() {
    let mut r = Lcg(0xA11CE);
    for _ in 0..500 {
        check2(&rrect(&mut r, 8), &rrect(&mut r, 8));
    }
}

#[test]
fn stress_point_and_line_vs_polygon() {
    let mut r = Lcg(0xB0B);
    for _ in 0..400 {
        let poly = rrect(&mut r, 8);
        check2(&poly, &rpt(&mut r, 9));
        check2(&poly, &rline(&mut r, 9));
    }
}

#[test]
fn stress_line_line_and_misc() {
    let mut r = Lcg(0xCAFE);
    for _ in 0..400 {
        check2(&rline(&mut r, 6), &rline(&mut r, 6));
    }
    let mut r = Lcg(0xD00D);
    for _ in 0..200 {
        let a = ls(&[(r.n(8), r.n(8)), (r.n(8), r.n(8)), (r.n(8), r.n(8))]);
        check2(&a, &rrect(&mut r, 8));
    }
}

fn rtri(r: &mut Lcg, g: i64) -> Option<Geometry> {
    let p = [(r.n(g), r.n(g)), (r.n(g), r.n(g)), (r.n(g), r.n(g))];
    let area2 = (p[1].0 - p[0].0) * (p[2].1 - p[0].1) - (p[2].0 - p[0].0) * (p[1].1 - p[0].1);
    if area2.abs() < 1.0 {
        return None;
    }
    Some(poly(&[p[0], p[1], p[2], p[0]]))
}

#[test]
fn stress_triangles_vs_areas() {
    let mut r = Lcg(0x7711);
    let mut n = 0;
    while n < 400 {
        let Some(t) = rtri(&mut r, 9) else { continue };
        check2(&t, &rrect(&mut r, 9));
        if let Some(t2) = rtri(&mut r, 9) {
            check2(&t, &t2);
        }
        n += 1;
    }
}
