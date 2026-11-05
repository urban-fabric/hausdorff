use super::*;

pub(crate) fn segments(g: &Geometry) -> Vec<(Coord, Coord)> {
    let mut out = Vec::new();
    push_segments(g, &mut out);
    out
}

pub(super) fn push_segments(g: &Geometry, out: &mut Vec<(Coord, Coord)>) {
    match g {
        Geometry::Point(p) => {
            let c = Coord::new(p.x(), p.y());
            out.push((c, c));
        }
        Geometry::Line(l) => out.push((l.start, l.end)),
        Geometry::LineString(ls) => push_line_string(ls, out),
        Geometry::Polygon(poly) => push_polygon(poly, out),
        Geometry::MultiPoint(mp) => {
            for p in mp.iter() {
                let c = Coord::new(p.x(), p.y());
                out.push((c, c));
            }
        }
        Geometry::MultiLineString(mls) => {
            for ls in mls.iter() {
                push_line_string(ls, out);
            }
        }
        Geometry::MultiPolygon(mp) => {
            for poly in mp.iter() {
                push_polygon(poly, out);
            }
        }
        Geometry::Rect(r) => push_rect(r, out),
        Geometry::Triangle(t) => push_triangle(t, out),
        Geometry::GeometryCollection(c) => {
            for g in c.iter() {
                push_segments(g, out);
            }
        }
    }
}

fn push_line_string(ls: &LineString, out: &mut Vec<(Coord, Coord)>) {
    match ls.0.as_slice() {
        [] => {}
        [only] => out.push((*only, *only)),
        cs => {
            for w in cs.windows(2) {
                out.push((w[0], w[1]));
            }
        }
    }
}

fn push_polygon(poly: &Polygon, out: &mut Vec<(Coord, Coord)>) {
    push_line_string(poly.exterior(), out);
    for hole in poly.interiors() {
        push_line_string(hole, out);
    }
}

fn push_rect(r: &Rect, out: &mut Vec<(Coord, Coord)>) {
    let (lo, hi) = (r.min(), r.max());
    let corners = [
        Coord::new(lo.x, lo.y),
        Coord::new(hi.x, lo.y),
        Coord::new(hi.x, hi.y),
        Coord::new(lo.x, hi.y),
    ];
    for i in 0..4 {
        out.push((corners[i], corners[(i + 1) % 4]));
    }
}

fn push_triangle(t: &Triangle, out: &mut Vec<(Coord, Coord)>) {
    out.push((t.0, t.1));
    out.push((t.1, t.2));
    out.push((t.2, t.0));
}
