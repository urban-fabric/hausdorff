use super::*;

#[derive(Clone, Copy, Debug)]
struct RefAffine {
    a: f64,
    b: f64,
    xoff: f64,
    d: f64,
    e: f64,
    yoff: f64,
}

impl RefAffine {
    fn identity() -> Self {
        Self {
            a: 1.0,
            b: 0.0,
            xoff: 0.0,
            d: 0.0,
            e: 1.0,
            yoff: 0.0,
        }
    }
    fn translate(x: f64, y: f64) -> Self {
        Self {
            a: 1.0,
            b: 0.0,
            xoff: x,
            d: 0.0,
            e: 1.0,
            yoff: y,
        }
    }
    fn scale(xf: f64, yf: f64, origin: (f64, f64)) -> Self {
        let (x0, y0) = origin;
        Self {
            a: xf,
            b: 0.0,
            xoff: x0 - x0 * xf,
            d: 0.0,
            e: yf,
            yoff: y0 - y0 * yf,
        }
    }
    fn rotate(degrees: f64, origin: (f64, f64)) -> Self {
        let (s, c) = degrees.to_radians().sin_cos();
        let (x0, y0) = origin;
        Self {
            a: c,
            b: -s,
            xoff: x0 - x0 * c + y0 * s,
            d: s,
            e: c,
            yoff: y0 - x0 * s - y0 * c,
        }
    }
    fn skew(xs: f64, ys: f64, origin: (f64, f64)) -> Self {
        let (x0, y0) = origin;
        let mut tanx = xs.to_radians().tan();
        let mut tany = ys.to_radians().tan();
        if tanx.abs() < 2.5e-16 {
            tanx = 0.0;
        }
        if tany.abs() < 2.5e-16 {
            tany = 0.0;
        }
        Self {
            a: 1.0,
            b: tanx,
            xoff: -y0 * tanx,
            d: tany,
            e: 1.0,
            yoff: -x0 * tany,
        }
    }
    /// `self` then `other`: the matrix product `other · self`.
    fn then(self, other: RefAffine) -> RefAffine {
        RefAffine {
            a: other.a * self.a + other.b * self.d,
            b: other.a * self.b + other.b * self.e,
            xoff: other.a * self.xoff + other.b * self.yoff + other.xoff,
            d: other.d * self.a + other.e * self.d,
            e: other.d * self.b + other.e * self.e,
            yoff: other.d * self.xoff + other.e * self.yoff + other.yoff,
        }
    }
    fn apply(&self, p: (f64, f64)) -> (f64, f64) {
        (
            self.a * p.0 + self.b * p.1 + self.xoff,
            self.d * p.0 + self.e * p.1 + self.yoff,
        )
    }
    fn params(&self) -> [f64; 6] {
        [self.a, self.b, self.xoff, self.d, self.e, self.yoff]
    }
}

fn assert_line_strings_eq_ref(h: &LineString, r: &[(f64, f64)], ctx: &str) {
    assert_eq!(h.0.len(), r.len(), "{ctx}: coord count");
    for (hc, &(rx, ry)) in h.0.iter().zip(r.iter()) {
        assert!(
            (hc.x - rx).abs() < EPS && (hc.y - ry).abs() < EPS,
            "{ctx}: ({}, {}) vs ref ({rx}, {ry})",
            hc.x,
            hc.y,
        );
    }
}

fn assert_polygons_eq_ref(h: &Polygon, ext: &[(f64, f64)], holes: &[Vec<(f64, f64)>], ctx: &str) {
    assert_line_strings_eq_ref(h.exterior(), ext, ctx);
    assert_eq!(h.interiors().len(), holes.len(), "{ctx}: holes");
    for (hr, gr) in h.interiors().iter().zip(holes) {
        assert_line_strings_eq_ref(hr, gr, ctx);
    }
}

/// Hausdorff stores rings closed; produce the closed form a path takes once
/// inside a `Polygon`, then map each coord through `t`.
fn apply_ring(t: &RefAffine, ring: &[(f64, f64)]) -> Vec<(f64, f64)> {
    closed(ring).iter().map(|&c| t.apply(c)).collect()
}

fn apply_path(t: &RefAffine, path: &[(f64, f64)]) -> Vec<(f64, f64)> {
    path.iter().map(|&c| t.apply(c)).collect()
}

// The fixtures used across the affine tests, kept identical to the originals.
fn poly_fixture() -> PolyFixture {
    let exterior = vec![(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0)];
    let hole = vec![(3.0, 3.0), (5.0, 3.0), (5.0, 5.0), (3.0, 5.0)];
    (
        h_polygon(&exterior, std::slice::from_ref(&hole)),
        exterior,
        vec![hole],
    )
}

fn ls_fixture() -> (LineString, Vec<(f64, f64)>) {
    let path = vec![(0.0, 0.0), (5.0, 1.0), (7.0, 6.0), (-2.0, 4.0)];
    (h_line_string(&path), path)
}

#[test]
fn affine_transform_application_matches_reference() {
    let transforms: Vec<(AffineTransform, RefAffine)> = vec![
        (AffineTransform::identity(), RefAffine::identity()),
        (
            AffineTransform::translate(3.5, -2.25),
            RefAffine::translate(3.5, -2.25),
        ),
        (
            AffineTransform::scale(2.0, -0.5, Coord::new(1.0, 2.0)),
            RefAffine::scale(2.0, -0.5, (1.0, 2.0)),
        ),
        (
            AffineTransform::rotate(37.0, Coord::new(-1.0, 4.0)),
            RefAffine::rotate(37.0, (-1.0, 4.0)),
        ),
        (
            AffineTransform::skew(15.0, -30.0, Coord::new(2.0, 2.0)),
            RefAffine::skew(15.0, -30.0, (2.0, 2.0)),
        ),
        (
            AffineTransform::translate(1.0, 2.0)
                .rotated(45.0, Coord::new(0.0, 0.0))
                .scaled(3.0, 0.25, Coord::new(5.0, 5.0))
                .skewed(10.0, 20.0, Coord::new(-2.0, 3.0)),
            RefAffine::translate(1.0, 2.0)
                .then(RefAffine::rotate(45.0, (0.0, 0.0)))
                .then(RefAffine::scale(3.0, 0.25, (5.0, 5.0)))
                .then(RefAffine::skew(10.0, 20.0, (-2.0, 3.0))),
        ),
    ];
    let (h_poly, ext, holes) = poly_fixture();
    let (h_ls, ls_path) = ls_fixture();
    for (i, (ht, rt)) in transforms.iter().enumerate() {
        assert_eq!(
            [ht.a(), ht.b(), ht.xoff(), ht.d(), ht.e(), ht.yoff()],
            rt.params(),
            "transform {i}: matrix"
        );
        let ref_ext = apply_ring(rt, &ext);
        let ref_holes: Vec<Vec<(f64, f64)>> = holes.iter().map(|h| apply_ring(rt, h)).collect();
        assert_polygons_eq_ref(
            &h_poly.affine_transform(ht),
            &ref_ext,
            &ref_holes,
            &format!("transform {i} on polygon"),
        );
        assert_line_strings_eq_ref(
            &h_ls.affine_transform(ht),
            &apply_path(rt, &ls_path),
            &format!("transform {i} on linestring"),
        );
        let mut h_mut = h_ls.clone();
        h_mut.affine_transform_mut(ht);
        assert_eq!(h_mut, h_ls.affine_transform(ht), "transform {i} mut");
    }
}

/// Centroid of a hausdorff geometry, taken from its own verified `Centroid`
/// kernel — the anchor for the centroid-pivot convenience transforms.
fn centroid_of(g: &impl Centroid) -> (f64, f64) {
    let c = g.centroid().expect("non-empty fixture has a centroid");
    (c.x(), c.y())
}

/// Bounding-box centre, computed independently from the geometry's coords.
fn bbox_center(coords: &[(f64, f64)]) -> (f64, f64) {
    let (mut minx, mut miny) = (f64::INFINITY, f64::INFINITY);
    let (mut maxx, mut maxy) = (f64::NEG_INFINITY, f64::NEG_INFINITY);
    for &(x, y) in coords {
        minx = minx.min(x);
        miny = miny.min(y);
        maxx = maxx.max(x);
        maxy = maxy.max(y);
    }
    ((minx + maxx) / 2.0, (miny + maxy) / 2.0)
}

#[test]
fn translate_rotate_scale_skew_match_reference() {
    let (h_poly, ext, holes) = poly_fixture();
    let (h_ls, ls_path) = ls_fixture();

    // The polygon's bounding box is over its closed exterior ring; the
    // linestring's over its own coords. Centroids come from hausdorff's kernel.
    let poly_bbc = bbox_center(&closed(&ext));
    let poly_centroid = centroid_of(&h_poly);
    let ls_centroid = centroid_of(&h_ls);

    let apply_poly = |t: &RefAffine| -> Rings {
        (
            apply_ring(t, &ext),
            holes.iter().map(|h| apply_ring(t, h)).collect(),
        )
    };

    // translate
    {
        let t = RefAffine::translate(17.0, -3.0);
        let (re, rh) = apply_poly(&t);
        assert_polygons_eq_ref(&h_poly.translate(17.0, -3.0), &re, &rh, "translate polygon");
    }
    {
        let t = RefAffine::translate(-0.5, 9.0);
        assert_line_strings_eq_ref(
            &h_ls.translate(-0.5, 9.0),
            &apply_path(&t, &ls_path),
            "translate linestring",
        );
    }

    // rotate
    {
        let t = RefAffine::rotate(-45.0, (10.0, 0.0));
        let (re, rh) = apply_poly(&t);
        assert_polygons_eq_ref(
            &h_poly.rotate_around_point(-45.0, Point::new(10.0, 0.0)),
            &re,
            &rh,
            "rotate around point",
        );
    }
    {
        let t = RefAffine::rotate(33.0, poly_centroid);
        let (re, rh) = apply_poly(&t);
        assert_polygons_eq_ref(
            &h_poly.rotate_around_centroid(33.0),
            &re,
            &rh,
            "rotate around centroid",
        );
    }
    {
        let t = RefAffine::rotate(-120.0, poly_bbc);
        let (re, rh) = apply_poly(&t);
        assert_polygons_eq_ref(
            &h_poly.rotate_around_center(-120.0),
            &re,
            &rh,
            "rotate around center",
        );
    }
    {
        let t = RefAffine::rotate(75.0, ls_centroid);
        assert_line_strings_eq_ref(
            &h_ls.rotate_around_centroid(75.0),
            &apply_path(&t, &ls_path),
            "rotate linestring around centroid",
        );
    }

    // scale (about bounding-box centre)
    {
        let t = RefAffine::scale(2.5, 2.5, poly_bbc);
        let (re, rh) = apply_poly(&t);
        assert_polygons_eq_ref(&h_poly.scale(2.5), &re, &rh, "scale");
    }
    {
        let t = RefAffine::scale(2.0, -4.0, poly_bbc);
        let (re, rh) = apply_poly(&t);
        assert_polygons_eq_ref(&h_poly.scale_xy(2.0, -4.0), &re, &rh, "scale_xy");
    }
    {
        let t = RefAffine::scale(2.0, 4.0, (100.0, 100.0));
        assert_line_strings_eq_ref(
            &h_ls.scale_around_point(2.0, 4.0, Coord::new(100.0, 100.0)),
            &apply_path(&t, &ls_path),
            "scale_around_point",
        );
    }

    // skew (about bounding-box centre)
    {
        let t = RefAffine::skew(30.0, 30.0, poly_bbc);
        let (re, rh) = apply_poly(&t);
        assert_polygons_eq_ref(&h_poly.skew(30.0), &re, &rh, "skew");
    }
    {
        let t = RefAffine::skew(30.0, 12.0, poly_bbc);
        let (re, rh) = apply_poly(&t);
        assert_polygons_eq_ref(&h_poly.skew_xy(30.0, 12.0), &re, &rh, "skew_xy");
    }
    {
        let t = RefAffine::skew(45.0, 10.0, (2.0, 2.0));
        assert_line_strings_eq_ref(
            &h_ls.skew_around_point(45.0, 10.0, Coord::new(2.0, 2.0)),
            &apply_path(&t, &ls_path),
            "skew_around_point",
        );
    }

    // mut variants must agree with their immutable counterparts.
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

// ---------------------------------------------------------------------------
// Independent extremes reference: a strict-comparison scan over the exterior
// coordinate walk (interior rings excluded, first coord wins ties).
// ---------------------------------------------------------------------------
