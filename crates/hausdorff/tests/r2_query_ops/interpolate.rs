use super::*;

fn seg_len(a: (f64, f64), b: (f64, f64)) -> f64 {
    dist2(a, b).sqrt()
}

fn path_len(path: &[(f64, f64)]) -> f64 {
    path.windows(2).map(|w| seg_len(w[0], w[1])).sum()
}

/// Reference for `Line::line_interpolate_point` / `LineString` over a polyline.
/// Returns `None` for `NaN` fractions or non-finite results, clamps out-of-range
/// fractions to the endpoints, and returns `None` for empty paths in range.
fn ref_interpolate(path: &[(f64, f64)], fraction: f64) -> Option<(f64, f64)> {
    if fraction.is_nan() {
        return None;
    }
    let f = fraction.clamp(0.0, 1.0);
    let total = path_len(path);
    let target = total * f;
    let mut cum = 0.0;
    for w in path.windows(2) {
        let (a, b) = (w[0], w[1]);
        let len = seg_len(a, b);
        if cum + len >= target {
            // Mirror the kernel exactly: a zero-length segment yields a NaN
            // segment-fraction (0/0), and interpolating that gives None.
            let sf = (target - cum) / len;
            return ref_interpolate_seg(a, b, sf);
        }
        cum += len;
    }
    None
}

fn ref_interpolate_seg(a: (f64, f64), b: (f64, f64), fraction: f64) -> Option<(f64, f64)> {
    if fraction.is_nan() {
        return None;
    }
    let f = fraction.clamp(0.0, 1.0);
    let r = (a.0 + (b.0 - a.0) * f, a.1 + (b.1 - a.1) * f);
    (r.0.is_finite() && r.1.is_finite()).then_some(r)
}

/// Perpendicular distance from `p` to the bounded segment `a`–`b`.
fn point_seg_dist(p: (f64, f64), a: (f64, f64), b: (f64, f64)) -> f64 {
    let (abx, aby) = (b.0 - a.0, b.1 - a.1);
    let len2 = abx * abx + aby * aby;
    let proj = if len2 == 0.0 {
        a
    } else {
        let t = ((p.0 - a.0) * abx + (p.1 - a.1) * aby) / len2;
        let t = t.clamp(0.0, 1.0);
        (a.0 + abx * t, a.1 + aby * t)
    };
    seg_len(p, proj)
}

/// Reference for `Line::line_locate_point`: fraction along the segment of the
/// clamped projection; `None` on non-finite input, `Some(0)` for zero-length.
fn ref_locate_seg(a: (f64, f64), b: (f64, f64), p: (f64, f64)) -> Option<f64> {
    let (vx, vy) = (b.0 - a.0, b.1 - a.1);
    let v_sq = vx * vx + vy * vy;
    if v_sq == 0.0 {
        return Some(0.0);
    }
    let l = ((p.0 - a.0) * vx + (p.1 - a.1) * vy) / v_sq;
    l.is_finite().then(|| l.clamp(0.0, 1.0))
}

/// Reference for `LineString::line_locate_point`: the fraction of total length
/// of the closest location along the polyline.
fn ref_locate_path(path: &[(f64, f64)], p: (f64, f64)) -> Option<f64> {
    let total = path_len(path);
    if total == 0.0 {
        return Some(0.0);
    }
    let mut cum = 0.0;
    let mut closest = f64::INFINITY;
    let mut fraction = 0.0;
    for w in path.windows(2) {
        let (a, b) = (w[0], w[1]);
        let d = point_seg_dist(p, a, b);
        let len = seg_len(a, b);
        let sf = ref_locate_seg(a, b, p)?;
        if d < closest {
            closest = d;
            fraction = (cum + sf * len) / total;
        }
        cum += len;
    }
    Some(fraction)
}

#[test]
fn line_interpolate_point_matches_reference() {
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
    let seg = ((-1.0, 0.0), (3.0, 2.0));
    let h_seg = Line::new(Coord::new(seg.0.0, seg.0.1), Coord::new(seg.1.0, seg.1.1));
    let paths = [
        vec![(-1.0, 0.0), (0.0, 0.0), (0.0, 1.0)],
        vec![(0.0, 0.0), (3.0, 4.0), (3.0, 10.0), (-2.0, 10.0)],
        vec![(1.0, 1.0), (1.0, 1.0)],
        vec![],
    ];
    for &f in &fractions {
        let h = h_seg.line_interpolate_point(f);
        let r = ref_interpolate_seg(seg.0, seg.1, f);
        match (h, r) {
            (Some(hp), Some(rp)) => assert_points_eq(hp, rp, &format!("line interpolate {f}")),
            (None, None) => {}
            (h, r) => panic!("line interpolate {f}: {h:?} vs ref {r:?}"),
        }
        for path in &paths {
            let h = h_line_string(path).line_interpolate_point(f);
            // A two-point zero-length path yields a finite point; an empty path
            // yields None.
            let r = if path.len() < 2 {
                None
            } else {
                ref_interpolate(path, f)
            };
            match (h, r) {
                (Some(hp), Some(rp)) => {
                    assert_points_eq(hp, rp, &format!("ls interpolate {path:?} at {f}"))
                }
                (None, None) => {}
                (h, r) => panic!("ls interpolate {path:?} at {f}: {h:?} vs ref {r:?}"),
            }
        }
    }
}

#[test]
fn line_locate_point_matches_reference() {
    let seg = ((-1.0, 0.0), (3.0, 2.0));
    let h_seg = Line::new(Coord::new(seg.0.0, seg.0.1), Coord::new(seg.1.0, seg.1.1));
    let zero = ((1.0, 1.0), (1.0, 1.0));
    let h_zero = Line::new(
        Coord::new(zero.0.0, zero.0.1),
        Coord::new(zero.1.0, zero.1.1),
    );
    let path = vec![(0.0, 0.0), (3.0, 4.0), (3.0, 10.0), (-2.0, 10.0)];
    let h_ls = h_line_string(&path);

    let mut queries = query_points();
    queries.push((f64::NAN, 1.0));
    for (x, y) in queries {
        let q = (x, y);
        let p = Point::new(x, y);
        for (name, h, r) in [
            (
                "segment",
                h_seg.line_locate_point(&p),
                ref_locate_seg(seg.0, seg.1, q),
            ),
            (
                "zero",
                h_zero.line_locate_point(&p),
                ref_locate_seg(zero.0, zero.1, q),
            ),
            (
                "linestring",
                h_ls.line_locate_point(&p),
                ref_locate_path(&path, q),
            ),
        ] {
            match (h, r) {
                (Some(hf), Some(rf)) => assert!(
                    (hf - rf).abs() < EPS,
                    "{name} locate ({x}, {y}): {hf} vs ref {rf}"
                ),
                (None, None) => {}
                (h, r) => panic!("{name} locate ({x}, {y}): {h:?} vs ref {r:?}"),
            }
        }
    }

    // An empty (zero-length) line string locates everything at fraction 0.
    assert_eq!(
        h_line_string(&[]).line_locate_point(&Point::new(1.0, 1.0)),
        ref_locate_path(&[], (1.0, 1.0)),
    );
    assert_eq!(
        h_line_string(&[]).line_locate_point(&Point::new(1.0, 1.0)),
        Some(0.0),
    );
}

#[test]
fn locate_interpolate_round_trip_is_consistent() {
    let path = vec![(0.0, 0.0), (3.0, 4.0), (3.0, 10.0), (-2.0, 10.0)];
    let h_ls = h_line_string(&path);
    for i in 0..=20 {
        let f = i as f64 / 20.0;
        let hp = h_ls.line_interpolate_point(f).unwrap();
        let rp = ref_interpolate(&path, f).unwrap();
        assert_points_eq(hp, rp, &format!("round-trip interpolate {f}"));
        // Locating the interpolated point must recover the original fraction.
        let hf = h_ls.line_locate_point(&hp).unwrap();
        let rf = ref_locate_path(&path, (hp.x(), hp.y())).unwrap();
        assert!(
            (hf - rf).abs() < EPS && (hf - f).abs() < 1e-12,
            "round-trip locate {f}: {hf} vs ref {rf}"
        );
    }
}

// ---------------------------------------------------------------------------
// Independent affine reference oracle: a 3×3 (six-parameter) matrix built and
// applied with our own algebra, mirroring the documented constructors.
// ---------------------------------------------------------------------------
