//! The oriented envelope: rotation invariance, enclosure, and agreement with
//! enumeration over every direction a pair of input points defines.

use super::*;

#[test]
fn the_oriented_envelope_of_a_rotated_rectangle_is_that_rectangle() {
    // The sharpest test here. A 6x2 rectangle has area 12 at any rotation, and the
    // oriented envelope must report 12 every time — while the axis-aligned
    // bounding box grows by up to sqrt(2). That gap is the whole reason this
    // descriptor exists.
    let base = Polygon::new(
        LineString::from(vec![
            (0.0, 0.0),
            (6.0, 0.0),
            (6.0, 2.0),
            (0.0, 2.0),
            (0.0, 0.0),
        ]),
        vec![],
    );
    let mut turned_bigger = 0;
    let mut i = 0;
    while i < 24 {
        let theta = std::f64::consts::PI * f64::from(i) / 24.0;
        let r = rotate(&base, theta);
        let env = r.oriented_envelope().expect("a rectangle has an envelope");
        close(
            oriented_envelope_area(env),
            12.0,
            1e-9,
            &format!("envelope area at {theta} rad"),
        );
        close(env.length, 6.0, 1e-9, "longer side");
        close(env.width, 2.0, 1e-9, "shorter side");
        close(env.aspect_ratio(), 1.0 / 3.0, 1e-9, "aspect ratio");

        // The reported angle must be the rotation, modulo pi since an axis has no
        // direction.
        let expected = theta.rem_euclid(std::f64::consts::PI);
        let delta = (env.angle - expected).abs();
        let delta = delta.min(std::f64::consts::PI - delta);
        assert!(
            delta < 1e-6,
            "envelope angle {} should be {expected} at rotation {theta}",
            env.angle
        );

        // The axis-aligned box, for contrast.
        let bbox = r.bounding_rect().expect("has a bbox");
        let bbox_area = (bbox.max().x - bbox.min().x) * (bbox.max().y - bbox.min().y);
        if bbox_area > 12.0 * 1.01 {
            turned_bigger += 1;
        }
        i += 1;
    }
    assert!(
        turned_bigger > 0,
        "the axis-aligned box is supposed to be worse at some rotations; if it \
         never is, this test is not testing the thing it claims to"
    );
}

#[test]
fn the_oriented_envelope_encloses_everything_and_beats_the_bounding_box() {
    let shapes: Vec<Geometry> = vec![
        Geometry::Polygon(rotate(&square(3.0), 0.4)),
        Geometry::Polygon(Polygon::new(
            LineString::from(vec![(0.0, 0.0), (10.0, 10.0), (10.5, 9.5), (0.5, -0.5)]),
            vec![],
        )),
        Geometry::MultiPoint(MultiPoint::new(vec![
            Point(Coord::new(0.0, 0.0)),
            Point(Coord::new(5.0, 5.0)),
            Point(Coord::new(10.0, 10.0)),
            Point(Coord::new(1.0, 0.0)),
        ])),
    ];
    for g in &shapes {
        let env = g.oriented_envelope().expect("non-empty");
        let rect = env.to_polygon();
        for p in g.to_coords() {
            // Generous tolerance on the boundary: a hull vertex lies exactly on the
            // envelope edge by construction, which `contains` does not guarantee.
            let inside = rect.contains(&Point(p))
                || env
                    .corners
                    .iter()
                    .any(|c| (c.x - p.x).abs() < 1e-9 && (c.y - p.y).abs() < 1e-9)
                || rect.exterior().0.windows(2).any(|w| {
                    let d = {
                        let (dx, dy) = (w[1].x - w[0].x, w[1].y - w[0].y);
                        let l2 = dx * dx + dy * dy;
                        if l2 <= 0.0 {
                            ((p.x - w[0].x).powi(2) + (p.y - w[0].y).powi(2)).sqrt()
                        } else {
                            let t =
                                (((p.x - w[0].x) * dx + (p.y - w[0].y) * dy) / l2).clamp(0.0, 1.0);
                            ((p.x - w[0].x - t * dx).powi(2) + (p.y - w[0].y - t * dy).powi(2))
                                .sqrt()
                        }
                    };
                    d < 1e-9
                });
            assert!(inside, "{p:?} escaped its oriented envelope {env:?}");
        }
        // The oriented envelope can never be worse than the axis-aligned box, since
        // the axis-aligned orientation is one of the candidates it considers.
        let bbox = g.bounding_rect().expect("has a bbox");
        let bbox_area = (bbox.max().x - bbox.min().x) * (bbox.max().y - bbox.min().y);
        assert!(
            oriented_envelope_area(env) <= bbox_area * (1.0 + 1e-9),
            "an oriented envelope of {} must not exceed the axis-aligned {bbox_area}",
            oriented_envelope_area(env)
        );
    }
}

#[test]
fn the_oriented_envelope_agrees_with_exhaustive_enumeration() {
    // The minimum-area enclosing rectangle always has a side flush with a hull
    // edge, so every direction defined by a pair of input points is a superset of
    // the candidate orientations. That makes the enumeration exact rather than
    // sampled — a sampled oracle could not distinguish a genuinely worse answer
    // from one a fraction of a degree off.
    let brute = |pts: &[Coord]| -> Option<f64> {
        let mut best: Option<f64> = None;
        for a in pts {
            for b in pts {
                let (dx, dy) = (b.x - a.x, b.y - a.y);
                let len = (dx * dx + dy * dy).sqrt();
                if len < 1e-12 {
                    continue;
                }
                let (ux, uy) = (dx / len, dy / len);
                let (mut lo_u, mut hi_u) = (f64::INFINITY, f64::NEG_INFINITY);
                let (mut lo_v, mut hi_v) = (f64::INFINITY, f64::NEG_INFINITY);
                for p in pts {
                    let u = p.x * ux + p.y * uy;
                    let v = -p.x * uy + p.y * ux;
                    lo_u = lo_u.min(u);
                    hi_u = hi_u.max(u);
                    lo_v = lo_v.min(v);
                    hi_v = hi_v.max(v);
                }
                let area = (hi_u - lo_u) * (hi_v - lo_v);
                if best.is_none_or(|b| area < b) {
                    best = Some(area);
                }
            }
        }
        best
    };

    let mut state = 0xfeed_face_dead_beef_u64;
    for trial in 0..3000 {
        let pts = grid_points(&mut state, 8);
        let mp = MultiPoint::new(pts.iter().map(|c| Point(*c)).collect());
        let (Some(rect), Some(want)) = (mp.oriented_envelope(), brute(&pts)) else {
            continue;
        };

        // Enclosure, measured in the rectangle's own frame rather than with a
        // point-in-polygon test, so a degenerate (zero-width) rectangle is still
        // checked rather than silently passing.
        let origin = rect.corners[0];
        let edge = |i: usize| {
            let e = (rect.corners[i].x - origin.x, rect.corners[i].y - origin.y);
            let l = (e.0 * e.0 + e.1 * e.1).sqrt();
            (e, l)
        };
        let ((e1, l1), (e2, l2)) = (edge(1), edge(3));
        if l1 > 1e-12 && l2 > 1e-12 {
            for p in &pts {
                let s = ((p.x - origin.x) * e1.0 + (p.y - origin.y) * e1.1) / l1;
                let t = ((p.x - origin.x) * e2.0 + (p.y - origin.y) * e2.1) / l2;
                assert!(
                    s >= -1e-9 && s <= l1 + 1e-9 && t >= -1e-9 && t <= l2 + 1e-9,
                    "trial {trial}: {p:?} lies outside the oriented envelope of {pts:?}"
                );
            }
        }
        let got = oriented_envelope_area(rect);
        assert!(
            got <= want * (1.0 + 1e-9) + 1e-12,
            "trial {trial}: envelope area {got} beats the exhaustive minimum {want} \
             for {pts:?} — the calipers missed an orientation"
        );
    }
}
