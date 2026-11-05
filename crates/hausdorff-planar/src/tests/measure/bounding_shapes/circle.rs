//! The minimum bounding circle: closed forms, enclosure, minimality, and
//! agreement with exhaustive enumeration over every 1-, 2- and 3-point candidate.

use super::*;

#[test]
fn the_bounding_circle_of_a_square_is_its_circumcircle() {
    // A closed form: the circle through the four corners of a side-s square has
    // centre (s/2, s/2) and radius s/sqrt(2).
    for side in [1.0_f64, 2.0, 7.5, 1e4] {
        let c = square(side)
            .minimum_bounding_circle()
            .expect("a square has a bounding circle");
        close(c.radius, side / 2.0_f64.sqrt(), 1e-9, "circumradius");
        close(c.centre.x, side / 2.0, 1e-9, "centre x");
        close(c.centre.y, side / 2.0, 1e-9, "centre y");
    }
}

#[test]
fn the_bounding_circle_of_two_points_is_their_diameter() {
    // The two-point case is its own branch in Welzl's algorithm, so it gets its
    // own closed form: the circle on the segment as diameter.
    let a = Coord::new(-3.0, 4.0);
    let b = Coord::new(5.0, -2.0);
    let c = MultiPoint::new(vec![Point(a), Point(b)])
        .minimum_bounding_circle()
        .expect("two points have a bounding circle");
    let expected_radius = ((a.x - b.x).powi(2) + (a.y - b.y).powi(2)).sqrt() / 2.0;
    close(c.radius, expected_radius, 1e-12, "half the separation");
    close(c.centre.x, (a.x + b.x) / 2.0, 1e-12, "midpoint x");
    close(c.centre.y, (a.y + b.y) / 2.0, 1e-12, "midpoint y");
}

#[test]
fn the_bounding_circle_encloses_everything_and_is_minimal() {
    // Enclosure alone would be satisfied by any large circle, so minimality is
    // asserted too: shrinking by a hair must leave something outside.
    let shapes: Vec<Geometry> = vec![
        Geometry::Polygon(square(4.0)),
        Geometry::Polygon(rotate(&square(4.0), 0.3)),
        Geometry::Polygon(Polygon::new(
            LineString::from(vec![(0.0, 0.0), (10.0, 0.0), (5.0, 1.0), (0.0, 0.0)]),
            vec![],
        )),
        Geometry::MultiPoint(MultiPoint::new(vec![
            Point(Coord::new(0.0, 0.0)),
            Point(Coord::new(1.0, 7.0)),
            Point(Coord::new(-4.0, 2.0)),
            Point(Coord::new(3.0, -3.0)),
            Point(Coord::new(0.5, 0.5)),
        ])),
        Geometry::LineString(LineString::from(vec![
            (0.0, 0.0),
            (1.0, 0.1),
            (2.0, -0.1),
            (3.0, 0.0),
        ])),
    ];
    for g in &shapes {
        let c = g.minimum_bounding_circle().expect("non-empty");
        let coords = g.to_coords();
        for p in &coords {
            assert!(
                c.contains(*p),
                "{p:?} is outside its own bounding circle {c:?}"
            );
        }
        // Minimality: a circle 0.1% smaller must fail to enclose.
        let shrunk = Circle {
            centre: c.centre,
            radius: c.radius * 0.999,
        };
        assert!(
            coords.iter().any(|p| !shrunk.contains(*p)),
            "the circle {c:?} is not minimal — it still encloses everything 0.1% smaller"
        );
    }
}

#[test]
fn the_circle_polygon_approximation_is_inscribed() {
    // `to_polygon` is documented as inscribed, so its area is below the true
    // circle's and approaches it as segments grow. Pinned because the opposite
    // convention (circumscribed) is equally common and a caller comparing areas
    // would be wrong by a predictable factor.
    let c = Circle {
        centre: Coord::new(1.0, -2.0),
        radius: 3.0,
    };
    let exact = std::f64::consts::PI * 9.0;
    let mut previous = 0.0;
    for segments in [8usize, 16, 64, 256] {
        let area = c.to_polygon(segments).unsigned_area();
        assert!(
            area < exact,
            "an inscribed polygon of {segments} sides must be smaller than the circle"
        );
        assert!(
            area > previous,
            "more segments must approximate better: {area} after {previous}"
        );
        previous = area;
    }
    // The deficit of an inscribed regular n-gon has a closed form:
    //
    //     1 - (n / 2pi) * sin(2pi / n)   ~   (2pi)^2 / (6 n^2)
    //
    // which is 1.004e-4 at n = 256. Comparing against the formula rather than a
    // round number avoids asserting a bound the geometry cannot meet — the first
    // version of this demanded under 1e-4 and was impossible by 0.4%.
    let n = 256.0_f64;
    let deficit = 1.0 - (n / std::f64::consts::TAU) * (std::f64::consts::TAU / n).sin();
    let measured = (exact - previous) / exact;
    assert!(
        (measured - deficit).abs() < 1e-9,
        "the 256-gon deficit should match the closed form {deficit:e}, got {measured:e}"
    );
}

// ---- adversarial inputs ----
//
// The descriptors above have closed forms on clean input. These check the cases
// where a closed form is not available, against oracles that are independent of
// the implementations: exhaustive enumeration for the two hull-based descriptors,
// and a hand-computed optimum for the inscribed circle.

#[test]
fn the_bounding_circle_agrees_with_exhaustive_enumeration() {
    // The minimum enclosing circle is always determined by one, two, or three of
    // the input points, so enumerating every such candidate and keeping the
    // smallest that encloses everything *is* the answer — O(n⁴) and obviously
    // correct, which is exactly what an oracle should be. Welzl's algorithm must
    // agree with it on every input, and the interesting inputs are the degenerate
    // ones: `from_three` returns `None` for a collinear triple, and the inner loop
    // then keeps a circle that does not contain the point it was called for.
    let brute = |pts: &[Coord]| -> Option<f64> {
        let encloses = |c: Coord, r: f64| {
            pts.iter().all(|p| {
                let (dx, dy) = (p.x - c.x, p.y - c.y);
                (dx * dx + dy * dy).sqrt() <= r * (1.0 + 1e-9) + 1e-9
            })
        };
        let mut best: Option<f64> = None;
        let mut offer = |c: Coord, r: f64| {
            if r.is_finite() && encloses(c, r) && best.is_none_or(|b| r < b) {
                best = Some(r);
            }
        };
        for p in pts {
            offer(*p, 0.0);
        }
        for a in pts {
            for b in pts {
                let c = Coord::new(f64::midpoint(a.x, b.x), f64::midpoint(a.y, b.y));
                offer(c, ((a.x - c.x).powi(2) + (a.y - c.y).powi(2)).sqrt());
            }
        }
        for a in pts {
            for b in pts {
                for d in pts {
                    let det = 2.0 * (a.x * (b.y - d.y) + b.x * (d.y - a.y) + d.x * (a.y - b.y));
                    if det.abs() < 1e-12 {
                        continue;
                    }
                    let (a2, b2, d2) = (
                        a.x * a.x + a.y * a.y,
                        b.x * b.x + b.y * b.y,
                        d.x * d.x + d.y * d.y,
                    );
                    let c = Coord::new(
                        (a2 * (b.y - d.y) + b2 * (d.y - a.y) + d2 * (a.y - b.y)) / det,
                        (a2 * (d.x - b.x) + b2 * (a.x - d.x) + d2 * (b.x - a.x)) / det,
                    );
                    offer(c, ((a.x - c.x).powi(2) + (a.y - c.y).powi(2)).sqrt());
                }
            }
        }
        best
    };

    let mut state = 0x1234_5678_9abc_def0_u64;
    for trial in 0..4000 {
        let pts = grid_points(&mut state, 8);
        let mp = MultiPoint::new(pts.iter().map(|c| Point(*c)).collect());
        let got = mp.minimum_bounding_circle().expect("non-empty input");
        let want = brute(&pts).expect("non-empty input");

        for p in &pts {
            let (dx, dy) = (p.x - got.centre.x, p.y - got.centre.y);
            assert!(
                (dx * dx + dy * dy).sqrt() <= got.radius * (1.0 + 1e-9) + 1e-9,
                "trial {trial}: {p:?} is outside the bounding circle of {pts:?} \
                 (centre {:?}, radius {})",
                got.centre,
                got.radius
            );
        }
        assert!(
            got.radius <= want * (1.0 + 1e-9) + 1e-9,
            "trial {trial}: radius {} exceeds the exhaustive minimum {want} for {pts:?}",
            got.radius
        );
    }
}
