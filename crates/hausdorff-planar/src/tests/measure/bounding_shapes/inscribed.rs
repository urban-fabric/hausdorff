//! The maximum inscribed circle, which is the approximate one of the three.
//!
//! Includes the dumbbell that showed the old coarse-to-fine search settling in a
//! local maximum and ignoring `tolerance` entirely.

use super::*;

#[test]
fn the_inscribed_circle_of_a_square_has_half_its_side() {
    // A closed form again: the largest circle in a side-s square has radius s/2,
    // centred at the middle.
    for side in [2.0_f64, 10.0] {
        let c = square(side)
            .maximum_inscribed_circle(side * 1e-4)
            .expect("a square has an interior");
        close(c.radius, side / 2.0, 1e-2, "inradius");
        close(c.centre.x, side / 2.0, 1e-2, "centre x");
        close(c.centre.y, side / 2.0, 1e-2, "centre y");
    }
}

#[test]
fn the_inscribed_circle_stays_inside_and_respects_holes() {
    // The hole is the case that matters: the largest circle in an annulus is
    // limited by the inner ring as much as the outer, and ignoring holes would
    // place a label in a lake.
    let annulus = Polygon::new(
        LineString::from(vec![
            (0.0, 0.0),
            (20.0, 0.0),
            (20.0, 20.0),
            (0.0, 20.0),
            (0.0, 0.0),
        ]),
        vec![LineString::from(vec![
            (6.0, 6.0),
            (14.0, 6.0),
            (14.0, 14.0),
            (6.0, 14.0),
            (6.0, 6.0),
        ])],
    );
    let c = annulus
        .maximum_inscribed_circle(0.05)
        .expect("the ring has an interior");
    // The centre must be in the ring, not the hole.
    assert!(
        annulus.contains(&Point(c.centre)),
        "the centre {:?} must be inside the annulus, not in its hole",
        c.centre
    );
    // The straight part of the ring is 6 units wide, which would cap the radius at
    // 3 — but the *corners* of the ring are wider than its sides, and that is where
    // the largest circle actually sits. On the diagonal at `(t, t)` the clearance is
    // `min(t, sqrt(2)(6 - t))`, maximised where the two are equal:
    //
    //     t = 6*sqrt(2) / (1 + sqrt(2)) = 3.5147...
    //
    // Asserting 3 here would have been asserting that a square ring is as tight at
    // its corners as along its edges, which it is not.
    let optimum = 6.0 * 2.0_f64.sqrt() / (1.0 + 2.0_f64.sqrt());
    assert!(
        (c.radius - optimum).abs() < 0.05,
        "the largest circle sits in a corner of the ring, radius {optimum}, got {}",
        c.radius
    );
    assert!(
        c.radius > 3.0,
        "the corner admits a circle larger than the 3 the straight side would allow, got {}",
        c.radius
    );
    // And the circle must not reach into the hole.
    let hole_centre = Coord::new(10.0, 10.0);
    let to_hole =
        ((c.centre.x - hole_centre.x).powi(2) + (c.centre.y - hole_centre.y).powi(2)).sqrt();
    assert!(
        to_hole > c.radius - 0.2,
        "the circle overlaps the hole: centre {to_hole} away, radius {}",
        c.radius
    );
}

#[test]
fn the_inscribed_circle_is_not_trapped_in_a_local_maximum() {
    // Clearance-to-boundary is not unimodal, and that is what broke the previous
    // coarse-to-fine search: it gridded the extent, kept the single best point,
    // then re-gridded *around that point*, which discards every basin but one. A
    // dumbbell has a local maximum in each lobe, so favouring the wrong lobe on
    // the coarse pass was unrecoverable — the error was set by the initial grid
    // spacing and did not shrink when `tolerance` did.
    //
    // Lobe A is 2.0 wide and lobe B is 2.2, joined by a corridor too thin to hold
    // either circle, so the exact answer is known: the largest circle inside
    // `[10, 12.2] × [0, 2.2]`, centred at (11.1, 1.1) with radius 1.1. The old
    // search returned 1.0 — lobe A's optimum — for every tolerance from 1e-3 down
    // to 1e-9.
    let dumbbell = Polygon::new(
        LineString::from(vec![
            (0.0, 0.0),
            (2.0, 0.0),
            (2.0, 0.9),
            (10.0, 0.9),
            (10.0, 0.0),
            (12.2, 0.0),
            (12.2, 2.2),
            (10.0, 2.2),
            (10.0, 1.1),
            (2.0, 1.1),
            (2.0, 2.0),
            (0.0, 2.0),
            (0.0, 0.0),
        ]),
        vec![],
    );

    for tol in [1e-2_f64, 1e-4, 1e-6, 1e-9] {
        let c = dumbbell
            .maximum_inscribed_circle(tol)
            .expect("the dumbbell has an interior");
        assert!(
            (c.radius - 1.1).abs() <= tol.max(1e-9),
            "tolerance {tol:.0e}: radius {} is not within {tol:.0e} of the true 1.1 \
             (centre {:?}) — 1.0 means the search settled in the smaller lobe",
            c.radius,
            c.centre
        );
        assert!(
            c.centre.x > 10.0,
            "tolerance {tol:.0e}: the optimum is in the right-hand lobe, got centre {:?}",
            c.centre
        );
        // The circle must actually fit: no boundary vertex inside it.
        for v in &dumbbell.exterior().0 {
            let (dx, dy) = (v.x - c.centre.x, v.y - c.centre.y);
            assert!(
                (dx * dx + dy * dy).sqrt() >= c.radius * (1.0 - 1e-9) - 1e-9,
                "tolerance {tol:.0e}: vertex {v:?} is inside the inscribed circle"
            );
        }
    }

    // Tightening the tolerance must not make the answer worse — monotone
    // refinement is the property a caller relies on when choosing one.
    let coarse = dumbbell.maximum_inscribed_circle(1e-2).unwrap().radius;
    let fine = dumbbell.maximum_inscribed_circle(1e-9).unwrap().radius;
    assert!(
        fine >= coarse - 1e-2,
        "a tighter tolerance returned a smaller circle: {fine} < {coarse}"
    );
}
