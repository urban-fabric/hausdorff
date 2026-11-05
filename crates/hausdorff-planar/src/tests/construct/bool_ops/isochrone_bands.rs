//! Regression: the boolean sweep must converge on near-coincident vertices.
//!
//! These two fixtures are the outer and inner ring of one isochrone band —
//! concave hulls of nested reachable sets, so they share long runs of collinear
//! edges. Three near-vertical edges meet at one vertex with the incident split
//! points a single ULP apart. Before the split guard in `divide_segment`, the
//! sweep split that sliver, put the two-ULP remainder back on the queue, found
//! it to "intersect" its neighbour an ULP further along, and split again: two
//! events per split, the intersection creeping up the edge, converging on
//! nothing. It passed three million splits and 114 GB of resident events before
//! the kernel killed the process that called it.
//!
//! The failure mode is non-termination, so the work runs on a thread and the
//! test fails on timeout rather than hanging a CI run forever.

use std::sync::mpsc;
use std::time::Duration;

use super::fixture;
use crate::space::r2::{Area, BooleanOps};

/// A generous ceiling: the same difference completes in well under a
/// millisecond once the sweep converges. Anything near this bound means it has
/// stopped converging again.
const BUDGET: Duration = Duration::from_secs(30);

#[test]
fn isochrone_band_difference_converges() {
    let (tx, rx) = mpsc::channel();
    std::thread::spawn(move || {
        let outer = fixture("isochrone_band_outer");
        let inner = fixture("isochrone_band_inner");
        let difference = outer.difference(&inner);
        let intersection = outer.intersection(&inner);
        let _ = tx.send((
            outer.unsigned_area(),
            difference.unsigned_area(),
            intersection.unsigned_area(),
            difference.0.len(),
        ));
    });

    let Ok((outer_area, difference_area, intersection_area, parts)) = rx.recv_timeout(BUDGET)
    else {
        panic!(
            "band difference did not converge within {BUDGET:?} — the sweep is subdividing without bound again"
        );
    };

    // Terminating is what this test exists for. The bounds below are what a
    // terminating answer must satisfy no matter how the sweep resolves the
    // vertex; they would also catch a "fix" that terminated by dropping real
    // splits and returning an empty or oversized difference.
    assert!(parts > 0, "the difference should not be empty");
    assert!(
        difference_area > 0.0 && difference_area <= outer_area,
        "0 < A\\B ≤ A: {difference_area} vs {outer_area}"
    );
    assert!(
        intersection_area > 0.0 && intersection_area <= outer_area,
        "0 < A∩B ≤ A: {intersection_area} vs {outer_area}"
    );

    // A\B = A - A∩B. This did not hold when the guard landed — the sweep was
    // mislabelling regions on these hulls, separately from the non-termination
    // this test exists for — and is asserted now that it does.
    assert!(
        (difference_area - (outer_area - intersection_area)).abs() < 1e-9,
        "A\\B = A-A∩B: {difference_area} vs {}",
        outer_area - intersection_area
    );
}
