//! Hierarchical span profile of the instrumented kernels, via huygens.
//!
//! Run: `cargo run --release --example profile_kernels --features profiling`
//!
//! Prints a `PerfReport` tree: each instrumented kernel and, for the boolean
//! overlay, its internal phases (fill / subdivide / connect_edges).

use std::hint::black_box;

use hausdorff::huygens::{PerfReport, ProfileMode, set_mode};
use hausdorff::space::r2::relate::relate;
use hausdorff::space::r2::{BooleanOps, Euclidean, frechet_distance, hausdorff_distance};
use hausdorff::{Geometry, LineString, MultiPolygon, Polygon};

/// A closed `n`-gon ring sampled on a circle.
fn ring(cx: f64, cy: f64, r: f64, n: usize) -> Vec<(f64, f64)> {
    (0..=n)
        .map(|i| {
            let t = i as f64 / n as f64 * std::f64::consts::TAU;
            (cx + r * t.cos(), cy + r * t.sin())
        })
        .collect()
}

fn main() {
    set_mode(ProfileMode::Heavy);

    let poly = |c: Vec<(f64, f64)>| Polygon::new(LineString::from(c), vec![]);
    // Two overlapping 512-gon discs, so the overlay/relate take their full
    // noding path (not the disjoint-bbox short-circuit).
    let a_mp = MultiPolygon::new(vec![poly(ring(0.0, 0.0, 10.0, 512))]);
    let b_mp = MultiPolygon::new(vec![poly(ring(8.0, 0.0, 10.0, 512))]);
    let a_geom = Geometry::Polygon(poly(ring(0.0, 0.0, 10.0, 512)));
    let b_geom = Geometry::Polygon(poly(ring(8.0, 0.0, 10.0, 512)));
    let la = Geometry::LineString(LineString::from(ring(0.0, 0.0, 10.0, 512)));
    let lb = Geometry::LineString(LineString::from(ring(8.0, 1.0, 9.0, 512)));

    let mut report = PerfReport::new();
    {
        let _root = report.span("kernels");
        black_box(a_mp.intersection(&b_mp));
        black_box(relate(&a_geom, &b_geom));
        black_box(hausdorff_distance(&Euclidean, &la, &lb));
        black_box(frechet_distance(&Euclidean, &la, &lb));
    }
    print!("{}", report.render_tree());
}
