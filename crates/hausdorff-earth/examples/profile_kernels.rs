//! A speed profile of the crate's kernels.
//!
//! Run: `cargo run -p hausdorff-earth --release --example profile_kernels`
//!
//! Prints a flat per-call timing table over every kernel, sorted slowest-first.

use std::hint::black_box;
use std::time::Instant;

use hausdorff_earth::earth::batch::LonLatBatch;
use hausdorff_earth::earth::{GeoPoint, haversine_m};
use hausdorff_index::index::geohash;
use hausdorff_index::index::h3::H3Cell;
use hausdorff_planar::format::{geojson, wkb, wkt};
use hausdorff_planar::space::r2::relate;
use hausdorff_planar::space::r2::{
    Area, BooleanOps, BoundingRect, Centroid, ConcaveHull, Contains, ConvexHull, Destination,
    Euclidean, Haversine, Intersects, PointBatch, Simplify, StrTree, euclidean_distance,
    frechet_distance, hausdorff_distance,
};
use hausdorff_planar::space::r2::{
    Coord, Geometry, LineString, MultiPoint, MultiPolygon, Point, Polygon, Rect,
};

/// Deterministic pseudo-random coordinate in `[-half, half]`.
struct Rng(u64);
impl Rng {
    fn next(&mut self) -> u64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        self.0 >> 33
    }
    fn f(&mut self, half: f64) -> f64 {
        ((self.next() & 0xFFFF_FFFF) as f64 / 4_294_967_295.0 * 2.0 - 1.0) * half
    }
}

/// A closed `n`-gon ring on a circle.
fn ring(cx: f64, cy: f64, r: f64, n: usize) -> Vec<(f64, f64)> {
    (0..=n)
        .map(|i| {
            let t = i as f64 / n as f64 * std::f64::consts::TAU;
            (cx + r * t.cos(), cy + r * t.sin())
        })
        .collect()
}

fn poly(c: Vec<(f64, f64)>) -> Polygon {
    Polygon::new(LineString::from(c), vec![])
}

struct Measurement {
    name: String,
    mean_ns: f64,
    min_ns: f64,
    max_ns: f64,
}

impl Measurement {
    fn mean_ns(&self) -> f64 {
        self.mean_ns
    }
}

/// Time one kernel: `iters` total calls, split into 50 samples (so the per-call
/// distribution amortizes the clock read for ns-scale kernels).
fn prof(out: &mut Vec<Measurement>, name: &str, iters: u64, mut f: impl FnMut()) {
    let samples = 50usize;
    let batch = (iters / samples as u64).max(1);
    for _ in 0..batch {
        f();
    }

    let mut sample_ns = Vec::with_capacity(samples);
    for _ in 0..samples {
        let start = Instant::now();
        for _ in 0..batch {
            f();
        }
        sample_ns.push(start.elapsed().as_nanos() as f64 / batch as f64);
    }

    let mean_ns = sample_ns.iter().sum::<f64>() / sample_ns.len() as f64;
    let min_ns = sample_ns.iter().copied().fold(f64::INFINITY, f64::min);
    let max_ns = sample_ns.iter().copied().fold(0.0, f64::max);
    out.push(Measurement {
        name: name.to_owned(),
        mean_ns,
        min_ns,
        max_ns,
    });
}

fn render_table(rows: &[Measurement]) -> String {
    let name_width = rows
        .iter()
        .map(|m| m.name.len())
        .max()
        .unwrap_or("kernel".len())
        .max("kernel".len());
    let mut out = format!(
        "{:<name_width$} {:>14} {:>14} {:>14}\n",
        "kernel",
        "mean ns",
        "min ns",
        "max ns",
        name_width = name_width
    );
    out.push_str(&format!(
        "{:-<name_width$} {:-<14} {:-<14} {:-<14}\n",
        "",
        "",
        "",
        "",
        name_width = name_width
    ));
    for m in rows {
        out.push_str(&format!(
            "{:<name_width$} {:>14.1} {:>14.1} {:>14.1}\n",
            m.name,
            m.mean_ns,
            m.min_ns,
            m.max_ns,
            name_width = name_width
        ));
    }
    out
}

fn main() {
    let mut rng = Rng(0x00C0_FFEE_1234_5678);

    let pa = poly(ring(0.0, 0.0, 5.0, 64));
    let ga = Geometry::Polygon(pa.clone());
    let gb = Geometry::Polygon(poly(ring(20.0, 3.0, 4.0, 64)));
    let ov_a = Geometry::Polygon(poly(ring(0.0, 0.0, 10.0, 64)));
    let ov_b = Geometry::Polygon(poly(ring(8.0, 0.0, 10.0, 64)));
    let pt_in = Coord::new(0.0, 0.0);
    let line_a = Geometry::LineString(LineString::from(ring(0.0, 0.0, 10.0, 64)));
    let line_b = Geometry::LineString(LineString::from(ring(8.0, 1.0, 9.0, 64)));

    let hull_pts: Vec<(f64, f64)> = (0..500).map(|_| (rng.f(100.0), rng.f(100.0))).collect();
    let mp_hull = MultiPoint::from(hull_pts.clone());
    let wiggly: Vec<(f64, f64)> = (0..1000)
        .map(|i| (i as f64 * 0.1, (i as f64 * 0.3).sin() * 5.0))
        .collect();
    let wig_ls = LineString::from(wiggly);

    let cloud_a = Geometry::MultiPoint(MultiPoint::from(
        (0..1000)
            .map(|_| (rng.f(100.0), rng.f(100.0)))
            .collect::<Vec<_>>(),
    ));
    let cloud_b = Geometry::MultiPoint(MultiPoint::from(
        (0..1000)
            .map(|_| (rng.f(100.0), rng.f(100.0)))
            .collect::<Vec<_>>(),
    ));

    let coords: Vec<Coord> = (0..1024)
        .map(|_| Coord::new(rng.f(1000.0), rng.f(1000.0)))
        .collect();
    let batch = PointBatch::from_coords(&coords);
    let geo_pts: Vec<GeoPoint> = (0..1024)
        .map(|_| GeoPoint::from_lat_lon(rng.f(80.0), rng.f(170.0)))
        .collect();
    let lonlat = LonLatBatch::from_points(&geo_pts);

    let (g1, g2) = (
        GeoPoint::from_lat_lon(52.52, 13.405),
        GeoPoint::from_lat_lon(48.857, 2.352),
    );
    let mp_sq = |x0, y0, x1, y1| {
        MultiPolygon::new(vec![poly(vec![
            (x0, y0),
            (x1, y0),
            (x1, y1),
            (x0, y1),
            (x0, y0),
        ])])
    };
    let (m1, m2) = (mp_sq(0.0, 0.0, 4.0, 4.0), mp_sq(2.0, 2.0, 6.0, 6.0));

    let idx_items: Vec<(Rect, usize)> = (0..5000)
        .map(|i| {
            let (x, y) = (rng.f(1000.0), rng.f(1000.0));
            (Rect::new(Coord::new(x, y), Coord::new(x + 2.0, y + 2.0)), i)
        })
        .collect();
    let tree = StrTree::new(idx_items.clone());
    let q = Rect::new(Coord::new(0.0, 0.0), Coord::new(50.0, 50.0));
    let hash = geohash::encode(13.405, 52.52, 9);
    let cell = H3Cell::for_point(52.52, 13.405, 9).unwrap();
    let wkb_bytes = wkb::to_wkb(&ga);
    let wkt_str = wkt::to_wkt(&ga);
    let gj_str = geojson::to_geojson(&ga);

    let mut r: Vec<Measurement> = Vec::new();

    prof(&mut r, "haversine_m", 1_000_000, || {
        black_box(haversine_m(g1, g2));
    });
    prof(&mut r, "destination (geodesic)", 1_000_000, || {
        black_box(Haversine.destination(Point::new(13.405, 52.52), 45.0, 1000.0));
    });
    prof(&mut r, "area (64-gon)", 1_000_000, || {
        black_box(pa.unsigned_area());
    });
    prof(&mut r, "centroid (64-gon)", 1_000_000, || {
        black_box(pa.centroid());
    });
    prof(&mut r, "bounding_rect (64-gon)", 1_000_000, || {
        black_box(pa.bounding_rect());
    });
    prof(&mut r, "point_in_polygon", 1_000_000, || {
        black_box(pa.contains(black_box(&pt_in)));
    });

    prof(&mut r, "geohash_encode (p9)", 1_000_000, || {
        black_box(geohash::encode(13.405, 52.52, 9));
    });
    prof(&mut r, "geohash_decode (p9)", 1_000_000, || {
        black_box(geohash::decode(&hash));
    });
    prof(&mut r, "h3_for_point (res9)", 500_000, || {
        black_box(H3Cell::for_point(52.52, 13.405, 9));
    });
    prof(&mut r, "h3_boundary", 500_000, || {
        black_box(cell.boundary_lat_lng());
    });

    prof(&mut r, "intersects (disjoint)", 500_000, || {
        black_box(ga.intersects(&gb));
    });
    prof(&mut r, "euclidean_distance (64-gon)", 200_000, || {
        black_box(euclidean_distance(&ga, &gb));
    });
    prof(&mut r, "relate (disjoint)", 200_000, || {
        black_box(relate(&ga, &gb).is_intersects());
    });
    prof(&mut r, "relate (overlap 64-gon)", 20_000, || {
        black_box(relate(&ov_a, &ov_b).is_intersects());
    });
    prof(&mut r, "min_dist clouds x1000", 2_000, || {
        black_box(euclidean_distance(&cloud_a, &cloud_b));
    });
    prof(&mut r, "hausdorff_distance (64-gon)", 100_000, || {
        black_box(hausdorff_distance(&Euclidean, &ga, &gb));
    });
    prof(&mut r, "frechet_distance (64-pt)", 100_000, || {
        black_box(frechet_distance(&Euclidean, &line_a, &line_b));
    });

    prof(&mut r, "convex_hull x500", 20_000, || {
        black_box(mp_hull.convex_hull());
    });
    prof(&mut r, "concave_hull x500", 5_000, || {
        black_box(mp_hull.concave_hull(2.0));
    });
    prof(&mut r, "simplify_rdp x1000", 50_000, || {
        black_box(wig_ls.simplify(0.5));
    });
    prof(&mut r, "bool_intersection", 100_000, || {
        black_box(m1.intersection(&m2));
    });

    prof(&mut r, "dist_fan x1024 (SIMD)", 50_000, || {
        black_box(batch.distances_to(Coord::new(3.0, -7.0)));
    });
    prof(&mut r, "haversine_batch x1024", 20_000, || {
        black_box(lonlat.haversine_to(g1));
    });

    prof(&mut r, "strtree_build x5000", 5_000, || {
        black_box(StrTree::new(idx_items.clone()));
    });
    prof(&mut r, "strtree_query x5000", 200_000, || {
        black_box(tree.query(q));
    });

    prof(&mut r, "wkb_roundtrip (64-gon)", 100_000, || {
        black_box(wkb::from_wkb(black_box(&wkb_bytes)).unwrap());
    });
    prof(&mut r, "wkt_parse (64-gon)", 50_000, || {
        black_box(wkt::from_wkt(black_box(&wkt_str)).unwrap());
    });
    prof(&mut r, "geojson_parse (64-gon)", 50_000, || {
        black_box(geojson::from_geojson(black_box(&gj_str)).unwrap());
    });

    r.sort_by(|a, b| b.mean_ns().total_cmp(&a.mean_ns()));
    println!("\nper-call kernel timings (slowest first)\n");
    print!("{}", render_table(&r));

    let big_a = MultiPolygon::new(vec![poly(ring(0.0, 0.0, 10.0, 512))]);
    let big_b = MultiPolygon::new(vec![poly(ring(8.0, 0.0, 10.0, 512))]);
    let gbig_a = Geometry::Polygon(poly(ring(0.0, 0.0, 10.0, 512)));
    let gbig_b = Geometry::Polygon(poly(ring(8.0, 0.0, 10.0, 512)));
    let start = Instant::now();
    black_box(big_a.intersection(&big_b));
    black_box(relate(&gbig_a, &gbig_b));
    black_box(hausdorff_distance(&Euclidean, &gbig_a, &gbig_b));
    println!(
        "\noverlap 512-gon combined pass: {:.3} ms\n",
        start.elapsed().as_secs_f64() * 1_000.0
    );
}
