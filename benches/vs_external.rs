//! Wall-clock comparison of hausdorff's kernels against the external reference
//! libraries (geo, h3o, geohash), so we know where hausdorff is faster or slower.
//!
//! Run with `cargo bench --bench vs_external` (release, `harness = false`). This
//! is a simple comparative harness — a deliberate, low-dependency stand-in for a
//! statistical tool like criterion — that times the same workload through both
//! implementations and prints a per-kernel table.

use std::hint::black_box;

use huygens::bench::{Bench, compare, fmt_ns};

use geo::{
    Area as _, BooleanOps as _, Contains as _, ConvexHull as _, Distance as _, Euclidean,
    Haversine, Relate as _, Simplify as _,
};
use geo_traits::to_geo::ToGeoGeometry as _;
use wkt::{ToWkt as _, TryFromWkt as _};

use hausdorff::earth::{GeoPoint, haversine_m};
use hausdorff::index::h3::H3Cell;
use hausdorff::space::r2::relate::relate;
use hausdorff::space::r2::{
    Area, BooleanOps, Contains, ConvexHull, PointBatch, Rect, Simplify, StrTree, euclidean_distance,
};
use hausdorff::{Coord, Geometry, LineString, MultiPoint, MultiPolygon, Polygon};

/// Timed samples collected per workload. Each sample times `iters / SAMPLES`
/// calls (the `batch`), so the total timed work matches the call budget the old
/// per-call-mean harness used, while yielding a real per-call distribution.
const SAMPLES: usize = 100;

/// Time both implementations of one kernel with the huygens benchmark harness
/// (TSC fast-clock, statistical per-call distribution) and print a comparison
/// row. `iters` is the old total call budget, split into [`SAMPLES`] samples.
fn row(name: &str, ext_name: &str, iters: u64, haus: impl FnMut(), ext: impl FnMut()) {
    let batch = (iters / SAMPLES as u64).max(1);
    let bench = Bench::new().warmup(batch).samples(SAMPLES).batch(batch);
    let h = bench.run("hausdorff", haus);
    let e = bench.run(ext_name, ext);
    // hausdorff is the candidate, the external library the baseline, so the
    // verdict ("Nx faster" / "Nx slower" / "even") reads from hausdorff's side.
    let verdict = compare(&e, &h).verdict();
    println!(
        "  {name:<26} {:>14}   {ext_name:<8} {:>14}   hausdorff {verdict}",
        fmt_ns(h.mean_ns()),
        fmt_ns(e.mean_ns()),
    );
}

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
        let r = (self.next() & 0xFFFF_FFFF) as f64 / 4_294_967_295.0;
        (r * 2.0 - 1.0) * half
    }
}

fn geo_ring(coords: &[(f64, f64)]) -> geo::LineString<f64> {
    geo::LineString(coords.iter().map(|&(x, y)| geo::Coord { x, y }).collect())
}

fn main() {
    println!(
        "\nhausdorff vs external reference libraries (per-call mean via huygens; lower = faster)\n"
    );
    println!(
        "  {:<26} {:>14}   {:<8} {:>14}",
        "kernel", "hausdorff", "external", ""
    );
    println!("  {}", "-".repeat(88));

    let mut rng = Rng(0x00C0_FFEE_1234_5678);

    // ---- haversine great-circle distance (vs geo) --------------------------
    let (a, b) = (
        GeoPoint::from_lat_lon(52.52, 13.405),
        GeoPoint::from_lat_lon(48.857, 2.352),
    );
    let (ga, gb) = (
        geo::Point::new(13.405, 52.52),
        geo::Point::new(2.352, 48.857),
    );
    row(
        "haversine_distance",
        "geo",
        2_000_000,
        || {
            black_box(haversine_m(black_box(a), black_box(b)));
        },
        || {
            black_box(Haversine.distance(black_box(ga), black_box(gb)));
        },
    );

    // ---- ST_Distance: polygon-polygon minimum distance (vs geo) ------------
    let ring = |cx: f64, cy: f64, r: f64, n: usize| -> Vec<(f64, f64)> {
        (0..=n)
            .map(|i| {
                let t = i as f64 / n as f64 * std::f64::consts::TAU;
                (cx + r * t.cos(), cy + r * t.sin())
            })
            .collect()
    };
    let (pa, pb) = (ring(0.0, 0.0, 5.0, 64), ring(20.0, 3.0, 4.0, 64));
    let ha = Geometry::Polygon(Polygon::new(LineString::from(pa.clone()), vec![]));
    let hb = Geometry::Polygon(Polygon::new(LineString::from(pb.clone()), vec![]));
    let gga = geo::Geometry::Polygon(geo::Polygon::new(geo_ring(&pa), vec![]));
    let ggb = geo::Geometry::Polygon(geo::Polygon::new(geo_ring(&pb), vec![]));
    row(
        "min_distance (ST_Distance)",
        "geo",
        200_000,
        || {
            black_box(euclidean_distance(black_box(&ha), black_box(&hb)));
        },
        || {
            black_box(Euclidean.distance(black_box(&gga), black_box(&ggb)));
        },
    );

    // ---- ST_Distance at scale: 512-gons vs geo (brute + bbox prune) -------
    let (la, lb) = (ring(0.0, 0.0, 5.0, 512), ring(20.0, 3.0, 4.0, 512));
    let hla = Geometry::Polygon(Polygon::new(LineString::from(la.clone()), vec![]));
    let hlb = Geometry::Polygon(Polygon::new(LineString::from(lb.clone()), vec![]));
    let ggla = geo::Geometry::Polygon(geo::Polygon::new(geo_ring(&la), vec![]));
    let gglb = geo::Geometry::Polygon(geo::Polygon::new(geo_ring(&lb), vec![]));
    row(
        "min_distance x512",
        "geo",
        10_000,
        || {
            black_box(euclidean_distance(black_box(&hla), black_box(&hlb)));
        },
        || {
            black_box(Euclidean.distance(black_box(&ggla), black_box(&gglb)));
        },
    );

    // ---- ST_Distance on overlapping point clouds: cost model picks the index
    //      (overlapping boxes + large), where geo's brute O(n·m) loses ---------
    let cloud_a: Vec<(f64, f64)> = (0..1000).map(|_| (rng.f(100.0), rng.f(100.0))).collect();
    let cloud_b: Vec<(f64, f64)> = (0..1000).map(|_| (rng.f(100.0), rng.f(100.0))).collect();
    let hca = Geometry::MultiPoint(MultiPoint::from(cloud_a.clone()));
    let hcb = Geometry::MultiPoint(MultiPoint::from(cloud_b.clone()));
    let gca = geo::Geometry::MultiPoint(
        cloud_a
            .iter()
            .map(|&(x, y)| geo::Point::new(x, y))
            .collect(),
    );
    let gcb = geo::Geometry::MultiPoint(
        cloud_b
            .iter()
            .map(|&(x, y)| geo::Point::new(x, y))
            .collect(),
    );
    row(
        "min_dist clouds x1000",
        "geo",
        2_000,
        || {
            black_box(euclidean_distance(black_box(&hca), black_box(&hcb)));
        },
        || {
            black_box(Euclidean.distance(black_box(&gca), black_box(&gcb)));
        },
    );

    // ---- batch distance fan: SIMD vs a geo per-point loop ------------------
    let coords: Vec<Coord> = (0..1024)
        .map(|_| Coord::new(rng.f(1000.0), rng.f(1000.0)))
        .collect();
    let batch = PointBatch::from_coords(&coords);
    let gpts: Vec<geo::Point<f64>> = coords.iter().map(|c| geo::Point::new(c.x, c.y)).collect();
    let (tgt, gtgt) = (Coord::new(3.0, -7.0), geo::Point::new(3.0, -7.0));
    row(
        "dist_fan x1024 (SIMD)",
        "geo",
        50_000,
        || {
            black_box(batch.distances_to(black_box(tgt)));
        },
        || {
            let v: Vec<f64> = gpts.iter().map(|p| Euclidean.distance(*p, gtgt)).collect();
            black_box(v);
        },
    );

    // ---- convex hull over 500 random points (vs geo) ----------------------
    let hull_pts: Vec<(f64, f64)> = (0..500).map(|_| (rng.f(100.0), rng.f(100.0))).collect();
    let hmp = MultiPoint::from(hull_pts.clone());
    let gmp: geo::MultiPoint<f64> = hull_pts
        .iter()
        .map(|&(x, y)| geo::Point::new(x, y))
        .collect();
    row(
        "convex_hull x500",
        "geo",
        20_000,
        || {
            black_box(hmp.convex_hull());
        },
        || {
            black_box(gmp.convex_hull());
        },
    );

    // ---- RDP simplify of a 1000-vertex wiggly line (vs geo) ---------------
    let wiggly: Vec<(f64, f64)> = (0..1000)
        .map(|i| (i as f64 * 0.1, (i as f64 * 0.3).sin() * 5.0))
        .collect();
    let hls = LineString::from(wiggly.clone());
    let gls = geo_ring(&wiggly);
    row(
        "simplify_rdp x1000",
        "geo",
        50_000,
        || {
            black_box(hls.simplify(black_box(0.5)));
        },
        || {
            black_box(gls.simplify(black_box(0.5)));
        },
    );

    // ---- H3 point indexing at res 9 (vs h3o) ------------------------------
    let (lat, lon) = (52.52, 13.405);
    let res9 = h3o::Resolution::try_from(9).unwrap();
    row(
        "h3_for_point (res9)",
        "h3o",
        1_000_000,
        || {
            black_box(H3Cell::for_point(black_box(lat), black_box(lon), 9));
        },
        || {
            black_box(h3o::LatLng::new(lat, lon).unwrap().to_cell(res9));
        },
    );

    // ---- geohash encode at precision 9 (vs the geohash crate) -------------
    row(
        "geohash_encode (p9)",
        "geohash",
        1_000_000,
        || {
            black_box(hausdorff::index::geohash::encode(
                black_box(lon),
                black_box(lat),
                9,
            ));
        },
        || {
            black_box(geohash::encode(geohash::Coord { x: lon, y: lat }, 9).unwrap());
        },
    );

    // ---- WKB roundtrip (encode + decode) vs the `wkb` crate ---------------
    row(
        "wkb_roundtrip (64-gon)",
        "wkb",
        100_000,
        || {
            let bytes = hausdorff::format::wkb::to_wkb(&ha);
            black_box(hausdorff::format::wkb::from_wkb(&bytes).unwrap());
        },
        || {
            let mut buf = Vec::new();
            wkb::writer::write_geometry(&mut buf, &gga, &wkb::writer::WriteOptions::default())
                .unwrap();
            black_box(wkb::reader::read_wkb(&buf).unwrap().to_geometry());
        },
    );

    // ---- WKT roundtrip (encode + decode) vs the `wkt` crate ---------------
    row(
        "wkt_roundtrip (64-gon)",
        "wkt",
        50_000,
        || {
            let s = hausdorff::format::wkt::to_wkt(&ha);
            black_box(hausdorff::format::wkt::from_wkt(&s).unwrap());
        },
        || {
            let s = gga.wkt_string();
            black_box(geo_types::Geometry::<f64>::try_from_wkt_str(&s).unwrap());
        },
    );

    // ---- polygon area (vs geo) --------------------------------------------
    let hpoly = Polygon::new(LineString::from(pa.clone()), vec![]);
    let gpoly = geo::Polygon::new(geo_ring(&pa), vec![]);
    row(
        "area (64-gon)",
        "geo",
        1_000_000,
        || {
            black_box(hpoly.unsigned_area());
        },
        || {
            black_box(gpoly.unsigned_area());
        },
    );

    // ---- point-in-polygon (vs geo) ----------------------------------------
    let (hp, gp) = (Coord::new(2.0, 1.0), geo::Point::new(2.0, 1.0));
    row(
        "point_in_polygon",
        "geo",
        1_000_000,
        || {
            black_box(hpoly.contains(black_box(&hp)));
        },
        || {
            black_box(gpoly.contains(black_box(&gp)));
        },
    );

    // ---- DE-9IM relate, polygon-polygon (vs geo) --------------------------
    row(
        "relate (poly-poly)",
        "geo",
        100_000,
        || {
            black_box(relate(black_box(&ha), black_box(&hb)).is_intersects());
        },
        || {
            black_box(gga.relate(&ggb).is_intersects());
        },
    );

    // ---- boolean intersection of overlapping squares (vs geo) -------------
    let mp = |x0: f64, y0: f64, x1: f64, y1: f64| {
        MultiPolygon::new(vec![Polygon::new(
            LineString::from(vec![(x0, y0), (x1, y0), (x1, y1), (x0, y1), (x0, y0)]),
            vec![],
        )])
    };
    let gmp_sq = |x0: f64, y0: f64, x1: f64, y1: f64| {
        geo::MultiPolygon(vec![geo::Polygon::new(
            geo_ring(&[(x0, y0), (x1, y0), (x1, y1), (x0, y1), (x0, y0)]),
            vec![],
        )])
    };
    let (hm1, hm2) = (mp(0.0, 0.0, 4.0, 4.0), mp(2.0, 2.0, 6.0, 6.0));
    let (gm1, gm2) = (gmp_sq(0.0, 0.0, 4.0, 4.0), gmp_sq(2.0, 2.0, 6.0, 6.0));
    row(
        "bool_intersection",
        "geo",
        100_000,
        || {
            black_box(hm1.intersection(black_box(&hm2)));
        },
        || {
            black_box(gm1.intersection(black_box(&gm2)));
        },
    );

    // ---- spatial range query: STR-tree vs brute-force scan ----------------
    let idx_items: Vec<(Rect, usize)> = (0..5000)
        .map(|i| {
            let (x, y) = (rng.f(1000.0), rng.f(1000.0));
            (Rect::new(Coord::new(x, y), Coord::new(x + 2.0, y + 2.0)), i)
        })
        .collect();
    let tree = StrTree::new(idx_items.clone());
    let q = Rect::new(Coord::new(0.0, 0.0), Coord::new(50.0, 50.0));
    row(
        "range_query x5000",
        "brute",
        100_000,
        || {
            black_box(tree.query(black_box(q)));
        },
        || {
            let hits: Vec<usize> = idx_items
                .iter()
                .filter(|(r, _)| {
                    r.min().x <= q.max().x
                        && q.min().x <= r.max().x
                        && r.min().y <= q.max().y
                        && q.min().y <= r.max().y
                })
                .map(|(_, i)| *i)
                .collect();
            black_box(hits);
        },
    );

    println!(
        "\n  (single-threaded; SIMD auto-detected; per-call mean over {SAMPLES} huygens samples, --release)\n"
    );
}
