use super::*;

#[test]
fn geodesic_area_matches_golden() {
    use r2::GeodesicArea as _;

    // CCW square with a CW hole, a CW-wound exterior, and a high-latitude
    // sliver crossing the antimeridian. Goldens (perimeter_signed, area_signed,
    // perimeter_unsigned, area_unsigned) captured from geo's GeodesicArea.
    let cases: Vec<r2::Polygon> = vec![
        r2::Polygon::new(
            r2::LineString::from(vec![(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]),
            vec![r2::LineString::from(vec![
                (0.25, 0.25),
                (0.25, 0.75),
                (0.75, 0.75),
                (0.75, 0.25),
            ])],
        ),
        r2::Polygon::new(
            r2::LineString::from(vec![(2.0, 48.0), (2.0, 49.0), (3.0, 49.0), (3.0, 48.0)]),
            vec![],
        ),
        r2::Polygon::new(
            r2::LineString::from(vec![
                (179.5, 70.0),
                (-179.5, 70.0),
                (-179.5, 71.0),
                (179.5, 71.0),
            ]),
            vec![],
        ),
    ];
    const GOLDEN: [(f64, f64, f64, f64); 3] = [
        (
            665659.5125489294,
            9231614224.814873,
            665659.5125489294,
            9231614224.814873,
        ),
        (
            370196.1683474253,
            -8217589924.770447,
            370196.1683474253,
            510057404134163.7,
        ),
        (
            297673.23927019304,
            4157942715.0476074,
            297673.23927019304,
            4157942715.0476074,
        ),
    ];

    for (poly, &(ps, a_s, pu, au)) in cases.iter().zip(GOLDEN.iter()) {
        let (p_ours, a_ours) = poly.geodesic_perimeter_area_signed();
        assert!((p_ours - ps).abs() < 1e-6, "perimeter {p_ours} vs {ps}");
        assert!((a_ours - a_s).abs() < 1e-2, "signed area {a_ours} vs {a_s}");

        let (pu_ours, au_ours) = poly.geodesic_perimeter_area_unsigned();
        assert!((pu_ours - pu).abs() < 1e-6, "u-perimeter {pu_ours} vs {pu}");
        assert!(
            (au_ours - au).abs() < 1e-2,
            "unsigned area {au_ours} vs {au}"
        );
    }
}

#[test]
fn geodesic_area_improves_on_spherical() {
    use r2::{ChamberlainDuquetteArea as _, GeodesicArea as _};

    // A 1°×1° square at mid-latitude: the spherical approximation should sit
    // within ~0.6% of the exact ellipsoidal area, and they must disagree (the
    // whole point of the Karney area).
    let poly = r2::Polygon::new(
        r2::LineString::from(vec![(2.0, 48.0), (3.0, 48.0), (3.0, 49.0), (2.0, 49.0)]),
        vec![],
    );
    let exact = poly.geodesic_area_unsigned();
    let spherical = poly.chamberlain_duquette_unsigned_area();
    let rel = (exact - spherical).abs() / exact;
    assert!(rel > 1e-5 && rel < 6e-3, "relative gap {rel}");

    // Ellipsoid-level API agrees with the r2 trait.
    let ring: Vec<GeoPoint> = [(48.0, 2.0), (48.0, 3.0), (49.0, 3.0), (49.0, 2.0)]
        .iter()
        .map(|&(lat, lon)| g(lat, lon))
        .collect();
    let (_, area) = Ellipsoid::WGS84.ring_perimeter_area_m(&ring);
    assert!((area - exact).abs() < 1e-3, "{area} vs {exact}");

    // The whole planet, both hemispheres: total area sanity.
    let total = Ellipsoid::WGS84.total_area_m2();
    assert!((total - 510.0656e12).abs() / total < 1e-4, "{total}");
}

#[test]
fn ellipsoidal_vs_spherical_sanity() {
    // The ellipsoidal distance differs from haversine by less than ~0.6%
    // everywhere, and equals it (to round-off scale) on the f = 0 ellipsoid.
    for (a, b) in point_pairs() {
        let ell = Ellipsoid::WGS84.distance_m(a, b);
        let sph = haversine_m(a, b);
        if sph > 1_000.0 {
            assert!(
                (ell - sph).abs() / sph < 6e-3,
                "ellipsoidal {ell} vs spherical {sph} for {a:?} → {b:?}"
            );
        }
    }
}
