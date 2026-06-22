use super::*;

#[test]
fn local_ij_to_cell_round_trips_and_matches_golden() {
    let mut f = Fnv::new();
    for raw in small_sample().into_iter().take(300) {
        let origin = H3Cell::try_new(raw).expect("origin");

        // Self-consistency: to_local_ij then to_cell recovers the original.
        for cell in origin.grid_disk(2) {
            if let Some(ij) = cell.to_local_ij(origin) {
                assert_eq!(ij.to_cell(), Some(cell), "round trip {raw:#x}");
                assert_eq!(H3Cell::from_local_ij(ij), Some(cell));
            }
        }

        let base = origin.to_local_ij(origin).expect("origin local IJ");
        for di in -3i32..=3 {
            for dj in -3i32..=3 {
                let ij = H3LocalIJ {
                    anchor: origin,
                    i: base.i + di,
                    j: base.j + dj,
                };
                f.opt_u64(ij.to_cell().map(H3Cell::as_u64));
            }
        }
    }
    assert_eq!(
        f.finish(),
        golden::LOCAL_IJ_TO_CELL,
        "local_ij_to_cell golden"
    );
}

#[test]
fn great_circle_distance_matches_golden() {
    let mut f = Fnv::new();
    let cells = small_sample();
    for pair in cells.windows(2).take(500) {
        let a = H3Cell::try_new(pair[0]).expect("a");
        let b = H3Cell::try_new(pair[1]).expect("b");

        // Self-consistency: typed accessor agrees with the m accessor, distance
        // to self is zero, and the metric is symmetric.
        assert_eq!(
            a.great_circle_distance(b).value(),
            a.great_circle_distance_m(b)
        );
        assert_eq!(a.great_circle_distance_rads(a), 0.0);
        assert_eq!(
            a.great_circle_distance_rads(b),
            b.great_circle_distance_rads(a)
        );

        f.f64q(a.great_circle_distance_rads(b));
        f.f64q(a.great_circle_distance_km(b));
        f.f64q(a.great_circle_distance_m(b));
    }
    assert_eq!(f.finish(), golden::GREAT_CIRCLE, "great circle golden");
}

#[test]
fn for_points_batch_matches_scalar() {
    let lats = [0.0, 12.3, -45.6, 80.0, 91.0, f64::NAN, 51.5];
    let lngs = [0.0, 34.5, -120.0, 179.9, 10.0, 0.0, -0.12];
    for res in [0u8, 5, 9, 15] {
        let batch = H3Cell::for_points(&lats, &lngs, res);
        let scalar: Vec<_> = lats
            .iter()
            .zip(&lngs)
            .map(|(&la, &lo)| H3Cell::for_point(la, lo, res))
            .collect();
        assert_eq!(batch, scalar, "batch != scalar at res {res}");
    }
    assert_eq!(H3Cell::for_points(&lats, &lngs, 16), vec![None; lats.len()]);
    assert_eq!(H3Cell::for_points(&lats, &lngs[..2], 9).len(), 2);
}

#[test]
fn for_points_approx_is_exact_or_edge_neighbour() {
    let mut lats = Vec::new();
    let mut lngs = Vec::new();
    let mut lat = -85.0;
    while lat <= 85.0 {
        let mut lng = -175.0;
        while lng <= 175.0 {
            lats.push(lat);
            lngs.push(lng);
            lng += 7.3;
        }
        lat += 5.1;
    }

    for res in [0u8, 5, 9, 12, 15] {
        let approx = H3Cell::for_points_approx(&lats, &lngs, res);
        let exact = H3Cell::for_points(&lats, &lngs, res);
        for ((a, e), (&la, &lo)) in approx.iter().zip(&exact).zip(lats.iter().zip(&lngs)) {
            match (a, e) {
                (Some(a), Some(e)) => assert!(
                    a == e || a.is_neighbor_with(*e),
                    "approx {a} vs exact {e} at ({lo},{la}) res {res}: not equal or adjacent"
                ),
                (None, None) => {}
                _ => panic!("approx/exact disagree on validity at ({lo},{la}) res {res}"),
            }
        }
    }

    assert_eq!(
        H3Cell::for_points_approx(&lats, &lngs, 16),
        vec![None; lats.len()]
    );
}
