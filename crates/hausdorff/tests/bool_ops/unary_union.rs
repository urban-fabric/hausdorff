use super::*;

mod unary_union_conformance {
    use super::*;
    use crate::basic::Lcg;
    use crate::random_rects::rand_rect;
    use hausdorff::space::r2::bool_ops::unary_union;

    /// Independent oracle for the cascaded `unary_union`: a linear left-fold of
    /// pairwise unions (a different code path through the same clipper).
    fn linear_fold(polys: &[Polygon]) -> MultiPolygon {
        polys.iter().fold(MultiPolygon::default(), |acc, p| {
            acc.union(&MultiPolygon::new(vec![p.clone()]))
        })
    }

    fn check_unary(polys: &[Vec<(f64, f64)>]) -> MultiPolygon {
        let ours: Vec<Polygon> = polys.iter().map(|p| hpoly(&[p])).collect();
        let our_union = unary_union(&ours);
        let folded = linear_fold(&ours);
        assert!(
            (ha(&our_union) - ha(&folded)).abs() < EPS,
            "unary_union area {} vs linear fold {}",
            ha(&our_union),
            ha(&folded)
        );
        our_union
    }

    #[test]
    fn overlapping_grid_of_squares() {
        let mut polys = Vec::new();
        for i in 0..4 {
            for j in 0..4 {
                let (x, y) = (i as f64, j as f64);
                polys.push(sq(x, y, x + 1.5, y + 1.5));
            }
        }
        let out = check_unary(&polys);
        assert_eq!(out.len(), 1);
        assert!((ha(&out) - 4.5 * 4.5).abs() < EPS);
    }

    #[test]
    fn touching_grid_of_squares() {
        let mut polys = Vec::new();
        for i in 0..3 {
            for j in 0..3 {
                let (x, y) = (i as f64, j as f64);
                polys.push(sq(x, y, x + 1.0, y + 1.0));
            }
        }
        let out = check_unary(&polys);
        assert_eq!(out.len(), 1);
        assert!((ha(&out) - 9.0).abs() < EPS);
    }

    #[test]
    fn disjoint_squares() {
        let polys: Vec<_> = (0..6)
            .map(|i| {
                let x = 3.0 * i as f64;
                sq(x, 0.0, x + 1.0, 1.0)
            })
            .collect();
        let out = check_unary(&polys);
        assert_eq!(out.len(), 6);
        assert!((ha(&out) - 6.0).abs() < EPS);
    }

    #[test]
    fn nested_rings() {
        let donut = |lo: f64, hi: f64| {
            let outer = sq(lo, lo, hi, hi);
            let inner: Vec<(f64, f64)> = sq(lo + 1.0, lo + 1.0, hi - 1.0, hi - 1.0)
                .into_iter()
                .rev()
                .collect();
            hpoly(&[&outer, &inner])
        };
        let polys = vec![
            donut(0.0, 12.0),
            donut(2.0, 10.0),
            hpoly(&[&sq(5.0, 5.0, 7.0, 7.0)]),
        ];
        let our_union = unary_union(&polys);
        assert!(
            (ha(&our_union) - ha(&linear_fold(&polys))).abs() < EPS,
            "nested rings area {} vs linear fold {}",
            ha(&our_union),
            ha(&linear_fold(&polys))
        );
    }

    #[test]
    fn multipolygon_operands() {
        let polys = vec![
            hpoly(&[&sq(0.0, 0.0, 2.0, 2.0)]),
            hpoly(&[&sq(1.0, 1.0, 3.0, 3.0)]),
            hpoly(&[&sq(10.0, 0.0, 11.0, 1.0)]),
        ];
        let our_union = unary_union(&polys);
        assert!((ha(&our_union) - ha(&linear_fold(&polys))).abs() < EPS);
    }

    #[test]
    fn area_is_conserved_on_known_fixture() {
        let a = sq(0.0, 0.0, 4.0, 4.0);
        let b = sq(2.0, 2.0, 6.0, 6.0);
        let c = sq(10.0, 10.0, 12.0, 12.0);
        let polys = vec![hpoly(&[&a]), hpoly(&[&b]), hpoly(&[&c])];
        let out = unary_union(&polys);
        let expected = 16.0 + 16.0 - 4.0 + 4.0;
        assert!(
            (ha(&out) - expected).abs() < EPS,
            "area {} vs expected {expected}",
            ha(&out)
        );
    }

    #[test]
    fn matches_linear_pairwise_fold() {
        let mut rng = Lcg(0x5151_5151_5151_5151);
        for _ in 0..20 {
            let polys: Vec<Polygon> = (0..8)
                .map(|_| hpoly(&[&sq2(&rand_rect(&mut rng, 10))]))
                .collect();
            let cascaded = unary_union(&polys);
            assert!(
                (ha(&cascaded) - ha(&linear_fold(&polys))).abs() < EPS,
                "cascaded {} vs linear {}",
                ha(&cascaded),
                ha(&linear_fold(&polys))
            );
        }
    }

    fn sq2(r: &(f64, f64, f64, f64)) -> Vec<(f64, f64)> {
        sq(r.0, r.1, r.2, r.3)
    }
}
