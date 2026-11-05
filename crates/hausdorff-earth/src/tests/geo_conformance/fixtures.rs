use super::*;

/// Distance from `p` to the bounded segment `a`–`b` (RDP's metric).
pub(super) fn seg_dist(p: (f64, f64), a: (f64, f64), b: (f64, f64)) -> f64 {
    let (abx, aby) = (b.0 - a.0, b.1 - a.1);
    let len2 = abx * abx + aby * aby;
    let t = if len2 == 0.0 {
        0.0
    } else {
        (((p.0 - a.0) * abx + (p.1 - a.1) * aby) / len2).clamp(0.0, 1.0)
    };
    let (qx, qy) = (a.0 + abx * t, a.1 + aby * t);
    ((p.0 - qx).powi(2) + (p.1 - qy).powi(2)).sqrt()
}

#[test]
fn real_wkt_fixtures_conformance() {
    use hausdorff_planar::format::wkt::{from_wkt, to_wkt};

    // Golden convex-hull areas captured from the validated implementation.
    for (name, wkt, hull_golden) in [
        (
            "poly1",
            include_str!("../fixtures/poly1.wkt"),
            154873.47276695012_f64,
        ),
        (
            "poly2",
            include_str!("../fixtures/poly2.wkt"),
            481365.4849054594_f64,
        ),
    ] {
        let g = from_wkt(wkt).unwrap_or_else(|e| panic!("parse {name}: {e:?}"));
        let coords: Vec<(f64, f64)> = match &g {
            hausdorff_planar::space::r2::Geometry::LineString(ls) => {
                ls.0.iter().map(|c| (c.x, c.y)).collect()
            }
            other => panic!("{name}: expected LineString, got {other:?}"),
        };
        assert!(
            coords.len() > 50,
            "{name}: real shape ({} pts)",
            coords.len()
        );

        match from_wkt(&to_wkt(&g)).expect("re-parse") {
            hausdorff_planar::space::r2::Geometry::LineString(ls) => {
                assert_eq!(ls.0.len(), coords.len())
            }
            other => panic!("{name}: round-trip changed type to {other:?}"),
        }

        let ha = hausdorff_planar::space::r2::MultiPoint::from(coords.clone())
            .convex_hull()
            .unsigned_area();
        assert!(
            (ha - hull_golden).abs() <= hull_golden.abs() * 1e-9,
            "{name}: hull area {ha} vs {hull_golden}"
        );

        let hls = hausdorff_planar::space::r2::LineString::from(coords.clone());
        let mut prev_len = coords.len();
        for eps in [0.5_f64, 2.0, 10.0] {
            let s: Vec<(f64, f64)> = hls.simplify(eps).0.iter().map(|c| (c.x, c.y)).collect();
            assert!(
                s.len() >= 2 && s.len() <= prev_len,
                "{name}: monotone at eps={eps}"
            );
            prev_len = s.len();
            assert_eq!(s[0], coords[0], "{name}: first endpoint kept");
            assert_eq!(
                *s.last().unwrap(),
                *coords.last().unwrap(),
                "{name}: last endpoint kept"
            );
            for &c in &s {
                assert!(
                    coords.contains(&c),
                    "{name}: kept point {c:?} is an original"
                );
            }
            let max_dev = coords
                .iter()
                .map(|&p| {
                    s.windows(2)
                        .map(|w| seg_dist(p, w[0], w[1]))
                        .fold(f64::INFINITY, f64::min)
                })
                .fold(0.0_f64, f64::max);
            assert!(
                max_dev <= eps + 1e-9,
                "{name}: simplify deviation {max_dev} > eps {eps}"
            );
        }
    }
}
