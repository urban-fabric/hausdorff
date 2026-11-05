use super::fixtures::seg_dist;
use super::*;

/// Distance from `(x, y)` to the nearest polygon edge is within `tol`.
pub(super) fn is_near_boundary(coords: &[(f64, f64)], x: f64, y: f64, tol: f64) -> bool {
    let n = coords.len();
    for i in 0..n {
        let (ax, ay) = coords[i];
        let (bx, by) = coords[(i + 1) % n];
        let (dx, dy) = (bx - ax, by - ay);
        let len2 = dx * dx + dy * dy;
        let t = if len2 == 0.0 {
            0.0
        } else {
            (((x - ax) * dx + (y - ay) * dy) / len2).clamp(0.0, 1.0)
        };
        let (px, py) = (ax + t * dx, ay + t * dy);
        if ((x - px).powi(2) + (y - py).powi(2)).sqrt() < tol {
            return true;
        }
    }
    false
}

pub(super) fn h_line(coords: &[(f64, f64)]) -> hausdorff_planar::space::r2::LineString {
    hausdorff_planar::space::r2::LineString::new(
        coords
            .iter()
            .map(|&(x, y)| hausdorff_planar::space::r2::Coord::new(x, y))
            .collect(),
    )
}

#[test]
fn coords_count_closes_the_ring() {
    // hausdorff's polygon iterator yields the closing coordinate (len + 1).
    for f in fixtures() {
        assert_eq!(
            h_poly(&f).coords_count(),
            f.len() + 1,
            "coords_count on {f:?}"
        );
    }
}

#[test]
fn simplify_satisfies_rdp_postconditions() {
    let wiggly = vec![
        (0.0, 0.0),
        (1.0, 0.1),
        (2.0, -0.1),
        (3.0, 5.0),
        (4.0, 0.05),
        (5.0, -0.02),
        (6.0, 0.0),
        (7.0, 8.0),
        (8.0, 0.0),
    ];
    let mut prev_len = wiggly.len() + 1;
    for eps in [0.2_f64, 1.0, 3.0] {
        let s: Vec<(f64, f64)> = h_line(&wiggly)
            .simplify(eps)
            .0
            .iter()
            .map(|c| (c.x, c.y))
            .collect();
        // Endpoints preserved, subset of the original, monotone in eps.
        assert_eq!(s[0], wiggly[0], "first endpoint kept at eps={eps}");
        assert_eq!(
            *s.last().unwrap(),
            *wiggly.last().unwrap(),
            "last endpoint kept at eps={eps}"
        );
        assert!(
            s.len() >= 2 && s.len() <= prev_len,
            "monotone vertex count at eps={eps}"
        );
        prev_len = s.len();
        for &c in &s {
            assert!(
                wiggly.contains(&c),
                "kept point {c:?} is an original at eps={eps}"
            );
        }
        // Every original point is within eps of the simplified polyline.
        let max_dev = wiggly
            .iter()
            .map(|&p| {
                s.windows(2)
                    .map(|w| seg_dist(p, w[0], w[1]))
                    .fold(f64::INFINITY, f64::min)
            })
            .fold(0.0_f64, f64::max);
        assert!(max_dev <= eps + 1e-9, "deviation {max_dev} > eps {eps}");
    }
}
