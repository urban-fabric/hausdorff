use super::*;
use crate::basic::Lcg;

/// A random axis-aligned rectangle on an integer grid.
pub(super) fn rand_rect(rng: &mut Lcg, grid: i64) -> (f64, f64, f64, f64) {
    let x0 = rng.coord(grid);
    let y0 = rng.coord(grid);
    let w = 1.0 + rng.coord(grid / 2);
    let h = 1.0 + rng.coord(grid / 2);
    (x0, y0, x0 + w, y0 + h)
}

#[test]
fn stress_random_rectangles() {
    let mut rng = Lcg(0x1234_5678_9abc_def0);
    for _ in 0..400 {
        let a = rand_rect(&mut rng, 12);
        let b = rand_rect(&mut rng, 12);
        check_rects(&[a], &[b]);
    }
}

#[test]
fn stress_random_multirects() {
    let mut rng = Lcg(0x0fed_cba9_8765_4321);
    for _ in 0..200 {
        let a1 = rand_rect(&mut rng, 10);
        let r2 = rand_rect(&mut rng, 10);
        let a2 = (r2.0 + 20.0, r2.1, r2.2 + 20.0, r2.3); // shifted clear of a1 (x < 10)
        let (bx, by) = (rng.coord(6), rng.coord(8));
        let b1 = (bx, by, bx + 18.0 + rng.coord(14), by + 2.0 + rng.coord(8));
        check_rects(&[a1, a2], &[b1]);
    }
}
