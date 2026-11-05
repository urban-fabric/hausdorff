use super::*;
use crate::tests::support::Lcg;

fn rand_tri(rng: &mut Lcg, grid: i64) -> Vec<(f64, f64)> {
    vec![
        (rng.modulo(grid), rng.modulo(grid)),
        (rng.modulo(grid), rng.modulo(grid)),
        (rng.modulo(grid), rng.modulo(grid)),
    ]
}

fn tri_area2(t: &[(f64, f64)]) -> f64 {
    let (a, b, c) = (t[0], t[1], t[2]);
    ((b.0 - a.0) * (c.1 - a.1) - (c.0 - a.0) * (b.1 - a.1)).abs()
}

#[test]
fn stress_random_triangles() {
    let mut rng = Lcg(0xdead_beef_cafe_0001);
    let mut tested = 0;
    while tested < 300 {
        let a = rand_tri(&mut rng, 16);
        let b = rand_tri(&mut rng, 16);
        if tri_area2(&a) < 2.0 || tri_area2(&b) < 2.0 {
            continue;
        }
        // Sutherland–Hodgman gives the independent intersection-area oracle.
        let expected = convex_intersection_area(&a, &b);
        check_ops(&hmp(&[&[&a]]), &hmp(&[&[&b]]), expected);
        tested += 1;
    }
}
