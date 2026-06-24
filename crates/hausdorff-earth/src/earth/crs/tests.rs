use super::*;

fn approx(a: f64, b: f64, eps: f64) -> bool {
    (a - b).abs() < eps
}

mod core;

mod datum;

mod registry_3d;

mod serialization;
