//! Internal: scalar-first SIMD dispatch for batch kernels, after RustFFT's
//! runtime feature selection.
//!
//! Each kernel has a portable **scalar reference** — the source of truth — plus,
//! where the CPU supports it, a vectorized path chosen at runtime via feature
//! detection. The vector path applies the *same* IEEE-754 operations, four `f64`
//! lanes at a time, so it is **bit-identical** to the scalar reference; the choice
//! is invisible to callers. Non-x86-64 targets (and x86-64 CPUs without AVX) use
//! the scalar path. The `unsafe` needed for the intrinsics is confined to this
//! module.
#![allow(unsafe_code)]

/// Euclidean distance from each `(xs[i], ys[i])` to the fixed point `(tx, ty)`,
/// written into `out`. `xs`, `ys`, and `out` must all have length `out.len()`.
///
/// Runs the AVX kernel on x86-64 CPUs that support it, else the scalar reference.
pub(crate) fn euclidean_distances(xs: &[f64], ys: &[f64], tx: f64, ty: f64, out: &mut [f64]) {
    debug_assert_eq!(xs.len(), out.len());
    debug_assert_eq!(ys.len(), out.len());

    #[cfg(target_arch = "x86_64")]
    {
        if std::is_x86_feature_detected!("avx") {
            // SAFETY: the `avx` target feature was just confirmed at runtime, and
            // the kernel only does in-bounds loads/stores (all slices length n).
            unsafe {
                x86::euclidean_distances_avx(xs, ys, tx, ty, out);
            }
            return;
        }
    }
    scalar::euclidean_distances(xs, ys, tx, ty, out);
}

mod scalar {
    /// Portable reference: one independent distance per element. Uses separate
    /// `mul`/`add` (not a fused `mul_add`) so the vector path can match it bit
    /// for bit.
    pub(super) fn euclidean_distances(xs: &[f64], ys: &[f64], tx: f64, ty: f64, out: &mut [f64]) {
        for ((o, &x), &y) in out.iter_mut().zip(xs).zip(ys) {
            let (dx, dy) = (x - tx, y - ty);
            *o = (dx * dx + dy * dy).sqrt();
        }
    }
}

#[cfg(target_arch = "x86_64")]
mod x86 {
    use std::arch::x86_64::{
        _mm256_add_pd, _mm256_loadu_pd, _mm256_mul_pd, _mm256_set1_pd, _mm256_sqrt_pd,
        _mm256_storeu_pd, _mm256_sub_pd,
    };

    /// AVX kernel: four `f64` lanes per iteration, with a scalar tail. Must only
    /// be called when the `avx` feature is available.
    #[target_feature(enable = "avx")]
    pub(super) unsafe fn euclidean_distances_avx(
        xs: &[f64],
        ys: &[f64],
        tx: f64,
        ty: f64,
        out: &mut [f64],
    ) {
        let n = out.len();
        let chunks = n / 4;
        // SAFETY: every load/store covers `[i, i+4)` with `i + 4 <= chunks*4 <= n`,
        // and xs/ys/out all have length n, so all accesses are in bounds.
        unsafe {
            let txv = _mm256_set1_pd(tx);
            let tyv = _mm256_set1_pd(ty);
            for c in 0..chunks {
                let i = c * 4;
                let dx = _mm256_sub_pd(_mm256_loadu_pd(xs.as_ptr().add(i)), txv);
                let dy = _mm256_sub_pd(_mm256_loadu_pd(ys.as_ptr().add(i)), tyv);
                let d2 = _mm256_add_pd(_mm256_mul_pd(dx, dx), _mm256_mul_pd(dy, dy));
                _mm256_storeu_pd(out.as_mut_ptr().add(i), _mm256_sqrt_pd(d2));
            }
        }
        for i in (chunks * 4)..n {
            let (dx, dy) = (xs[i] - tx, ys[i] - ty);
            out[i] = (dx * dx + dy * dy).sqrt();
        }
    }
}
