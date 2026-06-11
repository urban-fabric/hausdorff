//! `space` — this crate's slice of the unit-less geometry spaces: the
//! unit-sphere [`s2`]. The `s1`/`r3` substrates live in `hausdorff-base`, the
//! planar `r2` in `hausdorff-planar`; the `hausdorff` façade re-assembles all
//! four under one `space` module.

pub mod s2;
