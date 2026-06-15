//! `space` — the foundational unit-less spaces, organized after Google's S2
//! into dimensional pieces: [`r1`] for the real line (intervals), [`s1`] for the
//! unit circle (angles), and [`r3`] for the ℝ³ vector substrate the unit sphere
//! is built on. The planar `r2` and spherical `s2` geometries live in
//! `hausdorff-planar` / `hausdorff-sphere`; the `hausdorff` façade re-assembles
//! them all under one `space` module. All carry no physical unit (metres enter
//! only in [`crate::earth`]); angles live in [`s1`], not `euclid`.

pub mod r1;
pub mod r3;
pub mod s1;
pub mod sphere;
pub mod surface;
