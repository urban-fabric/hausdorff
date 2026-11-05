//! `space` — the foundational unit-less spaces, organized after Google's S2
//! into dimensional pieces: [`r1`] for the real line (intervals), [`s1`] for the
//! unit circle (angles), and [`r3`] for the ℝ³ vector substrate the unit sphere
//! is built on. The planar `r2` and spherical `s2` geometries live in
//! `hausdorff-planar` / `hausdorff-sphere`, each under its own crate's `space`
//! module; nothing re-assembles the five into one namespace, so a consumer
//! imports the spaces it uses. All carry no physical unit (metres enter only in
//! [`crate::earth`]); angles live in [`s1`], not `euclid-units`.

pub mod r1;
pub mod r3;
pub mod s1;
pub mod sphere;
pub mod surface;
