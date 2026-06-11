//! `space` — the unit-less geometry spaces, organized after Google's S2 into
//! dimensional pieces: [`r2`]/[`r3`] for Euclidean geometry and [`s1`]/[`s2`]
//! for unit-circle/unit-sphere spherical geometry. All carry no physical unit
//! (metres enter only in [`crate::earth`]); angles live in [`s1`], not `euclid`.

pub mod r2;
pub mod r3;
pub mod s1;
pub mod s2;
