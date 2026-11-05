//! The ℝ³ vector alias, shared by every group of `s2`.
//!
//! Stays at the top level of `s2` deliberately: `pub(super)` here means visible
//! throughout `s2`, including inside the groups, so nothing has to be widened to
//! `pub(crate)` for it.
//!
//! The ℝ³ vector backing the S2 cube-sphere projection. The canonical type now
//! lives in [`hausdorff_base::space::r3`]; this alias keeps the cell-indexing internals reading
//! as `R3Vector`.

pub(super) use hausdorff_base::space::r3::{Axis, Vector as R3Vector};
