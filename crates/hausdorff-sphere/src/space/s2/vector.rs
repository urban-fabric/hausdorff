//! The ℝ³ vector backing the S2 cube-sphere projection. The canonical type now
//! lives in [`hausdorff_base::space::r3`]; this alias keeps the cell-indexing internals reading
//! as `R3Vector`.

pub(super) use hausdorff_base::space::r3::{Axis, Vector as R3Vector};
