//! Claims that span many algorithms at once, and so belong to no single one.
//!
//! This group exists because these suites have nowhere else to go without being
//! arbitrarily filed under one of their subjects. [`scale_invariance`] checks
//! area, winding, centroid, polygonize, repair, simplification and the boolean ops
//! together, because the hazard it covers — a constant that should have been
//! relative — was one mistake repeated in five places. Filing it under `measure`
//! would have hidden that.
//!
//! [`geometry_dispatch`] is the same shape: it checks that every `Geometry` enum
//! variant reaches the right arm, across every trait implemented for the enum.
//! [`algebra`] checks the generic metric and vector-space kernels that the
//! concrete geometry then specialises.
//!
//! A suite belongs here when it would be *weakened* by being filed under one
//! subject, not merely when it touches several.

mod algebra;
mod geometry_dispatch;
mod scale_invariance;
