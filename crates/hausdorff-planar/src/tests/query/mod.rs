//! Finding geometries among many.
//!
//! An index is only ever a *prune*, so the property that matters is one-sided: it
//! may return extra candidates but must never drop a true one. Both suites check
//! against a nested-loop scan over the same input, which is the definition of the
//! answer rather than an approximation of it.

mod index;
mod join;
