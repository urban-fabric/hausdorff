//! Finding geometries among many: the spatial index and the joins built on it.
//!
//! Grouped together because they share a cost model rather than a shape. Every
//! other group here answers a question about *one* geometry (or a fixed pair);
//! these answer questions about a collection, where the work is dominated by how
//! many candidate pairs survive a prune rather than by the per-pair predicate.
//! [`join`] is the clearest case — it exists to call a predicate from the
//! `predicate` group as few times as possible.
//!
//! Modules are `pub(crate)`: the public surface is [`space::r2`](super)'s facade,
//! so this grouping can change without reaching a caller.

pub(crate) mod batch;
pub(crate) mod index;
pub(crate) mod join;
