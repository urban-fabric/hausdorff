//! The S2 cell hierarchy: the cube-to-sphere projection, cell ids, cell sets,
//! and covering.
//!
//! This is the half of `s2` that is an *index* rather than a geometry. [`coords`]
//! is the bottom of it — the face/uv/st quadratic projection that turns a unit
//! vector into a position on one of six cube faces — and [`cellid`] is the 64-bit
//! Hilbert-ordered encoding built on it. [`cellunion`] is a normalised set of
//! those ids, [`metric`] gives the cell-size constants (`AVG_EDGE`, `MIN_WIDTH`)
//! that let a caller choose a level, and [`cover`] approximates a geometry by
//! cells.
//!
//! [`coords`] carries most of this group's internal API, and three of its items
//! are `pub(crate)` rather than group-private because `shape` and `ops` need the
//! projection too — see the facade's cross-group block.

pub(crate) mod cellid;
pub(crate) mod cellunion;
pub(crate) mod coords;
pub(crate) mod cover;
pub(crate) mod metric;
