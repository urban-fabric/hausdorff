//! Operations that take shapes and produce answers about, or out of, them.
//!
//! The top of the dependency graph: [`bool_ops`] is the native spherical boolean
//! engine (an arrangement, not a sweep — a sphere has no total order to sweep
//! along), and [`relate`] is the DE-9IM classifier that builds on it and on
//! [`cover`](super::cell::cover).
//!
//! Both are here rather than beside the shapes they operate on because they are
//! the only modules that compose *several* shapes at once, and because they are
//! where the delicate arithmetic concentrates — `bool_ops` carries this module's
//! stated scope limits for operands larger than a hemisphere.

pub(crate) mod bool_ops;
pub(crate) mod relate;
