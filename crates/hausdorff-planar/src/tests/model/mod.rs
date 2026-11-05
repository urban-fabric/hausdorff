//! The geometry types themselves: construction, accessors, and the invariants of
//! the data structures.
//!
//! The group everything else rests on, and the only one with no algorithm in it.
//! A failure here means a shape cannot be built or read back correctly, which
//! makes every other suite's result meaningless — so these are the tests to look
//! at first when a lot of things fail at once.
//!
//! These were the type half of a directory called `r2`, whose own header
//! described it as "public-API tests for `crate::space::r2`, relocated out of
//! `src/`" — a migration landing zone rather than a category. Its algorithm half
//! went to the groups that own those algorithms.

mod coord;
mod dimensions;
mod geometry;
mod line_string;
mod macros;
mod polygon;
mod primitives;
