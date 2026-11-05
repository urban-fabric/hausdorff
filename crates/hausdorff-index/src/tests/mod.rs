//! Test suites for the discrete spatial indexes, organized by index system.
//!
//! These live here, inside the crate that owns the code, so that a grid failure
//! is reproducible without building the Earth layer above it, and so a suite can
//! reach `pub(crate)` internals when the public surface is not the thing under
//! test.
//!
//! | suite | question it answers |
//! |---|---|
//! | [`h3`] | does the in-house H3 port still agree with the reference it replaced? |
//! | [`geohash`] | does the base-32 grid round-trip, and are its cell bounds right? |
//!
//! Within `h3`, `pentagons` is the suite to read first: the twelve icosahedron
//! vertices are where a hexagonal grid cannot close, so every traversal and
//! hierarchy operation special-cases them, and that is where a port is most likely
//! to be wrong.
//!
//! `h3` is the large one, and it is a **golden-corpus** suite: it was developed
//! as a differential test against the `h3o` crate, and the oracle's answers are
//! now frozen into FNV-1a checksums and exact value lists. That keeps the full
//! differential coverage with no reference crate in the build graph. The
//! consequence to know before editing it: a golden constant is *evidence*, not a
//! preference. If one starts failing, the implementation changed — find out why
//! before touching the number.

mod geohash;
mod h3;
