//! Cross-cutting suites for the kernel crate — the ones that are about the
//! *crate's* promise rather than about any one formula.
//!
//! | suite | question it answers |
//! |---|---|
//! | [`f32_parity`] | is the crate actually carrier-generic, or only generic on paper? |
//!
//! Formula-specific tests live beside their implementations (`earth/tests.rs`,
//! `planar/tests.rs`, `projection/tests.rs`), where they can reach the private
//! series and crossover constants. Anything reachable through the public surface
//! belongs here.

mod f32_parity;
