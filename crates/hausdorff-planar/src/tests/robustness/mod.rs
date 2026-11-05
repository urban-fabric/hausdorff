//! Input that must not crash the process, as distinct from input that must fail.
//!
//! A separate group because the failure mode is different in kind. Unbounded
//! recursion and unvalidated allocation counts **abort** rather than panic, so they
//! cannot be caught by a caller or reported as a test failure — a regression here
//! takes the whole test binary down, which is worth recognising before debugging a
//! mysteriously truncated run.
//!
//! The assertion is therefore never "this returns the right answer" but "this
//! returns an error, and the process is still alive to report it".

mod codec_hostile;
mod hostile_geometry;
