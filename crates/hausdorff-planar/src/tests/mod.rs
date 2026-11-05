//! Test suites for the planar geometry, grouped by **the kind of claim they
//! make** — which is also what determines the oracle, and what you want nearby
//! when one fails.
//!
//! This is the largest suite in the workspace, because `space::r2` is the largest
//! surface and the one ported most directly from georust `geo`. These live here,
//! inside the crate that owns the code, so a planar failure is reproducible
//! without compiling the spherical, index, or Earth layers that sit above it.
//!
//! ## The taxonomy
//!
//! The first six groups mirror `space::r2`'s own, so a module and its tests are
//! always at the same coordinates — if the code is in `measure/`, the test is in
//! `measure/`. That mirroring costs nothing in expressiveness because the source
//! groups are themselves cut by question ("geometry in, number out"), not by
//! implementation.
//!
//! | group | the claim it makes | oracle it leans on |
//! |---|---|---|
//! | [`model`] | a shape can be built and read back | the type's own contract |
//! | [`measure`] | geometry in, number out | closed forms |
//! | [`predicate`] | geometry in, verdict out | relations agreeing with each other |
//! | [`construct`] | geometry in, new geometry out | set algebra, reimplementation |
//! | [`transform`] | the same geometry, moved | what did *not* change |
//! | [`query`] | finding geometries among many | a nested-loop scan |
//!
//! The last three have no counterpart in the source tree, which is exactly why
//! they are named rather than scattered:
//!
//! | group | the claim it makes |
//! |---|---|
//! | [`format`] | the codecs round-trip |
//! | [`invariant`] | a property holds across many algorithms at once |
//! | [`robustness`] | hostile input fails without taking the process down |
//!
//! Plus [`support`] for scaffolding two or more suites share, and [`tuning`] for
//! the config layer, which is not geometry.
//!
//! ## What this replaced, and why
//!
//! 31 modules sat at this level, 27 of them loose files, and the index above them
//! was a table of 18 rows — six of which bundled two to four files because they
//! answered one question, and one of which was literally "the remaining
//! algorithms". The taxonomy existed in that table; the directory did not
//! implement it. Three consequences had already followed:
//!
//! - **A bucket that was not a category.** `r2/` held fifteen files and described
//!   itself as tests "relocated out of `src/`". Half were type tests and half were
//!   algorithm tests, so neither half could be found from the name.
//! - **Duplicated coverage.** There were two `geometry_dispatch.rs` files, one at
//!   this level with eight tests and one inside `r2/` with a ninth, each with its
//!   own `square()` fixture. `Contains` was covered in two unrelated files.
//! - **A directory formed by scaffolding rather than subject.** `r2_query_ops/`
//!   held closest-point, interpolation, affine and extremes tests — three
//!   different groups' worth — solely because they shared a harness. Those helpers
//!   are now in [`support`], and the four suites are filed where their subjects
//!   are.
//!
//! ## Five house rules
//!
//! **A test lives here if it goes through the public API, and next to the code if
//! it does not.** The two locations are a deliberate split, not drift: the suites
//! in this directory reach their subject with `crate::space::r2::…`, so they
//! exercise it exactly as a consumer would and break when the *surface* breaks.
//! The eight `tests.rs` files that sit inside `space/r2` instead reach theirs with
//! `super::`, because what they check is private: a noding helper, a sweep event
//! ordering, an orientation kernel. Those cannot be tested from out here, and
//! moving them would mean widening visibility to suit the tests.
//!
//! So the question to ask of a new test is which import it needs. `super::` means
//! it belongs beside its module; `crate::` means it belongs in this directory —
//! and then the group is whichever of the nine above describes the claim.
//!
//! **Oracles are embedded, not depended on.** These suites were originally
//! differential tests against `geo`. That dependency is gone: the references are
//! now closed-form formulas (rectangle overlap), independent algorithms
//! (a Sutherland–Hodgman clipper for convex cases), set-algebra invariants that
//! tie the operations to each other, and golden constants where no closed form
//! exists. A golden is *evidence*, not a preference — if one starts failing, the
//! implementation changed, so find out why before editing the number.
//!
//! **An oracle used by one suite stays with that suite.** [`support`] is for
//! scaffolding two or more groups need. The reference closest-point implementation
//! lives beside the closest-point tests, not in `support`, because nothing else
//! can use it and the comparison is easier to audit when both halves are on one
//! screen. Getting this backwards is what produced `r2_query_ops`.
//!
//! **Large fixtures are embedded with `include_str!`**, from
//! `src/tests/fixtures/`, never read at run time. Their coordinates must survive
//! verbatim because the defects they provoke are topological; and embedding them
//! means a test cannot fail because of a working directory.
//!
//! **No file carrying tests goes over 400 lines.** Split by sub-topic into a
//! directory whose `mod.rs` header explains the split and holds the shared
//! fixtures, as [`construct::bool_ops`], [`construct::triangulate`] and
//! [`measure::bounding_shapes`] do.

// `pub(crate)` rather than private: the white-box `tests` modules that sit inside
// `space/r2` share the fixture generator too, so the scaffolding has to be
// reachable from there. It is still `#[cfg(test)]`, so none of it ships.
pub(crate) mod support;

mod construct;
mod format;
mod invariant;
mod measure;
mod model;
mod predicate;
mod query;
mod robustness;
mod transform;
mod tuning;
