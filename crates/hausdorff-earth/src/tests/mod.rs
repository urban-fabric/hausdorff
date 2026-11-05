//! Test suites for the integrated Earth layer — the parts of `earth` that need
//! the planar, spherical, or index machinery beneath them.
//!
//! This crate sits at the top of the stack, so its suites are the ones that
//! legitimately span layers. A test belongs here when its expectation cannot be
//! stated without something from two crates below: a projection needs `r2`
//! metres, a batch needs H3 cells, `Sphere` needs a radius to turn steradians
//! into m².
//!
//! | suite | question it answers |
//! |---|---|
//! | [`projection`] | do the map projections round-trip and preserve what they claim? |
//! | [`projection_proj_reference`] | do they match PROJ's published reference coordinates? |
//! | [`geo_conformance`] | does the in-house stack still agree with georust `geo`'s semantics? |
//! | [`batch`] | do the columnar lon/lat bulk kernels match their scalar peers? |
//! | [`convert`] | does the planar <-> geographic bridge preserve coordinates? |
//! | [`similarity`] | are the great-circle Hausdorff/Frechet measures metrics? |
//! | [`sphere`] | does the planet-flexible metric layer scale correctly? |
//!
//! [`projection_proj_reference`] is the one to understand before adding a
//! projection. It checks forward and inverse against coordinates produced by
//! PROJ itself, which is the only external reference in the suite and the reason
//! a new projection family can be trusted at all. A new projection should gain an
//! entry there, not just a round-trip test — a round trip confirms the inverse
//! undoes the forward, and says nothing about whether either is the projection
//! it claims to be.
//!
//! [`geo_conformance`] serves a different purpose: it pins the behaviours callers
//! migrating off `geo` depend on, so a divergence shows up as a named conformance
//! failure rather than as a surprise downstream. Its large polygon fixtures are
//! embedded from `src/tests/fixtures/` with `include_str!`.

mod batch;
mod batch_parity;
mod convert;
mod crs_round_trips;
mod geo_conformance;
mod projection;
mod projection_proj_reference;
mod similarity;
mod sphere;
