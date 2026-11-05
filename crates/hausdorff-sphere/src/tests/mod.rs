//! Test suites for the unit-sphere geometry, organized by the surface under
//! test.
//!
//! These live here, inside the crate that owns the code, so a cell-id or
//! edge-crossing failure is reproducible without compiling the Earth layer above
//! this crate.
//!
//! | suite | question it answers |
//! |---|---|
//! | [`s2`] | does the 64-bit cell-id codec and hierarchy match the Go reference? |
//! | [`s2_regions`] | do caps, loops, and region coverings cover what they claim? |
//! | [`spherical`] | do the robust measures and predicates match analytic values? |
//!
//! Everything here is **radius-free**: areas are steradians and lengths are
//! angles, computed from unit vectors with no projection and no antimeridian or
//! pole special cases. That is the property the suites are really defending, so
//! a test that needs a planet radius to state its expectation belongs in
//! `hausdorff-earth` instead — as `sphere_is_planet_flexible` does.
//!
//! Where a suite does name `EARTH_RADIUS_M`, it is only scaling a steradian
//! result into m² to make an assertion legible in familiar units; the geometry
//! under test never sees it.

mod cellid_seams;
mod s2;
mod s2_regions;
mod spherical;
