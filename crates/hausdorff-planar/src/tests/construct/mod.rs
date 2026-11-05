//! Geometry in, new geometry out — the group where the hard failures live.
//!
//! Producing geometry means deciding topology, so the characteristic defect here
//! is not a wrong number but a wrong *structure*: a dropped component, a hole
//! attached to the wrong shell, a ring traced inside out. Those are invisible to a
//! tolerance check and obvious to a structural one, which is why these suites
//! assert component and hole counts alongside areas.
//!
//! Two kinds of oracle carry the group. [`bool_ops`] uses **set algebra** —
//! inclusion–exclusion, idempotence, De Morgan — which ties the four operations to
//! each other and to `area` without any reference implementation.
//! [`triangulate`] and [`closest_point`] use **independent reimplementations**:
//! brute-force nearest-site scans and a closed-form point-to-segment projection,
//! written in the test and deliberately slower than the code they check.

mod bool_ops;
mod closest_point;
mod concave_hull;
mod convex_hull;
mod interior_point;
mod interpolate;
mod polygonize;
mod triangulate;
