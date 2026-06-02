//! `index` — discrete spatial indexes that bucket points/cells over the globe:
//! the [`h3`] hexagonal hierarchy and the [`geohash`] base-32 grid. (The S2 cell
//! index lives with the spherical geometry it is built from, in
//! [`crate::space::s2`].)

pub mod geohash;
pub mod h3;
