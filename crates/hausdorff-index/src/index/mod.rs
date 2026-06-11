//! `index` — discrete spatial indexes that bucket points/cells over the globe:
//! the [`h3`] hexagonal hierarchy, the [`geohash`] base-32 grid, and a uniform
//! lon/lat point [`grid`] for radius / k-nearest joins. (The S2 cell index lives
//! with the spherical geometry it is built from, in `hausdorff-sphere`'s `space::s2`.)

pub mod geohash;
pub mod grid;
pub mod h3;
