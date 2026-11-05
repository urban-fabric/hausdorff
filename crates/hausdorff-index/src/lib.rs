//! Discrete spatial indexes that bucket points/cells over the globe, under
//! [`index`]: the H3 hexagonal hierarchy (ported from Uber's H3; see the repo
//! `NOTICE`), the geohash base-32 grid, and a uniform lon/lat point grid for
//! radius / k-nearest joins.

pub mod index;

#[cfg(test)]
mod tests;
