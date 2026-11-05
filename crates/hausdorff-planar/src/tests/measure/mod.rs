//! Geometry in, number out: areas, centroids, extents, distances, metrics.
//!
//! The oracle is usually a **closed form**, which is what makes this group the
//! cheapest to trust — the area of a unit square is 1 and no reference
//! implementation is needed to say so. Where there is no closed form the fallback
//! is exhaustive enumeration ([`bounding_shapes`]) or a golden constant
//! ([`geodesic`]).
//!
//! [`geodesic`] and [`rhumb`] are here rather than in a geography group because
//! they are measures that happen to be spherical: they answer "how far", and the
//! only difference from [`min_distance`] is which metric trait is in play.

mod area;
mod bounding_rect;
mod bounding_shapes;
mod centroid;
mod extremes;
mod geodesic;
mod min_distance;
mod rhumb;
