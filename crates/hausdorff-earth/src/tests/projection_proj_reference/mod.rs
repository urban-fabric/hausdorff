// PROJ reference values are transcribed verbatim; grouping their decimals would
// only obscure the comparison with the generator output.
#![allow(clippy::inconsistent_digit_grouping)]

//! Cross-validation of every map projection against **PROJ 9.4.0** (the
//! authoritative reference). Each reference coordinate was produced by pyproj
//! transforming from the projection's *own-datum* geographic CRS (so PROJ applies
//! the projection only, with no datum shift) — isolating the projection formula.
//!
//! Generated with:
//! ```text
//! Transformer.from_crs(src, tgt, always_xy=True).transform(lon, lat)
//! ```
//! for (tgt, src) in {(3857,4326),(32630/32633/32701,4326),(5070,4269),
//! (2154,4171),(3995/3031,4326),(27700,4277)}.

use crate::earth::GeoPoint;
use crate::earth::crs::{Crs, CrsCoord, GeoDim, transform};
use crate::earth::projection::{
    AlbersEqualArea, Ellipsoid, LambertConformalConic, PolarStereographic, Projection,
    TransverseMercator, WebMercator,
};
use hausdorff_planar::space::r2::Coord;

/// `(lat, lon, proj_easting, proj_northing)` — PROJ ground truth.
type Ref = (f64, f64, f64, f64);

fn airy1830() -> Ellipsoid {
    Ellipsoid {
        a: 6_377_563.396,
        f: 1.0 / 299.324_964_6,
    }
}

/// Forward+inverse against PROJ, reporting the worst easting/northing error and
/// the worst inverse round-trip error (metres / degrees), asserting bounds.
fn check<P: Projection>(name: &str, p: &P, refs: &[Ref], fwd_tol_m: f64) {
    let (mut max_e, mut max_n, mut max_inv) = (0.0_f64, 0.0_f64, 0.0_f64);
    for &(lat, lon, e, n) in refs {
        let got = p.project(GeoPoint::from_lat_lon(lat, lon));
        max_e = max_e.max((got.x - e).abs());
        max_n = max_n.max((got.y - n).abs());
        // Inverse PROJ's own easting/northing and compare to the input lat/lon.
        let back = p.inverse(Coord::new(e, n));
        max_inv = max_inv
            .max((back.lat - lat).abs())
            .max((back.lon - lon).abs());
    }
    eprintln!(
        "{name:12} forward max |ΔE|={max_e:.6} m |ΔN|={max_n:.6} m   inverse max Δ={max_inv:.3e}°"
    );
    assert!(
        max_e < fwd_tol_m,
        "{name}: forward ΔE {max_e} m exceeds {fwd_tol_m} m"
    );
    assert!(
        max_n < fwd_tol_m,
        "{name}: forward ΔN {max_n} m exceeds {fwd_tol_m} m"
    );
    // Inverting PROJ's coordinate must recover the original lat/lon (our forward
    // and inverse are mutually consistent and agree with PROJ's forward).
    assert!(max_inv < 1e-6, "{name}: inverse Δ {max_inv}° too large");
}

mod core;

mod conic;

mod registry;
