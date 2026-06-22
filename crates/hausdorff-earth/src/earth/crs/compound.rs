use super::*;

/// A **vertical CRS** — a one-dimensional height system.
///
/// Only the ellipsoidal datum is modelled today: [`Ellipsoidal`](Self::Ellipsoidal)
/// height is the geometric height the ECEF pivot already carries, so no grid is
/// needed. **Orthometric (geoid) height** — EGM2008, NAVD88, … — requires a
/// geoid-model grid lookup and is deferred (see `plans/active/crs-follow-ons.md`
/// §2); it would join this enum as a `Geoid` variant once the grid path lands.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum VerticalCrs {
    /// Height in metres above the reference ellipsoid — the geometric height the
    /// [`transform`] pivot carries natively. The trivial vertical datum.
    Ellipsoidal { epsg: Option<u32> },
}

/// A **compound CRS** — a horizontal [`Crs`] paired with a [`VerticalCrs`], the
/// EPSG model for "where, plus how high". With an [`Ellipsoidal`](VerticalCrs::Ellipsoidal)
/// vertical and a *geographic* horizontal CRS this is exactly a geographic-3D CRS
/// ([`as_crs`](Self::as_crs) returns the equivalent [`Crs`], so it flows through
/// [`transform`]). Projected-horizontal compounds and geoid verticals are
/// deferred (`plans/active/crs-follow-ons.md` §2).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CompoundCrs {
    /// The horizontal (geographic or projected) CRS.
    pub horizontal: Crs,
    /// The vertical (height) CRS.
    pub vertical: VerticalCrs,
}

impl CompoundCrs {
    /// A compound CRS from a horizontal and a vertical component.
    #[must_use]
    pub fn new(horizontal: Crs, vertical: VerticalCrs) -> CompoundCrs {
        CompoundCrs {
            horizontal,
            vertical,
        }
    }

    /// The equivalent single [`Crs`], when one exists: a *geographic* horizontal
    /// CRS over an [`Ellipsoidal`](VerticalCrs::Ellipsoidal) vertical is the
    /// geographic-3D CRS on the same datum (height retained). Returns `None` for
    /// the not-yet-collapsible cases (projected horizontal — there is no
    /// projected-with-height [`CrsCoord`] yet), so callers can fall back or error.
    #[must_use]
    pub fn as_crs(self) -> Option<Crs> {
        match (self.horizontal, self.vertical) {
            (Crs::Geographic { datum, .. }, VerticalCrs::Ellipsoidal { .. }) => {
                Some(Crs::Geographic {
                    datum,
                    epsg: None,
                    dim: GeoDim::D3,
                })
            }
            _ => None,
        }
    }
}
