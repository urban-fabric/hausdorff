//! `datum` — geodetic [`Datum`]s and the [`Helmert`] / NTv2-grid [`DatumShift`]
//! that takes each to the WGS84 frame (the datum half of the [`super`] CRS
//! module): the reference ellipsoid plus the transform onto WGS84.

use crate::earth::gridshift::{self, GridShift};
use hausdorff_base::earth::{EcefPoint, Ellipsoid};

/// A **Helmert 7-parameter** (Bursa–Wolf) transform between geocentric (ECEF)
/// frames — the datum shift. Stored in SI units: translation in metres, rotation
/// in **radians**, and `scale` as the full factor (`1.0` = no scale). Construct
/// from the EPSG/PROJ `towgs84` units (metres, arc-seconds, ppm) with
/// [`bursa_wolf`](Self::bursa_wolf).
///
/// Position-vector convention: `p' = T + scale · R · p` with the small-rotation
/// matrix `R = I + [θ]ₓ`. [`inverse_apply`](Self::inverse_apply) is the **exact**
/// inverse (not the negated-parameter approximation), so `apply` round-trips.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Helmert {
    /// Translation (metres).
    pub tx: f64,
    pub ty: f64,
    pub tz: f64,
    /// Rotation (radians), position-vector convention.
    pub rx: f64,
    pub ry: f64,
    pub rz: f64,
    /// Scale factor (`1.0` for none; `1 + ppm·1e-6`).
    pub scale: f64,
}

impl Helmert {
    /// The identity transform (same frame).
    pub const IDENTITY: Helmert = Helmert {
        tx: 0.0,
        ty: 0.0,
        tz: 0.0,
        rx: 0.0,
        ry: 0.0,
        rz: 0.0,
        scale: 1.0,
    };

    /// From EPSG/PROJ `towgs84` units: translation metres, rotation **arc-seconds**,
    /// scale **ppm** (parts per million). E.g. OSGB36→WGS84 is
    /// `bursa_wolf(446.448, -125.157, 542.06, 0.15, 0.247, 0.842, -20.489)`.
    #[must_use]
    pub fn bursa_wolf(
        tx: f64,
        ty: f64,
        tz: f64,
        rx_sec: f64,
        ry_sec: f64,
        rz_sec: f64,
        ds_ppm: f64,
    ) -> Self {
        let sec = std::f64::consts::PI / (180.0 * 3600.0); // arc-seconds → radians
        Helmert {
            tx,
            ty,
            tz,
            rx: rx_sec * sec,
            ry: ry_sec * sec,
            rz: rz_sec * sec,
            scale: 1.0 + ds_ppm * 1e-6,
        }
    }

    /// True iff this is the identity (no shift) — lets transforms skip the work.
    #[must_use]
    pub fn is_identity(self) -> bool {
        self == Helmert::IDENTITY
    }

    /// Apply the transform: `p' = T + scale · R · p`.
    #[must_use]
    pub fn apply(self, p: EcefPoint) -> EcefPoint {
        if self.is_identity() {
            return p;
        }
        EcefPoint {
            x: self.tx + self.scale * (p.x - self.rz * p.y + self.ry * p.z),
            y: self.ty + self.scale * (self.rz * p.x + p.y - self.rx * p.z),
            z: self.tz + self.scale * (-self.ry * p.x + self.rx * p.y + p.z),
        }
    }

    /// The **exact** inverse of [`apply`](Self::apply): `p = R⁻¹·(p' − T)/scale`,
    /// using `R⁻¹ = (I − [θ]ₓ + θθᵀ)/(1 + |θ|²)` (the closed inverse of `I + [θ]ₓ`),
    /// so `inverse_apply(apply(p)) = p`.
    #[must_use]
    pub fn inverse_apply(self, p: EcefPoint) -> EcefPoint {
        if self.is_identity() {
            return p;
        }
        let (qx, qy, qz) = (
            (p.x - self.tx) / self.scale,
            (p.y - self.ty) / self.scale,
            (p.z - self.tz) / self.scale,
        );
        let (rx, ry, rz) = (self.rx, self.ry, self.rz);
        let d = 1.0 + rx * rx + ry * ry + rz * rz;
        EcefPoint {
            x: ((1.0 + rx * rx) * qx + (rz + rx * ry) * qy + (-ry + rx * rz) * qz) / d,
            y: ((-rz + rx * ry) * qx + (1.0 + ry * ry) * qy + (rx + ry * rz) * qz) / d,
            z: ((ry + rx * rz) * qx + (-rx + ry * rz) * qy + (1.0 + rz * rz) * qz) / d,
        }
    }
}

/// How a datum reaches WGS84: a single rigid [`Helmert`] transform (~1 m), or a
/// high-accuracy **NTv2 grid shift** (~cm) that models the residual distortion.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum DatumShift {
    /// A 7-parameter Helmert (Bursa–Wolf) ECEF transform to WGS84.
    Helmert(Helmert),
    /// An embedded NTv2 grid that shifts this datum's geographic coordinates to
    /// WGS84/ETRS89 directly (bilinear interpolation).
    Grid(&'static GridShift),
}

/// A **geodetic datum**: a reference [`Ellipsoid`] plus the [`DatumShift`] taking
/// its coordinates to the WGS84 frame. WGS84 itself has the identity Helmert.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Datum {
    /// The reference ellipsoid.
    pub ellipsoid: Ellipsoid,
    /// Transform from this datum to WGS84.
    pub to_wgs84: DatumShift,
}

impl Datum {
    /// This datum's native ECEF → WGS84 ECEF. A Helmert acts on ECEF directly; a
    /// grid shift is applied at the geographic level (ECEF → datum geodetic →
    /// grid → WGS84 geodetic → ECEF), the ellipsoidal height carried through.
    pub(super) fn to_wgs84_ecef(self, ecef: EcefPoint) -> EcefPoint {
        match self.to_wgs84 {
            DatumShift::Helmert(h) => h.apply(ecef),
            DatumShift::Grid(g) => {
                let (geo, h) = self.ellipsoid.ecef_to_geodetic(ecef);
                Ellipsoid::WGS84.geodetic_to_ecef(g.forward(geo), h)
            }
        }
    }

    /// WGS84 ECEF → this datum's native ECEF (the inverse of [`to_wgs84_ecef`]).
    pub(super) fn decode_wgs84_ecef(self, ecef: EcefPoint) -> EcefPoint {
        match self.to_wgs84 {
            DatumShift::Helmert(h) => h.inverse_apply(ecef),
            DatumShift::Grid(g) => {
                let (geo, h) = Ellipsoid::WGS84.ecef_to_geodetic(ecef);
                self.ellipsoid.geodetic_to_ecef(g.inverse(geo), h)
            }
        }
    }
}

impl Datum {
    /// WGS84 — the hub datum (identity shift).
    pub const WGS84: Datum = Datum {
        ellipsoid: Ellipsoid::WGS84,
        to_wgs84: DatumShift::Helmert(Helmert::IDENTITY),
    };

    /// GRS80 / NAD83 — within centimetres of WGS84 at the datum level; treated as
    /// a zero shift here (the ellipsoid differs only in the 12th significant
    /// figure of the flattening).
    pub const GRS80: Datum = Datum {
        ellipsoid: Ellipsoid::GRS80,
        to_wgs84: DatumShift::Helmert(Helmert::IDENTITY),
    };

    /// **DHDN** (Deutsches Hauptdreiecksnetz) on the Bessel 1841 ellipsoid, with
    /// the **BETA2007 NTv2 grid** shift to ETRS89/WGS84 — ~cm accuracy, vs ~1 m
    /// for a Helmert. The grid is embedded (no runtime files).
    #[must_use]
    pub fn dhdn() -> Datum {
        Datum {
            ellipsoid: bessel1841(),
            to_wgs84: DatumShift::Grid(gridshift::beta2007()),
        }
    }

    /// **NTF** (Nouvelle Triangulation Française) on the Clarke 1880 IGN
    /// ellipsoid, with the **ntf_r93 NTv2 grid** shift to RGF93/WGS84.
    #[must_use]
    pub fn ntf() -> Datum {
        Datum {
            ellipsoid: clarke1880ign(),
            to_wgs84: DatumShift::Grid(gridshift::ntf_r93()),
        }
    }

    /// **NZGD49** (New Zealand Geodetic Datum 1949) on the International 1924
    /// ellipsoid, with the **nzgd2kgrid0005 NTv2 grid** shift to NZGD2000/WGS84.
    #[must_use]
    pub fn nzgd49() -> Datum {
        Datum {
            ellipsoid: intl1924(),
            to_wgs84: DatumShift::Grid(gridshift::nzgd2000()),
        }
    }

    /// **OSGB36** (Ordnance Survey Great Britain 1936) on the Airy 1830 ellipsoid,
    /// with the official OSTN-mean 7-parameter shift to WGS84. The datum of the
    /// British National Grid (EPSG:27700).
    #[must_use]
    pub fn osgb36() -> Datum {
        Datum {
            ellipsoid: airy1830(),
            to_wgs84: DatumShift::Helmert(Helmert::bursa_wolf(
                446.448, -125.157, 542.06, 0.15, 0.247, 0.842, -20.489,
            )),
        }
    }

    /// **ED50** (European Datum 1950) on the International 1924 ellipsoid, with a
    /// mean western-Europe 3-parameter shift to WGS84.
    #[must_use]
    pub fn ed50() -> Datum {
        Datum {
            ellipsoid: intl1924(),
            to_wgs84: DatumShift::Helmert(Helmert::bursa_wolf(-87.0, -98.0, -121.0, 0.0, 0.0, 0.0, 0.0)),
        }
    }

    /// **NAD27** (North American Datum 1927) on the Clarke 1866 ellipsoid, with
    /// the mean CONUS 3-parameter shift to WGS84.
    #[must_use]
    pub fn nad27() -> Datum {
        Datum {
            ellipsoid: clarke1866(),
            to_wgs84: DatumShift::Helmert(Helmert::bursa_wolf(-8.0, 160.0, 176.0, 0.0, 0.0, 0.0, 0.0)),
        }
    }
}

/// Airy 1830 — the ellipsoid of OSGB36 / the British National Grid.
fn airy1830() -> Ellipsoid {
    Ellipsoid {
        a: 6_377_563.396,
        f: 1.0 / 299.324_964_6,
    }
}

/// International 1924 (Hayford) — the ellipsoid of ED50.
fn intl1924() -> Ellipsoid {
    Ellipsoid {
        a: 6_378_388.0,
        f: 1.0 / 297.0,
    }
}

/// Clarke 1866 — the ellipsoid of NAD27.
fn clarke1866() -> Ellipsoid {
    Ellipsoid {
        a: 6_378_206.4,
        f: 1.0 / 294.978_698_2,
    }
}

/// Bessel 1841 — the ellipsoid of DHDN (and many older European datums).
fn bessel1841() -> Ellipsoid {
    Ellipsoid {
        a: 6_377_397.155,
        f: 1.0 / 299.152_812_8,
    }
}

/// Clarke 1880 (IGN) — the ellipsoid of NTF (France).
fn clarke1880ign() -> Ellipsoid {
    Ellipsoid {
        a: 6_378_249.2,
        f: 1.0 / 293.466_021_293_627,
    }
}
