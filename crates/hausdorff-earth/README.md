# hausdorff-earth

The integrated Earth layer: projections, the EPSG-keyed CRS stack, columnar
batches, and the planar ↔ geographic bridge. Part of [hausdorff](../..).

```toml
[dependencies]
hausdorff-earth = { git = "ssh://git@github.com/urban-fabric/hausdorff.git", branch = "main" }
```

## The problem

Some Earth operations need more than one geometry layer at once, and so cannot
live in any single one of them. A map projection consumes a geodetic latitude and
produces a planar `r2` coordinate. A polygon measured natively on the curved Earth
needs the spherical excess kernels *and* a planet radius. A columnar batch of
points needs the geodesics, the planar coordinates, and the H3 encoder together.

Putting any of those in a lower crate would force that crate to depend upward, so
this is the top of the stack: the place where the planar, spherical, and index
layers are allowed to meet.

If you need only `GeoPoint`, `BBox`, and a distance, depend on **`hausdorff-base`**
instead and skip this crate and everything it pulls in.

## What is here

| module | contents |
|---|---|
| `earth::projection` | transverse Mercator, UTM, web Mercator, Lambert conformal conic, Albers equal-area, polar stereographic, local tangent plane, Tissot indicatrix |
| `earth::crs` | EPSG-keyed coordinate reference systems, datum transforms through WGS84 ECEF, and PROJ-string / WKT2 / PROJJSON parsing and emission |
| `earth::gridshift` | NTv2 grid-shift datum transforms, embedded at compile time |
| `earth::sphere` | `Sphere` and `AuthalicSphere` — the radius that turns steradians into m² |
| `earth::arc` | great-circle arcs over the S2 kernels, in lat/lon terms |
| `earth::polygon` | WGS84 polygons measured natively on the curved Earth — exact area, no projection |
| `earth::batch` | columnar lon/lat batches with bulk kernels |
| `convert` | the planar ↔ geographic bridge |

`earth` re-exports the whole `hausdorff-base` earth core, so a consumer of this
layer gets one coherent `earth` module rather than importing half of it from each
of two crates. That is the one deliberate glob re-export in the workspace.

## Transforms pivot through ECEF, not lat/lon

`crs::transform` converts between any two CRSes by going
**native → datum geodetic → datum ECEF → (Helmert or grid shift) → WGS84 ECEF**
and then back out. Routing through geocentric Cartesian rather than through
lat/lon is what lets a **datum shift** sit on the pivot: a Helmert transform is a
rigid motion in ECEF and has no clean lat/lon form. Skip the ECEF step and
cross-datum conversions are wrong by up to a few hundred metres — a plausible
error, which is the dangerous kind.

## Grid shifts are embedded, not loaded

A 7-parameter Helmert shift is one rigid transform for a whole datum, good to
about a metre. National agencies publish **NTv2 grids** that model the residual
distortion node by node, reaching centimetres.

Those `.gsb` files are committed and pulled in with `include_bytes!`, so they are
baked into the library at compile time and parsed in memory. There is **no runtime
file access** — no PROJ data directory, no paths, no mtimes, and nothing to
install or mis-configure. The committed artifact stays the verifiable official
grid, so you can diff it against the agency distribution or drop in a newer one.

The licensing consequence is real and is called out in [`NOTICE`](NOTICE): the
three embedded grids (BKG Germany, IGN France, LINZ New Zealand) remain under
their issuing agencies' terms, **not** this crate's Apache-2.0. Check those terms
before redistributing a build that embeds them.

## Equal-area work goes through the authalic sphere

For a polygon's true area on the spheroid there are two routes. The Karney series
is accurate but numerically validated rather than proved. The **authalic sphere**
is the proof-clean route: warp each latitude `φ → ξ` onto the sphere of equal total
surface area, then use the *proven* spherical-excess kernel. The two agree to
about 1e-4 %.

That path is also where a real defect lived. The authalic latitude is built on
Snyder 3-12, which carries a literal `1/(2e)` — a pole at zero eccentricity. Since
`Ellipsoid::SPHERE` is a public constant, every equal-area quantity on a spherical
datum came back `NaN`, **including the entire Albers projection**, with no unusual
input at all. The kernel is now written through `atanh`, which has no pole and is
also more accurate at small eccentricity.

## License

Licensed under the [Apache License, Version 2.0](LICENSE).
Copyright 2026 Raymond Mead. See [`NOTICE`](NOTICE) — note the embedded NTv2 grid
data, which carries its own terms.
