# hausdorff

Urban Fabric's geometry and spatial-indexing substrate. It follows the
[fleet glossary](https://github.com/urban-fabric/shannon/blob/main/GLOSSARY.md) and defines no domain vocabulary: it knows
nothing about census geography, OSM, or tiles, and every spatial claim the fleet
makes ultimately reduces to a predicate proved here.

It replaces the external geospatial stack outright — `geo`/`geo-types` for planar
geometry, golang/geo for the unit sphere, Uber H3 for the hexagonal index,
GeographicLib for ellipsoidal geodesics — so that numeric behaviour, tuning, and
correctness are ours to fix rather than to work around.

**There is no umbrella crate.** Depend on the crates you need; each one documents
itself.

## The crates

Each has its own README with the problem it solves, what is in it, and the
gotchas worth knowing before using it. Start with the one whose question matches
yours.

| crate | the question it answers |
|---|---|
| [`hausdorff-traits`](crates/hausdorff-traits) | *How is the same geometric argument shared across five different spaces?* The abstract trait surface — orientation, vector spaces, simplices, regions, intervals — and the algorithms written once over it. |
| [`hausdorff-kernels`](crates/hausdorff-kernels) | *How is one formula made to run on scalars, tensors, and query engines alike?* Carrier-generic geodesy and projection kernels with no data types attached. |
| [`hausdorff-base`](crates/hausdorff-base) | *Where do I get a point, a bbox, and a distance without compiling a geometry engine?* WGS84 primitives, spherical and ellipsoidal geodesy, and the unit-less `r1`/`s1`/`r3` substrates. The deliberate floor of the workspace. |
| [`hausdorff-planar`](crates/hausdorff-planar) | *How do planar algorithms survive real data?* The `r2` geometry model, boolean operations, DE-9IM relate, exact integer predicates, and the WKB/WKT/GeoJSON/shapefile codecs. |
| [`hausdorff-sphere`](crates/hausdorff-sphere) | *How do I get the right answer at the antimeridian and the poles?* Unit-vector spherical geometry and the S2 cell hierarchy — no projection, no seam. |
| [`hausdorff-index`](crates/hausdorff-index) | *How do I turn a position into a database key?* H3 hexagons, geohash, and a uniform lon/lat grid, for grouping, joining, and range scans. |
| [`hausdorff-earth`](crates/hausdorff-earth) | *What needs more than one layer at once?* Map projections, the EPSG-keyed CRS stack with NTv2 datum shifts, native curved-Earth polygon measures, and columnar batches. |

## How they fit together

The split exists so a consumer takes on only the surface it uses. A thin one
needing `GeoPoint`, `haversine_m`, and `BBox` depends on `hausdorff-base` alone
and never builds the planar, spherical, or index layers.

```text
{hausdorff-kernels, hausdorff-traits} → hausdorff-base → hausdorff-planar
  → {hausdorff-sphere, hausdorff-index} → hausdorff-earth
```

Reading it as a stack of responsibilities rather than a dependency graph:

- **`hausdorff-traits`** and **`hausdorff-kernels`** are abstraction and formula,
  with no data. Both are carrier-generic over
  [euclid](https://github.com/urban-fabric/euclid), so a formula written once
  runs on `f64`, `f32`, or a consumer's tensor and expression types.
- **`hausdorff-base`** is where metres become real, and the only crate most
  consumers need.
- **`hausdorff-planar`**, **`hausdorff-sphere`**, and **`hausdorff-index`** are
  three different models of the same Earth: flat, round, and bucketed. Which one
  is right depends on whether you care about the seam at ±180°, and how much
  precision you are willing to trade for a key you can `GROUP BY`.
- **`hausdorff-earth`** is where those three are allowed to meet.

The public surface is organized by the **space** it describes, and each space
lives in exactly one crate — so a module path tells you which crate to depend on:

| what you want | where it lives |
|---|---|
| `GeoPoint`, `BBox`, geodesy, ECEF | `hausdorff_base::earth` |
| `r1` intervals, `s1` angles, `r3` vectors | `hausdorff_base::space` |
| planar geometry, codecs, `algebra`, `config` | `hausdorff_planar` |
| unit-sphere geometry and the S2 cell index | `hausdorff_sphere::space::s2` |
| H3, geohash, the lon/lat grid | `hausdorff_index::index` |
| projections, CRS, batches, `Sphere` | `hausdorff_earth::earth` |
| carrier-generic formula kernels | `hausdorff_kernels` |

```rust
use hausdorff_base::earth::{GeoPoint, haversine_m};
use hausdorff_planar::space::r2::{Area, Polygon};
```

Crate names are hyphenated; in Rust they are imported with underscores. The
construction macros are at the crate root rather than under `space::r2`:
`hausdorff_planar::{coord, point, line_string, polygon}`.

## Dependencies

[euclid](https://github.com/urban-fabric/euclid) for numeric and unit contracts
(`euclid-numeric` and `euclid-units`, pinned to `^0.3.0`), and
`bytemuck`/`serde`/`serde_json` for POD layouts and the native JSON codec. There
is **no third-party geo crate in the graph**.

## Consumers

`galois` uses hausdorff for Arrow-native spatial kernels and S2 partitioning;
`hypatia` for OSM geometry and region coverings; `cantor` for point-in-polygon
resolution. Those repositories own the I/O, partition layout, and governance
around the geometry — hausdorff owns only whether the geometry is right.

## Contributing

[CONTRIBUTING.md](CONTRIBUTING.md) covers the `jj` workflow, the
branch-and-squash-merge flow, and the five local checks that stand in for CI
(there is none, by choice). In short, before pushing:

```sh
cargo fmt --all --check
cargo test --all-targets --workspace
cargo clippy --all-targets --workspace -- -D warnings
cargo test --doc --workspace
RUSTDOCFLAGS="-D warnings" cargo doc --workspace --no-deps --document-private-items
```

## License

Licensed under the [Apache License, Version 2.0](LICENSE).
Copyright 2026 Raymond Mead.

Ported code and embedded data carry their own attributions. Each crate has a
[`NOTICE`](NOTICE) scoped to what it actually contains, and the root `NOTICE` is
the union. The one to read before redistributing a binary is
[`hausdorff-earth`'s](crates/hausdorff-earth/NOTICE): the three embedded NTv2
grid files remain under their issuing agencies' terms rather than Apache-2.0.
