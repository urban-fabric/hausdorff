# hausdorff

Hausdorff is Urban Fabric's geometry and spatial-indexing substrate. It follows
the [fleet glossary](../glossary.md) and defines no domain vocabulary: it knows
nothing about census geography, OSM, or tiles, and every spatial claim the fleet
makes ultimately reduces to a predicate proved here.

It replaces the external geospatial stack outright — `geo`/`geo-types` for
planar geometry, golang/geo for the unit sphere, Uber H3 for the hexagonal
index — so that numeric behavior, tuning, and correctness are ours to fix rather
than to work around.

## Layers

The public surface is organized by the space it describes, not by the crate it
lives in:

| Layer | Covers |
|---|---|
| `algebra` | Cross-cutting metric- and vector-space abstractions. |
| `space` | Unit-less continuous geometry: Euclidean `r1`/`r2`/`r3` and spherical `s1`/`s2`. Angles live here, not in `euclid`. |
| `earth` | The Earth-specific layer — WGS84 lat/lon, great-circle metres, generalization to any planet radius, and map projections. |
| `kernels` | Backend-neutral numeric formulas over `euclid::Real`, expressible by scalar, tensor, or expression backends. |
| `index` | Discrete grids that bucket points over the globe: `h3`, `geohash`, and a uniform lon/lat grid. The S2 cell index lives with its geometry in `space::s2`. |

## Crates

The workspace is flat; every crate lives directly under `crates/`. The split
exists so a thin consumer that needs only `GeoPoint`, `haversine_m`, and `BBox`
can depend on `hausdorff-base` alone, while the `hausdorff` façade preserves the
original single-crate module paths unchanged.

| Crate | Role |
|---|---|
| `hausdorff-traits` | The abstract geometry-trait surface — `Orientation`, `VectorSpace`, `Simplex`, `Interpolate`, `PointRegion`/`Prunable` — plus the generic algorithms written once over them. Upstream of every space, so each space implements the traits beside its own point and region types. |
| `hausdorff-kernels` | Backend-neutral formula kernels over `euclid`'s lazy numeric traits, deliberately free of eager records, encodings, and indexes. |
| `hausdorff-base` | Geospatial foundations: `GeoPoint`/`BBox`, haversine and Karney WGS84 geodesy, ECEF, spheres, and the unit-less `s1` angle and `r3` vector substrates. |
| `hausdorff-planar` | Planar `f64` geometry (the georust `geo` port): `space::r2` types and algorithms, the `algebra` abstractions, machine-tunable kernel parameters, and native WKB/WKT/GeoJSON/shapefile codecs. |
| `hausdorff-sphere` | Unit-sphere geometry (the golang/geo S2 port): cell ids and hierarchy, caps, loops, polylines, edge crossings, orientation predicates, measures, and region coverings. |
| `hausdorff-index` | Discrete spatial indexes: the H3 hexagonal hierarchy (Uber H3 port), the geohash base-32 grid, and a uniform lon/lat point grid for radius and k-nearest joins. |
| `hausdorff-earth` | The integrated Earth layer: great-circle arcs over the S2 kernels, map projections, geographic Hausdorff/Fréchet similarity, columnar lon/lat batches with bulk kernels, and the planar ↔ geographic bridge. |
| `hausdorff` | The façade re-assembling every layer under the original single-crate public API. |

The dependency DAG is:

```text
{hausdorff-kernels, hausdorff-traits} → hausdorff-base → hausdorff-planar
  → {hausdorff-sphere, hausdorff-index} → hausdorff-earth → hausdorff
```

Ported code retains its upstream attribution in the repository `NOTICE`.

## Consumers

Hausdorff is consumed as a Git dependency, not a path dependency. `galois` uses
it for Arrow-native spatial kernels and S2 partitioning; `hypatia` uses it for
OSM geometry and region coverings; `cantor` uses it for point-in-polygon
resolution. Those repositories own the I/O, partition layout, and governance
around the geometry — hausdorff owns only whether the geometry is right.

## Build and test

```bash
cargo build
cargo test --workspace
```
