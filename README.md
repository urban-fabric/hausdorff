# hausdorff

Internal geospatial primitives and index implementations for UrbanFabric — an
in-house, dependency-light replacement for the `geo` / `h3o` / `s2` / `geohash`
stack, so UrbanFabric projects share one geometry core. Proprietary software;
see [`NOTICE`](NOTICE) for the open-source works it ports from.

The H3 core was the first piece, moved out of Kolmogorov so other projects can
use the same internal implementation.

## Modules

- **`index::h3`** — from-scratch Rust H3 hexagonal index. Cell
  construction/validation, lat/lon indexing, cell centers/boundaries, hierarchy
  operations, child positions, compaction/uncompaction, pentagon/base-cell
  inspection, grid disk/ring traversal, neighbor checks, grid distance,
  great-circle cell distance, local IJ coordinates (both directions), grid
  paths, the grid-path polyline cover (`H3Line`), directed edges and edge
  geometry, canonical vertex indexes, exact cell area, polygon-to-cells covers
  in all three containment modes (centroid / overlapping / full, including
  date-line-crossing rings), cell-set boundary dissolve, and the standard
  resolution metadata tables. Differential-tested against `h3o`.
- **`index::geohash`** — native geohash encode/decode of lon/lat points.
- **`space::r2`** — planar (`f64`) geometry: a replacement for `geo`/`geo-types`
  whose shapes (Coord, Point, Line, LineString, Polygon, MultiPolygon, …) and
  trait names mirror `geo` so callers migrate by changing the import path.
  Includes `Area`, `Centroid`, `BoundingRect`, `ConvexHull`, `ConcaveHull`,
  `Extremes`, `Simplify`, `Distance` / `min_distance`, `ClosestPoint`,
  `Intersects`, `Length` / `Perimeter` / `Densify` (under the `Euclidean`,
  `Haversine`, or rhumb-line `Rhumb` metric), linear referencing
  (`LineInterpolatePoint` / `LineLocatePoint`), affine transforms
  (`AffineTransform` with `Translate` / `Rotate` / `Scale` / `Skew`), `Orient`,
  `Validation` (OGC validity with typed reasons), `RemoveRepeatedPoints`,
  `InteriorPoint`, `BooleanOps` plus cascaded `unary_union`, `relate` (DE-9IM),
  the `StrTree` spatial index, `PreparedGeometry`, spatial joins, and set/curve
  distances (`hausdorff_distance` / `frechet_distance`).
- **`space::s2`** — spherical S2 primitives ported from golang/geo: cell ids,
  normalized cell-set unions (`CellUnion` with union / intersection /
  difference), caps, lat/lng rectangles (`Rect`, antimeridian-aware),
  polylines, edge crossings, orientation predicates, and measures.
  **`RegionCoverer` is experimental** (see *Stability*).
- **`space::r3` / `space::s1`** — supporting unit-less vector / angle types.
- **`earth`** — the WGS84 lat/lon layer: `GeoPoint`, haversine and geodesic
  distances, rhumb-line (loxodrome) distance / bearing / destination /
  interpolation (`earth::rhumb`), bounding boxes, `RadialGrid` (a
  concentric-ring sampler around a point), and the Hausdorff / Fréchet
  curve-similarity measures (`earth::similarity`).
- **`algebra`** — metric / vector-space traits (`MetricSpace`, `Norm`,
  `InnerProduct`, `VectorSpace`) and the generic Hausdorff / Fréchet kernels.
- **`format`** — native WKT, WKB, and GeoJSON codecs over `space::r2::Geometry`,
  plus a read-only ESRI Shapefile (`.shp`) geometry reader.
- **`convert`** — moves geometry between the planar (`r2`) and spherical
  (`earth`) layers.

## Stability

Pre-1.0; APIs may change. As a working classification (subject to revision —
treat this as guidance, not a guarantee):

- **Stable** — exercised by conformance suites and used downstream:
  `index::h3`, `index::geohash`, the `space::r2` geometry / algorithm surface,
  `format` (wkt / wkb / geojson, shapefile read), `earth` (points, geodesics,
  bbox, `RadialGrid`, `similarity`), and `algebra`. The `space::r2` trait names
  intentionally mirror `geo` as **migration shims** — they exist so callers can
  swap the import path, and may diverge from `geo` over time.
- **Experimental** — `space::s2::RegionCoverer`. Fixed-level covers, adaptive
  multi-level covers (normalized `CellUnion` output under a soft `max_cells`
  budget), and interior coverings exist, but the cell-intersection test is
  conservative — it counts boundary touches and either-way nesting, so it
  over-covers rather than dropping boundary-coincident cells — and it is not
  tuned for minimal covers. Good for bucketing/partitioning; not yet a
  minimal-cover engine.

### Semantics worth knowing

- **Empty-operand distances return `None`, not `0.0`.** `hausdorff_distance`,
  `frechet_distance`, and `directed_hausdorff_distance` in `space::r2`, and the
  `*_m` / typed forms in `earth::similarity`, return `Option` so an empty operand
  is not silently read as "identical". The generic `algebra` kernels remain total
  (`0.0` for empty slices).
- **`POINT EMPTY` is lossy in type.** `Point` wraps a single coordinate and has
  no empty form, so `POINT EMPTY` parses to an empty `MultiPoint` and re-emits as
  `MULTIPOINT EMPTY`. Every other `… EMPTY` round-trips with its own type
  preserved.

## Attribution

hausdorff ports code and algorithms from georust `geo`, Uber's H3, golang/geo
(S2), and others. See [`NOTICE`](NOTICE).
