# hausdorff

Internal geospatial primitives and index implementations for UrbanFabric — an
in-house, dependency-light replacement for the `geo` / `h3o` / `s2` / `geohash`
stack, so UrbanFabric projects share one geometry core. Proprietary software;
see [`NOTICE`](NOTICE) for the open-source works it ports from.

The H3 core was the first piece, moved out of Kolmogorov so other projects can
use the same internal implementation.

## Workspace layout

A cargo workspace (flat layout, every crate under `crates/`). The `hausdorff`
crate is a façade that re-exports the subcrates under the original single-crate
module surface, so consumers depending on package `hausdorff` are unaffected by
the split. Consumers that only need a slice can depend on a subcrate directly
(same git URL, different `package` name):

| crate | contents | depends on |
|---|---|---|
| `hausdorff-base` | `earth` core (`GeoPoint`, `BBox`, haversine/rhumb/Karney geodesy, `ecef`, `ellipsoid`, `RadialGrid`), `space::s1`, `space::r3` | — |
| `hausdorff-planar` | `space::r2`, `algebra`, `config` (+ the `profiling` feature hooks) | base |
| `hausdorff-sphere` | `space::s2` | base, planar |
| `hausdorff-index` | `index::h3`, `index::geohash`, `index::grid` | base, planar |
| `hausdorff-format` | `format` (WKB/WKT/GeoJSON/shapefile) | planar |
| `hausdorff-earth` | `earth::{arc, projection, sphere, similarity, batch}`, `convert`; re-exports the base earth core as the complete `earth` module | base, planar, sphere, index |
| `hausdorff` | the façade: re-exports everything under the original paths; hosts the integration test suite, benches, examples | all of the above |

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
  `Haversine`, ellipsoidal `Geodesic`, or rhumb-line `Rhumb` metric),
  `GeodesicArea` (exact ellipsoidal polygon area/perimeter) beside the
  spherical `ChamberlainDuquetteArea`, linear referencing
  (`LineInterpolatePoint` / `LineLocatePoint`), metric interpolation for
  Euclidean, Haversine, and rhumb-line paths, affine transforms
  (`AffineTransform` with `Translate` / `Rotate` / `Scale` / `Skew`),
  `MapCoords`, `Orient`, `Validation` (OGC validity with typed reasons),
  `RemoveRepeatedPoints`, `InteriorPoint`, `BooleanOps` plus cascaded
  `unary_union`, `relate` (DE-9IM), the `StrTree` spatial index,
  `PreparedGeometry`, spatial joins, exact fixed-grid predicates
  (`exact::orient`, `exact::segment_intersection`,
  `exact::point_in_triangle`, `exact::in_circle_position`), and set/curve
  distances (`hausdorff_distance` / `frechet_distance`).
- **`space::s2`** — spherical S2 primitives ported from golang/geo: cell ids,
  normalized cell-set unions (`CellUnion` with union / intersection /
  difference), caps, lat/lng rectangles (`Rect`, antimeridian-aware),
  polylines, edge crossings, orientation predicates, measures, and
  great-circle arc kernels (`s2::arc`: slerp, cross/along-track, latitude
  extremes, arc rect bounds). `Cap`, `Polyline`, `Loop`, `Polygon`, and
  `Region` expose conservative lat/lng rectangle bounds. **`RegionCoverer` is
  experimental** (see *Stability*).
- **`space::r3` / `space::s1`** — supporting unit-less vector / angle types.
- **`earth`** — the WGS84 lat/lon layer: `GeoPoint`, spherical (haversine) and
  ellipsoidal geodesics (Karney direct/inverse/area on any `Ellipsoid` — WGS84,
  GRS80, Mars, …; differential-tested against `geographiclib-rs`), geodesic
  lines for cheap waypoint sampling (`Ellipsoid::line`), ellipsoidal polygon
  perimeter/area (`Ellipsoid::ring_perimeter_area_m`), rhumb-line (loxodrome)
  distance / bearing / destination / interpolation (`earth::rhumb`),
  planet-flexible sphere geodesics (`sphere::Sphere`), great-circle arc
  operations (`earth::arc`: interpolation, cross/along-track distance, arc
  intersection, spherical point-in-ring, Clairaut-aware arc bboxes),
  geocentric/local frames (`earth::ecef`: ECEF ↔ geodetic, ENU tangent
  frames), map projections (`projection`: Web Mercator, Transverse Mercator /
  UTM, with whole-geometry `project_geometry` / `inverse_geometry`), bounding
  boxes, `RadialGrid` (a
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
  `format` (wkt / wkb / geojson, shapefile read), `earth` (points, geodesics —
  spherical and ellipsoidal — ellipsoids, projections, bbox, `RadialGrid`,
  `similarity`), and `algebra`. The `space::r2` trait names
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
- **Proof-backed APIs are algebraic contracts, not a blanket IEEE proof.**
  Riemann theorem names in docs identify the exact real/integer formulas the Rust
  code implements. Floating-point kernels still have ordinary `f64` roundoff;
  fixed-grid `space::r2::exact` predicates report overflow instead of guessing.

## Attribution

hausdorff ports code and algorithms from georust `geo`, Uber's H3, golang/geo
(S2), and others. See [`NOTICE`](NOTICE).
