# hausdorff-index

Discrete spatial indexes: the H3 hexagonal hierarchy, geohash, and a uniform
lon/lat grid. Part of [hausdorff](../..).

```toml
[dependencies]
hausdorff-index = { git = "ssh://git@github.com/urban-fabric/hausdorff.git", branch = "main" }
```

## The problem

Continuous geometry answers "where exactly is this". Indexing answers a different
question: "give me a **key** I can group, join, and range-scan on". Those are not
the same problem, and the second one is mostly about databases.

A discrete index turns a position into an integer or a short string, so that
points near each other usually get keys near each other. Once that holds you can
do spatial work with ordinary data machinery — `GROUP BY` a cell, join two tables
on a cell, shard by a key prefix — with no spatial index and no geometry engine in
the query path. That is why every one of these exists, and it is why they live in
a crate separate from the geometry.

| | shape | key | good at |
|---|---|---|---|
| **H3** | hexagons (+12 pentagons) | `u64` | aggregation and neighbourhoods; uniform adjacency |
| **geohash** | lat/lon rectangles | base-32 string | prefix queries; human-readable; works as a plain string key |
| **grid** | uniform lon/lat cells | cell index | radius and k-nearest joins where simplicity beats uniformity |

The S2 cell index is *not* here — it lives with its geometry in
`hausdorff-sphere`, because an `S2CellId` is as much a region as a key.

## Why hexagons

Hexagons have one property squares do not: **every neighbour shares an edge, and
every neighbour's centre is the same distance away**. A square grid has eight
neighbours at two different distances, so "within one step" is ambiguous and any
distance-weighted aggregation over a square grid is subtly anisotropic. For
movement, flow, and coverage analysis — the things H3 was built for — that
matters more than the convenience of squares.

## The twelve pentagons

You cannot tile a sphere with hexagons. H3 projects an icosahedron, an
icosahedron has twelve vertices, and a hexagonal grid cannot close around a vertex
of degree five — so at each of those twelve places one cell is a **pentagon** and
the grid is locally wrong by construction.

Every traversal, hierarchy, and area operation has to special-case them, which
makes pentagons by far the most likely place for an H3 port to be subtly wrong.
`src/tests/h3/pentagons.rs` is built around the counting identities that catch it:

| shape | ring at `k` | disk out to `k` | descendants at resolution `r` |
|---|---|---|---|
| hexagon | `6k` | `1 + 3k(k+1)` | `7^r` |
| pentagon | `5k` | `1 + 5k(k+1)/2` | `1 + 5(7^r − 1)/6` |

There is a trap in using them, and the suite states which regime each case is in.
**The ideal hexagon counts only hold where the disk contains no pentagon.** At
resolution 0 there are 122 cells and 12 are pentagons, so almost every base cell
has one within two steps and its disk is legitimately *short* of the ideal — 18
instead of 19 at `k = 2`. A test asserting the ideal count on a base cell would be
asserting that pentagons do not exist, and "fixing" the implementation to satisfy
it would break the grid.

## The H3 suite has no runtime oracle

The port was developed as a differential test against the `h3o` crate. That
dependency is **gone**: every input is deterministic (fixed lat/lon grids, the
native base-cell set, a seeded LCG), so the oracle's answers are frozen into
FNV-1a checksums and exact value lists. The differential coverage is intact with
no reference crate in the build graph.

The consequence for editing it: a golden constant is *evidence*, not a
preference. If one starts failing, the implementation changed — find out why
before touching the number.

## License

Licensed under the [Apache License, Version 2.0](LICENSE).
Copyright 2026 Raymond Mead. See [`NOTICE`](NOTICE) — `index::h3` is ported from
Uber's H3 (Apache-2.0) and its constant tables are lifted verbatim.
