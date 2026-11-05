# hausdorff-sphere

Unit-sphere geometry and the S2 cell hierarchy — the golang/geo S2 port. Part of
[hausdorff](../..).

```toml
[dependencies]
hausdorff-sphere = { git = "ssh://git@github.com/urban-fabric/hausdorff.git", branch = "main" }
```

## The problem

Planar geometry on lat/lon is wrong, and it is wrong in ways that pass tests.
Treating degrees as a Cartesian plane means the antimeridian is a wall, the poles
are a line, and a polygon that spans either is silently mangled. Project instead
and the distortion moves rather than disappearing: every projection is wrong about
*something*, and which thing depends on the projection.

The fix is to stop projecting. Represent a point on the sphere as a **unit vector
in ℝ³** and the pathologies go away, because a unit vector has no seam and no
singular point. A great-circle arc is then a plane through the origin, "which side
of this edge" is the sign of a triple product, and a spherical polygon's area is a
sum of triangle excesses. No projection, no antimeridian special case, no pole
special case.

That is what this crate is: the geometry you want when the answer has to be right
near ±180° or at a pole.

## Everything here is radius-free

Areas are **steradians** and lengths are **angles**. Nothing in this crate knows
a planet radius, which is deliberate — the geometry is the same on Earth, Mars,
and the unit sphere, and the radius is a separate concern. Multiply by one to get
metres, via `Sphere` in `hausdorff-earth`.

The practical test of that boundary: a claim that needs a radius to even *state*
belongs in `hausdorff-earth`, not here. There is exactly one such test
(`sphere_is_planet_flexible`) and it lives over there for that reason.

## What is here

| area | contents |
|---|---|
| cells | `S2CellId` — the 64-bit codec, validity, hex tokens, the cube-to-sphere quadratic projection |
| hierarchy | `parent` / `children` / child positions, `range_min`/`range_max` containment, Hilbert-order traversal |
| neighbours | edge, vertex, and all-neighbour queries |
| regions | `Cap`, `Loop`, `Polygon`, `Rect`, `CellUnion`, and region coverings |
| predicates | robust orientation, edge crossing, point-on-arc |
| measures | triangle area (Girard and van Oosterom), angles, perimeters |
| relate | all ten OGC spatial predicates, for points and for areal operands with holes |
| bool_ops | union, intersection, difference, symmetric difference on spherical polygons |

## The cell id is a space-filling curve, and that is the point

An `S2CellId` is not an opaque handle. The hierarchy is encoded in the bits, so
`parent` is a mask, containment is an integer range comparison, and cells sorted
numerically are in **Hilbert order** — spatially adjacent cells are usually
numerically adjacent.

That is what makes S2 a partitioning tool rather than just a geometry library: a
range scan over sorted cell ids is a spatial query, so a database that can do
range scans can do spatial indexing with no spatial index. `galois` uses it for
exactly that.

## Boolean operations, natively on the sphere

`bool_ops` gives union, intersection, difference, and symmetric difference on
spherical polygons with holes and multiple components — the operations a planar
engine cannot do honestly near the seams.

It is **not a sweep**. A Martinez–Rueda sweep needs a total order to advance
along, and a sphere has none: sweeping by longitude reintroduces exactly the
antimeridian and pole degeneracies the crate exists to avoid. So the architecture
is an *arrangement* — node every arc against the other operand into one shared
vertex pool, classify each resulting arc by a single point test, trace the
survivors into rings, then nest them by orientation and containment. The same
shape `hausdorff-planar`'s `make_valid` uses, and for the same reason.

The oracles are closed-form lune areas, set algebra in steradians, agreement with
the planar engine on a small flat patch, and rotation invariance in longitude —
the last being the one that shows the module was worth writing, since a planar
reading is not rotation-invariant across ±180°.

## Scope, honestly stated

`relate` now answers all ten OGC predicates for point-against-areal and for
areal-against-areal including holes. What still returns `None` — the signal to
fall back to the planar predicate — is any **linear** operand, and the
dimension-dependent predicates over *mixed* dimensions. Returning `None` rather
than guessing is what makes wiring it in safe: it cannot regress a case it does
not cover.

Two limits worth knowing about the boolean engine. Noding is `O(n·m)` with no
spatial index yet, so two large coverings are quadratic — the cell machinery here
is the obvious prune, and that is a performance change rather than a correctness
one. And input loops must be simple, which is `Loop`'s existing precondition;
self-intersecting input gets a best effort rather than a panic.

## License

Licensed under the [Apache License, Version 2.0](LICENSE).
Copyright 2026 Raymond Mead. See [`NOTICE`](NOTICE) — a port of golang/geo's S2
(Apache-2.0, Google Inc.).
