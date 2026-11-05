# hausdorff-traits

The abstract geometry trait surface, and the algorithms written once over it.
Part of [hausdorff](../..).

```toml
[dependencies]
hausdorff-traits = { git = "ssh://git@github.com/urban-fabric/hausdorff.git", branch = "main" }
```

## The problem

The same geometric question gets asked in five different spaces. "Is this point
inside that triangle" is a question about the plane, about the unit sphere, about
ℝ³, and about an ellipsoid — and the answer is the *same argument* every time: a
sign test against each edge, with "sign" meaning the planar cross product in one
space and the spherical triple product in another.

Written the obvious way, that argument gets implemented once per space. The
copies then drift: a boundary case is handled in the planar version and not the
spherical one, and the two spaces start disagreeing about whether a point on an
edge is inside. That divergence is invisible — both implementations look right in
isolation — and it is exactly the kind of inconsistency this workspace exists to
remove.

This crate holds the argument, once, as a trait plus a generic algorithm. A space
supplies its own orientation predicate; the containment logic above it is shared.

## Why it sits at the bottom

The traits are **upstream of every concrete space**, including `r1`, `s1`, and
`r3` in `hausdorff-base`. That placement is forced rather than stylistic: those
substrates need to implement the abstractions, and the abstractions used to live
in `hausdorff-planar`, which is *downstream* of base. A substrate cannot
implement a trait defined in a crate that depends on it, so the trait surface has
to be its own crate below all of them.

The practical consequence is the one worth knowing: an instance lives beside the
*type*, not beside the trait. `Orientation for S2Point` is in `hausdorff-sphere`,
next to `S2Point`.

## What is here

| module | contract | generic algorithms over it |
|---|---|---|
| `orientation` | `Orientation` — the signed turn predicate | `in_triangle`, `left_of_all`, `is_collinear`, `is_convex_ring`, `contains_convex`, `convex_winding` |
| `vector` | `VectorSpace` / `InnerProduct` / `Norm` | `norm`, `norm_squared`, `normalized`, `dot` |
| `simplex` | `Simplex` / `Barycentric` | `barycentric`, `strictly_inside`, `contains` |
| `interpolate` | `Interpolate` (re-exported from `euclid-numeric`) | `sample`, `sample_at`, midpoints |
| `region` | `PointRegion` / `Prunable` — bounding volumes | `containment_masks`, `filter_contained`, `count_contained`, `any_contained`, `disjoint` |
| `interval` | `ClosedInterval` | `contains`, `intersection`, `union`, `expanded`, hull |

`interval` is the clearest payoff: `r1::Interval` on the real line and
`s2::LatInterval` on the sphere are the same algebra, and sharing it means a
latitude interval cannot develop a different notion of `union` than a plain one.

## Scalars are carrier-generic

Bounds are on `euclid-numeric`'s traits, not on `f64`, so the generic algorithms
run over tensor and expression carriers as well as host scalars. Two rules follow
from that and both are load-bearing:

- **Comparisons return masks, not `bool`.** A predicate yields `Compare::Mask`,
  consumed by `select`. You cannot branch on a graph node, so nothing here does.
- **Values are `Clone`, not `Copy`.**

The region filters return mask vectors for the same reason — a spatial prune has
to be expressible as a columnar operation, not only as a loop with an `if`.

## License

Licensed under the [Apache License, Version 2.0](LICENSE).
Copyright 2026 Raymond Mead. See [`NOTICE`](NOTICE).
