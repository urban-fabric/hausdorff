# hausdorff-kernels

Backend-neutral geospatial formula kernels, written once over any numeric
carrier. Part of [hausdorff](../..).

```toml
[dependencies]
hausdorff-kernels = { git = "ssh://git@github.com/urban-fabric/hausdorff.git", branch = "main" }
```

## The problem

A haversine distance is needed in four places at once: on host `f64` at an API
edge, inside a dataframe expression pushed down to a query engine, elementwise
over a column of a million points, and staged into a graph for deferred
execution. Written the ordinary way that is four implementations, and they drift
— a coefficient is corrected in one, a `log(0)` guard exists in two of four, a
pole clamp is right in the batch path and wrong at the edge.

The drift is worse here than in most domains because the failures are *plausible*.
A distance that is wrong by a kilometre still looks like a distance; nothing
crashes, and the test that would have caught it exists only in the path that was
already correct.

This crate makes the formula the single artifact and the carrier a parameter. One
`haversine_on_m`, bound on `euclid-numeric`'s `Real`, runs on `f64`, on `f32`, and
on whatever a consumer defines.

## What "backend-neutral" costs

Three rules, all inherited from `euclid-numeric`, and all of them exist so that a
*deferred* carrier stays possible — you cannot branch on a graph node, `Copy` a
symbolic handle, or write a float literal into an expression tree:

- comparisons return masks, never a host `bool`
- values are `Clone`, never `Copy`
- literals are lifted through `Real::scalar`, never written as `2.0`

The visible consequence is that **iterations run a fixed count rather than
breaking on a tolerance**. `tauf` does exactly five Newton steps and
`lat_from_authalic_q` fifteen fixed-point steps, because "stop when the residual
is small" needs a scalar `bool` that a tensor backend does not have. Both counts
reach full `f64` precision well inside the bound; the comment at each site says
so.

The same rule is why singularities are removed algebraically instead of branched
around. `authalic_q` used to carry Snyder's literal `1/(2e)`, which is a pole at
zero eccentricity — so every equal-area quantity on a spherical datum came back
`NaN`. It is now written through `atanh`, which has no pole, and the one genuine
branch left uses mask `select` with the divisor sanitized *before* the division so
an elementwise carrier never forms the `0/0` at all.

## What is here

| module | kernels |
|---|---|
| `earth` | haversine and its radius-parameterized forms, forward/direct destination, initial bearing, the `Sphere<T>` carrier, lat/lon validity masks |
| `projection` | transverse Mercator, web Mercator, Lambert conformal conic, Albers equal-area, polar stereographic, plus the shared auxiliary latitudes (conformal, authalic, rectifying) |
| `planar` | orientation and turn predicates, segment and point distances, rectangle overlap, curve helpers |
| `interval` | closed-interval algebra on the real line |
| `spherical1` | angular normalization and sector algebra on the circle |

What is deliberately **absent** is as much the point as what is present: no eager
data records, no storage encodings, no spatial indexes, no geometry types. A
kernel takes scalars and returns scalars. `GeoPoint`-typed, unit-typed, and
indexed wrappers over these formulas live in `hausdorff-base` and
`hausdorff-earth`, which is what keeps this crate lowerable to a backend that has
no structs.

## License

Licensed under the [Apache License, Version 2.0](LICENSE).
Copyright 2026 Raymond Mead. See [`NOTICE`](NOTICE) — the projection formulas are
Snyder's published series, implemented rather than ported.
