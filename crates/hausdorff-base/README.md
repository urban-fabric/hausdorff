# hausdorff-base

WGS84 lat/lon primitives, geodesy, and the unit-less substrates — the thin slice
most consumers need. Part of [hausdorff](../..).

```toml
[dependencies]
hausdorff-base = { git = "ssh://git@github.com/urban-fabric/hausdorff.git", branch = "main" }
```

## The problem

Most code that touches geography needs a very small amount of it: a point, a
bounding box, and the distance between two points. Pulling a full geometry stack
in to get those three things means compiling a boolean-operations sweep, a DE-9IM
engine, an H3 port, and four codecs to call `haversine_m`.

So this crate is the deliberate floor of the workspace. **Depend on it alone and
the planar, spherical, index, and Earth layers are never built.** That is the
entire reason the workspace is split the way it is, and it is why there is no
umbrella crate: a single `hausdorff` dependency used to pull all seven in
regardless of what was used, which made the layering real on paper and absent
from everyone's build.

## Two things live here, and they are different

**`earth` — WGS84, where metres are real.** `GeoPoint`, `BBox`, `Bearing`,
spherical (haversine/rhumb) geodesics, the ellipsoidal Karney solver, ECEF
conversion, and the six-latitude system of an oblate spheroid.

**`space` — the unit-less substrates.** `r1` intervals, `s1` angles, `r3`
vectors. These carry **no radius and no physical unit**; they are the
dimensionless geometry the spherical layer is built from. Angles live here rather
than in `euclid-units` because the spherical code computes in them constantly and
wants them concrete.

The split matters when reading results. An `s1::Angle` has no length. A metre
appears only when something in `earth` supplies a radius.

## Two geodesics, two conventions

| | `haversine_m` (spherical) | `Ellipsoid::distance_m` (Karney) |
|---|---|---|
| model | sphere of mean radius | oblate spheroid |
| accuracy | ~0.5% | round-off, ~15 nm |
| antipodal pairs | fine | *this is why it exists* |
| latitude outside ±90° | finite, meaningless | `NaN` |

The last row is a genuine divergence and it is documented at both call sites
rather than left to be discovered. The spherical path does not validate latitude —
the check costs more than the formula on a hot path — while the ellipsoidal path
follows GeographicLib and signals. Code that switches metrics should not infer a
convention from whichever it tried first.

The Karney solver is the interesting one. Vincenty's method, the usual choice,
*fails to converge* on nearly-antipodal pairs, and carrying a more complicated
algorithm is the price of those inputs working. `src/tests/geodesy_extremes.rs`
is built around that: antipodes, both poles, sub-micrometre separations, and the
result worth not "correcting" — antipodal points on the equator are 20 003 931 m
apart, not half the equatorial circumference, because on an oblate ellipsoid the
shortest path between them goes over the pole.

## Units are a trait, not a struct

`LengthMeters` and `AreaSquareMeters` are this crate's own carriers implementing
`euclid-units`' dimension contracts. That is the division of labour euclid
intends: euclid defines what a length *is* and that `length × length = area`;
a consumer defines how one is stored. So `radius() * radius()` type-checks to an
area and `radius() * angle.radians()` to a length, and a unit confusion that
would otherwise compile does not.

Every geodesic has a raw `_m` form and a unit-typed form. Prefer the typed one;
the `_m` kernels exist for hot paths and for callers already carrying bare `f64`.

## License

Licensed under the [Apache License, Version 2.0](LICENSE).
Copyright 2026 Raymond Mead. See [`NOTICE`](NOTICE) — the Karney kernel is ported
from GeographicLib (MIT) and the substrates mirror golang/geo.
