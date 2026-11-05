# hausdorff-planar

Planar `f64` geometry, the metric abstractions over it, and the native codecs —
the in-house replacement for `geo`/`geo-types`. Part of [hausdorff](../..).

```toml
[dependencies]
hausdorff-planar = { git = "ssh://git@github.com/urban-fabric/hausdorff.git", branch = "main" }
```

This is the largest crate in the workspace, because the plane is where most
geometry actually happens.

## The problem

Planar geometry is a solved problem with unsolved edges. The algorithms are
decades old and the textbook versions are short; what is hard is that they are
*fragile at degeneracies* — collinear points, coincident vertices, a boundary
touched exactly, two rings meeting at a single point — and every real dataset is
full of them.

Depending on an external stack for this means that when a sweep fails to
terminate on a pinched isochrone band, the options are to work around it
downstream or to wait. Owning the implementation means numeric behaviour,
tuning, and correctness are ours to fix. That is the whole trade: more code in
exchange for the defects being reachable.

The repository history is the argument. A boolean union once returned *less area
than one of its operands* on pinch geometry; a sweep once failed to converge on
near-coincident vertices, passing three million splits and 114 GB of resident
events before the kernel killed the process. Both are now fixtures.

## What is here

| module | contents |
|---|---|
| `space::r2` | the geometry model (`Coord`, `Point`, `Line`, `LineString`, `Polygon`, `Multi*`, `Rect`, `Triangle`, `Geometry`) and the algorithms over it |
| `algebra` | metric- and vector-space abstractions: Hausdorff/Fréchet distances, generic over the metric |
| `format` | native WKB, WKT, GeoJSON, and shapefile codecs |
| `config` | machine-tunable kernel parameters |

Within `space::r2`, the pieces worth knowing about:

- **`bool_ops`** — a Martinez–Rueda sweep for union, intersection, difference,
  and xor.
- **`relate`** — the DE-9IM matrix, with the full matrix exposed so a caller can
  ask for a pattern this type does not name.
- **`exact`** — orientation and in-circle predicates over checked `i128`.
- **`make_valid`** — odd-even polygon repair by arrangement, not by ray casting.
- **`prepared`** / **`index`** — an STR-tree and a prepared-geometry fast path.

## Exact predicates, and why they are not optional

A floating-point determinant can return the **wrong sign** — not an imprecise
magnitude, a wrong sign — when three points are nearly collinear. That is
categorically worse than an imprecise number, because the sweep and the hull walk
*branch* on it: one bad sign and the algorithm takes a structurally impossible
path and either produces nonsense or loops forever.

`exact` answers in `i128` on an integer grid, with overflow **reported** rather
than silently wrapped — a wrapped sign would be indistinguishable from a real
answer, which is the worst failure available to a predicate whose entire value is
sign correctness. The two predicates have different ceilings, because `orient` is
degree 2 and `in_circle` degree 4: around 2⁴² a coordinate is comfortable for one
and already past the other. That asymmetry is pinned by a test rather than left
to be discovered.

Do not replace an exact predicate with a tolerance comparison because it looks
simpler.

## The codecs read bytes they did not write

`format` is the only place in the workspace that parses untrusted input, which
makes "what does this do on garbage" a correctness question. The bar is that
every input yields `Ok` or `Err` — never a panic, and in particular never an
*abort*, which a caller cannot catch. Two failure modes that did abort:

- nine bytes of WKB claiming `u32::MAX` members asked the allocator for tens of
  gigabytes, and
- nested collections recursed without a depth limit, so ~180 KB of WKB and ~1 MB
  of WKT each exhausted the stack.

Counts are now checked against the remaining bytes before anything is allocated,
so allocation stays proportional to input, and both parsers carry a depth limit.
`src/tests/codec_hostile.rs` is what keeps them bounded — and it fails by
*crashing* if a guard is removed, which is worth knowing before debugging a
mysteriously truncated test run.

## License

Licensed under the [Apache License, Version 2.0](LICENSE).
Copyright 2026 Raymond Mead. See [`NOTICE`](NOTICE) — `space::r2` is a port of
georust `geo` (MIT OR Apache-2.0).
