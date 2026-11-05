# Contributing to hausdorff

hausdorff is Urban Fabric's geometry and spatial-indexing substrate. It is a
**beta** repository: the public surface is still moving, and downstream projects
(`galois`, `hypatia`, `cantor`) are expected to track it rather than pin to a
stable contract.

What makes this repository different from the rest of the fleet is that it is
mostly *ported* code — georust `geo`, golang/geo's S2, Uber's H3, Karney's
GeographicLib — re-implemented so that numeric behavior and correctness are ours
to fix rather than to work around. That shapes the review bar: the question is
usually not "does this compile" but "does this still agree with the reference it
replaced". See [`NOTICE`](NOTICE) for what came from where.

## No Umbrella Crate

**Every type has exactly one public path, and consumers name the crates they
use.** There is no `hausdorff` crate re-exporting the other seven; one existed
until 0.2.0 and was removed.

This is worth stating as a rule because the pressure to re-add it is real — a
single `hausdorff` dependency reads as more convenient than four. What it
actually bought was worse than that:

- **It undid the split.** The workspace is cut so a consumer needing `GeoPoint`
  and `haversine_m` builds `hausdorff-base` and nothing else. Depending on the
  façade built all seven crates regardless of what was used, so the layering
  existed on paper and not in anyone's build.
- **It gave most types two paths.** `hausdorff::Point` and
  `hausdorff_planar::space::r2::Point` were the same type, so call sites drifted
  between them, imports in one file disagreed with the next, and "where does this
  live" stopped having one answer. That is also what made its own test file
  necessary: it existed only to prove the re-export paths still resolved, which
  is not a behavior worth a test.

**The same rule applies inside a crate, where it is easier to break by
accident.** A module is either private with its items re-exported by its parent,
or a namespace the parent re-exports nothing from. Declaring it `pub mod` *and*
re-exporting its contents gives every item in it two paths, which is the umbrella
crate's mistake at one level down. `hausdorff-planar`'s `space::r2` had drifted
into exactly that: 38 modules were public and re-exported, so 107 names had two
paths, while 26 other public items had only the long one — which is what made the
long path look load-bearing. Both halves are the same rule.

Pick the namespace form when the module has a vocabulary of its own that reads
better qualified, as `space::r2::validation` does — thirteen fault locators and
`Invalid*` enums, one of which (`GeometryIndex`) collides outright with an
unrelated type in the facade. A collision like that is a good sign the module
really is a separate vocabulary. Otherwise keep the module private and let the
parent's `pub use` be the whole surface: the file layout then stays an internal
matter, free to regroup without touching a caller.

So: do not add a crate whose `src/` is re-exports, and do not add a glob
re-export that republishes a sibling crate's surface. The one re-export that
remains deliberate is `hausdorff_earth::earth`, which republishes the
`hausdorff-base` earth core because `hausdorff-earth` genuinely *extends* that
same namespace with `arc`, `crs`, `projection`, `sphere`, `batch`, and
`polygon` — a consumer of the Earth layer would otherwise have to import half of
one module from each of two crates. If you find yourself adding a second case
like that, it probably means two crates want to be one.

## Version Control: Jujutsu

**This repository uses [Jujutsu](https://github.com/jj-vcs/jj) (`jj`), not git
directly.** It is a *colocated* repo — `.jj` and `.git` sit side by side — so git
tooling can still read it, but `jj` is the interface you should work through.

One consequence to know up front: in a colocated repo git's `HEAD` is normally
**detached**, because jj owns the working copy. That is expected, not a problem to
fix. Running `git commit` here will create a commit git can see but that sits
outside jj's view of history until the next `jj` invocation imports it, so it is
worth avoiding.

The everyday loop:

```sh
jj status                   # what has changed (jj auto-tracks new files)
jj diff                     # review it
jj commit -m "message"      # commit the working copy and start a fresh change
jj log                      # history
```

Useful beyond that:

```sh
jj describe -m "message"    # set/replace the current change's message
jj new                      # start a new empty change on top
jj edit <rev>               # move the working copy to an existing change
jj squash                   # fold the working copy into its parent
jj undo                     # undo the last jj operation
jj op log                   # every operation, so nothing is truly lost
```

## Publishing: A Branch and a Squash Merge, Never a Direct Push

**Do not push to `main`.** Push the work to a new branch, open a pull request, and
squash-merge that into `main`.

The reason is the section below: there is no CI here, so the pull request is the
only review gate the project has. A direct push to `main` does not skip a
formality — it skips the only point at which a change is read by someone other
than its author. And because the merge is a *squash*, the PR title and body are
not commentary on the commit, they **become** it: whatever you write there is what
lands on `main` and is all anyone reading `jj log` later will have.

Bookmarks (jj's branches) **do not advance on their own**, which works in your
favour here — after `jj commit`, nothing has moved and nothing has left your
machine. Publishing is a separate, deliberate step:

```sh
jj git push --named my-feature=@-   # create the bookmark and push it in one go
```

`--named` creates the remote branch and starts tracking it. Push follow-up commits
on the same branch with `jj bookmark set my-feature -r @-` and then
`jj git push --bookmark my-feature`. If you do not care what the branch is called,
`jj git push -c @-` will name one after the change id.

Then open the PR and squash-merge it on GitHub.

**After the squash merge, reconcile deliberately.** A squash merge puts a *new,
different* commit on `main` — it does not fast-forward to yours, so your local
commits are now redundant copies of work that has already landed under another
hash. Do not merge or rebase them forward:

```sh
jj git fetch                 # pick up the squashed commit on main
jj new main@origin           # start the next change on top of it
jj abandon <your-local-commits>   # drop the now-redundant originals
```

Skipping that is how the same change ends up in the history twice.

Pushing and merging both reach outside your machine, so treat them as publish
actions rather than part of the commit loop — and run the checks below first.

## Checks Run Locally, Not in CI

**There is no CI, by choice.** We do not want to pay for external compute to run
a suite that takes a couple of seconds on a laptop. The consequence is that
**the checks are your responsibility before you push** — nothing downstream will
catch a regression for you.

Run all five, from the workspace root:

```sh
cargo fmt --all --check
cargo test --all-targets --workspace
cargo clippy --all-targets --workspace -- -D warnings
cargo test --doc --workspace
RUSTDOCFLAGS="-D warnings" cargo doc --workspace --no-deps --document-private-items
```

All five must be clean.

The last one is not redundant with `cargo test --doc`: that runs doctests, while
this resolves **intra-doc links**. Nothing else catches a `[`foo`](path)` whose
path has gone stale, and a broken link is worse than a missing one — it silently
points somewhere plausible. This matters more here than in most repositories
because the crate split moved whole modules between crates: a link to
`earth::sphere` written when it lived in `hausdorff-base` now has to be prose,
since `hausdorff-base` sits *upstream* of the crate that owns it and cannot name
it. `--document-private-items` additionally catches public docs that link into
private items, which break for an outside reader.

One rustdoc trap specific to this layout: do **not** put a `///` description on a
`pub mod` whose file already has a `//!` header. rustdoc concatenates the two and
then resolves the whole block in the *parent* module's scope, so every link that
referred to an item inside the submodule silently stops resolving. Keep module
docs in the module's own file.

**A git `pre-push` hook will not help you here.** `jj` does not run git hooks, so
a `.git/hooks/pre-push` script is silently skipped on `jj git push` — exactly when
you would want it. This repository used to carry a `.githooks/pre-push` gate for
that reason; it has been removed, because a gate that does not run on the command
people actually use is worse than no gate, and it cannot be fixed from inside the
repository. Use a `jj` alias instead, in your own `~/.config/jj/config.toml`:

```toml
[aliases]
check = ["util", "exec", "--", "sh", "-c", """
  set -e
  cargo fmt --all --check
  cargo test --all-targets --workspace
  cargo clippy --all-targets --workspace -- -D warnings
  cargo test --doc --workspace
  RUSTDOCFLAGS="-D warnings" cargo doc --workspace --no-deps --document-private-items
"""]
```

Then `jj check` before `jj git push`. That lives in your user config rather than
the repo, so it stays a personal habit rather than a project mandate — and it
costs nothing but your own CPU.

Lints are declared once in `[workspace.lints]` at the root and inherited by every
crate, so `clippy -- -D warnings` is the same gate everywhere. `unsafe_code` is
denied outright.

Two lints euclid denies are **not** denied here yet: `unwrap_used` and
`expect_used`. There are a few hundred sites in `crates/*/src`, most of them
load-bearing in ported code where the upstream relied on the same invariant, so
turning them on is its own piece of work rather than a flag flip. Do not add new
ones in code you write; prefer returning an error or restructuring so the
invariant is visible.

## Tests Live Beside the Code They Test

Tests are `#[cfg(test)]` modules under each crate's `src/tests/`, not integration
tests in a `tests/` directory. The crate that owns a behavior owns its suite:

| crate | its `src/tests/` covers |
|---|---|
| `hausdorff-base` | the WGS84 earth core and `space::r1` |
| `hausdorff-planar` | `space::r2`, boolean ops, DE-9IM relate, the `Geodesic`/`Rhumb` metric traits, the codecs, `algebra`, `config` |
| `hausdorff-sphere` | `space::s2` cells, regions, and the robust spherical measures |
| `hausdorff-index` | `index::h3` and `index::geohash` |
| `hausdorff-earth` | projections, the CRS layer, `geo` conformance, batches, `convert`, similarity, `Sphere` |

Each crate's `src/tests/mod.rs` header is the authoritative index of its own
suites — what each one answers and the house rules it follows. Read it before
adding a suite, and update it when you do. It lives there rather than here
because a list in two places is a list that goes stale in one of them.

This layout is deliberate, and it is the opposite of where the suite started.
The reason is that these crates are a *stack*: `hausdorff-planar` has no idea
`hausdorff-earth` exists, so a failure in planar geometry should be reproducible
without compiling the Earth layer at all. The whole suite used to be integration
tests under a `hausdorff` umbrella crate, which made every failure a
whole-workspace failure and gave the bottom of the stack no suite of its own.
Tests beside the code also reach `pub(crate)` internals, which is where the
ported numerics actually need pinning.

There is now **no integration-test directory anywhere in the workspace**, and no
crate should grow one. If you find yourself wanting one, the thing you are
actually testing is a re-export, and a re-export is not a behavior — it is a
path. `cargo doc` and the compiler already check paths.

Where a suite lands is decided by what it *needs*, not by what it is about, and
the two can differ. The ellipsoidal-geodesic suite lives in `hausdorff-planar`
rather than in `hausdorff-base` where the Karney kernel lives, because most of
its assertions drive the `r2` `Geodesic` metric trait — and `hausdorff-base` sits
below the crate that defines it. Conversely `sphere_is_planet_flexible` sits in
`hausdorff-earth` and not with its `s2` siblings, because stating it needs a
planet radius, and `hausdorff-sphere` is deliberately radius-free. When a test
will not compile in the crate that "should" own it, that is usually the stack
telling you something true about the layering.

The two `examples/` are the only targets outside `src/`:
`hausdorff-planar/examples/bool_ops_corpus.rs` and
`hausdorff-earth/examples/profile_kernels.rs`. Each sits in the lowest crate that
can compile it, same rule as the suites. They are built by
`cargo test --all-targets`, so they are kept compiling even though nothing
asserts on them.

**A test belongs beside the implementation in its own file** (`foo/tests.rs` next
to `foo.rs`) when it needs private constants or helpers the `src/tests/` tree
cannot see. Anything reachable through the public surface belongs in `src/tests/`,
where it keeps passing when the algorithm behind it is replaced.

**No file carrying tests goes over 400 lines.** Split by sub-topic into a
directory whose `mod.rs` header explains the split, as `src/tests/r2/` and
`earth/crs/tests/serialization/` do. Pure implementation surface is exempt,
having no sub-topics to split along — `index::grid` and `space::s2::loops` are
both over 400 lines of implementation with a small test module at the end, and
that is fine. The rule is about the tests.

## What a Change Should Come With

- **A ported algorithm** — say which upstream revision or reference you ported
  from, and keep the differential or golden test that shows you agree with it.
  Where the port deliberately diverges (dropped capability masks in the Karney
  kernel, ignored Z/M ordinates in the WKB reader), say so in the module header.
  A silent divergence from the reference is the failure mode that matters most
  here, because the whole point of the port is that callers can trust it behaves
  like the thing it replaced.
- **A geometric predicate** — cover the degenerate cases, not just the general
  one: collinear points, zero-length segments, repeated vertices, a point exactly
  on a boundary, rings that touch at a single vertex, and the antimeridian and
  poles for anything geographic. These are where the ported implementations earn
  their keep, and where a plausible-looking rewrite breaks.
- **Anything the contract leaves open** — boundary behavior especially. "A point
  exactly on an edge is not guaranteed either way" is a fine contract; leaving it
  unsaid is not, because two call sites will assume opposite answers.
- **Docs are part of the deliverable.** Document panics, accepted error bounds,
  and the units a number is in — `_m` for metres, radians versus degrees, and
  steradians versus m². A unit confusion that type-checks is this repository's
  most expensive class of bug.

## Numerical Work

Numeric changes carry risk that ordinary review does not catch.

Prefer **identity and round-trip tests** over hand-copied reference digits where
you can: an inverse/direct geodesic round trip, a projection's forward/inverse
pair, `parent(child(c)) == c` on a cell hierarchy, or a polygon area computed two
independent ways validates itself and cannot be transcribed wrongly. Where a
reference value is needed, derive it from an exact closed form or an independent
implementation, and say in a comment where it came from.

Never "fix" a test by loosening a tolerance without understanding why the error
grew. Tolerances here are set from measured error against a reference; an
unexplained loosening hides the regression it was meant to catch.

Exact predicates are exact on purpose. The orientation and incircle tests in
`space::r2::exact` use adaptive error bounds so that a sign is never wrong, and
the boolean-ops sweep depends on that. Do not replace one with a tolerance
comparison because it looks simpler.

## Commit and PR Messages

State what changed and why it had to change. The subject is imperative and
describes the outcome, not the activity — "Make Country the organization-wide
country definition", not "update country enum". If a change fixes something
subtle, say what the failure was: a reader six months later needs the reason,
and the diff already shows the mechanics.

Avoid narrating the work. History records that a decision was made; the message
should record the decision.

## No CHANGELOG

**There is no CHANGELOG, by choice.** While hausdorff is beta, the version number
and the `jj log` are the record. Maintaining a hand-written changelog for a
surface that is still being reshaped costs more than it informs, and a stale
changelog is worse than none.

Instead:

- Put the explanation in the **pull request body**, which the squash merge turns
  into the commit message on `main`. That is where a changed predicate
  convention or a moved module gets described. Write it as the permanent record
  it becomes, not as a note to the reviewer.
- Bump the workspace version in the root `Cargo.toml` when the public surface
  moves. All seven crates share `[workspace.package]`, so they version together.
- Say so in the PR title when module paths or semantics change, since `galois`,
  `hypatia`, and `cantor` read titles to decide whether they need to react — and
  since under a squash merge the title is the one line that survives into
  `jj log`.

## Documentation

Repository documentation describes **implemented behavior, accepted limits, and
operational contracts**. Forward-looking work lives in `plan/<repo>/`, and change
history lives in version control rather than in prose — a document says what the
system is, not what it once was or when it changed.

Shared vocabulary is defined once in the fleet glossary. A repository document
may add genuinely local terms; it defers to the glossary for every shared noun,
and must not redefine one.

## License

By contributing you agree your contributions are licensed under the
[Apache License, Version 2.0](LICENSE), per section 5 of that license.

If you port code in, it must be under a licence compatible with Apache-2.0, and
you must attribute it in the same change — naming what was taken, where it now
lives, and the licence it arrived under. The same applies to embedded data: the
NTv2 grid files stay under their issuing agencies' terms, and anything similar
you add needs its own entry.

**There are eight `NOTICE` files and you need two of them.** Each crate carries
its own, scoped to what that crate actually contains, and the root one is the
union. Add the entry to *both* the owning crate's `NOTICE` and the root. The
per-crate copy is the one that matters legally: a consumer who takes
`hausdorff-planar` alone receives the georust attribution and not the H3 one,
which is only true if the attribution lives next to the code it describes. The
root copy is for someone reading the repository.

Every crate also carries its own `LICENSE` and `README.md`, which is what makes
each one self-contained in a published tarball. Keep documentation local: a
crate's README explains the problem *that crate* solves, and the root README is
an overview that points into them. Detail added to the root instead of to a crate
is detail that drifts out of sync with the code it describes.
