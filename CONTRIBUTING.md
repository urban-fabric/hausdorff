# Contributing to hausdorff

Geometry and spatial-index substrate. Keep changes focused and easy to review.

## Merging

**Squash and merge.** A pull request lands on `main` as one commit whose message
is the PR title and body, so the branch's intermediate commits never reach
`main`. Two consequences worth knowing:

- Write the PR title and body as the commit you want in the history. Local
  commit messages on the branch are working notes.
- After a merge, your local branch commits are *not* ancestors of `origin/main`.
  `git pull --ff-only` will refuse and leave you on a stale `main`; use
  `git checkout -B main origin/main` once the branch is merged. Check for
  uncommitted work first — a squash discards nothing on the remote, but a local
  reset discards whatever is in your working tree.

Branch names are `<kind>/<subject>`: `feat/`, `fix/`, `docs/`, `chore/`,
`refactor/`.

## Commit and PR messages

State what changed and why it had to change. The subject is imperative and
describes the outcome, not the activity — "Make Country the organization-wide
country definition", not "update country enum". If a change fixes something
subtle, say what the failure was: a reader six months later needs the reason,
and the diff already shows the mechanics.

Avoid narrating the work. History records that a decision was made; the message
should record the decision.

## Documentation

Repository documentation describes **implemented behavior, accepted limits, and
operational contracts**. Forward-looking work lives in `plan/<repo>/`, and change
history lives in Git rather than in prose — a document says what the system is,
not what it once was or when it changed.

Shared vocabulary is defined once in the fleet glossary. A repository document
may add genuinely local terms; it defers to the glossary for every shared noun,
and must not redefine one.

## Verification

Run these before opening a pull request:

```bash
cargo build
cargo test --workspace
```

