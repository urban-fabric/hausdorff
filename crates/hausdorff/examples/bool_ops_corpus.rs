//! Set-algebra invariants over a corpus of real difference operands.
//!
//! **The corpus is not in this repository, and is not expected to be.** It was
//! captured from kolmogorov's isochrone band builder — concave hulls of nested
//! reachable sets, pinched and sliver-ridden and self-touching in ways the fixed
//! fixtures in `tests/bool_ops.rs` are not — and it is neither redistributable
//! here nor small. This example is the harness; the operands are supplied by
//! whoever runs it.
//!
//! It is kept because the counts asserted elsewhere in the tree came from it and
//! are otherwise unreproducible: the thirty-two band pairs behind
//! [`hausdorff::space::r2::bool_ops`]'s tracing notes, the six that reached the
//! tracer with an odd-degree vertex, the three whose union was smaller than one
//! of their operands. Two of those pairs are checked in under `tests/fixtures`
//! and assert exact areas; a suite of two cannot say what a rate is, and this
//! can. Run it against a fresh capture before trusting a change to the sweep or
//! to contour tracing on real hulls.
//!
//! The input contract is deliberately dumb so a corpus can be rebuilt from any
//! source of polygon pairs, not just that one: `<dir>` holds one `.hex` file per
//! case, each two non-blank lines of WKB hex, `A` then `B`. Cases are reported
//! by file stem, so name them for what they are.
//!
//! Usage: `cargo run -p hausdorff --release --example bool_ops_corpus -- <dir>`

use hausdorff::format::wkb::from_wkb;
use hausdorff::{Area, BooleanOps, Geometry, MultiPolygon};

const EPS: f64 = 1e-6;

fn unhex(s: &str) -> Vec<u8> {
    let s = s.trim();
    (0..s.len() / 2)
        .map(|i| u8::from_str_radix(&s[2 * i..2 * i + 2], 16).expect("hex"))
        .collect()
}

fn as_multipolygon(g: Geometry) -> MultiPolygon {
    match g {
        Geometry::MultiPolygon(mp) => mp,
        Geometry::Polygon(p) => MultiPolygon::new(vec![p]),
        other => panic!("expected a polygonal geometry, got {other:?}"),
    }
}

/// Every invariant that ties the four ops to one another and to `unsigned_area`.
/// Returns `(union_violations, difference_violations, worst_identity_residual)`.
fn check(a: &MultiPolygon, b: &MultiPolygon) -> (usize, usize, f64, Vec<String>) {
    let inter = a.intersection(b).unsigned_area();
    let uni = a.union(b).unsigned_area();
    let diff = a.difference(b).unsigned_area();
    let xor = a.xor(b).unsigned_area();
    let (aa, ba) = (a.unsigned_area(), b.unsigned_area());

    let mut notes = Vec::new();
    let mut union_bad = 0;
    let mut diff_bad = 0;

    // Monotone bounds — these are the ones that catch an invented chord
    // directly, because a double-counted ring inflates a result past its bound.
    if uni < aa.max(ba) - EPS {
        union_bad += 1;
        notes.push(format!("union {uni:.9} < max(A,B) {:.9}", aa.max(ba)));
    }
    if diff > aa + EPS {
        diff_bad += 1;
        notes.push(format!("A\\B {diff:.9} > A {aa:.9}"));
    }
    if inter > aa.min(ba) + EPS {
        notes.push(format!("A∩B {inter:.9} > min(A,B) {:.9}", aa.min(ba)));
    }

    // Inclusion–exclusion and the derived identities.
    let residuals = [
        ("A∪B = A+B-A∩B", (uni - (aa + ba - inter)).abs()),
        ("A\\B = A-A∩B", (diff - (aa - inter)).abs()),
        ("A⊕B = A∪B-A∩B", (xor - (uni - inter)).abs()),
    ];
    let mut worst = 0.0f64;
    for (name, r) in residuals {
        worst = worst.max(r);
        if r >= EPS {
            notes.push(format!("{name}: residual {r:.6e}"));
        }
    }

    (union_bad, diff_bad, worst, notes)
}

fn main() {
    let Some(dir) = std::env::args().nth(1) else {
        eprintln!(
            "usage: cargo run -p hausdorff --release --example bool_ops_corpus -- <dir>\n\
             \n\
             <dir> holds one file per case, two lines of WKB hex each (A then B).\n\
             The band corpus these were captured from is not in the repository;\n\
             the two pairs worth keeping are fixtures under tests/fixtures."
        );
        std::process::exit(2);
    };

    let mut files: Vec<_> = std::fs::read_dir(&dir)
        .unwrap_or_else(|e| panic!("read {dir}: {e}"))
        .filter_map(|e| e.ok())
        .map(|e| e.path())
        .filter(|p| p.extension().is_some_and(|x| x == "hex"))
        .collect();
    files.sort();

    let (mut union_bad, mut diff_bad, mut ident_bad, mut pairs) = (0, 0, 0, 0);
    let mut worst_overall = 0.0f64;
    let mut worst_case = String::new();

    for path in &files {
        let text = std::fs::read_to_string(path).expect("read case");
        let lines: Vec<&str> = text.lines().filter(|l| !l.trim().is_empty()).collect();
        assert_eq!(lines.len(), 2, "{}: expected 2 WKB lines", path.display());

        let a = as_multipolygon(from_wkb(&unhex(lines[0])).expect("A"));
        let b = as_multipolygon(from_wkb(&unhex(lines[1])).expect("B"));
        let name = path.file_stem().unwrap().to_string_lossy().to_string();

        pairs += 1;
        let (u, d, worst, notes) = check(&a, &b);
        union_bad += u;
        diff_bad += d;
        if notes.iter().any(|n| n.contains("residual")) {
            ident_bad += 1;
        }
        if worst > worst_overall {
            worst_overall = worst;
            worst_case = name.clone();
        }
        if !notes.is_empty() {
            println!("{name}  [A={} B={}]", a.0.len(), b.0.len());
            for n in notes {
                println!("    {n}");
            }
        }
    }

    println!("\n--- summary over {pairs} pairs ---");
    println!("union violations      (A∪B < max(A,B)) : {union_bad}");
    println!("difference violations (A\\B > A)        : {diff_bad}");
    println!("identity violations   (residual ≥ {EPS:.0e}): {ident_bad}");
    println!("worst identity residual                : {worst_overall:.6e}  ({worst_case})");
}
