//! The matrix's own reads: rendering, the pattern alphabet, and agreement
//! between the named predicates and the patterns their documentation claims.
//!
//! That last one is the point of the file. Each `is_*` predicate is
//! hand-written against a handful of cells and carries its OGC pattern in a
//! doc comment. Nothing checked that the code and the comment agreed until the
//! pattern matcher existed; now the comment is executable, and a predicate
//! rewritten to read a different cell fails here rather than silently
//! disagreeing with the standard it cites.

use super::*;
use crate::space::r2::{LineString, Point, Polygon};

fn square(x0: f64, y0: f64, x1: f64, y1: f64) -> Geometry {
    Geometry::Polygon(Polygon::new(
        LineString::from(vec![(x0, y0), (x1, y0), (x1, y1), (x0, y1), (x0, y0)]),
        Vec::new(),
    ))
}

fn point(x: f64, y: f64) -> Geometry {
    Geometry::Point(Point::new(x, y))
}

fn line(coords: &[(f64, f64)]) -> Geometry {
    Geometry::LineString(LineString::from(coords.to_vec()))
}

#[test]
fn disjoint_squares_render_the_canonical_disjoint_matrix() {
    let im = relate(&square(0.0, 0.0, 1.0, 1.0), &square(5.0, 5.0, 6.0, 6.0));
    // Two disjoint areal geometries: no interior or boundary meeting, and each
    // one's parts all lie in the other's exterior.
    assert_eq!(im.de9im(), "FF2FF1212");
    assert!(im.is_disjoint());
    assert!(!im.is_intersects());
}

#[test]
fn identical_squares_render_the_canonical_equal_matrix() {
    let im = relate(&square(0.0, 0.0, 1.0, 1.0), &square(0.0, 0.0, 1.0, 1.0));
    assert_eq!(im.de9im(), "2FFF1FFF2");
    assert!(im.is_equal_topo());
}

#[test]
fn a_rendered_matrix_never_contains_the_pattern_character() {
    // `T` means "any dimension" and belongs to a pattern, not to a matrix that
    // already knows the dimension. A renderer emitting it would be throwing
    // away what it computed.
    for (a, b) in [
        (square(0.0, 0.0, 1.0, 1.0), square(0.5, 0.5, 2.0, 2.0)),
        (point(0.0, 0.0), square(0.0, 0.0, 1.0, 1.0)),
        (line(&[(0.0, 0.0), (1.0, 1.0)]), square(0.0, 0.0, 1.0, 1.0)),
    ] {
        let rendered = relate(&a, &b).de9im();
        assert_eq!(rendered.len(), 9, "{rendered}");
        assert!(
            !rendered.contains('T') && !rendered.contains('t'),
            "a rendered matrix must spell dimensions, not `T`: {rendered}"
        );
        assert!(
            rendered.chars().all(|c| matches!(c, 'F' | '0' | '1' | '2')),
            "unexpected character in {rendered}"
        );
    }
}

#[test]
fn cells_agree_with_the_rendered_string() {
    let im = relate(&square(0.0, 0.0, 1.0, 1.0), &square(0.5, 0.5, 2.0, 2.0));
    let rendered: Vec<char> = im.de9im().chars().collect();
    for (row_index, row) in im.cells().iter().enumerate() {
        for (column, dim) in row.iter().enumerate() {
            let expected = match dim {
                0 => '0',
                1 => '1',
                2 => '2',
                _ => 'F',
            };
            assert_eq!(rendered[row_index * 3 + column], expected);
        }
    }
}

#[test]
fn the_pattern_alphabet_reads_each_cell() {
    let im = relate(&square(0.0, 0.0, 1.0, 1.0), &square(0.5, 0.5, 2.0, 2.0));
    let exact = im.de9im();
    // A matrix always satisfies its own rendering, and always satisfies the
    // wholly unconstrained pattern.
    assert_eq!(im.matches(&exact), Some(true));
    assert_eq!(im.matches("*********"), Some(true));
    // `T` is satisfied by any dimension, so replacing every digit with `T`
    // must still match, while replacing them with `F` must not.
    let anything: String = exact
        .chars()
        .map(|c| if c == 'F' { 'F' } else { 'T' })
        .collect();
    assert_eq!(im.matches(&anything), Some(true));
    let nothing: String = exact.chars().map(|_| 'F').collect();
    assert_eq!(im.matches(&nothing), Some(false));
    // Lowercase is accepted.
    assert_eq!(im.matches(&anything.to_lowercase()), Some(true));
}

#[test]
fn a_malformed_pattern_is_rejected_rather_than_unsatisfied() {
    let im = relate(&square(0.0, 0.0, 1.0, 1.0), &square(0.5, 0.5, 2.0, 2.0));
    // Eight characters, ten characters, and an unknown character are all
    // `None`. A caller that got `Some(false)` for these could not tell a typo
    // from a geometry that genuinely fails the pattern.
    assert_eq!(im.matches("TTTTTTTT"), None);
    assert_eq!(im.matches("TTTTTTTTTT"), None);
    assert_eq!(im.matches("TTTTTTTTx"), None);
    assert_eq!(im.matches(""), None);
}

/// Each named predicate against the OGC pattern(s) its doc comment cites.
#[test]
fn the_named_predicates_agree_with_their_documented_patterns() {
    let cases = [
        (square(0.0, 0.0, 1.0, 1.0), square(5.0, 5.0, 6.0, 6.0)),
        (square(0.0, 0.0, 1.0, 1.0), square(0.5, 0.5, 2.0, 2.0)),
        (square(0.0, 0.0, 1.0, 1.0), square(0.0, 0.0, 1.0, 1.0)),
        (square(0.0, 0.0, 4.0, 4.0), square(1.0, 1.0, 2.0, 2.0)),
        (square(0.0, 0.0, 1.0, 1.0), square(1.0, 0.0, 2.0, 1.0)),
        (point(0.5, 0.5), square(0.0, 0.0, 1.0, 1.0)),
        (point(0.0, 0.0), square(0.0, 0.0, 1.0, 1.0)),
        (line(&[(0.0, 0.0), (2.0, 2.0)]), square(0.0, 0.0, 1.0, 1.0)),
    ];
    for (a, b) in &cases {
        let im = relate(a, b);
        let rendered = im.de9im();
        let matches = |pattern: &str| {
            im.matches(pattern)
                .unwrap_or_else(|| panic!("pattern {pattern} rejected"))
        };

        // `disjoint` is `FF*FF****`; `intersects` is its negation.
        assert_eq!(
            im.is_disjoint(),
            matches("FF*FF****"),
            "disjoint {rendered}"
        );
        assert_eq!(
            im.is_intersects(),
            !im.is_disjoint(),
            "intersects {rendered}"
        );

        // `within` is `T*F**F***`, and `contains` is its transpose.
        assert_eq!(im.is_within(), matches("T*F**F***"), "within {rendered}");
        assert_eq!(
            im.is_contains(),
            matches("T*****FF*"),
            "contains {rendered}"
        );

        // `coveredBy` is the union of four patterns; `covers` is its transpose.
        assert_eq!(
            im.is_coveredby(),
            matches("T*F**F***")
                || matches("*TF**F***")
                || matches("**FT*F***")
                || matches("**F*TF***"),
            "coveredby {rendered}"
        );
        assert_eq!(
            im.is_covers(),
            matches("T*****FF*")
                || matches("*T****FF*")
                || matches("***T**FF*")
                || matches("****T*FF*"),
            "covers {rendered}"
        );

        // `equals` is `T*F**FFF*`.
        assert_eq!(
            im.is_equal_topo(),
            matches("T*F**FFF*"),
            "equals {rendered}"
        );

        // `touches` is `FT*******` | `F**T*****` | `F***T****`.
        assert_eq!(
            im.is_touches(),
            matches("FT*******") || matches("F**T*****") || matches("F***T****"),
            "touches {rendered}"
        );
    }
}

#[test]
fn transpose_swaps_the_operands() {
    let a = square(0.0, 0.0, 4.0, 4.0);
    let b = square(1.0, 1.0, 2.0, 2.0);
    let forward = relate(&a, &b);
    let reverse = relate(&b, &a);
    // `contains` and `within` are each other's transpose, and the rendered
    // matrices are transposes of one another.
    assert!(forward.is_contains());
    assert!(reverse.is_within());
    let forward_cells = forward.cells();
    let reverse_cells = reverse.cells();
    for row in 0..3 {
        for column in 0..3 {
            assert_eq!(
                forward_cells[row][column], reverse_cells[column][row],
                "cell ({row},{column}) is not the transpose"
            );
        }
    }
}
