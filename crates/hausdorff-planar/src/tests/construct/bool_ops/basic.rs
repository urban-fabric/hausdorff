use super::*;

#[test]
fn disjoint() {
    check_rects(&[(0.0, 0.0, 1.0, 1.0)], &[(5.0, 5.0, 6.0, 6.0)]);
}

#[test]
fn nested() {
    check_rects(&[(0.0, 0.0, 10.0, 10.0)], &[(2.0, 2.0, 4.0, 4.0)]);
}

#[test]
fn crossing_squares() {
    check_rects(&[(0.0, 0.0, 4.0, 4.0)], &[(2.0, 2.0, 6.0, 6.0)]);
}

#[test]
fn identical() {
    check_rects(&[(0.0, 0.0, 3.0, 3.0)], &[(0.0, 0.0, 3.0, 3.0)]);
}

#[test]
fn shared_edge() {
    check_rects(&[(0.0, 0.0, 2.0, 2.0)], &[(2.0, 0.0, 4.0, 2.0)]);
}

#[test]
fn partial_edge_overlap() {
    check_rects(&[(0.0, 0.0, 4.0, 2.0)], &[(2.0, 0.0, 6.0, 2.0)]);
}

#[test]
fn t_vertex_touching() {
    check_rects(&[(0.0, 0.0, 2.0, 2.0)], &[(2.0, 0.0, 4.0, 1.0)]);
}

#[test]
fn t_vertex_grid_tiles() {
    // left rectangle vs two stacked right rectangles (all touching along x = 4).
    check_rects(
        &[(0.0, 0.0, 4.0, 2.0)],
        &[(4.0, 0.0, 6.0, 1.0), (4.0, 1.0, 6.0, 2.0)],
    );
}

#[test]
fn t_vertex_overlap() {
    check_rects(&[(0.0, 0.0, 4.0, 4.0)], &[(2.0, 1.0, 6.0, 3.0)]);
}

#[test]
fn multipolygon_inputs() {
    check_rects(
        &[(0.0, 0.0, 2.0, 2.0), (4.0, 4.0, 6.0, 6.0)],
        &[(1.0, 1.0, 5.0, 5.0)],
    );
}

// ---- concave / holed fixtures: golden intersection areas captured from the
// validated implementation (HOLE and STAIR additionally cross-checked by hand). ----

#[test]
fn t_vertex_staircase_union() {
    let big = sq(0.0, 0.0, 6.0, 6.0);
    let comb = vec![
        (3.0, 0.0),
        (6.0, 0.0),
        (6.0, 2.0),
        (4.0, 2.0),
        (4.0, 4.0),
        (6.0, 4.0),
        (6.0, 6.0),
        (3.0, 6.0),
    ];
    check_ops(&hmp(&[&[&big]]), &hmp(&[&[&comb]]), 14.0);
}

#[test]
fn concave_l_shapes() {
    let l1 = [
        (0.0, 0.0),
        (6.0, 0.0),
        (6.0, 2.0),
        (2.0, 2.0),
        (2.0, 6.0),
        (0.0, 6.0),
    ];
    let l2 = [
        (1.0, 1.0),
        (5.0, 1.0),
        (5.0, 5.0),
        (3.0, 5.0),
        (3.0, 3.0),
        (1.0, 3.0),
    ];
    check_ops(&hmp(&[&[&l1]]), &hmp(&[&[&l2]]), 5.0);
}

#[test]
fn star_into_square() {
    let star = [(0.0, 0.0), (6.0, 0.0), (6.0, 6.0), (3.0, 3.0), (0.0, 6.0)];
    check_ops(&hmp(&[&[&star]]), &hmp(&[&[&sq(1.0, 1.0, 5.0, 5.0)]]), 12.0);
}

#[test]
fn subject_with_hole() {
    let outer = sq(0.0, 0.0, 10.0, 10.0);
    let hole = sq(3.0, 3.0, 7.0, 7.0);
    check_ops(
        &hmp(&[&[&outer, &hole]]),
        &hmp(&[&[&sq(5.0, 5.0, 12.0, 12.0)]]),
        21.0,
    );
}
