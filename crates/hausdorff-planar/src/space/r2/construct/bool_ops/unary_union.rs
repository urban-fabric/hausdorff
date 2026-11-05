//! n-way [`union`](BooleanOps::union) of many polygons at once, mirroring
//! `geo`'s `unary_union`. The inputs' component polygons are sorted by bounding
//! box and merged pairwise up a balanced binary tree (a cascaded union), so the
//! pairwise Martinez–Rueda kernel sees spatially adjacent operands of similar
//! size — more robust and asymptotically cheaper than a linear left fold.
//! Validated against `geo`'s `unary_union` in `tests/bool_ops.rs`.

use super::BooleanOps;
use crate::space::r2::{BoundingRect, MultiPolygon, Polygon};

/// An operand of [`unary_union`]: any geometry viewed as a set of component
/// polygons.
pub trait UnaryUnionOperand {
    /// The component polygons.
    fn component_polygons(&self) -> impl Iterator<Item = &Polygon>;
}

impl UnaryUnionOperand for Polygon {
    fn component_polygons(&self) -> impl Iterator<Item = &Polygon> {
        std::iter::once(self)
    }
}

impl UnaryUnionOperand for MultiPolygon {
    fn component_polygons(&self) -> impl Iterator<Item = &Polygon> {
        self.iter()
    }
}

/// The union of every input geometry, computed as a cascaded (balanced binary
/// tree) sequence of pairwise unions over spatially sorted components. Returns
/// an empty [`MultiPolygon`] for an empty input.
pub fn unary_union<'a, B: UnaryUnionOperand + 'a>(
    inputs: impl IntoIterator<Item = &'a B>,
) -> MultiPolygon {
    let mut layer: Vec<MultiPolygon> = inputs
        .into_iter()
        .flat_map(UnaryUnionOperand::component_polygons)
        .map(|p| MultiPolygon::new(vec![p.clone()]))
        .collect();

    layer.sort_by(|a, b| {
        let key = |mp: &MultiPolygon| {
            mp.bounding_rect()
                .map(|r| (r.min().x, r.min().y))
                .unwrap_or((f64::INFINITY, f64::INFINITY))
        };
        let (ka, kb) = (key(a), key(b));
        ka.0.total_cmp(&kb.0).then(ka.1.total_cmp(&kb.1))
    });

    while layer.len() > 1 {
        let mut next = Vec::with_capacity(layer.len().div_ceil(2));
        for pair in layer.chunks(2) {
            next.push(match pair {
                [a, b] => a.union(b),
                [a] => a.clone(),
                _ => unreachable!("chunks(2) yields 1 or 2 items"),
            });
        }
        layer = next;
    }
    layer.pop().unwrap_or_default()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::space::r2::{Area, Coord, LineString};

    fn square(x0: f64, y0: f64, size: f64) -> Polygon {
        Polygon::new(
            LineString::new(vec![
                Coord::new(x0, y0),
                Coord::new(x0 + size, y0),
                Coord::new(x0 + size, y0 + size),
                Coord::new(x0, y0 + size),
            ]),
            vec![],
        )
    }

    #[test]
    fn empty_input_yields_empty_multipolygon() {
        let polygons: Vec<Polygon> = Vec::new();
        assert!(unary_union(&polygons).is_empty());
    }

    #[test]
    fn single_polygon_round_trips() {
        let p = square(0.0, 0.0, 2.0);
        let out = unary_union(std::iter::once(&p));
        assert_eq!(out.len(), 1);
        assert!((out.unsigned_area() - 4.0).abs() < 1e-9);
    }

    #[test]
    fn disjoint_squares_keep_their_total_area() {
        let polygons = vec![square(0.0, 0.0, 1.0), square(3.0, 0.0, 1.0)];
        let out = unary_union(&polygons);
        assert_eq!(out.len(), 2);
        assert!((out.unsigned_area() - 2.0).abs() < 1e-9);
    }

    #[test]
    fn overlapping_squares_subtract_the_overlap_once() {
        let polygons = vec![square(0.0, 0.0, 2.0), square(1.0, 0.0, 2.0)];
        let out = unary_union(&polygons);
        assert_eq!(out.len(), 1);
        assert!((out.unsigned_area() - 6.0).abs() < 1e-9);
    }

    #[test]
    fn multipolygon_operands_are_flattened() {
        let mps = vec![
            MultiPolygon::new(vec![square(0.0, 0.0, 1.0), square(1.0, 0.0, 1.0)]),
            MultiPolygon::new(vec![square(2.0, 0.0, 1.0)]),
        ];
        let out = unary_union(&mps);
        assert_eq!(out.len(), 1);
        assert!((out.unsigned_area() - 3.0).abs() < 1e-9);
    }
}
