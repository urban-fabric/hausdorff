//! Extreme (minimum / maximum x and y) coordinates of a geometry with their
//! positions — mirror of `geo`'s `Extremes` trait and its `Outcome` /
//! `Extreme` result types.
//!
//! Indices count along [`CoordsIter::exterior_coords_iter`], so interior rings
//! never produce an extreme. Strict comparisons mean the *first* coordinate at
//! an extreme value wins ties, matching `geo`.

use super::{Coord, CoordsIter};

/// One extreme coordinate and its index along the exterior coordinate walk.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Extreme {
    pub index: usize,
    pub coord: Coord,
}

/// The four extremes of a geometry.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Outcome {
    pub x_min: Extreme,
    pub y_min: Extreme,
    pub x_max: Extreme,
    pub y_max: Extreme,
}

/// Find the extreme coordinates and indices of a geometry, or `None` when it
/// has no coordinates.
///
/// Proof-backed by `riemann` `Hausdorff/Extremes`: the running min/max fold
/// returns a genuine extreme — its key bounds every candidate (`argmin_le_mem` /
/// `argmax_ge_mem`) and is **attained** by an input coordinate (`argmin_mem` /
/// `argmax_mem`), with the strict `<`/`>` keeping the first coordinate on ties.
pub trait Extremes {
    fn extremes(&self) -> Option<Outcome>;
}

impl<G: CoordsIter> Extremes for G {
    fn extremes(&self) -> Option<Outcome> {
        let mut iter = self.exterior_coords_iter().enumerate();

        let mut outcome = iter.next().map(|(index, coord)| Outcome {
            x_min: Extreme { index, coord },
            y_min: Extreme { index, coord },
            x_max: Extreme { index, coord },
            y_max: Extreme { index, coord },
        })?;

        for (index, coord) in iter {
            if coord.x < outcome.x_min.coord.x {
                outcome.x_min = Extreme { index, coord };
            }
            if coord.y < outcome.y_min.coord.y {
                outcome.y_min = Extreme { index, coord };
            }
            if coord.x > outcome.x_max.coord.x {
                outcome.x_max = Extreme { index, coord };
            }
            if coord.y > outcome.y_max.coord.y {
                outcome.y_max = Extreme { index, coord };
            }
        }

        Some(outcome)
    }
}

#[cfg(test)]
mod tests {
    use super::super::{LineString, MultiPoint, Polygon};
    use super::*;

    #[test]
    fn diamond_polygon() {
        let polygon = Polygon::new(
            LineString::from(vec![
                (1.0, 0.0),
                (2.0, 1.0),
                (1.0, 2.0),
                (0.0, 1.0),
                (1.0, 0.0),
            ]),
            vec![],
        );

        let actual = polygon.extremes();

        assert_eq!(
            Some(Outcome {
                x_min: Extreme {
                    index: 3,
                    coord: Coord::new(0.0, 1.0)
                },
                y_min: Extreme {
                    index: 0,
                    coord: Coord::new(1.0, 0.0)
                },
                x_max: Extreme {
                    index: 1,
                    coord: Coord::new(2.0, 1.0)
                },
                y_max: Extreme {
                    index: 2,
                    coord: Coord::new(1.0, 2.0)
                },
            }),
            actual
        );
    }

    #[test]
    fn empty_is_none() {
        assert_eq!(MultiPoint::new(vec![]).extremes(), None);
        assert_eq!(LineString::new(vec![]).extremes(), None);
    }

    #[test]
    fn first_wins_ties() {
        let mp = MultiPoint::from(vec![(0.0, 0.0), (1.0, 0.0), (0.0, 1.0), (1.0, 1.0)]);
        let extremes = mp.extremes().unwrap();
        assert_eq!(extremes.x_min.index, 0);
        assert_eq!(extremes.y_min.index, 0);
        assert_eq!(extremes.x_max.index, 1);
        assert_eq!(extremes.y_max.index, 2);
    }

    #[test]
    fn interior_rings_are_ignored() {
        let polygon = Polygon::new(
            LineString::from(vec![(0.0, 0.0), (4.0, 0.0), (4.0, 4.0), (0.0, 4.0)]),
            vec![LineString::from(vec![
                (1.0, 1.0),
                (3.0, 1.0),
                (3.0, 3.0),
                (1.0, 3.0),
            ])],
        );
        let extremes = polygon.extremes().unwrap();
        assert_eq!(extremes.x_max.coord, Coord::new(4.0, 0.0));
        assert!(extremes.x_max.index < polygon.exterior().coords_count());
    }
}
