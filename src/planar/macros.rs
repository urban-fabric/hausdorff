//! Construction macros mirroring `geo`'s `coord!` / `point!` / `line_string!` /
//! `polygon!`, so call sites migrate by swapping the import path only.

/// Build a [`Coord`](crate::planar::Coord): `coord! { x: 1.0, y: 2.0 }`.
#[macro_export]
macro_rules! coord {
    (x: $x:expr, y: $y:expr $(,)?) => {
        $crate::planar::Coord { x: $x, y: $y }
    };
}

/// Build a [`Point`](crate::planar::Point): `point! { x: 1.0, y: 2.0 }`.
#[macro_export]
macro_rules! point {
    (x: $x:expr, y: $y:expr $(,)?) => {
        $crate::planar::Point::new($x, $y)
    };
}

/// Build a [`LineString`](crate::planar::LineString) from `(x:, y:)` vertices:
/// `line_string![(x: 0., y: 0.), (x: 1., y: 1.)]`.
#[macro_export]
macro_rules! line_string {
    () => {
        $crate::planar::LineString::new(::std::vec::Vec::new())
    };
    ($((x: $x:expr, y: $y:expr $(,)?)),* $(,)?) => {
        $crate::planar::LineString::new(::std::vec![
            $( $crate::planar::Coord { x: $x, y: $y } ),*
        ])
    };
    ($($coord:expr),* $(,)?) => {
        $crate::planar::LineString::new(::std::vec![
            $( $crate::planar::Coord::from($coord) ),*
        ])
    };
}

/// Build a [`Polygon`](crate::planar::Polygon).
///
/// Exterior-only from `(x:, y:)` vertices:
/// `polygon![(x: 0., y: 0.), (x: 1., y: 0.), (x: 0., y: 1.)]`,
/// or with holes:
/// `polygon!(exterior: [(x:,y:), ...], interiors: [[(x:,y:), ...], ...])`.
#[macro_export]
macro_rules! polygon {
    () => {
        $crate::planar::Polygon::new($crate::line_string![], ::std::vec::Vec::new())
    };
    (
        exterior: [ $((x: $ex:expr, y: $ey:expr $(,)?)),* $(,)? ],
        interiors: [ $([ $((x: $ix:expr, y: $iy:expr $(,)?)),* $(,)? ]),* $(,)? ]
        $(,)?
    ) => {
        $crate::planar::Polygon::new(
            $crate::line_string![ $((x: $ex, y: $ey)),* ],
            ::std::vec![
                $( $crate::line_string![ $((x: $ix, y: $iy)),* ] ),*
            ],
        )
    };
    ($((x: $x:expr, y: $y:expr $(,)?)),* $(,)?) => {
        $crate::planar::Polygon::new(
            $crate::line_string![ $((x: $x, y: $y)),* ],
            ::std::vec::Vec::new(),
        )
    };
}

#[cfg(test)]
mod tests {
    use crate::planar::{Coord, Point};

    #[test]
    fn coord_and_point_macros() {
        assert_eq!(coord! { x: 1.0, y: 2.0 }, Coord::new(1.0, 2.0));
        assert_eq!(point! { x: 1.0, y: 2.0 }, Point::new(1.0, 2.0));
    }

    #[test]
    fn line_string_macro_tuples() {
        let ls = line_string![(x: 0.0, y: 0.0), (x: 1.0, y: 1.0)];
        assert_eq!(ls.coords_count(), 2);
        assert_eq!(ls.0[1], Coord::new(1.0, 1.0));
    }

    #[test]
    fn polygon_macro_exterior_only() {
        let p = polygon![(x: 0.0, y: 0.0), (x: 1.0, y: 0.0), (x: 0.0, y: 1.0)];
        assert!(p.exterior().is_closed());
        assert_eq!(p.exterior().coords_count(), 4); // closed
        assert!(p.interiors().is_empty());
    }

    #[test]
    fn polygon_macro_with_hole() {
        let p = polygon!(
            exterior: [(x: 0.0, y: 0.0), (x: 4.0, y: 0.0), (x: 4.0, y: 4.0), (x: 0.0, y: 4.0)],
            interiors: [[(x: 1.0, y: 1.0), (x: 2.0, y: 1.0), (x: 2.0, y: 2.0), (x: 1.0, y: 2.0)]]
        );
        assert_eq!(p.interiors().len(), 1);
        assert!(p.interiors()[0].is_closed());
    }
}
