//! Construction macros mirroring `geo`'s `coord!` / `point!` / `line_string!` /
//! `polygon!`, so call sites migrate by swapping the import path only.

/// Build a [`Coord`](crate::r2::Coord): `coord! { x: 1.0, y: 2.0 }`.
#[macro_export]
macro_rules! coord {
    (x: $x:expr, y: $y:expr $(,)?) => {
        $crate::r2::Coord { x: $x, y: $y }
    };
}

/// Build a [`Point`](crate::r2::Point): `point! { x: 1.0, y: 2.0 }`.
#[macro_export]
macro_rules! point {
    (x: $x:expr, y: $y:expr $(,)?) => {
        $crate::r2::Point::new($x, $y)
    };
}

/// Build a [`LineString`](crate::r2::LineString) from `(x:, y:)` vertices:
/// `line_string![(x: 0., y: 0.), (x: 1., y: 1.)]`.
#[macro_export]
macro_rules! line_string {
    () => {
        $crate::r2::LineString::new(::std::vec::Vec::new())
    };
    ($((x: $x:expr, y: $y:expr $(,)?)),* $(,)?) => {
        $crate::r2::LineString::new(::std::vec![
            $( $crate::r2::Coord { x: $x, y: $y } ),*
        ])
    };
    ($($coord:expr),* $(,)?) => {
        $crate::r2::LineString::new(::std::vec![
            $( $crate::r2::Coord::from($coord) ),*
        ])
    };
}

/// Build a [`Polygon`](crate::r2::Polygon).
///
/// Exterior-only from `(x:, y:)` vertices:
/// `polygon![(x: 0., y: 0.), (x: 1., y: 0.), (x: 0., y: 1.)]`,
/// or with holes:
/// `polygon!(exterior: [(x:,y:), ...], interiors: [[(x:,y:), ...], ...])`.
#[macro_export]
macro_rules! polygon {
    () => {
        $crate::r2::Polygon::new($crate::line_string![], ::std::vec::Vec::new())
    };
    (
        exterior: [ $((x: $ex:expr, y: $ey:expr $(,)?)),* $(,)? ],
        interiors: [ $([ $((x: $ix:expr, y: $iy:expr $(,)?)),* $(,)? ]),* $(,)? ]
        $(,)?
    ) => {
        $crate::r2::Polygon::new(
            $crate::line_string![ $((x: $ex, y: $ey)),* ],
            ::std::vec![
                $( $crate::line_string![ $((x: $ix, y: $iy)),* ] ),*
            ],
        )
    };
    ($((x: $x:expr, y: $y:expr $(,)?)),* $(,)?) => {
        $crate::r2::Polygon::new(
            $crate::line_string![ $((x: $x, y: $y)),* ],
            ::std::vec::Vec::new(),
        )
    };
}
