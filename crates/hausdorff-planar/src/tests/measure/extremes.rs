use crate::space::r2::Extremes;
use crate::space::r2::{Geometry, MultiPoint};
use crate::tests::support::{EPS, closed, line_string, polygon};

#[derive(Clone, Copy, Debug, PartialEq)]
struct RefExtreme {
    index: usize,
    coord: (f64, f64),
}

#[derive(Clone, Copy, Debug)]
struct RefExtremes {
    x_min: RefExtreme,
    y_min: RefExtreme,
    x_max: RefExtreme,
    y_max: RefExtreme,
}

/// Compute extremes over a flat coordinate walk, exactly mirroring the kernel:
/// the *first* coordinate attaining an extreme value wins (strict `<` / `>`).
fn ref_extremes(coords: &[(f64, f64)]) -> Option<RefExtremes> {
    let &first = coords.first()?;
    let mut out = RefExtremes {
        x_min: RefExtreme {
            index: 0,
            coord: first,
        },
        y_min: RefExtreme {
            index: 0,
            coord: first,
        },
        x_max: RefExtreme {
            index: 0,
            coord: first,
        },
        y_max: RefExtreme {
            index: 0,
            coord: first,
        },
    };
    for (index, &coord) in coords.iter().enumerate().skip(1) {
        if coord.0 < out.x_min.coord.0 {
            out.x_min = RefExtreme { index, coord };
        }
        if coord.1 < out.y_min.coord.1 {
            out.y_min = RefExtreme { index, coord };
        }
        if coord.0 > out.x_max.coord.0 {
            out.x_max = RefExtreme { index, coord };
        }
        if coord.1 > out.y_max.coord.1 {
            out.y_max = RefExtreme { index, coord };
        }
    }
    Some(out)
}

#[test]
fn extremes_match_reference() {
    let diamond = vec![(1.0, 0.0), (2.0, 1.0), (1.0, 2.0), (0.0, 1.0), (1.0, 0.0)];
    let irregular = vec![
        (1.0, -3.1),
        (3.7, 2.7),
        (0.9, 7.6),
        (-4.8, 6.7),
        (-7.5, 0.9),
    ];

    // The exterior coordinate walk a geometry exposes to `Extremes`:
    //  - a polygon: its closed exterior ring;
    //  - a linestring: its coords as given;
    //  - a multipoint: its points in order.
    let h_cases: Vec<Geometry> = vec![
        Geometry::Polygon(polygon(&diamond, &[])),
        Geometry::Polygon(polygon(&irregular, &[])),
        Geometry::LineString(line_string(&irregular)),
        Geometry::MultiPoint(MultiPoint::from(irregular.clone())),
    ];
    let r_walks: Vec<Vec<(f64, f64)>> = vec![
        closed(&diamond),
        closed(&irregular),
        irregular.clone(),
        irregular.clone(),
    ];

    for (h, walk) in h_cases.iter().zip(&r_walks) {
        let ho = h.extremes().unwrap();
        let ro = ref_extremes(walk).unwrap();
        for (name, he, re) in [
            ("x_min", ho.x_min, ro.x_min),
            ("y_min", ho.y_min, ro.y_min),
            ("x_max", ho.x_max, ro.x_max),
            ("y_max", ho.y_max, ro.y_max),
        ] {
            assert_eq!(he.index, re.index, "{name} index on {h:?}");
            assert!(
                (he.coord.x - re.coord.0).abs() < EPS && (he.coord.y - re.coord.1).abs() < EPS,
                "{name} coord on {h:?}: {:?} vs ref {:?}",
                he.coord,
                re.coord
            );
        }
    }

    assert!(MultiPoint::new(vec![]).extremes().is_none());
}

#[test]
fn extremes_ignore_interior_rings() {
    let exterior = vec![(0.0, 0.0), (4.0, 0.0), (4.0, 4.0), (0.0, 4.0)];
    let hole = vec![(1.0, 1.0), (3.0, 1.0), (3.0, 3.0), (1.0, 3.0)];
    let h = polygon(&exterior, std::slice::from_ref(&hole));
    // The reference walks only the (closed) exterior ring, never the hole.
    let ro = ref_extremes(&closed(&exterior)).unwrap();
    let ho = h.extremes().unwrap();
    assert_eq!(ho.x_max.index, ro.x_max.index);
    assert_eq!(ho.x_max.coord.x, ro.x_max.coord.0);
    assert!(ho.x_max.index < h.exterior().coords_count());
}
