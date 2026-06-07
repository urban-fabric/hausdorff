//! Cover an arbitrary [`Geometry`] with H3 cells — the geometry-enum dispatch
//! over [`H3Polygon`] and [`H3Line`], so a caller indexing a heterogeneous
//! geometry column (e.g. galois's `wkb_h3_cell_ids`) does not re-derive the
//! per-variant logic. The `Series`/WKB plumbing stays in the caller; this is the
//! pure geometry → cells core.

use std::collections::BTreeSet;

use crate::space::r2::{Geometry, LineString, Polygon};

use super::{H3Cell, H3Line, H3Polygon};

/// H3 cells *approximately* covering `geom` at `resolution`:
/// - **points** → the containing cell;
/// - **lines** → the cells along the H3 grid path between consecutive vertex
///   cells ([`H3Line`]);
/// - **areal** → the interior centroid cells ([`H3Polygon::to_cells`]) plus the
///   grid-path cover of the boundary rings, which pulls in many edge-straddling
///   cells the centroid fill omits;
/// - **collections** recurse.
///
/// **This is an approximation, not a guaranteed superset.** The boundary is
/// traced by [`H3Line`], which walks the *shortest H3 grid path* between
/// consecutive vertex cells — not an exact segment↔cell intersection — so a cell
/// the true edge clips can still be missed (and the centroid interior fill can
/// miss thin features). It is suitable for bucketing / coarse partitioning; do
/// **not** rely on it as a no-false-negative spatial index without first
/// densifying toward a true boundary-intersection cover.
///
/// Returns `None` when no cover can be produced at all: a point is unmappable, a
/// ring is invalid (degenerate, out-of-range coordinate, or transmeridian), a
/// grid path is undefined (vertices too far apart or split by a pentagon), or
/// `resolution` is out of range.
#[must_use]
pub fn cover_geometry(geom: &Geometry, resolution: u8) -> Option<Vec<H3Cell>> {
    let mut cells = BTreeSet::new();
    collect(geom, resolution, &mut cells)?;
    Some(cells.into_iter().collect())
}

fn collect(geom: &Geometry, res: u8, out: &mut BTreeSet<H3Cell>) -> Option<()> {
    match geom {
        Geometry::Point(p) => {
            out.insert(H3Cell::for_point(p.y(), p.x(), res)?);
        }
        Geometry::MultiPoint(mp) => {
            for p in mp.iter() {
                out.insert(H3Cell::for_point(p.y(), p.x(), res)?);
            }
        }
        Geometry::Line(l) => {
            cover_line(vec![(l.start.y, l.start.x), (l.end.y, l.end.x)], res, out)?;
        }
        Geometry::LineString(ls) => cover_line(lat_lng(ls), res, out)?,
        Geometry::MultiLineString(mls) => {
            for ls in mls.iter() {
                cover_line(lat_lng(ls), res, out)?;
            }
        }
        Geometry::Polygon(poly) => cover_polygon(poly, res, out)?,
        Geometry::MultiPolygon(mp) => {
            for poly in mp.iter() {
                cover_polygon(poly, res, out)?;
            }
        }
        Geometry::Rect(r) => cover_polygon(&r.to_polygon(), res, out)?,
        Geometry::Triangle(t) => cover_polygon(&t.to_polygon(), res, out)?,
        Geometry::GeometryCollection(gs) => {
            for g in gs {
                collect(g, res, out)?;
            }
        }
    }
    Some(())
}

/// `(latitude, longitude)` vertices of a ring — the convention `H3Polygon`/
/// `H3Line` take (geo stores `x = lon`, `y = lat`).
fn lat_lng(ring: &LineString) -> Vec<(f64, f64)> {
    ring.coords().map(|c| (c.y, c.x)).collect()
}

fn cover_line(points: Vec<(f64, f64)>, res: u8, out: &mut BTreeSet<H3Cell>) -> Option<()> {
    out.extend(H3Line::new(points).ok()?.to_cells(res).ok()?);
    Some(())
}

fn cover_polygon(poly: &Polygon, res: u8, out: &mut BTreeSet<H3Cell>) -> Option<()> {
    let exterior = lat_lng(poly.exterior());
    let holes: Vec<Vec<(f64, f64)>> = poly.interiors().iter().map(lat_lng).collect();

    out.extend(
        H3Polygon::new(exterior.clone(), holes.clone())
            .ok()?
            .to_cells(res)
            .ok()?,
    );

    for ring in std::iter::once(exterior).chain(holes) {
        out.extend(H3Line::new(ring).ok()?.to_cells(res).ok()?);
    }
    Some(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::space::r2::{LineString, Point, Polygon};

    #[test]
    fn point_covers_its_cell() {
        let g = Geometry::Point(Point::new(10.0, 50.0));
        let cells = cover_geometry(&g, 8).unwrap();
        assert_eq!(cells.len(), 1);
        assert_eq!(cells[0], H3Cell::for_point(50.0, 10.0, 8).unwrap());
    }

    #[test]
    fn line_is_covered_via_grid_path() {
        let g = Geometry::LineString(LineString::from(vec![(10.0, 50.0), (10.5, 50.0)]));
        let cells = cover_geometry(&g, 8).unwrap();
        assert!(cells.len() >= 2);
        assert!(cells.contains(&H3Cell::for_point(50.0, 10.0, 8).unwrap()));
        assert!(cells.contains(&H3Cell::for_point(50.0, 10.5, 8).unwrap()));
    }

    #[test]
    fn polygon_cover_includes_boundary_vertex_cells() {
        let verts = vec![
            (9.5, 49.5),
            (10.5, 49.5),
            (10.5, 50.5),
            (9.5, 50.5),
            (9.5, 49.5),
        ];
        let g = Geometry::Polygon(Polygon::new(LineString::from(verts.clone()), vec![]));
        let cells: std::collections::BTreeSet<H3Cell> =
            cover_geometry(&g, 5).unwrap().into_iter().collect();
        assert!(!cells.is_empty());
        for &(lon, lat) in &verts[..4] {
            assert!(
                cells.contains(&H3Cell::for_point(lat, lon, 5).unwrap()),
                "boundary vertex ({lon},{lat}) cell missing"
            );
        }
    }
}
