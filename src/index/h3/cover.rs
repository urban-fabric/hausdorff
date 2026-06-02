//! Cover an arbitrary [`Geometry`] with H3 cells — the geometry-enum dispatch
//! over [`H3Polygon`] and [`H3Line`], so a caller indexing a heterogeneous
//! geometry column (e.g. galois's `wkb_h3_cell_ids`) does not re-derive the
//! per-variant logic. The `Series`/WKB plumbing stays in the caller; this is the
//! pure geometry → cells core.

use std::collections::BTreeSet;

use crate::space::r2::{Geometry, LineString, Polygon};

use super::{H3Cell, H3Line, H3Polygon};

/// H3 cells at `resolution` covering `geom`:
/// - **points** → the containing cell;
/// - **lines** → the grid-path cover ([`H3Line`]), every cell the path crosses;
/// - **areal** → the interior centroid cells ([`H3Polygon::to_cells`]) **plus**
///   the boundary rings' grid-path cover, so edge-overlapping cells (centre
///   outside) are included — a boundary-complete cover, not a centroid one;
/// - **collections** recurse.
///
/// Returns `None` when a point is unmappable, a polygon is invalid/transmeridian,
/// a grid path is undefined (vertices too far apart or split by a pentagon), or
/// `resolution` is out of range — i.e. when no honest cover can be produced.
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

    // Interior fill: cells whose centre is inside the polygon.
    out.extend(H3Polygon::new(exterior.clone(), holes.clone()).ok()?.to_cells(res).ok()?);

    // Boundary completion: the grid-path cover of every ring catches the
    // edge-straddling cells the centroid fill omits. The r2 rings are already
    // closed (first == last), so every edge — including the closing one — is
    // walked.
    for ring in std::iter::once(exterior).chain(holes) {
        if let Ok(line) = H3Line::new(ring) {
            if let Ok(cells) = line.to_cells(res) {
                out.extend(cells);
            }
        }
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
        assert!(cells.len() >= 2, "a spanning line covers several cells");
        // Endpoints' cells are in the cover.
        assert!(cells.contains(&H3Cell::for_point(50.0, 10.0, 8).unwrap()));
        assert!(cells.contains(&H3Cell::for_point(50.0, 10.5, 8).unwrap()));
    }

    #[test]
    fn polygon_cover_is_boundary_complete() {
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
        // Every boundary vertex's cell is covered (the centroid-only fill misses
        // these when the vertex cell's centre falls outside the polygon).
        for &(lon, lat) in &verts[..4] {
            assert!(
                cells.contains(&H3Cell::for_point(lat, lon, 5).unwrap()),
                "boundary vertex ({lon},{lat}) cell missing"
            );
        }
    }
}
