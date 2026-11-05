use crate::format::shapefile::{ShapeType, ShapefileError, ShapefileReader, read_geometries};
use crate::space::r2::{Geometry, MultiLineString, Point};

fn i32_be(out: &mut Vec<u8>, v: i32) {
    out.extend_from_slice(&v.to_be_bytes());
}

fn i32_le(out: &mut Vec<u8>, v: i32) {
    out.extend_from_slice(&v.to_le_bytes());
}

fn f64_le(out: &mut Vec<u8>, v: f64) {
    out.extend_from_slice(&v.to_le_bytes());
}

fn record(number: i32, payload: Vec<u8>) -> Vec<u8> {
    assert_eq!(payload.len() % 2, 0);
    let mut out = Vec::new();
    i32_be(&mut out, number);
    i32_be(&mut out, (payload.len() / 2) as i32);
    out.extend(payload);
    out
}

fn shapefile(shape_type: ShapeType, records: Vec<Vec<u8>>) -> Vec<u8> {
    let len = 100 + records.iter().map(Vec::len).sum::<usize>();
    assert_eq!(len % 2, 0);

    let mut out = Vec::new();
    i32_be(&mut out, 9994);
    out.extend_from_slice(&[0; 20]);
    i32_be(&mut out, (len / 2) as i32);
    i32_le(&mut out, 1000);
    i32_le(&mut out, shape_type.code());
    for _ in 0..8 {
        f64_le(&mut out, 0.0);
    }
    assert_eq!(out.len(), 100);
    for record in records {
        out.extend(record);
    }
    out
}

fn point_payload(shape_type: ShapeType, x: f64, y: f64) -> Vec<u8> {
    let mut out = Vec::new();
    i32_le(&mut out, shape_type.code());
    f64_le(&mut out, x);
    f64_le(&mut out, y);
    out
}

fn point_z_payload(x: f64, y: f64, z: f64, m: f64) -> Vec<u8> {
    let mut out = point_payload(ShapeType::PointZ, x, y);
    f64_le(&mut out, z);
    f64_le(&mut out, m);
    out
}

fn multipart_payload(shape_type: ShapeType, parts: &[&[(f64, f64)]]) -> Vec<u8> {
    let mut out = Vec::new();
    i32_le(&mut out, shape_type.code());

    let mut all = parts.iter().flat_map(|part| part.iter().copied());
    let first = all.next().unwrap_or((0.0, 0.0));
    let (mut min_x, mut min_y, mut max_x, mut max_y) = (first.0, first.1, first.0, first.1);
    for (x, y) in all {
        min_x = min_x.min(x);
        min_y = min_y.min(y);
        max_x = max_x.max(x);
        max_y = max_y.max(y);
    }
    for v in [min_x, min_y, max_x, max_y] {
        f64_le(&mut out, v);
    }

    i32_le(&mut out, parts.len() as i32);
    i32_le(
        &mut out,
        parts.iter().map(|part| part.len()).sum::<usize>() as i32,
    );
    let mut start = 0;
    for part in parts {
        i32_le(&mut out, start);
        start += part.len() as i32;
    }
    for part in parts {
        for &(x, y) in *part {
            f64_le(&mut out, x);
            f64_le(&mut out, y);
        }
    }
    out
}

#[test]
fn reads_pointz_as_xy_point() {
    let bytes = shapefile(
        ShapeType::PointZ,
        vec![record(1, point_z_payload(1.25, -2.5, 99.0, 7.0))],
    );
    let geometries = read_geometries(bytes.as_slice()).unwrap();
    assert_eq!(
        geometries,
        vec![Some(Geometry::Point(Point::new(1.25, -2.5)))]
    );
}

#[test]
fn reads_multipart_polyline() {
    let a = &[(0.0, 0.0), (1.0, 1.0)][..];
    let b = &[(5.0, 5.0), (6.0, 7.0), (8.0, 9.0)][..];
    let bytes = shapefile(
        ShapeType::PolyLine,
        vec![record(7, multipart_payload(ShapeType::PolyLine, &[a, b]))],
    );
    let geometries = read_geometries(bytes.as_slice()).unwrap();
    let Some(Geometry::MultiLineString(MultiLineString(lines))) = &geometries[0] else {
        panic!("expected multipart line");
    };
    assert_eq!(lines.len(), 2);
    assert_eq!(lines[0].coords_count(), 2);
    assert_eq!(lines[1].coords_count(), 3);
}

#[test]
fn reads_polygon_outer_with_inner_ring() {
    let outer_clockwise = &[(0.0, 0.0), (0.0, 4.0), (4.0, 4.0), (4.0, 0.0), (0.0, 0.0)][..];
    let hole_counter_clockwise = &[(1.0, 1.0), (2.0, 1.0), (2.0, 2.0), (1.0, 2.0), (1.0, 1.0)][..];
    let bytes = shapefile(
        ShapeType::Polygon,
        vec![record(
            3,
            multipart_payload(
                ShapeType::Polygon,
                &[outer_clockwise, hole_counter_clockwise],
            ),
        )],
    );
    let geometries = read_geometries(bytes.as_slice()).unwrap();
    let Some(Geometry::Polygon(poly)) = &geometries[0] else {
        panic!("expected polygon");
    };
    assert_eq!(poly.exterior().coords_count(), 5);
    assert_eq!(poly.interiors().len(), 1);
    assert_eq!(poly.interiors()[0].coords_count(), 5);
}

#[test]
fn rejects_wrong_file_code() {
    let mut bytes = shapefile(
        ShapeType::Point,
        vec![record(1, point_payload(ShapeType::Point, 0.0, 0.0))],
    );
    bytes[..4].copy_from_slice(&42_i32.to_be_bytes());
    match ShapefileReader::new(bytes.as_slice()) {
        Err(ShapefileError::InvalidFileCode(42)) => {}
        Err(err) => panic!("expected invalid file code, got {err}"),
        Ok(_) => panic!("expected invalid file code, got reader"),
    }
}
