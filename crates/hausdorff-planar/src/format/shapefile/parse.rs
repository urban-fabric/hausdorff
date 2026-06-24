use super::*;

pub(super) fn read_header<R: Read>(source: &mut R) -> Result<ShapefileHeader, ShapefileError> {
    let file_code = read_i32_be(source)?;
    if file_code != FILE_CODE {
        return Err(ShapefileError::InvalidFileCode(file_code));
    }

    skip_exact(source, 20)?;
    let file_length_bytes = words_to_bytes(read_i32_be(source)?, "file length")?;
    if file_length_bytes < HEADER_SIZE_BYTES {
        return Err(ShapefileError::InvalidHeader(
            "declared file length is smaller than header",
        ));
    }

    let version = read_i32_le(source)?;
    let shape_type = ShapeType::from_code(read_i32_le(source)?)?;
    let min_x = read_f64_le(source)?;
    let min_y = read_f64_le(source)?;
    let max_x = read_f64_le(source)?;
    let max_y = read_f64_le(source)?;
    let min_z = read_f64_le(source)?;
    let max_z = read_f64_le(source)?;
    let min_m = read_f64_le(source)?;
    let max_m = read_f64_le(source)?;

    Ok(ShapefileHeader {
        file_length_bytes,
        version,
        shape_type,
        bbox: Rect::new(Coord::new(min_x, min_y), Coord::new(max_x, max_y)),
        z_range: [min_z, max_z],
        m_range: [min_m, max_m],
    })
}

pub(super) fn parse_record(
    number: i32,
    payload: &[u8],
    header_shape_type: ShapeType,
) -> Result<ShapeRecord, ShapefileError> {
    let mut r = SliceReader::new(payload);
    let shape_type = ShapeType::from_code(r.i32_le()?)?;
    if shape_type != ShapeType::NullShape && shape_type != header_shape_type {
        return Err(ShapefileError::MismatchedRecordShapeType {
            header: header_shape_type,
            record: shape_type,
        });
    }

    let geometry = match shape_type {
        ShapeType::NullShape => None,
        ShapeType::Point | ShapeType::PointM | ShapeType::PointZ => Some(read_point(&mut r)?),
        ShapeType::MultiPoint | ShapeType::MultiPointM | ShapeType::MultiPointZ => {
            Some(read_multipoint(&mut r)?)
        }
        ShapeType::PolyLine | ShapeType::PolyLineM | ShapeType::PolyLineZ => {
            Some(read_polyline(&mut r)?)
        }
        ShapeType::Polygon | ShapeType::PolygonM | ShapeType::PolygonZ => {
            Some(read_polygon(&mut r)?)
        }
        ShapeType::MultiPatch => return Err(ShapefileError::UnsupportedGeometry(shape_type)),
    };

    Ok(ShapeRecord {
        number,
        shape_type,
        geometry,
    })
}

fn read_point(r: &mut SliceReader<'_>) -> Result<Geometry, ShapefileError> {
    let x = r.f64_le()?;
    let y = r.f64_le()?;
    Ok(Geometry::Point(Point::new(x, y)))
}

fn read_multipoint(r: &mut SliceReader<'_>) -> Result<Geometry, ShapefileError> {
    skip_bbox(r)?;
    let num_points = nonnegative_usize(r.i32_le()?, "num_points")?;
    let points = read_points(r, num_points)?;
    Ok(Geometry::MultiPoint(MultiPoint::new(
        points.into_iter().map(Point).collect(),
    )))
}

fn read_polyline(r: &mut SliceReader<'_>) -> Result<Geometry, ShapefileError> {
    let parts = read_parts_shape(r)?;
    let lines = parts.into_iter().map(LineString::new).collect::<Vec<_>>();
    if lines.len() == 1 {
        Ok(Geometry::LineString(lines.into_iter().next().unwrap()))
    } else {
        Ok(Geometry::MultiLineString(MultiLineString::new(lines)))
    }
}

fn read_polygon(r: &mut SliceReader<'_>) -> Result<Geometry, ShapefileError> {
    let rings = read_parts_shape(r)?;
    let mut polygons: Vec<Polygon> = Vec::new();
    let mut current: Option<(LineString, Vec<LineString>)> = None;

    for ring in rings {
        let line = LineString::new(ring);
        if ring_is_inner(&line) {
            let Some((_, holes)) = current.as_mut() else {
                return Err(ShapefileError::OrphanedInnerRing);
            };
            holes.push(line);
        } else {
            if let Some((exterior, holes)) = current.take() {
                polygons.push(Polygon::new(exterior, holes));
            }
            current = Some((line, Vec::new()));
        }
    }

    if let Some((exterior, holes)) = current.take() {
        polygons.push(Polygon::new(exterior, holes));
    }

    if polygons.len() == 1 {
        Ok(Geometry::Polygon(polygons.into_iter().next().unwrap()))
    } else {
        Ok(Geometry::MultiPolygon(MultiPolygon::new(polygons)))
    }
}

fn read_parts_shape(r: &mut SliceReader<'_>) -> Result<Vec<Vec<Coord>>, ShapefileError> {
    skip_bbox(r)?;
    let num_parts = nonnegative_usize(r.i32_le()?, "num_parts")?;
    let num_points = nonnegative_usize(r.i32_le()?, "num_points")?;
    let mut starts = Vec::with_capacity(num_parts);
    for _ in 0..num_parts {
        let start = nonnegative_usize(r.i32_le()?, "part_start")?;
        if start > num_points {
            return Err(ShapefileError::InvalidRecord(
                "part start is greater than point count",
            ));
        }
        starts.push(start);
    }

    let points = read_points(r, num_points)?;
    let mut parts = Vec::with_capacity(num_parts);
    for i in 0..starts.len() {
        let start = starts[i];
        let end = starts.get(i + 1).copied().unwrap_or(num_points);
        if end < start || end > num_points {
            return Err(ShapefileError::InvalidRecord("invalid part range"));
        }
        parts.push(points[start..end].to_vec());
    }
    Ok(parts)
}

fn read_points(r: &mut SliceReader<'_>, n: usize) -> Result<Vec<Coord>, ShapefileError> {
    let mut points = Vec::with_capacity(n);
    for _ in 0..n {
        points.push(Coord::new(r.f64_le()?, r.f64_le()?));
    }
    Ok(points)
}

fn skip_bbox(r: &mut SliceReader<'_>) -> Result<(), ShapefileError> {
    r.skip(4 * size_of::<f64>())
}

fn ring_is_inner(ring: &LineString) -> bool {
    ring.0
        .windows(2)
        .map(|pts| (pts[1].x - pts[0].x) * (pts[1].y + pts[0].y))
        .sum::<f64>()
        / 2.0
        < 0.0
}

fn nonnegative_usize(value: i32, field: &'static str) -> Result<usize, ShapefileError> {
    if value < 0 {
        return Err(ShapefileError::NegativeCount { field, value });
    }
    Ok(value as usize)
}

pub(super) fn words_to_bytes(words: i32, field: &'static str) -> Result<usize, ShapefileError> {
    let words = nonnegative_usize(words, field)?;
    words
        .checked_mul(2)
        .ok_or(ShapefileError::InvalidRecord("byte length overflow"))
}

fn skip_exact<R: Read>(source: &mut R, n: usize) -> Result<(), ShapefileError> {
    let mut buf = [0; 32];
    let mut remaining = n;
    while remaining > 0 {
        let take = remaining.min(buf.len());
        source.read_exact(&mut buf[..take])?;
        remaining -= take;
    }
    Ok(())
}

pub(super) fn read_i32_be<R: Read>(source: &mut R) -> Result<i32, ShapefileError> {
    let mut b = [0; 4];
    source.read_exact(&mut b)?;
    Ok(i32::from_be_bytes(b))
}

fn read_i32_le<R: Read>(source: &mut R) -> Result<i32, ShapefileError> {
    let mut b = [0; 4];
    source.read_exact(&mut b)?;
    Ok(i32::from_le_bytes(b))
}

fn read_f64_le<R: Read>(source: &mut R) -> Result<f64, ShapefileError> {
    let mut b = [0; 8];
    source.read_exact(&mut b)?;
    Ok(f64::from_le_bytes(b))
}

struct SliceReader<'a> {
    bytes: &'a [u8],
    pos: usize,
}

impl<'a> SliceReader<'a> {
    fn new(bytes: &'a [u8]) -> Self {
        Self { bytes, pos: 0 }
    }

    fn take(&mut self, n: usize) -> Result<&'a [u8], ShapefileError> {
        let end = self
            .pos
            .checked_add(n)
            .ok_or(ShapefileError::InvalidRecord("record cursor overflow"))?;
        if end > self.bytes.len() {
            return Err(ShapefileError::RecordTooShort {
                needed: n,
                remaining: self.bytes.len().saturating_sub(self.pos),
            });
        }
        let out = &self.bytes[self.pos..end];
        self.pos = end;
        Ok(out)
    }

    fn skip(&mut self, n: usize) -> Result<(), ShapefileError> {
        self.take(n).map(|_| ())
    }

    fn i32_le(&mut self) -> Result<i32, ShapefileError> {
        let b = self.take(4)?;
        Ok(i32::from_le_bytes([b[0], b[1], b[2], b[3]]))
    }

    fn f64_le(&mut self) -> Result<f64, ShapefileError> {
        let b = self.take(8)?;
        Ok(f64::from_le_bytes([
            b[0], b[1], b[2], b[3], b[4], b[5], b[6], b[7],
        ]))
    }
}
