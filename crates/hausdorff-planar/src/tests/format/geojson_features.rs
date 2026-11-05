//! `Feature` / `FeatureCollection` round trips, and the spec details parsers get
//! wrong.
//!
//! The geometry codec was already tested; what is new here is the *document*
//! layer, and almost every case below is a thing RFC 7946 permits that a naive
//! reader rejects or silently destroys: a null geometry, a numeric `id`, a foreign
//! member, a `bbox` that disagrees with the geometry it describes.
//!
//! The through-line is **losslessness**. A codec that reads and writes the same
//! format is used to move data between tools, so anything it drops is data a user
//! loses by passing it through — and dropping is silent by nature. So the round
//! trips here assert equality of the whole structure, not just of the geometry.

use crate::format::geojson::{Feature, FeatureCollection, GeoJson, from_geojson_any};
use crate::space::r2::{Coord, Geometry, LineString, Point, Polygon};
use serde_json::{Map, Value, json};

fn point(x: f64, y: f64) -> Geometry {
    Geometry::Point(Point(Coord::new(x, y)))
}

#[test]
fn a_feature_collection_round_trips_losslessly() {
    let text = r#"{
      "type": "FeatureCollection",
      "features": [
        {"type":"Feature","geometry":{"type":"Point","coordinates":[1,2]},
         "properties":{"name":"a","pop":42}},
        {"type":"Feature","geometry":{"type":"LineString","coordinates":[[0,0],[1,1]]},
         "properties":{"name":"b"},"id":7}
      ]
    }"#;
    let fc = FeatureCollection::from_geojson(text).expect("a real-world collection parses");
    assert_eq!(fc.features.len(), 2);
    assert_eq!(fc.features[0].properties["name"], json!("a"));
    assert_eq!(fc.features[0].properties["pop"], json!(42));
    assert_eq!(fc.features[1].id, Some(json!(7)));
    assert_eq!(fc.geometries().count(), 2, "both features carry geometry");

    // Re-parsing our own output must give the same structure back.
    let again = FeatureCollection::from_geojson(&fc.to_geojson()).expect("our output parses");
    assert_eq!(again, fc, "the collection must survive a round trip intact");
}

#[test]
fn a_null_geometry_is_legal_and_survives() {
    // RFC 7946 section 3.2: a feature with attributes and no location is valid.
    // A reader that errors here fails on perfectly ordinary exports.
    let text = r#"{"type":"Feature","geometry":null,"properties":{"name":"nowhere"}}"#;
    let f = Feature::from_geojson(text).expect("a null geometry is legal");
    assert!(f.geometry.is_none());
    assert_eq!(f.properties["name"], json!("nowhere"));

    let again = Feature::from_geojson(&f.to_geojson()).expect("round trips");
    assert_eq!(again, f, "null must stay null rather than becoming empty");
    assert!(
        f.to_geojson().contains("\"geometry\":null"),
        "the member is required, so it is emitted as an explicit null: {}",
        f.to_geojson()
    );
    // And a collection of them has no geometries without that being an error.
    let fc = FeatureCollection::new(vec![f]);
    assert_eq!(fc.geometries().count(), 0);
}

#[test]
fn an_id_may_be_a_string_or_a_number_and_is_not_coerced() {
    // Coercing to String would make 1 and "1" the same feature; coercing to a
    // number would drop UUIDs. So the type is preserved as written.
    let numeric =
        Feature::from_geojson(r#"{"type":"Feature","geometry":null,"properties":{},"id":1}"#)
            .expect("numeric id");
    let string =
        Feature::from_geojson(r#"{"type":"Feature","geometry":null,"properties":{},"id":"1"}"#)
            .expect("string id");
    assert_eq!(numeric.id, Some(json!(1)));
    assert_eq!(string.id, Some(json!("1")));
    assert_ne!(
        numeric.id, string.id,
        "1 and \"1\" are different ids and must stay different"
    );
    for f in [&numeric, &string] {
        assert_eq!(
            &Feature::from_geojson(&f.to_geojson()).expect("round trips"),
            f,
            "the id type must survive a round trip"
        );
    }
    // A non-scalar id is not legal.
    assert!(
        Feature::from_geojson(r#"{"type":"Feature","geometry":null,"properties":{},"id":[1]}"#)
            .is_err(),
        "an array id is not a string or a number"
    );
}

#[test]
fn foreign_members_are_preserved_rather_than_dropped() {
    // Section 6.1 permits them, and tools attach CRS hints, styling, and vendor
    // metadata. Dropping unknown keys silently destroys a user's data.
    let text = r#"{
      "type":"FeatureCollection",
      "name":"my layer",
      "crs":{"type":"name","properties":{"name":"EPSG:4326"}},
      "features":[
        {"type":"Feature","geometry":{"type":"Point","coordinates":[1,2]},
         "properties":{},"style":{"color":"red"}}
      ]
    }"#;
    let fc = FeatureCollection::from_geojson(text).expect("foreign members are legal");
    assert_eq!(fc.foreign["name"], json!("my layer"));
    assert!(fc.foreign.contains_key("crs"), "the crs hint must be kept");
    assert_eq!(fc.features[0].foreign["style"], json!({"color":"red"}));
    // Known members must not leak into `foreign`.
    for key in ["type", "features", "bbox"] {
        assert!(
            !fc.foreign.contains_key(key),
            "{key} is a known member and must not be duplicated into foreign"
        );
    }
    let again = FeatureCollection::from_geojson(&fc.to_geojson()).expect("round trips");
    assert_eq!(again, fc, "foreign members must survive a round trip");
}

#[test]
fn bbox_is_carried_verbatim_and_never_recomputed() {
    // `bbox` is advisory, may disagree with the geometry, and is 4 or 6 numbers.
    // Recomputing it on write would emit a different document than the one read;
    // validating it would reject files other tools accept. So it is carried.
    let wrong = r#"{"type":"Feature","geometry":{"type":"Point","coordinates":[1,2]},
                    "properties":{},"bbox":[-180,-90,180,90]}"#;
    let f = Feature::from_geojson(wrong).expect("a disagreeing bbox is still legal");
    assert_eq!(f.bbox, Some(vec![-180.0, -90.0, 180.0, 90.0]));
    assert!(
        f.to_geojson().contains("[-180.0,-90.0,180.0,90.0]")
            || f.to_geojson().contains("[-180,-90,180,90]"),
        "the bbox must be written back as given: {}",
        f.to_geojson()
    );
    // Six numbers, for a 3-D bbox, is equally legal even though the geometry
    // model is 2-D.
    let three_d = r#"{"type":"Feature","geometry":null,"properties":{},
                      "bbox":[0,0,0,1,1,1]}"#;
    assert_eq!(
        Feature::from_geojson(three_d)
            .expect("a 6-element bbox parses")
            .bbox,
        Some(vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0])
    );
    // A malformed bbox is an error rather than a silent None: it means the
    // document is not what it claims to be.
    for bad in [
        r#"{"type":"Feature","geometry":null,"properties":{},"bbox":["a","b"]}"#,
        r#"{"type":"Feature","geometry":null,"properties":{},"bbox":3}"#,
        r#"{"type":"Feature","geometry":null,"properties":{},"bbox":{"a":1}}"#,
    ] {
        assert!(Feature::from_geojson(bad).is_err(), "malformed bbox: {bad}");
    }
    // Absent and explicitly null both read as None.
    for ok in [
        r#"{"type":"Feature","geometry":null,"properties":{}}"#,
        r#"{"type":"Feature","geometry":null,"properties":{},"bbox":null}"#,
    ] {
        assert_eq!(Feature::from_geojson(ok).expect("parses").bbox, None);
    }
}

#[test]
fn absent_or_null_properties_read_as_empty() {
    // The member is required by the spec but is routinely absent or null in real
    // files, and an empty object is the only useful reading.
    for text in [
        r#"{"type":"Feature","geometry":null,"properties":{}}"#,
        r#"{"type":"Feature","geometry":null,"properties":null}"#,
        r#"{"type":"Feature","geometry":null}"#,
    ] {
        let f = Feature::from_geojson(text).expect("parses");
        assert!(f.properties.is_empty(), "{text}");
    }
    // But a non-object properties is a genuine error.
    assert!(Feature::from_geojson(r#"{"type":"Feature","geometry":null,"properties":3}"#).is_err());
}

#[test]
fn from_geojson_any_dispatches_on_the_type_member() {
    // The entry point for a file boundary, where the document kind is not known
    // until it has been read.
    /// A document alongside the test for which `GeoJson` variant it must become.
    type Case = (&'static str, fn(&GeoJson) -> bool);
    let cases: Vec<Case> = vec![
        (r#"{"type":"Point","coordinates":[1,2]}"#, |g| {
            matches!(g, GeoJson::Geometry(_))
        }),
        (
            r#"{"type":"GeometryCollection","geometries":[{"type":"Point","coordinates":[0,0]}]}"#,
            |g| matches!(g, GeoJson::Geometry(_)),
        ),
        (
            r#"{"type":"Feature","geometry":null,"properties":{}}"#,
            |g| matches!(g, GeoJson::Feature(_)),
        ),
        (r#"{"type":"FeatureCollection","features":[]}"#, |g| {
            matches!(g, GeoJson::FeatureCollection(_))
        }),
    ];
    for (text, is_right_kind) in cases {
        let parsed = from_geojson_any(text).unwrap_or_else(|e| panic!("{text}: {e:?}"));
        assert!(is_right_kind(&parsed), "wrong kind for {text}");
        // And whatever the kind, re-parsing the output gives the same thing.
        assert_eq!(
            from_geojson_any(&parsed.to_geojson()).expect("round trips"),
            parsed
        );
    }
    // `geometries()` works without matching on the kind.
    assert_eq!(
        from_geojson_any(r#"{"type":"Point","coordinates":[1,2]}"#)
            .expect("parses")
            .geometries()
            .len(),
        1
    );
    assert_eq!(
        from_geojson_any(r#"{"type":"Feature","geometry":null,"properties":{}}"#)
            .expect("parses")
            .geometries()
            .len(),
        0,
        "a null geometry contributes nothing"
    );
}

#[test]
fn malformed_documents_are_rejected_not_panicked() {
    for bad in [
        "",
        "{}",
        "null",
        "[]",
        r#"{"type":"Nonsense"}"#,
        r#"{"type":3}"#,
        r#"{"features":[]}"#,
        r#"{"type":"FeatureCollection"}"#,
        r#"{"type":"FeatureCollection","features":{}}"#,
        r#"{"type":"FeatureCollection","features":[{"type":"Point","coordinates":[0,0]}]}"#,
        r#"{"type":"Feature","geometry":{"type":"Point"},"properties":{}}"#,
        r#"{"type":"Feature","geometry":{"type":"Bogus","coordinates":[0,0]},"properties":{}}"#,
        r#"{"type":"Point","coordinates":[1,2]"#,
    ] {
        let any = from_geojson_any(bad);
        let feat = Feature::from_geojson(bad);
        let coll = FeatureCollection::from_geojson(bad);
        assert!(
            any.is_err() || feat.is_err() || coll.is_err(),
            "{bad:?} should not parse as every kind at once"
        );
    }
    // A collection whose members are geometries rather than features is a common
    // malformation and must be refused rather than half-read.
    assert!(
        FeatureCollection::from_geojson(
            r#"{"type":"FeatureCollection","features":[{"type":"Point","coordinates":[0,0]}]}"#
        )
        .is_err(),
        "a geometry is not a Feature"
    );
}

#[test]
fn an_empty_collection_is_legal_and_round_trips() {
    let fc = FeatureCollection::from_geojson(r#"{"type":"FeatureCollection","features":[]}"#)
        .expect("an empty collection is legal");
    assert!(fc.features.is_empty());
    assert_eq!(fc.geometries().count(), 0);
    assert_eq!(
        FeatureCollection::from_geojson(&fc.to_geojson()).expect("round trips"),
        fc
    );
}

#[test]
fn constructed_features_serialize_to_parseable_geojson() {
    // The writer must produce what the reader accepts — a writer and parser that
    // disagree is how a format quietly becomes a private dialect.
    let mut properties = Map::new();
    properties.insert("name".into(), Value::from("somewhere"));
    properties.insert("nested".into(), json!({"a":[1,2,3]}));
    let feature = Feature {
        geometry: Some(Geometry::Polygon(Polygon::new(
            LineString::from(vec![(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 0.0)]),
            vec![],
        ))),
        properties,
        id: Some(Value::from("abc")),
        bbox: Some(vec![0.0, 0.0, 1.0, 1.0]),
        foreign: Map::new(),
    };
    let fc = FeatureCollection::new(vec![
        feature.clone(),
        Feature::from_geometry(point(5.0, 6.0)),
    ]);
    let round = FeatureCollection::from_geojson(&fc.to_geojson()).expect("our writer is parseable");
    assert_eq!(round, fc);
    assert_eq!(round.geometries().count(), 2);
    assert_eq!(round.features[0].properties["nested"], json!({"a":[1,2,3]}));
}
