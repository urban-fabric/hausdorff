//! The tuning config: bundled defaults parse, partial overrides fill from
//! defaults, and unknown keys are rejected.

use crate::config::{BundledConfig, Config, config};

#[test]
fn bundled_defaults_parse() {
    let c = Config::bundled().expect("bundled config must parse");
    assert_eq!(c.cost_model.min_distance.index_build_weight, 6.0);
    assert_eq!(c.cost_model.min_distance.index_query_weight, 6.0);
    assert_eq!(c.index.node_capacity, 16);
    assert_eq!(config().index.node_capacity, 16);
    assert_eq!(config().cost_model.min_distance.index_build_weight, 6.0);
}

#[test]
fn partial_override_fills_remaining_from_defaults() {
    let c = Config::from_toml_str("[cost_model.min_distance]\nindex_build_weight = 99.0\n")
        .expect("override parses");
    assert_eq!(c.cost_model.min_distance.index_build_weight, 99.0);
    assert_eq!(c.cost_model.min_distance.index_query_weight, 6.0);
    assert_eq!(c.index.node_capacity, 16);
}

#[test]
fn unknown_field_is_rejected() {
    assert!(Config::from_toml_str("[index]\nbogus = 1\n").is_err());
}
