# hausdorff

Internal geospatial index implementations for UrbanFabric.

The first implementation is `hausdorff::h3`, a from-scratch Rust H3 cell
implementation moved out of Kolmogorov so other UrbanFabric projects can use
the same internal H3 core.

Current H3 coverage includes cell construction/validation, lat/lon indexing,
cell centers and boundaries, hierarchy operations, child positions,
compaction/uncompaction, pentagon/base-cell inspection, grid disk/ring
traversal, neighbor checks, grid distance, local IJ coordinates, grid paths,
directed edges and edge geometry, canonical vertex indexes, exact cell area,
centroid polygon cover, cell-set boundary dissolve, and the standard H3
resolution metadata tables.

The `geo` module also owns reusable geodesic sampling helpers such as
`RadialGrid`, a concentric ring sampler around a WGS84 point.
