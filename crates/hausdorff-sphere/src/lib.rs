//! Unit-sphere geometry — the golang/geo S2 port (see the repo `NOTICE`):
//! cell ids and the cell hierarchy, caps, loops, polylines, edge crossings,
//! orientation predicates, point/edge measures, and region coverings, under
//! [`space::s2`]. The `s1` angle and `r3` vector substrates come from
//! `hausdorff-base`; the lat/lon boundary (`S2Point ↔ GeoPoint`) uses
//! `hausdorff-base`'s `earth` core.

pub mod space;
