//! [`H3Cell`] edges and vertices: the directed-edge constructors
//! ([`edge`](H3Cell::edge) / [`edges`](H3Cell::edges)) and the topological
//! vertex machinery ([`vertex`](H3Cell::vertex) / [`vertexes`](H3Cell::vertexes)
//! and the direction↔vertex rotation helpers behind them).

use super::*;

const NUM_HEX_VERTS: u8 = 6;
const NUM_PENT_VERTS: u8 = 5;
const DIRECTIONS_CCW: [u8; 6] = [2, 3, 1, 5, 4, 6];
const REV_NEIGHBOR_DIRECTIONS_HEX: [u8; 6] = [5, 3, 4, 1, 0, 2];
const VERTEX_TO_DIRECTION_HEX: [u8; 6] = [6, 2, 3, 1, 5, 4];
const VERTEX_TO_DIRECTION_PENT: [u8; 5] = [6, 2, 3, 5, 4];
const DIRECTION_TO_VERTEX_HEX: [u8; 6] = [3, 1, 2, 5, 4, 0];
const DIRECTION_TO_VERTEX_PENT: [u8; 5] = [1, 2, 4, 3, 0];
const HEXAGON_EDGE_ORDER: [u8; 6] = [1, 5, 4, 6, 2, 3];
const PENTAGON_EDGE_ORDER: [u8; 5] = [2, 3, 5, 4, 6];

impl H3Cell {
    /// Directed edge from this cell to `destination`, if the cells are
    /// neighbors.
    pub fn edge(self, destination: Self) -> Option<H3DirectedEdge> {
        let direction = grid::direction_for_neighbor(self.0, destination.0)?;
        H3DirectedEdge::from_origin_and_direction(self, direction)
    }

    /// All directed edges originating at this cell.
    pub fn edges(self) -> Vec<H3DirectedEdge> {
        let directions: &[u8] = if self.is_pentagon() {
            &PENTAGON_EDGE_ORDER
        } else {
            &HEXAGON_EDGE_ORDER
        };
        directions
            .iter()
            .copied()
            .filter_map(|direction| H3DirectedEdge::from_origin_and_direction(self, direction))
            .collect()
    }

    /// Number of topological vertices on this cell.
    pub fn vertex_count(self) -> u8 {
        if self.is_pentagon() {
            NUM_PENT_VERTS
        } else {
            NUM_HEX_VERTS
        }
    }

    /// Canonical vertex index for `vertex` (`0..5` for hexagons, `0..4` for
    /// pentagons).
    pub fn vertex(self, vertex: u8) -> Option<H3Vertex> {
        let vertex_count = self.vertex_count();
        if vertex >= vertex_count {
            return None;
        }

        let mut owner = self;
        let mut owner_vertex = vertex;

        if self.resolution() == 0 || get_direction(self.0, self.resolution()) != 0 {
            let left = self.direction_for_vertex(vertex)?;
            let (left_bits, left_rotation) = grid::neighbor_rotations(self.0, left, 0)?;
            let left_cell = Self(left_bits);
            if left_cell < owner {
                owner = left_cell;
            }

            if self.resolution() == 0 || get_direction(left_bits, self.resolution()) != 0 {
                let right_vertex = (vertex + vertex_count - 1) % vertex_count;
                let right = self.direction_for_vertex(right_vertex)?;
                let (right_bits, right_rotation) = grid::neighbor_rotations(self.0, right, 0)?;
                let right_cell = Self(right_bits);

                if right_cell < owner {
                    owner = right_cell;
                    let direction = if owner.is_pentagon() {
                        grid::direction_for_neighbor(owner.0, self.0)?
                    } else {
                        let offset = (REV_NEIGHBOR_DIRECTIONS_HEX[usize::from(right) - 1]
                            + right_rotation)
                            % NUM_HEX_VERTS;
                        DIRECTIONS_CCW[usize::from(offset)]
                    };
                    owner_vertex = owner.vertex_for_direction(direction)?;
                }
            }

            if owner == left_cell {
                let direction = if owner.is_pentagon() {
                    grid::direction_for_neighbor(owner.0, self.0)?
                } else {
                    let offset = (REV_NEIGHBOR_DIRECTIONS_HEX[usize::from(left) - 1]
                        + left_rotation)
                        % NUM_HEX_VERTS;
                    DIRECTIONS_CCW[usize::from(offset)]
                };
                owner_vertex = (owner.vertex_for_direction(direction)? + 1) % owner.vertex_count();
            }
        }

        Some(H3Vertex(set_reserved_bits(
            set_mode(owner.0, VERTEX_MODE),
            owner_vertex,
        )))
    }

    /// All canonical vertex indexes for this cell.
    pub fn vertexes(self) -> Vec<H3Vertex> {
        (0..self.vertex_count())
            .filter_map(|vertex| self.vertex(vertex))
            .collect()
    }

    fn vertex_rotations(self) -> u8 {
        let fijk =
            faceijk::FaceIJK::from_cell(self.0, self.base_cell(), self.resolution() as usize);
        let mut rotations =
            faceijk::FaceIJK::base_cell_rotation_count(self.base_cell(), fijk.face).unwrap_or(0);

        if base_cells::is_pentagon(self.base_cell()) {
            let direction_faces = base_cells::pentagon_direction_faces(self.base_cell());
            let home_face = base_cells::home_face(self.base_cell());

            if fijk.face != home_face
                && (base_cells::is_polar_pentagon(self.base_cell())
                    || fijk.face == direction_faces[3])
            {
                rotations = (rotations + 1) % NUM_HEX_VERTS;
            }

            if first_axe(self.0) == Some(3) && fijk.face == direction_faces[3] {
                rotations = (rotations + 5) % NUM_HEX_VERTS;
            } else if first_axe(self.0) == Some(5) && fijk.face == direction_faces[1] {
                rotations = (rotations + 1) % NUM_HEX_VERTS;
            }
        }

        rotations
    }

    pub(super) fn direction_for_vertex(self, vertex: u8) -> Option<u8> {
        if vertex >= self.vertex_count() {
            return None;
        }
        let rotations = self.vertex_rotations();
        if self.is_pentagon() {
            let index = (vertex + rotations) % NUM_PENT_VERTS;
            Some(VERTEX_TO_DIRECTION_PENT[usize::from(index)])
        } else {
            let index = (vertex + rotations) % NUM_HEX_VERTS;
            Some(VERTEX_TO_DIRECTION_HEX[usize::from(index)])
        }
    }

    pub(super) fn vertex_for_direction(self, direction: u8) -> Option<u8> {
        if direction == 0 || direction > 6 || (self.is_pentagon() && direction == K_AXIS_DIGIT) {
            return None;
        }
        let rotations = self.vertex_rotations();
        if self.is_pentagon() {
            let vertex = DIRECTION_TO_VERTEX_PENT[usize::from(direction - 2)];
            Some((vertex + NUM_PENT_VERTS - rotations) % NUM_PENT_VERTS)
        } else {
            let vertex = DIRECTION_TO_VERTEX_HEX[usize::from(direction - 1)];
            Some((vertex + NUM_HEX_VERTS - rotations) % NUM_HEX_VERTS)
        }
    }
}
