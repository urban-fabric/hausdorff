//! Per-face icosahedron data: centers, axes, neighbor orientations.
//!
//! Standard H3 data (20 icosahedron faces), lifted from the reference tables.

use super::coord::{CoordIJK, Vec3d};

pub(super) const NUM_ICOSA_FACES: usize = 20;

/// Icosahedron face centers as unit-sphere x/y/z (n-vectors).
#[rustfmt::skip]
pub(super) static CENTER_POINT: [Vec3d; NUM_ICOSA_FACES] = [
    Vec3d::new( 0.2199307791404606,  0.6583691780274996,  0.7198475378926182),
    Vec3d::new(-0.2139234834501421,  0.1478171829550703,  0.9656017935214205),
    Vec3d::new( 0.1092625278784797, -0.481195157287321,   0.8697775121287253),
    Vec3d::new( 0.7428567301586791, -0.3593941678278028,  0.5648005936517033),
    Vec3d::new( 0.8112534709140969,  0.3448953237639384,  0.472138773641393),
    Vec3d::new(-0.1055498149613921,  0.9794457296411413,  0.1718874610009365),
    Vec3d::new(-0.8075407579970092,  0.1533552485898818,  0.5695261994882688),
    Vec3d::new(-0.2846148069787907, -0.8644080972654206,  0.4144792552473539),
    Vec3d::new( 0.7405621473854482, -0.6673299564565524, -0.0789837646326737),
    Vec3d::new( 0.8512303986474293,  0.4722343788582681, -0.2289137388687808),
    Vec3d::new(-0.7405621473854481,  0.6673299564565524,  0.0789837646326737),
    Vec3d::new(-0.8512303986474292, -0.4722343788582682,  0.2289137388687808),
    Vec3d::new( 0.1055498149613919, -0.9794457296411413, -0.1718874610009365),
    Vec3d::new( 0.8075407579970092, -0.1533552485898819, -0.5695261994882688),
    Vec3d::new( 0.2846148069787908,  0.8644080972654204, -0.4144792552473539),
    Vec3d::new(-0.7428567301586791,  0.3593941678278027, -0.5648005936517033),
    Vec3d::new(-0.811253470914097,  -0.3448953237639382, -0.472138773641393),
    Vec3d::new(-0.2199307791404607, -0.6583691780274996, -0.7198475378926182),
    Vec3d::new( 0.213923483450142,  -0.1478171829550704, -0.9656017935214205),
    Vec3d::new(-0.1092625278784796,  0.481195157287321,  -0.8697775121287253),
];

/// Face `ijk` axes as azimuth (radians) from face center to vertex 0.
pub(super) static AXES_AZ_RADS_CII: [f64; NUM_ICOSA_FACES] = [
    5.6199582685239395,
    5.7603390817141875,
    0.78021365439343,
    0.4304693639799999,
    6.130269123335111,
    2.692877706530643,
    2.982963003477244,
    3.532912002790141,
    3.494305004259568,
    3.0032141694995382,
    5.930472956509812,
    0.13837848409025486,
    0.4487149470591504,
    0.15862965011254937,
    5.891865957979238,
    2.711123289609793,
    3.294508837434268,
    3.80481969224544,
    3.6644388790551923,
    2.361378999196363,
];

/// How to transform into an adjacent face's IJK system.
#[derive(Debug, Clone, Copy)]
pub(super) struct FaceOrientIJK {
    /// Adjacent face number.
    pub(super) face: u8,
    /// Resolution-0 translation relative to the primary face.
    pub(super) translate: CoordIJK,
    /// Number of 60° CCW rotations relative to the primary face.
    pub(super) ccw_rot60: u8,
}

/// IJ quadrant index into [`NEIGHBORS`].
pub(super) const IJ: usize = 1;
/// KI quadrant index into [`NEIGHBORS`].
pub(super) const KI: usize = 2;
/// JK quadrant index into [`NEIGHBORS`].
pub(super) const JK: usize = 3;

macro_rules! foijk {
    ($face:literal, ($i:literal, $j:literal, $k:literal), $rot:literal) => {
        FaceOrientIJK {
            face: $face,
            translate: CoordIJK::new($i, $j, $k),
            ccw_rot60: $rot,
        }
    };
}

/// Which faces neighbor each face (central, IJ, KI, JK quadrants).
#[rustfmt::skip]
pub(super) static NEIGHBORS: [[FaceOrientIJK; 4]; NUM_ICOSA_FACES] = [
    [foijk!(0,(0,0,0),0), foijk!(4,(2,0,2),1), foijk!(1,(2,2,0),5), foijk!(5,(0,2,2),3)],
    [foijk!(1,(0,0,0),0), foijk!(0,(2,0,2),1), foijk!(2,(2,2,0),5), foijk!(6,(0,2,2),3)],
    [foijk!(2,(0,0,0),0), foijk!(1,(2,0,2),1), foijk!(3,(2,2,0),5), foijk!(7,(0,2,2),3)],
    [foijk!(3,(0,0,0),0), foijk!(2,(2,0,2),1), foijk!(4,(2,2,0),5), foijk!(8,(0,2,2),3)],
    [foijk!(4,(0,0,0),0), foijk!(3,(2,0,2),1), foijk!(0,(2,2,0),5), foijk!(9,(0,2,2),3)],
    [foijk!(5,(0,0,0),0), foijk!(10,(2,2,0),3), foijk!(14,(2,0,2),3), foijk!(0,(0,2,2),3)],
    [foijk!(6,(0,0,0),0), foijk!(11,(2,2,0),3), foijk!(10,(2,0,2),3), foijk!(1,(0,2,2),3)],
    [foijk!(7,(0,0,0),0), foijk!(12,(2,2,0),3), foijk!(11,(2,0,2),3), foijk!(2,(0,2,2),3)],
    [foijk!(8,(0,0,0),0), foijk!(13,(2,2,0),3), foijk!(12,(2,0,2),3), foijk!(3,(0,2,2),3)],
    [foijk!(9,(0,0,0),0), foijk!(14,(2,2,0),3), foijk!(13,(2,0,2),3), foijk!(4,(0,2,2),3)],
    [foijk!(10,(0,0,0),0), foijk!(5,(2,2,0),3), foijk!(6,(2,0,2),3), foijk!(15,(0,2,2),3)],
    [foijk!(11,(0,0,0),0), foijk!(6,(2,2,0),3), foijk!(7,(2,0,2),3), foijk!(16,(0,2,2),3)],
    [foijk!(12,(0,0,0),0), foijk!(7,(2,2,0),3), foijk!(8,(2,0,2),3), foijk!(17,(0,2,2),3)],
    [foijk!(13,(0,0,0),0), foijk!(8,(2,2,0),3), foijk!(9,(2,0,2),3), foijk!(18,(0,2,2),3)],
    [foijk!(14,(0,0,0),0), foijk!(9,(2,2,0),3), foijk!(5,(2,0,2),3), foijk!(19,(0,2,2),3)],
    [foijk!(15,(0,0,0),0), foijk!(16,(2,0,2),1), foijk!(19,(2,2,0),5), foijk!(10,(0,2,2),3)],
    [foijk!(16,(0,0,0),0), foijk!(17,(2,0,2),1), foijk!(15,(2,2,0),5), foijk!(11,(0,2,2),3)],
    [foijk!(17,(0,0,0),0), foijk!(18,(2,0,2),1), foijk!(16,(2,2,0),5), foijk!(12,(0,2,2),3)],
    [foijk!(18,(0,0,0),0), foijk!(19,(2,0,2),1), foijk!(17,(2,2,0),5), foijk!(13,(0,2,2),3)],
    [foijk!(19,(0,0,0),0), foijk!(15,(2,0,2),1), foijk!(18,(2,2,0),5), foijk!(14,(0,2,2),3)],
];

/// Direction (quadrant) from origin face `i` to destination face `j`, encoded
/// 3 bits per destination: `001`=IJ, `010`=KI, `011`=JK, `000`=central,
/// `111`=invalid.
pub(super) fn adjacent_face_dir(i: u8, j: u8) -> u8 {
    ((ADJACENT_FACE_DIR[i as usize] >> (j as usize * 3)) & 0b111) as u8
}

macro_rules! adj {
    (central = $c:literal, ij = $ij:literal, ki = $ki:literal, jk = $jk:literal) => {
        !(0u64
            | (0b111 << (3 * $c))
            | (0b110 << (3 * $ij))
            | (0b101 << (3 * $ki))
            | (0b100 << (3 * $jk)))
    };
}

#[rustfmt::skip]
static ADJACENT_FACE_DIR: [u64; NUM_ICOSA_FACES] = [
    adj!(central=0,  ij=4,  ki=1,  jk=5),
    adj!(central=1,  ij=0,  ki=2,  jk=6),
    adj!(central=2,  ij=1,  ki=3,  jk=7),
    adj!(central=3,  ij=2,  ki=4,  jk=8),
    adj!(central=4,  ij=3,  ki=0,  jk=9),
    adj!(central=5,  ij=10, ki=14, jk=0),
    adj!(central=6,  ij=11, ki=10, jk=1),
    adj!(central=7,  ij=12, ki=11, jk=2),
    adj!(central=8,  ij=13, ki=12, jk=3),
    adj!(central=9,  ij=14, ki=13, jk=4),
    adj!(central=10, ij=5,  ki=6,  jk=15),
    adj!(central=11, ij=6,  ki=7,  jk=16),
    adj!(central=12, ij=7,  ki=8,  jk=17),
    adj!(central=13, ij=8,  ki=9,  jk=18),
    adj!(central=14, ij=9,  ki=5,  jk=19),
    adj!(central=15, ij=16, ki=19, jk=10),
    adj!(central=16, ij=17, ki=15, jk=11),
    adj!(central=17, ij=18, ki=16, jk=12),
    adj!(central=18, ij=19, ki=17, jk=13),
    adj!(central=19, ij=15, ki=18, jk=14),
];
