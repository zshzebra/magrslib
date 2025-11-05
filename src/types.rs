//! Core types for magrslib
//!
//! This module defines the fundamental types used throughout the library:
//! - Vec3: 3D vectors
//! - Position: Single position in 3D space
//! - Path: Sequence of positions (for animation/time evolution)
//! - Orientation: 3D rotation using quaternions
//! - Observers: Collection of observer positions

use nalgebra::{Unit, UnitQuaternion, Vector3};
use ndarray::{Array1, Array2};

/// Efficient 3D vector type
pub type Vec3 = [f64; 3];

/// Single position in 3D space (in meters)
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Position {
    pub coords: Vec3,
}

impl Position {
    /// Create a new position
    pub fn new(x: f64, y: f64, z: f64) -> Self {
        Self { coords: [x, y, z] }
    }

    /// Create position at origin
    pub fn origin() -> Self {
        Self {
            coords: [0.0, 0.0, 0.0],
        }
    }

    /// Get x coordinate
    pub fn x(&self) -> f64 {
        self.coords[0]
    }

    /// Get y coordinate
    pub fn y(&self) -> f64 {
        self.coords[1]
    }

    /// Get z coordinate
    pub fn z(&self) -> f64 {
        self.coords[2]
    }

    /// Convert to nalgebra Vector3
    pub fn to_vector3(&self) -> Vector3<f64> {
        Vector3::new(self.coords[0], self.coords[1], self.coords[2])
    }
}

impl From<Vec3> for Position {
    fn from(coords: Vec3) -> Self {
        Self { coords }
    }
}

impl From<(f64, f64, f64)> for Position {
    fn from((x, y, z): (f64, f64, f64)) -> Self {
        Self::new(x, y, z)
    }
}

/// Path representing a sequence of positions (for animation/time evolution)
///
/// In magpylib, sources and observers can have time-varying positions.
/// This is represented as a "path" - a sequence of positions.
#[derive(Debug, Clone, PartialEq)]
pub struct Path {
    positions: Vec<Position>,
}

impl Path {
    /// Create a new path from a vector of positions
    pub fn new(positions: Vec<Position>) -> Self {
        assert!(!positions.is_empty(), "Path must have at least one position");
        Self { positions }
    }

    /// Create a path with a single position
    pub fn from_position(pos: Position) -> Self {
        Self {
            positions: vec![pos],
        }
    }

    /// Get the number of positions in the path
    pub fn len(&self) -> usize {
        self.positions.len()
    }

    /// Check if path is empty (should never be true due to invariant)
    pub fn is_empty(&self) -> bool {
        self.positions.is_empty()
    }

    /// Get a position at a specific index
    pub fn get(&self, index: usize) -> Option<&Position> {
        self.positions.get(index)
    }

    /// Iterate over positions
    pub fn iter(&self) -> impl Iterator<Item = &Position> {
        self.positions.iter()
    }

    /// Get reference to all positions
    pub fn positions(&self) -> &[Position] {
        &self.positions
    }

    /// Convert to ndarray Array2 with shape (path_length, 3)
    pub fn to_array2(&self) -> Array2<f64> {
        let mut arr = Array2::zeros((self.len(), 3));
        for (i, pos) in self.positions.iter().enumerate() {
            arr[[i, 0]] = pos.coords[0];
            arr[[i, 1]] = pos.coords[1];
            arr[[i, 2]] = pos.coords[2];
        }
        arr
    }
}

impl From<Position> for Path {
    fn from(pos: Position) -> Self {
        Self::from_position(pos)
    }
}

impl From<Vec<Position>> for Path {
    fn from(positions: Vec<Position>) -> Self {
        Self::new(positions)
    }
}

impl From<Vec<Vec3>> for Path {
    fn from(coords: Vec<Vec3>) -> Self {
        Self::new(coords.into_iter().map(Position::from).collect())
    }
}

/// 3D orientation using unit quaternions
///
/// Provides a robust representation of rotations without gimbal lock.
/// Internally uses nalgebra's UnitQuaternion.
#[derive(Debug, Clone, PartialEq)]
pub struct Orientation {
    pub(crate) quat: UnitQuaternion<f64>,
}

impl Orientation {
    /// Create identity orientation (no rotation)
    pub fn identity() -> Self {
        Self {
            quat: UnitQuaternion::identity(),
        }
    }

    /// Create orientation from unit quaternion
    pub fn from_quaternion(quat: UnitQuaternion<f64>) -> Self {
        Self { quat }
    }

    /// Create orientation from axis-angle representation
    ///
    /// # Arguments
    /// * `axis` - Rotation axis (will be normalized)
    /// * `angle` - Rotation angle in radians
    pub fn from_axis_angle(axis: Vec3, angle: f64) -> Self {
        let axis_vec = Vector3::new(axis[0], axis[1], axis[2]);
        let axis_unit = Unit::new_normalize(axis_vec);
        Self {
            quat: UnitQuaternion::from_axis_angle(&axis_unit, angle),
        }
    }

    /// Create orientation from Euler angles (ZYX convention)
    ///
    /// # Arguments
    /// * `roll` - Rotation around x-axis (radians)
    /// * `pitch` - Rotation around y-axis (radians)
    /// * `yaw` - Rotation around z-axis (radians)
    pub fn from_euler_angles(roll: f64, pitch: f64, yaw: f64) -> Self {
        Self {
            quat: UnitQuaternion::from_euler_angles(roll, pitch, yaw),
        }
    }

    /// Create orientation from rotation around x-axis
    pub fn from_x_rotation(angle: f64) -> Self {
        Self::from_axis_angle([1.0, 0.0, 0.0], angle)
    }

    /// Create orientation from rotation around y-axis
    pub fn from_y_rotation(angle: f64) -> Self {
        Self::from_axis_angle([0.0, 1.0, 0.0], angle)
    }

    /// Create orientation from rotation around z-axis
    pub fn from_z_rotation(angle: f64) -> Self {
        Self::from_axis_angle([0.0, 0.0, 1.0], angle)
    }

    /// Get the inverse orientation
    pub fn inverse(&self) -> Self {
        Self {
            quat: self.quat.inverse(),
        }
    }

    /// Rotate a 3D vector
    pub fn rotate_vector(&self, v: Vec3) -> Vec3 {
        let vec = Vector3::new(v[0], v[1], v[2]);
        let rotated = self.quat * vec;
        [rotated.x, rotated.y, rotated.z]
    }

    /// Rotate a 3D vector using inverse orientation
    pub fn rotate_vector_inverse(&self, v: Vec3) -> Vec3 {
        let vec = Vector3::new(v[0], v[1], v[2]);
        let rotated = self.quat.inverse() * vec;
        [rotated.x, rotated.y, rotated.z]
    }

    /// Get access to internal quaternion
    pub fn quaternion(&self) -> &UnitQuaternion<f64> {
        &self.quat
    }
}

impl Default for Orientation {
    fn default() -> Self {
        Self::identity()
    }
}

/// Collection of observer positions where fields are evaluated
///
/// This is a lightweight wrapper around ndarray Array2 with shape (n_observers, 3).
/// In the Python version, observers can be:
/// - Sensor objects
/// - Collections
/// - Raw position arrays
///
/// For the minimal port, we'll just use position arrays directly.
#[derive(Debug, Clone)]
pub struct Observers {
    pub(crate) positions: Array2<f64>,
}

impl Observers {
    /// Create observers from an array with shape (n, 3)
    pub fn from_array(positions: Array2<f64>) -> Self {
        assert_eq!(
            positions.ncols(),
            3,
            "Observers array must have 3 columns (x, y, z)"
        );
        Self { positions }
    }

    /// Create observers from a vector of positions
    pub fn from_positions(positions: Vec<Position>) -> Self {
        let n = positions.len();
        let mut arr = Array2::zeros((n, 3));
        for (i, pos) in positions.iter().enumerate() {
            arr[[i, 0]] = pos.coords[0];
            arr[[i, 1]] = pos.coords[1];
            arr[[i, 2]] = pos.coords[2];
        }
        Self { positions: arr }
    }

    /// Create observers from a single position
    pub fn from_position(pos: Position) -> Self {
        Self::from_positions(vec![pos])
    }

    /// Create observers from Vec<Vec3>
    pub fn from_vec3s(coords: Vec<Vec3>) -> Self {
        Self::from_positions(coords.into_iter().map(Position::from).collect())
    }

    /// Get the number of observers
    pub fn len(&self) -> usize {
        self.positions.nrows()
    }

    /// Check if there are no observers
    pub fn is_empty(&self) -> bool {
        self.positions.nrows() == 0
    }

    /// Get a column view (x=0, y=1, z=2)
    pub fn column(&self, index: usize) -> Array1<f64> {
        self.positions.column(index).to_owned()
    }

    /// Get reference to underlying array
    pub fn as_array(&self) -> &Array2<f64> {
        &self.positions
    }

    /// Create a 2D grid of observers in a specified plane
    ///
    /// This is useful for visualization, similar to numpy's mgrid.
    ///
    /// # Arguments
    /// * `x_range` - (start, end, num_points) for first axis
    /// * `y_range` - (start, end, num_points) for second axis
    /// * `plane` - Which plane to use: "xy", "xz", or "yz"
    /// * `fixed_coord` - Value for the fixed coordinate (z for xy, y for xz, x for yz)
    ///
    /// # Returns
    /// Observers arranged in a grid, ordered row-major (iterating y, then x)
    ///
    /// # Example
    /// ```
    /// use magrslib::types::Observers;
    ///
    /// // Create 10x10 grid in xy-plane at z=0, from -0.05 to 0.05 meters
    /// let observers = Observers::from_grid_2d(
    ///     (-0.05, 0.05, 10),
    ///     (-0.05, 0.05, 10),
    ///     "xy",
    ///     0.0
    /// );
    /// assert_eq!(observers.len(), 100); // 10 * 10
    /// ```
    pub fn from_grid_2d(
        x_range: (f64, f64, usize),
        y_range: (f64, f64, usize),
        plane: &str,
        fixed_coord: f64,
    ) -> Self {
        let (x_start, x_end, nx) = x_range;
        let (y_start, y_end, ny) = y_range;

        let mut positions = Vec::with_capacity(nx * ny);

        // Generate grid points
        for j in 0..ny {
            let y = if ny == 1 {
                y_start
            } else {
                y_start + (y_end - y_start) * (j as f64) / ((ny - 1) as f64)
            };

            for i in 0..nx {
                let x = if nx == 1 {
                    x_start
                } else {
                    x_start + (x_end - x_start) * (i as f64) / ((nx - 1) as f64)
                };

                // Map to appropriate 3D coordinates based on plane
                let pos = match plane {
                    "xy" => Position::new(x, y, fixed_coord),
                    "xz" => Position::new(x, fixed_coord, y),
                    "yz" => Position::new(fixed_coord, x, y),
                    _ => panic!("Invalid plane: {}. Use 'xy', 'xz', or 'yz'", plane),
                };
                positions.push(pos);
            }
        }

        Self::from_positions(positions)
    }

    /// Get grid dimensions if this was created from a grid
    ///
    /// This is a helper for visualization - it returns None if the number
    /// of observers doesn't match a perfect grid.
    pub fn infer_grid_shape(&self) -> Option<(usize, usize)> {
        let n = self.len();
        // Try to find factors
        for nx in 1..=n {
            if n % nx == 0 {
                let ny = n / nx;
                if nx * ny == n {
                    return Some((nx, ny));
                }
            }
        }
        None
    }

    /// Create a 3D grid of observers
    ///
    /// This generates a structured 3D grid of observation points, useful for
    /// volume visualization and 3D field analysis.
    ///
    /// # Arguments
    /// * `x_range` - (start, end, num_points) for x-axis
    /// * `y_range` - (start, end, num_points) for y-axis
    /// * `z_range` - (start, end, num_points) for z-axis
    ///
    /// # Returns
    /// Observers arranged in a 3D grid, ordered: iterate z, then y, then x (VTK standard)
    ///
    /// # Example
    /// ```
    /// use magrslib::types::Observers;
    ///
    /// // Create 10x10x10 grid from -0.05 to 0.05 meters in all dimensions
    /// let observers = Observers::from_grid_3d(
    ///     (-0.05, 0.05, 10),
    ///     (-0.05, 0.05, 10),
    ///     (-0.05, 0.05, 10)
    /// );
    /// assert_eq!(observers.len(), 1000); // 10 * 10 * 10
    /// ```
    pub fn from_grid_3d(
        x_range: (f64, f64, usize),
        y_range: (f64, f64, usize),
        z_range: (f64, f64, usize),
    ) -> Self {
        let (x_start, x_end, nx) = x_range;
        let (y_start, y_end, ny) = y_range;
        let (z_start, z_end, nz) = z_range;

        let mut positions = Vec::with_capacity(nx * ny * nz);

        // VTK ordering: iterate z (slowest), then y, then x (fastest)
        // This matches PyVista ImageData expectations
        for k in 0..nz {
            let z = if nz == 1 {
                z_start
            } else {
                z_start + (z_end - z_start) * (k as f64) / ((nz - 1) as f64)
            };

            for j in 0..ny {
                let y = if ny == 1 {
                    y_start
                } else {
                    y_start + (y_end - y_start) * (j as f64) / ((ny - 1) as f64)
                };

                for i in 0..nx {
                    let x = if nx == 1 {
                        x_start
                    } else {
                        x_start + (x_end - x_start) * (i as f64) / ((nx - 1) as f64)
                    };

                    positions.push(Position::new(x, y, z));
                }
            }
        }

        Self::from_positions(positions)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_position() {
        let pos = Position::new(1.0, 2.0, 3.0);
        assert_eq!(pos.x(), 1.0);
        assert_eq!(pos.y(), 2.0);
        assert_eq!(pos.z(), 3.0);
    }

    #[test]
    fn test_path() {
        let positions = vec![
            Position::new(0.0, 0.0, 0.0),
            Position::new(1.0, 0.0, 0.0),
            Position::new(2.0, 0.0, 0.0),
        ];
        let path = Path::new(positions);
        assert_eq!(path.len(), 3);

        let arr = path.to_array2();
        assert_eq!(arr.shape(), &[3, 3]);
        assert_relative_eq!(arr[[1, 0]], 1.0);
    }

    #[test]
    fn test_orientation() {
        // Identity rotation
        let orient = Orientation::identity();
        let v = [1.0, 0.0, 0.0];
        let rotated = orient.rotate_vector(v);
        assert_relative_eq!(rotated[0], 1.0);
        assert_relative_eq!(rotated[1], 0.0);

        // 90 degree rotation around z-axis
        let orient = Orientation::from_z_rotation(std::f64::consts::PI / 2.0);
        let v = [1.0, 0.0, 0.0];
        let rotated = orient.rotate_vector(v);
        assert_relative_eq!(rotated[0], 0.0, epsilon = 1e-10);
        assert_relative_eq!(rotated[1], 1.0, epsilon = 1e-10);
        assert_relative_eq!(rotated[2], 0.0, epsilon = 1e-10);
    }

    #[test]
    fn test_observers() {
        let observers = Observers::from_vec3s(vec![[0.0, 0.0, 0.0], [1.0, 1.0, 1.0]]);
        assert_eq!(observers.len(), 2);

        let x_col = observers.column(0);
        assert_relative_eq!(x_col[0], 0.0);
        assert_relative_eq!(x_col[1], 1.0);
    }
}
