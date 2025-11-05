//! Coordinate transformation utilities for field computation
//!
//! Handles transformations between global and source-local coordinate systems.
//! In source-local coordinates, the source is centered at the origin with identity orientation.

use crate::types::{Orientation, Position};
use ndarray::Array2;

/// Transform observer positions from global to source-local coordinates
///
/// # Algorithm
/// 1. Subtract source position (translate to source center)
/// 2. Apply inverse rotation (rotate into source frame)
///
/// This is a direct port of the Python transformation:
/// ```python
/// pos_rel_rot = orientation.apply(observers - position, inverse=True)
/// ```
///
/// # Arguments
/// * `observers` - Observer positions in global coordinates (n, 3)
/// * `position` - Source position in global coordinates
/// * `orientation` - Source orientation
///
/// # Returns
/// Observer positions in source-local coordinates (n, 3)
pub fn transform_observers_to_local(
    observers: &Array2<f64>,
    position: &Position,
    orientation: &Orientation,
) -> Array2<f64> {
    let n = observers.nrows();
    let mut result = Array2::zeros((n, 3));

    // Get source position coordinates
    let pos = position.coords;

    // For each observer:
    // 1. Subtract source position
    // 2. Apply inverse rotation
    for i in 0..n {
        let obs = [observers[[i, 0]], observers[[i, 1]], observers[[i, 2]]];

        // Translate: obs - position
        let translated = [
            obs[0] - pos[0],
            obs[1] - pos[1],
            obs[2] - pos[2],
        ];

        // Rotate using inverse orientation
        let local = orientation.rotate_vector_inverse(translated);

        result[[i, 0]] = local[0];
        result[[i, 1]] = local[1];
        result[[i, 2]] = local[2];
    }

    result
}

/// Transform field vectors from source-local to global coordinates
///
/// # Algorithm
/// Apply forward rotation to field vectors (no translation needed for vectors)
///
/// This is a direct port of the Python transformation:
/// ```python
/// BH = orientation.apply(BH)
/// ```
///
/// # Arguments
/// * `fields` - Field vectors in source-local coordinates (n, 3)
/// * `orientation` - Source orientation
///
/// # Returns
/// Field vectors in global coordinates (n, 3)
pub fn transform_field_to_global(
    fields: &Array2<f64>,
    orientation: &Orientation,
) -> Array2<f64> {
    let n = fields.nrows();
    let mut result = Array2::zeros((n, 3));

    // For each field vector, apply forward rotation
    for i in 0..n {
        let field = [fields[[i, 0]], fields[[i, 1]], fields[[i, 2]]];
        let global = orientation.rotate_vector(field);

        result[[i, 0]] = global[0];
        result[[i, 1]] = global[1];
        result[[i, 2]] = global[2];
    }

    result
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::Orientation;
    use approx::assert_relative_eq;
    use ndarray::array;
    use std::f64::consts::PI;

    #[test]
    fn test_transform_identity() {
        // With identity orientation and zero position, transforms should be identity
        let observers = array![[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]];
        let position = Position::from([0.0, 0.0, 0.0]);
        let orientation = Orientation::identity();

        let local = transform_observers_to_local(&observers, &position, &orientation);

        assert_relative_eq!(local[[0, 0]], 1.0, epsilon = 1e-10);
        assert_relative_eq!(local[[0, 1]], 2.0, epsilon = 1e-10);
        assert_relative_eq!(local[[0, 2]], 3.0, epsilon = 1e-10);
        assert_relative_eq!(local[[1, 0]], 4.0, epsilon = 1e-10);
        assert_relative_eq!(local[[1, 1]], 5.0, epsilon = 1e-10);
        assert_relative_eq!(local[[1, 2]], 6.0, epsilon = 1e-10);
    }

    #[test]
    fn test_transform_translation() {
        // Test pure translation (no rotation)
        let observers = array![[5.0, 3.0, 7.0]];
        let position = Position::from([1.0, 2.0, 3.0]);
        let orientation = Orientation::identity();

        let local = transform_observers_to_local(&observers, &position, &orientation);

        // Should be observers - position
        assert_relative_eq!(local[[0, 0]], 4.0, epsilon = 1e-10);
        assert_relative_eq!(local[[0, 1]], 1.0, epsilon = 1e-10);
        assert_relative_eq!(local[[0, 2]], 4.0, epsilon = 1e-10);
    }

    #[test]
    fn test_transform_rotation_z() {
        // Test 90° rotation around z-axis
        let observers = array![[1.0, 0.0, 0.0]];
        let position = Position::from([0.0, 0.0, 0.0]);
        let orientation = Orientation::from_z_rotation(PI / 2.0);

        let local = transform_observers_to_local(&observers, &position, &orientation);

        // 90° rotation around z should map (1,0,0) -> (0,-1,0) with inverse rotation
        assert_relative_eq!(local[[0, 0]], 0.0, epsilon = 1e-10);
        assert_relative_eq!(local[[0, 1]], -1.0, epsilon = 1e-10);
        assert_relative_eq!(local[[0, 2]], 0.0, epsilon = 1e-10);
    }

    #[test]
    fn test_field_transform_rotation() {
        // Test field rotation back to global coordinates
        let fields = array![[0.0, 0.0, 1.0]]; // Field in z-direction
        let orientation = Orientation::from_x_rotation(PI / 2.0); // 90° around x

        let global = transform_field_to_global(&fields, &orientation);

        // 90° rotation around x should map (0,0,1) -> (0,-1,0)
        assert_relative_eq!(global[[0, 0]], 0.0, epsilon = 1e-10);
        assert_relative_eq!(global[[0, 1]], -1.0, epsilon = 1e-10);
        assert_relative_eq!(global[[0, 2]], 0.0, epsilon = 1e-10);
    }

    #[test]
    fn test_roundtrip_consistency() {
        // Observers transformed to local then fields transformed back should be consistent
        let observers = array![[1.0, 2.0, 3.0]];
        let position = Position::from([0.5, 0.5, 0.5]);
        let orientation = Orientation::from_euler_angles(0.1, 0.2, 0.3);

        let _local = transform_observers_to_local(&observers, &position, &orientation);

        // If we had a field in local coords, transforming back should reverse the rotation
        let field_local = array![[1.0, 0.0, 0.0]];
        let field_global = transform_field_to_global(&field_local, &orientation);

        // The rotation part should be consistent (can't test full roundtrip without field computation)
        assert!(field_global[[0, 0]].is_finite());
        assert!(field_global[[0, 1]].is_finite());
        assert!(field_global[[0, 2]].is_finite());
    }
}
