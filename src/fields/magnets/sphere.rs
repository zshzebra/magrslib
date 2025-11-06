//! B-field computation for homogeneously magnetized sphere magnets
//!
//! Algorithm from Python magpylib field_BH_sphere.py
//!
//! Inside the sphere: B = (2/3) * J (polarization)
//! Outside the sphere: Dipole field B = [(3(J·r)r - J*r²) / r⁵] * (r_sphere³ / 3)

use ndarray::{Array1, Array2};

/// Compute B-field for a sphere magnet
///
/// # Arguments
/// * `observers` - Observer positions (n_points, 3) relative to sphere center
/// * `diameter` - Sphere diameter (meters)
/// * `polarization` - Magnetic polarization vector [Jx, Jy, Jz] (Tesla)
///
/// # Returns
/// * B-field vectors (n_points, 3) in Tesla
///
/// # Algorithm
/// - Inside (r <= r_sphere): B = (2/3) * J
/// - Outside (r > r_sphere): Dipole field
///
/// # Reference
/// "Theoretical Physics", Bertelmann
pub fn compute_sphere_field(
    observers: &Array2<f64>,
    diameter: f64,
    polarization: &[f64; 3],
) -> Array2<f64> {
    let n_points = observers.nrows();
    let mut b_field = Array2::zeros((n_points, 3));

    let r_sphere = diameter.abs() / 2.0;
    let r_sphere_cubed = r_sphere * r_sphere * r_sphere;

    // Precompute polarization magnitudes
    let pol_x = polarization[0];
    let pol_y = polarization[1];
    let pol_z = polarization[2];

    for i in 0..n_points {
        let x = observers[[i, 0]];
        let y = observers[[i, 1]];
        let z = observers[[i, 2]];

        // Compute distance from sphere center
        let r = (x * x + y * y + z * z).sqrt();

        if r <= r_sphere {
            // Inside sphere: B = (2/3) * polarization
            b_field[[i, 0]] = (2.0 / 3.0) * pol_x;
            b_field[[i, 1]] = (2.0 / 3.0) * pol_y;
            b_field[[i, 2]] = (2.0 / 3.0) * pol_z;
        } else {
            // Outside sphere: Dipole field
            // B = [(3(J·r)r - J*r²) / r⁵] * (r_sphere³ / 3)

            let r_squared = r * r;
            let r_to_5 = r_squared * r_squared * r;

            // Compute J·r (dot product)
            let j_dot_r = pol_x * x + pol_y * y + pol_z * z;

            // Compute 3(J·r)r
            let three_jdotr_r_x = 3.0 * j_dot_r * x;
            let three_jdotr_r_y = 3.0 * j_dot_r * y;
            let three_jdotr_r_z = 3.0 * j_dot_r * z;

            // Compute J*r²
            let j_rsq_x = pol_x * r_squared;
            let j_rsq_y = pol_y * r_squared;
            let j_rsq_z = pol_z * r_squared;

            // Numerator: 3(J·r)r - J*r²
            let num_x = three_jdotr_r_x - j_rsq_x;
            let num_y = three_jdotr_r_y - j_rsq_y;
            let num_z = three_jdotr_r_z - j_rsq_z;

            // Full formula: numerator / r⁵ * (r_sphere³ / 3)
            let scale = r_sphere_cubed / (3.0 * r_to_5);

            b_field[[i, 0]] = num_x * scale;
            b_field[[i, 1]] = num_y * scale;
            b_field[[i, 2]] = num_z * scale;
        }
    }

    b_field
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use ndarray::array;

    #[test]
    fn test_sphere_inside() {
        // Test point inside sphere: should be (2/3) * polarization
        let observers = array![[0.0, 0.0, 0.0]];  // Center
        let diameter = 2.0;  // 2m diameter (1m radius)
        let polarization = [1.0, 0.0, 0.0];  // 1T in x direction

        let b_field = compute_sphere_field(&observers, diameter, &polarization);

        assert_relative_eq!(b_field[[0, 0]], 2.0 / 3.0, epsilon = 1e-10);
        assert_relative_eq!(b_field[[0, 1]], 0.0, epsilon = 1e-10);
        assert_relative_eq!(b_field[[0, 2]], 0.0, epsilon = 1e-10);
    }

    #[test]
    fn test_sphere_at_boundary() {
        // Test point exactly at sphere surface
        let observers = array![[1.0, 0.0, 0.0]];  // At surface (r = radius)
        let diameter = 2.0;  // 2m diameter (1m radius)
        let polarization = [1.0, 0.0, 0.0];

        let b_field = compute_sphere_field(&observers, diameter, &polarization);

        // At boundary, should use inside formula
        assert_relative_eq!(b_field[[0, 0]], 2.0 / 3.0, epsilon = 1e-10);
    }

    #[test]
    fn test_sphere_outside_on_axis() {
        // Test point outside on polarization axis
        let observers = array![[2.0, 0.0, 0.0]];  // 2m from center on x-axis
        let diameter = 1.0;  // 1m diameter (0.5m radius)
        let polarization = [1.0, 0.0, 0.0];  // 1T in x direction

        let b_field = compute_sphere_field(&observers, diameter, &polarization);

        // Outside: Dipole field
        // For point on axis parallel to polarization:
        // B_x = (r_sphere³ / 3) * (3 * J * x² - J * r²) / r⁵
        //     = (r_sphere³ / 3) * J * (3x² - r²) / r⁵
        //     = (r_sphere³ / 3) * J * (3*4 - 4) / 32
        //     = (r_sphere³ / 3) * J * 8 / 32
        //     = (r_sphere³ / 3) * J / 4
        let r_sphere = 0.5;
        let r_sphere_cubed = r_sphere * r_sphere * r_sphere;
        let expected_bx = r_sphere_cubed / (3.0 * 4.0);  // J=1, factor=1/4

        assert_relative_eq!(b_field[[0, 0]], expected_bx, epsilon = 1e-10);
        assert_relative_eq!(b_field[[0, 1]], 0.0, epsilon = 1e-10);
        assert_relative_eq!(b_field[[0, 2]], 0.0, epsilon = 1e-10);
    }

    #[test]
    fn test_sphere_outside_off_axis() {
        // Test point outside off-axis
        let observers = array![[1.0, 1.0, 0.0]];
        let diameter = 0.5;
        let polarization = [0.0, 0.0, 1.0];  // Polarized in z

        let b_field = compute_sphere_field(&observers, diameter, &polarization);

        // Should have non-zero x, y components from dipole field
        // (because observer is off-axis)
        assert!(b_field[[0, 0]].abs() < 0.1);  // Small but not exactly zero
        assert!(b_field[[0, 1]].abs() < 0.1);
        // z-component should be negative (dipole field pattern)
        assert!(b_field[[0, 2]] < 0.0);
    }

    #[test]
    fn test_sphere_multiple_points() {
        // Test multiple observation points
        let observers = array![
            [0.0, 0.0, 0.0],  // Inside (center)
            [2.0, 0.0, 0.0],  // Outside
            [0.0, 2.0, 0.0],  // Outside
        ];
        let diameter = 1.0;
        let polarization = [1.0, 0.0, 0.0];

        let b_field = compute_sphere_field(&observers, diameter, &polarization);

        // First point (inside): (2/3) * polarization
        assert_relative_eq!(b_field[[0, 0]], 2.0 / 3.0, epsilon = 1e-10);

        // Second and third points (outside): should be dipole field
        assert!(b_field[[1, 0]].abs() > 0.0);
        assert!(b_field[[2, 0]].abs() > 0.0);
    }
}
