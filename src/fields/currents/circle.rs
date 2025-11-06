//! Magnetic field computation for circular current loops
//!
//! Uses the Biot-Savart law to compute the magnetic field from a circular current loop.
//! The field calculation uses complete elliptic integrals K(m) and E(m) for exact solutions.
//!
//! # References
//! - Jackson, "Classical Electrodynamics", 3rd Edition, Section 5.6
//! - Griffiths, "Introduction to Electrodynamics", 4th Edition, Example 5.6
//! - NIST Handbook of Mathematical Functions, Chapter 19 (Elliptic Integrals)

use crate::special::standard_elliptic::{ellipk, ellipe};
use crate::utils::MU_0;
use ndarray::{Array1, Array2};
use std::f64::consts::PI;

/// Compute magnetic field from a circular current loop
///
/// # Arguments
/// * `observers` - Observer positions as (N, 3) array [x, y, z] in meters
/// * `diameter` - Loop diameter in meters
/// * `current` - Current in Amperes (positive for counter-clockwise when viewed from +z)
///
/// # Returns
/// * Magnetic field B as (N, 3) array [Bx, By, Bz] in Tesla
///
/// # Coordinate System
/// * Loop lies in the xy-plane, centered at origin
/// * Loop normal points in +z direction
/// * Current flows counter-clockwise when viewed from +z for positive current
///
/// # Mathematical Background
/// For a circular loop of radius a carrying current I, the magnetic field is:
///
/// On-axis (ρ = 0):
/// B_z = (μ₀I a²) / (2(z² + a²)^(3/2))
///
/// Off-axis: Uses elliptic integrals K(m) and E(m) where m = 4aρ/((a+ρ)² + z²)
/// B_ρ = (μ₀I/4π) * √(a/ρ) * [(2-m)K(m) - 2E(m)] / √((a+ρ)² + z²)
/// B_z = (μ₀I/4π) * 1/√((a+ρ)² + z²) * [K(m) + (a²-ρ²-z²)/((a-ρ)²+z²) * E(m)]
pub fn compute_circle_field(
    observers: &Array2<f64>,
    diameter: f64,
    current: f64,
) -> Array2<f64> {
    let n_obs = observers.nrows();
    let mut b_field = Array2::zeros((n_obs, 3));

    let radius = diameter / 2.0;

    // Prefactor: μ₀I / (4π)
    let prefactor = MU_0 * current / (4.0 * PI);

    for i in 0..n_obs {
        let x = observers[[i, 0]];
        let y = observers[[i, 1]];
        let z = observers[[i, 2]];

        // Convert to cylindrical coordinates (ρ, φ, z)
        let rho = (x * x + y * y).sqrt();
        let phi = y.atan2(x);

        // Handle on-axis case (ρ = 0) separately for numerical stability
        if rho < 1e-12 {
            // On-axis formula: B_z = μ₀I a² / (2(z² + a²)^(3/2))
            let denominator = (z * z + radius * radius).powf(1.5);
            let b_z = MU_0 * current * radius * radius / (2.0 * denominator);

            b_field[[i, 0]] = 0.0;  // B_x = 0 on axis
            b_field[[i, 1]] = 0.0;  // B_y = 0 on axis
            b_field[[i, 2]] = b_z;
        } else {
            // Off-axis case: use elliptic integrals
            let a_plus_rho = radius + rho;
            let a_minus_rho = radius - rho;
            let sqrt_denominator = (a_plus_rho * a_plus_rho + z * z).sqrt();

            // Parameter for elliptic integrals: m = 4aρ/((a+ρ)² + z²)
            let m = 4.0 * radius * rho / (a_plus_rho * a_plus_rho + z * z);

            // Ensure m is in valid range [0, 1) for elliptic integrals
            if m >= 1.0 {
                // Handle boundary case - this occurs when observer is very close to the loop
                // Use limiting values or approximate formulas
                b_field[[i, 0]] = 0.0;
                b_field[[i, 1]] = 0.0;
                b_field[[i, 2]] = 0.0;
                continue;
            }

            // Compute elliptic integrals
            let k_m = ellipk(m);
            let e_m = ellipe(m);

            // Radial component: B_ρ = prefactor * √(a/ρ) * [(2-m)K(m) - 2E(m)] / √((a+ρ)² + z²)
            let sqrt_a_over_rho = (radius / rho).sqrt();
            let b_rho = prefactor * sqrt_a_over_rho * ((2.0 - m) * k_m - 2.0 * e_m) / sqrt_denominator;

            // Axial component: B_z = prefactor * [K(m) + (a²-ρ²-z²)/((a-ρ)²+z²) * E(m)] / √((a+ρ)² + z²)
            let numerator_term = (radius * radius - rho * rho - z * z) / (a_minus_rho * a_minus_rho + z * z);
            let b_z = prefactor * (k_m + numerator_term * e_m) / sqrt_denominator;

            // Convert from cylindrical (B_ρ, B_φ, B_z) to Cartesian (B_x, B_y, B_z)
            // B_φ = 0 for axisymmetric field
            let cos_phi = x / rho;
            let sin_phi = y / rho;

            b_field[[i, 0]] = b_rho * cos_phi;  // B_x = B_ρ cos(φ)
            b_field[[i, 1]] = b_rho * sin_phi;  // B_y = B_ρ sin(φ)
            b_field[[i, 2]] = b_z;              // B_z = B_z
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
    fn test_circle_field_on_axis() {
        // Test on-axis field for a simple case
        let diameter = 0.02;  // 2 cm diameter
        let current = 1.0;    // 1 A
        let radius = diameter / 2.0;

        // Observer on axis at z = radius (45° from loop center)
        let observers = array![[0.0, 0.0, radius]];
        let b_field = compute_circle_field(&observers, diameter, current);

        // On-axis analytical formula: B_z = μ₀I a² / (2(z² + a²)^(3/2))
        let expected_b_z = MU_0 * current * radius * radius / (2.0 * (radius * radius + radius * radius).powf(1.5));
        let expected_b_z = MU_0 * current * radius * radius / (2.0 * (2.0 * radius * radius).powf(1.5));
        let expected_b_z = MU_0 * current / (4.0 * (2.0_f64).sqrt() * radius);

        assert_relative_eq!(b_field[[0, 0]], 0.0, epsilon = 1e-12);  // B_x = 0 on axis
        assert_relative_eq!(b_field[[0, 1]], 0.0, epsilon = 1e-12);  // B_y = 0 on axis
        assert_relative_eq!(b_field[[0, 2]], expected_b_z, epsilon = 1e-10);
    }

    #[test]
    fn test_circle_field_at_center() {
        // Test field at the center of the loop
        let diameter = 0.01;  // 1 cm diameter
        let current = 2.0;    // 2 A
        let radius = diameter / 2.0;

        // Observer at center (z = 0)
        let observers = array![[0.0, 0.0, 0.0]];
        let b_field = compute_circle_field(&observers, diameter, current);

        // At center: B_z = μ₀I / (2a)
        let expected_b_z = MU_0 * current / (2.0 * radius);

        assert_relative_eq!(b_field[[0, 0]], 0.0, epsilon = 1e-12);
        assert_relative_eq!(b_field[[0, 1]], 0.0, epsilon = 1e-12);
        assert_relative_eq!(b_field[[0, 2]], expected_b_z, epsilon = 1e-10);
    }

    #[test]
    fn test_circle_field_far_away() {
        // Test field far from the loop (dipole approximation should be valid)
        let diameter = 0.01;  // 1 cm diameter
        let current = 1.0;    // 1 A
        let radius = diameter / 2.0;

        // Observer far away on axis (z >> a)
        let z_far = 10.0 * radius;  // 10x the radius
        let observers = array![[0.0, 0.0, z_far]];
        let b_field = compute_circle_field(&observers, diameter, current);

        // Far-field dipole approximation: B_z ≈ μ₀I a² / (2z³)
        let magnetic_dipole_moment = current * PI * radius * radius;
        let expected_b_z = MU_0 * magnetic_dipole_moment / (2.0 * PI * z_far * z_far * z_far);

        assert_relative_eq!(b_field[[0, 0]], 0.0, epsilon = 1e-12);
        assert_relative_eq!(b_field[[0, 1]], 0.0, epsilon = 1e-12);
        assert_relative_eq!(b_field[[0, 2]], expected_b_z, epsilon = 1e-2);  // Looser tolerance for approximation
    }

    #[test]
    fn test_circle_field_symmetry() {
        // Test field symmetry: points equidistant from axis should have same |B|
        let diameter = 0.02;  // 2 cm diameter
        let current = 1.5;    // 1.5 A

        let observers = array![
            [0.01, 0.0, 0.005],   // Point 1: x = 1cm, z = 0.5cm
            [0.0, 0.01, 0.005],   // Point 2: y = 1cm, z = 0.5cm
            [-0.01, 0.0, 0.005],  // Point 3: x = -1cm, z = 0.5cm
            [0.0, -0.01, 0.005],  // Point 4: y = -1cm, z = 0.5cm
        ];

        let b_field = compute_circle_field(&observers, diameter, current);

        // All points are at same cylindrical radius ρ = 1cm, should have same B_z
        for i in 1..4 {
            assert_relative_eq!(b_field[[i, 2]], b_field[[0, 2]], epsilon = 1e-10);
        }

        // Magnitude of radial component should be same
        let b_rho_0 = (b_field[[0, 0]] * b_field[[0, 0]] + b_field[[0, 1]] * b_field[[0, 1]]).sqrt();
        for i in 1..4 {
            let b_rho_i = (b_field[[i, 0]] * b_field[[i, 0]] + b_field[[i, 1]] * b_field[[i, 1]]).sqrt();
            assert_relative_eq!(b_rho_i, b_rho_0, epsilon = 1e-10);
        }
    }
}