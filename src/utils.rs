//! Utility functions and physical constants

use crate::types::Vec3;
use std::f64::consts::PI;

/// Permeability of free space (μ₀) in T·m/A or H/m
/// From scipy.constants.mu_0
pub const MU_0: f64 = 1.25663706212e-6;

/// 4π constant (used frequently in field calculations)
pub const FOUR_PI: f64 = 4.0 * PI;

/// Conversion factor from H-field to B-field in vacuum: B = μ₀ * H
pub const H_TO_B: f64 = MU_0;

/// Tolerance for numerical comparisons
pub const EPSILON: f64 = 1e-10;

/// Check if a value is effectively zero
pub fn is_zero(x: f64) -> bool {
    x.abs() < EPSILON
}

/// Safe division that returns 0 if denominator is zero
pub fn safe_div(numerator: f64, denominator: f64) -> f64 {
    if is_zero(denominator) {
        0.0
    } else {
        numerator / denominator
    }
}

/// Compute the norm (magnitude) of a 3D vector
pub fn norm(v: Vec3) -> f64 {
    (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
}

/// Normalize a 3D vector (returns zero vector if input is zero)
pub fn normalize(v: Vec3) -> Vec3 {
    let n = norm(v);
    if is_zero(n) {
        [0.0, 0.0, 0.0]
    } else {
        [v[0] / n, v[1] / n, v[2] / n]
    }
}

/// Dot product of two 3D vectors
pub fn dot(a: Vec3, b: Vec3) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Cross product of two 3D vectors
pub fn cross(a: Vec3, b: Vec3) -> Vec3 {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// Vector subtraction
pub fn sub(a: Vec3, b: Vec3) -> Vec3 {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

/// Vector addition
pub fn add(a: Vec3, b: Vec3) -> Vec3 {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

/// Scalar multiplication
pub fn scale(v: Vec3, s: f64) -> Vec3 {
    [v[0] * s, v[1] * s, v[2] * s]
}

/// Convert Cartesian coordinates to cylindrical coordinates
///
/// # Arguments
/// * `x`, `y`, `z` - Cartesian coordinates
///
/// # Returns
/// * `(r, phi, z)` - Cylindrical coordinates where:
///   - `r` is the radial distance from the z-axis (always >= 0)
///   - `phi` is the azimuthal angle in radians (atan2(y, x), range [-π, π])
///   - `z` is the height (same as input)
pub fn cart_to_cyl(x: f64, y: f64, z: f64) -> (f64, f64, f64) {
    let r = (x * x + y * y).sqrt();
    let phi = y.atan2(x);
    (r, phi, z)
}

/// Convert cylindrical coordinates to Cartesian coordinates
///
/// # Arguments
/// * `r` - Radial distance from z-axis (should be >= 0)
/// * `phi` - Azimuthal angle in radians
/// * `z` - Height
///
/// # Returns
/// * `(x, y, z)` - Cartesian coordinates
pub fn cyl_to_cart(r: f64, phi: f64, z: f64) -> (f64, f64, f64) {
    let x = r * phi.cos();
    let y = r * phi.sin();
    (x, y, z)
}

/// Transform a vector field from cylindrical to Cartesian components
///
/// Given a vector field in cylindrical coordinates (B_r, B_phi, B_z) at a point
/// with azimuthal angle phi, compute the Cartesian components (B_x, B_y, B_z).
///
/// # Arguments
/// * `b_r` - Radial component (points away from z-axis)
/// * `b_phi` - Azimuthal component (tangent to circle around z-axis)
/// * `b_z` - Axial component (along z-axis)
/// * `phi` - Azimuthal angle of the point (radians)
///
/// # Returns
/// * `(b_x, b_y, b_z)` - Cartesian components
///
/// # Formula
/// ```text
/// B_x = B_r * cos(phi) - B_phi * sin(phi)
/// B_y = B_r * sin(phi) + B_phi * cos(phi)
/// B_z = B_z (unchanged)
/// ```
pub fn cyl_field_to_cart(b_r: f64, b_phi: f64, b_z: f64, phi: f64) -> (f64, f64, f64) {
    let cos_phi = phi.cos();
    let sin_phi = phi.sin();

    let b_x = b_r * cos_phi - b_phi * sin_phi;
    let b_y = b_r * sin_phi + b_phi * cos_phi;

    (b_x, b_y, b_z)
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_norm() {
        let v = [3.0, 4.0, 0.0];
        assert_relative_eq!(norm(v), 5.0);
    }

    #[test]
    fn test_normalize() {
        let v = [3.0, 4.0, 0.0];
        let n = normalize(v);
        assert_relative_eq!(n[0], 0.6);
        assert_relative_eq!(n[1], 0.8);
        assert_relative_eq!(norm(n), 1.0);
    }

    #[test]
    fn test_dot() {
        let a = [1.0, 2.0, 3.0];
        let b = [4.0, 5.0, 6.0];
        assert_relative_eq!(dot(a, b), 32.0);
    }

    #[test]
    fn test_cross() {
        let a = [1.0, 0.0, 0.0];
        let b = [0.0, 1.0, 0.0];
        let c = cross(a, b);
        assert_relative_eq!(c[0], 0.0);
        assert_relative_eq!(c[1], 0.0);
        assert_relative_eq!(c[2], 1.0);
    }

    #[test]
    fn test_safe_div() {
        assert_relative_eq!(safe_div(10.0, 2.0), 5.0);
        assert_relative_eq!(safe_div(10.0, 0.0), 0.0);
    }

    #[test]
    fn test_cart_to_cyl() {
        // Point on positive x-axis
        let (r, phi, z) = cart_to_cyl(1.0, 0.0, 5.0);
        assert_relative_eq!(r, 1.0);
        assert_relative_eq!(phi, 0.0);
        assert_relative_eq!(z, 5.0);

        // Point on positive y-axis
        let (r, phi, z) = cart_to_cyl(0.0, 1.0, 5.0);
        assert_relative_eq!(r, 1.0);
        assert_relative_eq!(phi, PI / 2.0);
        assert_relative_eq!(z, 5.0);

        // Point on negative x-axis
        let (r, phi, z) = cart_to_cyl(-1.0, 0.0, 5.0);
        assert_relative_eq!(r, 1.0);
        assert_relative_eq!(phi, PI, epsilon = 1e-10);
        assert_relative_eq!(z, 5.0);

        // Point at 45 degrees in first quadrant
        let (r, phi, z) = cart_to_cyl(1.0, 1.0, 2.0);
        assert_relative_eq!(r, std::f64::consts::SQRT_2);
        assert_relative_eq!(phi, PI / 4.0);
        assert_relative_eq!(z, 2.0);

        // Origin
        let (r, phi, z) = cart_to_cyl(0.0, 0.0, 0.0);
        assert_relative_eq!(r, 0.0);
        assert_relative_eq!(z, 0.0);
        // phi is undefined at origin but atan2(0,0) returns 0
        assert_relative_eq!(phi, 0.0);
    }

    #[test]
    fn test_cyl_to_cart() {
        // Point at phi=0 (positive x-axis)
        let (x, y, z) = cyl_to_cart(1.0, 0.0, 5.0);
        assert_relative_eq!(x, 1.0);
        assert_relative_eq!(y, 0.0, epsilon = 1e-10);
        assert_relative_eq!(z, 5.0);

        // Point at phi=π/2 (positive y-axis)
        let (x, y, z) = cyl_to_cart(1.0, PI / 2.0, 5.0);
        assert_relative_eq!(x, 0.0, epsilon = 1e-10);
        assert_relative_eq!(y, 1.0);
        assert_relative_eq!(z, 5.0);

        // Point at phi=π (negative x-axis)
        let (x, y, z) = cyl_to_cart(1.0, PI, 5.0);
        assert_relative_eq!(x, -1.0, epsilon = 1e-10);
        assert_relative_eq!(y, 0.0, epsilon = 1e-10);
        assert_relative_eq!(z, 5.0);

        // Origin
        let (x, y, z) = cyl_to_cart(0.0, 0.0, 0.0);
        assert_relative_eq!(x, 0.0);
        assert_relative_eq!(y, 0.0);
        assert_relative_eq!(z, 0.0);
    }

    #[test]
    fn test_cart_cyl_roundtrip() {
        // Test that cart -> cyl -> cart gives original point
        let original = (3.0, 4.0, 5.0);
        let (r, phi, z_cyl) = cart_to_cyl(original.0, original.1, original.2);
        let (x, y, z) = cyl_to_cart(r, phi, z_cyl);
        assert_relative_eq!(x, original.0, epsilon = 1e-10);
        assert_relative_eq!(y, original.1, epsilon = 1e-10);
        assert_relative_eq!(z, original.2, epsilon = 1e-10);
    }

    #[test]
    fn test_cyl_field_to_cart() {
        // Pure radial field at phi=0 (positive x-axis)
        // B_r = 1, B_phi = 0, B_z = 0 at phi=0
        // Should give B_x = 1, B_y = 0, B_z = 0
        let (b_x, b_y, b_z) = cyl_field_to_cart(1.0, 0.0, 0.0, 0.0);
        assert_relative_eq!(b_x, 1.0);
        assert_relative_eq!(b_y, 0.0, epsilon = 1e-10);
        assert_relative_eq!(b_z, 0.0);

        // Pure azimuthal field at phi=0
        // B_r = 0, B_phi = 1, B_z = 0 at phi=0
        // Should give B_x = 0, B_y = 1, B_z = 0
        let (b_x, b_y, b_z) = cyl_field_to_cart(0.0, 1.0, 0.0, 0.0);
        assert_relative_eq!(b_x, 0.0, epsilon = 1e-10);
        assert_relative_eq!(b_y, 1.0);
        assert_relative_eq!(b_z, 0.0);

        // Pure radial field at phi=π/2 (positive y-axis)
        // B_r = 1, B_phi = 0, B_z = 0 at phi=π/2
        // Should give B_x = 0, B_y = 1, B_z = 0
        let (b_x, b_y, b_z) = cyl_field_to_cart(1.0, 0.0, 0.0, PI / 2.0);
        assert_relative_eq!(b_x, 0.0, epsilon = 1e-10);
        assert_relative_eq!(b_y, 1.0);
        assert_relative_eq!(b_z, 0.0);

        // Pure azimuthal field at phi=π/2
        // B_r = 0, B_phi = 1, B_z = 0 at phi=π/2
        // Should give B_x = -1, B_y = 0, B_z = 0
        let (b_x, b_y, b_z) = cyl_field_to_cart(0.0, 1.0, 0.0, PI / 2.0);
        assert_relative_eq!(b_x, -1.0, epsilon = 1e-10);
        assert_relative_eq!(b_y, 0.0, epsilon = 1e-10);
        assert_relative_eq!(b_z, 0.0);

        // Mixed field at phi=π/4
        let (b_x, b_y, b_z) = cyl_field_to_cart(2.0, 3.0, 5.0, PI / 4.0);
        let cos_45 = (PI / 4.0).cos();
        let sin_45 = (PI / 4.0).sin();
        assert_relative_eq!(b_x, 2.0 * cos_45 - 3.0 * sin_45, epsilon = 1e-10);
        assert_relative_eq!(b_y, 2.0 * sin_45 + 3.0 * cos_45, epsilon = 1e-10);
        assert_relative_eq!(b_z, 5.0);
    }
}
