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
}
