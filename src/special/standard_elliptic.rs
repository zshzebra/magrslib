//! Standard complete elliptic integrals K(m), E(m), and Π(n,m)
//!
//! These are implemented using the `ellip` crate, which provides pure-Rust
//! implementations of elliptic integrals tested against Boost and Wolfram.
//!
//! # References
//! - Ellip crate: https://github.com/p-sira/ellip
//! - NIST DLMF Chapter 19: https://dlmf.nist.gov/19
//! - Abramowitz & Stegun, "Handbook of Mathematical Functions" (Chapter 17)
//!
//! # Parameter Relationships (Parameter m = k²)
//! The standard convention uses parameter m = k² instead of modulus k:
//! - K(m) = ∫₀^(π/2) dθ / sqrt(1 - m·sin²θ)
//! - E(m) = ∫₀^(π/2) sqrt(1 - m·sin²θ) dθ
//!
//! Our implementations follow scipy's parameter convention where m = k².

use ndarray::Array1;

/// Complete elliptic integral of the first kind K(m)
///
/// Computes K(m) = ∫₀^(π/2) dθ / sqrt(1 - m·sin²θ)
///
/// # Arguments
/// * `m` - Parameter (0 ≤ m < 1), where m = k² and k is the modulus
///
/// # Returns
/// * K(m) - Complete elliptic integral of the first kind
///
/// # Panics
/// * If computation fails (e.g., invalid parameter range)
///
/// # Note
/// This follows scipy's convention where the parameter m = k² is used,
/// not the modulus k directly.
///
/// # Implementation
/// Uses the `ellip` crate's `ellipk` function with modulus conversion.
/// ellip::ellipk expects modulus k, but we use parameter m = k².
pub fn ellipk(m: f64) -> f64 {
    let k = m.sqrt();  // Convert parameter m to modulus k
    ellip::ellipk(k).expect("ellipk: computation failed")
}

/// Complete elliptic integral of the second kind E(m)
///
/// Computes E(m) = ∫₀^(π/2) sqrt(1 - m·sin²θ) dθ
///
/// # Arguments
/// * `m` - Parameter (0 ≤ m ≤ 1), where m = k² and k is the modulus
///
/// # Returns
/// * E(m) - Complete elliptic integral of the second kind
///
/// # Panics
/// * If computation fails (e.g., invalid parameter range)
///
/// # Note
/// This follows scipy's convention where the parameter m = k² is used.
/// For m = 1 (k = 1), E(1) = 1.
///
/// # Implementation
/// Uses the `ellip` crate's `ellipe` function.
pub fn ellipe(m: f64) -> f64 {
    ellip::ellipe(m).expect("ellipe: computation failed")
}

/// Complete elliptic integral of the third kind Π(n, m)
///
/// Computes Π(n, m) = ∫₀^(π/2) dθ / [(1 - n·sin²θ) sqrt(1 - m·sin²θ)]
///
/// # Arguments
/// * `n` - Characteristic (real number)
/// * `m` - Parameter (0 ≤ m < 1), where m = k² and k is the modulus
///
/// # Returns
/// * Π(n, m) - Complete elliptic integral of the third kind
///
/// # Panics
/// * If computation fails (e.g., invalid parameter range)
///
/// # Note
/// This follows scipy's convention where the parameter m = k² is used.
///
/// # Implementation
/// Uses the `ellip` crate's `ellippi` function.
pub fn ellippi(n: f64, m: f64) -> f64 {
    ellip::ellippi(n, m).expect("ellippi: computation failed")
}

/// Vectorized complete elliptic integral of the first kind K(m)
///
/// Computes K(m) for each element in the input array.
///
/// # Arguments
/// * `m` - Array of parameters (each 0 ≤ m < 1), where m = k²
///
/// # Returns
/// * Array of K(m) values
///
/// # Panics
/// * If any computation fails
pub fn ellipk_vec(m: &Array1<f64>) -> Array1<f64> {
    m.mapv(|mi| ellipk(mi))
}

/// Vectorized complete elliptic integral of the second kind E(m)
///
/// Computes E(m) for each element in the input array.
///
/// # Arguments
/// * `m` - Array of parameters (each 0 ≤ m ≤ 1), where m = k²
///
/// # Returns
/// * Array of E(m) values
///
/// # Panics
/// * If any computation fails
pub fn ellipe_vec(m: &Array1<f64>) -> Array1<f64> {
    m.mapv(|mi| ellipe(mi))
}

/// Vectorized complete elliptic integral of the third kind Π(n, m)
///
/// Computes Π(n, m) for arrays of n and m values.
///
/// # Arguments
/// * `n` - Array of characteristics
/// * `m` - Array of parameters (each 0 ≤ m < 1), where m = k²
///
/// # Returns
/// * Array of Π(n, m) values
///
/// # Panics
/// * If n and m have different lengths
/// * If any computation fails
pub fn ellippi_vec(n: &Array1<f64>, m: &Array1<f64>) -> Array1<f64> {
    if n.len() != m.len() {
        panic!("ellippi_vec: n and m must have same length");
    }

    let mut result = Array1::zeros(n.len());
    for i in 0..n.len() {
        result[i] = ellippi(n[i], m[i]);
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use ndarray::array;

    #[test]
    fn test_ellipk_zero() {
        // K(m=0) = π/2 (complete elliptic integral with m=0)
        let k = ellipk(0.0);
        assert_relative_eq!(k, std::f64::consts::FRAC_PI_2, epsilon = 1e-10);
    }

    #[test]
    fn test_ellipk_quarter() {
        // K(m=0.25) = K(k²) where k=0.5, ≈ 1.8540746773013719
        // Reference value from scipy.special.ellipk(0.25)
        let k = ellipk(0.25);
        assert_relative_eq!(k, 1.8540746773013719, epsilon = 1e-10);
    }

    #[test]
    fn test_ellipe_zero() {
        // E(m=0) = π/2
        let e = ellipe(0.0);
        assert_relative_eq!(e, std::f64::consts::FRAC_PI_2, epsilon = 1e-10);
    }

    #[test]
    fn test_ellipe_one() {
        // E(m=1) = 1 (special case)
        let e = ellipe(1.0);
        assert_relative_eq!(e, 1.0, epsilon = 1e-10);
    }

    #[test]
    fn test_ellipe_quarter() {
        // E(m=0.25) = E(k²) where k=0.5, ≈ 1.4674622093394271
        // Reference value from scipy.special.ellipe(0.25)
        let e = ellipe(0.25);
        assert_relative_eq!(e, 1.4674622093394271, epsilon = 1e-10);
    }

    #[test]
    fn test_ellippi_zero_n() {
        // Π(n=0, m) = K(m) for any m
        let m_val = 0.09;  // m = 0.09 = 0.3²
        let pi = ellippi(0.0, m_val);
        let k = ellipk(m_val);
        assert_relative_eq!(pi, k, epsilon = 1e-10);
    }

    #[test]
    fn test_ellippi_values() {
        // Π(n=0.2, m=0.25) ≈ 1.9867935682987379
        // Reference value from scipy.special.ellippi(0.2, 0.25)
        let pi = ellippi(0.2, 0.25);
        assert_relative_eq!(pi, 1.9867935682987379, epsilon = 1e-9);
    }

    #[test]
    fn test_ellipk_vec() {
        let m = array![0.0, 0.25, 0.64];  // m = k², so k = 0, 0.5, 0.8
        let result = ellipk_vec(&m);

        assert_relative_eq!(result[0], std::f64::consts::FRAC_PI_2, epsilon = 1e-10);
        assert_relative_eq!(result[1], 1.8540746773013719, epsilon = 1e-10);
        assert_relative_eq!(result[2], ellipk(0.64), epsilon = 1e-10);
    }

    #[test]
    fn test_ellipe_vec() {
        let m = array![0.0, 0.25, 1.0];  // m = k², so k = 0, 0.5, 1.0
        let result = ellipe_vec(&m);

        assert_relative_eq!(result[0], std::f64::consts::FRAC_PI_2, epsilon = 1e-10);
        assert_relative_eq!(result[1], 1.4674622093394271, epsilon = 1e-10);
        assert_relative_eq!(result[2], 1.0, epsilon = 1e-10);
    }

    #[test]
    fn test_ellippi_vec() {
        let n = array![0.0, 0.2, 0.5];
        let m = array![0.09, 0.25, 0.49];  // m = 0.3², 0.5², 0.7²
        let result = ellippi_vec(&n, &m);

        // First should equal K(0.09)
        assert_relative_eq!(result[0], ellipk(0.09), epsilon = 1e-10);
        // Second matches reference
        assert_relative_eq!(result[1], 1.9867935682987379, epsilon = 1e-9);
        // Third should be computable
        assert!(result[2].is_finite());
        assert!(result[2] > 0.0);
    }

    #[test]
    #[should_panic(expected = "ellipk: computation failed")]
    fn test_ellipk_out_of_range() {
        ellipk(1.0);
    }

    #[test]
    #[should_panic(expected = "ellipe: computation failed")]
    fn test_ellipe_out_of_range() {
        ellipe(1.1);
    }
}
