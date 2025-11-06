//! Complete elliptic integrals using Bulirsch's CEL algorithm
//!
//! Direct port from Python magpylib special_cel.py
//! Used for computing magnetic fields from circular current loops and cylinders.
//!
//! Reference:
//! Bulirsch, R. (1965). Numerical calculation of elliptic integrals and elliptic functions.
//! Numerische Mathematik, 7(1), 78-90.

use ndarray::Array1;
use std::f64::consts::PI;

/// Half pi constant (π/2)
const HALF_PI: f64 = PI / 2.0;

/// Error tolerance for convergence
const ERR_TOL: f64 = 0.000001;

/// Maximum iterations to prevent infinite loops
const MAX_ITER: usize = 100;

/// Complete Elliptic Integral (scalar version)
///
/// Computes the complete elliptic integral using Bulirsch's CEL algorithm.
/// Direct port of Python `_cel0` function.
///
/// # Arguments
/// * `kc` - Complementary modulus (kc ≠ 0)
/// * `p` - Characteristic parameter
/// * `c` - C parameter
/// * `s` - S parameter
///
/// # Returns
/// The complete elliptic integral value
///
/// # Panics
/// Panics if kc == 0 (not allowed in CEL algorithm)
///
/// # Example
/// ```
/// use magrslib::special::elliptic::cel0;
///
/// let result = cel0(0.5, 1.0, 1.0, 1.0);
/// assert!(result.is_finite());
/// ```
pub fn cel0(kc: f64, p: f64, c: f64, s: f64) -> f64 {
    if kc == 0.0 {
        panic!("FAIL cel: kc=0 not allowed.");
    }

    let k = kc.abs();
    let mut pp = p;
    let mut cc = c;
    let mut ss = s;
    let mut em = 1.0;

    if p > 0.0 {
        pp = p.sqrt();
        ss = s / pp;
    } else {
        let mut f = kc * kc;
        let mut q = 1.0 - f;
        let g = 1.0 - pp;
        f = f - pp;
        q = q * (ss - c * pp);
        pp = (f / g).sqrt();
        cc = (c - ss) / g;
        ss = -q / (g * g * pp) + cc * pp;
    }

    let mut f = cc;
    cc = cc + ss / pp;
    let mut g = k / pp;
    ss = 2.0 * (ss + f * g);
    pp = g + pp;
    g = em;
    em = k + em;
    let mut kk = k;
    let mut k_iter = k;  // Separate variable for iteration

    let mut iter_count = 0;
    while (g - k_iter).abs() > g * ERR_TOL && iter_count < MAX_ITER {
        k_iter = 2.0 * kk.sqrt();
        kk = k_iter * em;
        f = cc;
        cc = cc + ss / pp;
        g = kk / pp;
        ss = 2.0 * (ss + f * g);
        pp = g + pp;
        g = em;
        em = k_iter + em;
        iter_count += 1;
    }

    HALF_PI * (ss + cc * em) / (em * (em + pp))
}

/// Complete Elliptic Integral (vectorized version)
///
/// Vectorized implementation of CEL for multiple values simultaneously.
/// Direct port of Python `_celv` function.
///
/// # Arguments
/// * `kc` - Complementary moduli array
/// * `p` - Characteristic parameters array
/// * `c` - C parameters array
/// * `s` - S parameters array
///
/// # Returns
/// Array of complete elliptic integral values
///
/// # Note
/// If any kc[i] == 0, the corresponding result will be NaN
pub fn celv(kc: Array1<f64>, p: Array1<f64>, c: Array1<f64>, s: Array1<f64>) -> Array1<f64> {
    let n = kc.len();

    let k = kc.mapv(f64::abs);
    let mut em = Array1::ones(n);

    let mut cc = c.clone();
    let mut pp = p.clone();
    let mut ss = s.clone();

    // Apply mask for p > 0 and p <= 0 cases
    for i in 0..n {
        if p[i] > 0.0 {
            pp[i] = p[i].sqrt();
            ss[i] = s[i] / pp[i];
        } else {
            let f = kc[i] * kc[i];
            let q = 1.0 - f;
            let g = 1.0 - pp[i];
            let f = f - pp[i];
            let q = q * (ss[i] - c[i] * pp[i]);
            pp[i] = (f / g).sqrt();
            cc[i] = (c[i] - ss[i]) / g;
            ss[i] = -q / (g * g * pp[i]) + cc[i] * pp[i];
        }
    }

    let mut f = cc.clone();
    cc = &cc + &ss / &pp;
    let mut g = &k / &pp;
    ss = 2.0 * (&ss + &f * &g);
    pp = &g + &pp;
    g = em.clone();
    em = &k + &em;
    let mut kk = k.clone();
    let mut k_iter = k.clone();  // Separate variable for iteration

    // Iterate with dynamic mask for convergence
    let mut mask = Array1::from_elem(n, true);
    let mut iter_count = 0;

    while mask.iter().any(|&m| m) && iter_count < MAX_ITER {
        for i in 0..n {
            if mask[i] {
                k_iter[i] = 2.0 * kk[i].sqrt();
            }
        }

        for i in 0..n {
            if mask[i] {
                kk[i] = k_iter[i] * em[i];
            }
        }

        for i in 0..n {
            if mask[i] {
                f[i] = cc[i];
            }
        }

        for i in 0..n {
            if mask[i] {
                cc[i] = cc[i] + ss[i] / pp[i];
            }
        }

        for i in 0..n {
            if mask[i] {
                g[i] = kk[i] / pp[i];
            }
        }

        for i in 0..n {
            if mask[i] {
                ss[i] = 2.0 * (ss[i] + f[i] * g[i]);
            }
        }

        for i in 0..n {
            if mask[i] {
                pp[i] = g[i] + pp[i];
            }
        }

        for i in 0..n {
            if mask[i] {
                g[i] = em[i];
            }
        }

        for i in 0..n {
            if mask[i] {
                em[i] = k_iter[i] + em[i];
            }
        }

        // Redefine mask - only continue where not converged
        for i in 0..n {
            mask[i] = (g[i] - k_iter[i]).abs() > g[i] * ERR_TOL;
        }

        iter_count += 1;
    }

    HALF_PI * (&ss + &cc * &em) / (&em * (&em + &pp))
}

/// Complete Elliptic Integral (adaptive dispatcher)
///
/// Automatically chooses between scalar loop and vectorized implementation
/// based on input size for optimal performance.
/// Direct port of Python `_cel` function.
///
/// # Arguments
/// * `kc` - Complementary moduli array
/// * `p` - Characteristic parameters array
/// * `c` - C parameters array
/// * `s` - S parameters array
///
/// # Returns
/// Array of complete elliptic integral values
///
/// # Note
/// For n < 10, uses scalar loop. For n >= 10, uses vectorized implementation.
pub fn cel(kc: Array1<f64>, p: Array1<f64>, c: Array1<f64>, s: Array1<f64>) -> Array1<f64> {
    let n_input = kc.len();

    if n_input < 10 {
        // Use scalar loop for small arrays
        let mut result = Array1::zeros(n_input);
        for i in 0..n_input {
            result[i] = cel0(kc[i], p[i], c[i], s[i]);
        }
        result
    } else {
        // Use vectorized version for larger arrays
        celv(kc, p, c, s)
    }
}

/// Iterative part of Bulirsch CEL algorithm (vectorized)
///
/// This is a direct port of the Python `_cel_iter` function.
/// Computes complete elliptic integrals using iterative refinement.
///
/// # Arguments
/// * `qc` - Modified modulus parameter
/// * `p` - Characteristic parameter
/// * `g` - Temporary variable for iteration
/// * `cc` - C coefficient
/// * `ss` - S coefficient
/// * `em` - M coefficient
/// * `kk` - K coefficient
///
/// # Returns
/// Complete elliptic integral value
pub fn cel_iter(
    mut qc: Array1<f64>,
    mut p: Array1<f64>,
    mut g: Array1<f64>,
    mut cc: Array1<f64>,
    mut ss: Array1<f64>,
    mut em: Array1<f64>,
    mut kk: Array1<f64>,
) -> Array1<f64> {
    let n = qc.len();

    // For small arrays, could use scalar loop (optimization)
    // For now, use vectorized version for all sizes

    let mut f = cc.clone();
    let mut iter_count = 0;

    // Iterate until convergence
    while iter_count < MAX_ITER {
        // Check convergence: |g - qc| < qc * tolerance
        let mut all_converged = true;
        for i in 0..n {
            if (g[i] - qc[i]).abs() >= qc[i] * ERR_TOL {
                all_converged = false;
                break;
            }
        }

        if all_converged {
            break;
        }

        // Update qc
        for i in 0..n {
            qc[i] = 2.0 * kk[i].sqrt();
        }

        // kk = qc * em
        for i in 0..n {
            kk[i] = qc[i] * em[i];
        }

        // f = cc (copy current cc values)
        f.assign(&cc);

        // cc = cc + ss / p
        for i in 0..n {
            cc[i] = cc[i] + ss[i] / p[i];
        }

        // g = kk / p
        for i in 0..n {
            g[i] = kk[i] / p[i];
        }

        // ss = 2 * (ss + f * g)
        for i in 0..n {
            ss[i] = 2.0 * (ss[i] + f[i] * g[i]);
        }

        // p = p + g
        for i in 0..n {
            p[i] = p[i] + g[i];
        }

        // g = em (copy current em values)
        g.assign(&em);

        // em = em + qc
        for i in 0..n {
            em[i] = em[i] + qc[i];
        }

        iter_count += 1;
    }

    // Compute final result: (π/2) * (ss + cc * em) / (em * (em + p))
    let mut result = Array1::zeros(n);
    for i in 0..n {
        result[i] = HALF_PI * (ss[i] + cc[i] * em[i]) / (em[i] * (em[i] + p[i]));
    }

    result
}

/// Scalar version of CEL iteration (for small arrays or single values)
///
/// More efficient for single evaluations or very small arrays.
pub fn cel_iter_scalar(
    mut qc: f64,
    mut p: f64,
    mut g: f64,
    mut cc: f64,
    mut ss: f64,
    mut em: f64,
    mut kk: f64,
) -> f64 {
    let mut iter_count = 0;

    while (g - qc).abs() >= qc * ERR_TOL && iter_count < MAX_ITER {
        qc = 2.0 * kk.sqrt();
        kk = qc * em;
        let f = cc;
        cc = cc + ss / p;
        g = kk / p;
        ss = 2.0 * (ss + f * g);
        p = p + g;
        g = em;
        em = em + qc;
        iter_count += 1;
    }

    HALF_PI * (ss + cc * em) / (em * (em + p))
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use ndarray::array;

    #[test]
    fn test_cel0_basic() {
        // Test basic CEL computation
        // Values from Python: cel0(0.5, 1.0, 1.0, 1.0)
        let result = cel0(0.5, 1.0, 1.0, 1.0);
        assert!(result.is_finite());
        assert!(result > 0.0);
    }

    #[test]
    #[should_panic(expected = "FAIL cel: kc=0 not allowed")]
    fn test_cel0_zero_kc_panics() {
        // CEL should panic when kc == 0
        cel0(0.0, 1.0, 1.0, 1.0);
    }

    #[test]
    fn test_cel0_p_positive() {
        // Test case where p > 0
        let result = cel0(0.7, 1.5, 1.0, 1.0);
        assert!(result.is_finite());
    }

    #[test]
    fn test_cel0_p_negative() {
        // Test case where p <= 0
        let result = cel0(0.7, -0.5, 1.0, 1.0);
        assert!(result.is_finite());
    }

    #[test]
    fn test_cel_vs_cel0() {
        // Test that cel dispatcher gives same result as cel0 for single value
        let kc = 0.5;
        let p = 1.0;
        let c = 1.0;
        let s = 1.0;

        let result_cel0 = cel0(kc, p, c, s);
        let result_cel = cel(array![kc], array![p], array![c], array![s]);

        assert_relative_eq!(result_cel0, result_cel[0], epsilon = 1e-10);
    }

    #[test]
    fn test_celv_multiple_values() {
        // Test vectorized CEL with multiple values
        let kc = array![0.5, 0.6, 0.7];
        let p = array![1.0, 1.2, 0.8];
        let c = array![1.0, 1.0, 1.0];
        let s = array![1.0, 1.0, 1.0];

        let result = celv(kc.clone(), p.clone(), c.clone(), s.clone());

        assert_eq!(result.len(), 3);
        for i in 0..3 {
            assert!(result[i].is_finite());
            // Check consistency with scalar version
            let scalar_result = cel0(kc[i], p[i], c[i], s[i]);
            assert_relative_eq!(result[i], scalar_result, epsilon = 1e-10);
        }
    }

    #[test]
    fn test_cel_small_array() {
        // For n < 10, cel should use scalar loop
        let kc = array![0.5, 0.6, 0.7, 0.8];
        let p = array![1.0, 1.2, 0.8, 1.5];
        let c = array![1.0, 1.0, 1.0, 1.0];
        let s = array![1.0, 1.0, 1.0, 1.0];

        let result = cel(kc.clone(), p.clone(), c.clone(), s.clone());

        assert_eq!(result.len(), 4);
        for i in 0..4 {
            assert!(result[i].is_finite());
        }
    }

    #[test]
    fn test_cel_large_array() {
        // For n >= 10, cel should use vectorized version
        let n = 15;
        let kc = Array1::from_vec(vec![0.5; n]);
        let p = Array1::from_vec(vec![1.0; n]);
        let c = Array1::from_vec(vec![1.0; n]);
        let s = Array1::from_vec(vec![1.0; n]);

        let result = cel(kc, p, c, s);

        assert_eq!(result.len(), n);
        for i in 0..n {
            assert!(result[i].is_finite());
        }
    }

    #[test]
    fn test_cel_consistency_scalar_vs_vectorized() {
        // Test that scalar and vectorized versions give same results
        let n = 20;
        let mut kc_vec = Vec::new();
        let mut p_vec = Vec::new();
        let mut c_vec = Vec::new();
        let mut s_vec = Vec::new();

        // Generate test values
        for i in 0..n {
            kc_vec.push(0.3 + i as f64 * 0.03);
            p_vec.push(0.5 + i as f64 * 0.05);
            c_vec.push(1.0);
            s_vec.push(1.0);
        }

        let kc = Array1::from_vec(kc_vec.clone());
        let p = Array1::from_vec(p_vec.clone());
        let c = Array1::from_vec(c_vec.clone());
        let s = Array1::from_vec(s_vec.clone());

        // Compute using dispatcher (will use vectorized for n=20)
        let result_vec = cel(kc, p, c, s);

        // Compute using scalar version
        for i in 0..n {
            let scalar_result = cel0(kc_vec[i], p_vec[i], c_vec[i], s_vec[i]);
            assert_relative_eq!(result_vec[i], scalar_result, epsilon = 1e-10);
        }
    }

    #[test]
    fn test_cel_iter_scalar() {
        // Test values derived from Python implementation
        // These are approximate values for validation
        let qc = 0.5;
        let p = 1.0;
        let g = 1.0;
        let cc = 1.0;
        let ss = 1.0;
        let em = 1.0;
        let kk = 0.25;

        let result = cel_iter_scalar(qc, p, g, cc, ss, em, kk);

        // Result should be finite and reasonable
        assert!(result.is_finite());
        assert!(result > 0.0);
    }

    #[test]
    fn test_cel_iter_vectorized() {
        // Test with small array
        let qc = array![0.5, 0.6];
        let p = array![1.0, 1.0];
        let g = array![1.0, 1.0];
        let cc = array![1.0, 1.0];
        let ss = array![1.0, 1.0];
        let em = array![1.0, 1.0];
        let kk = array![0.25, 0.36];

        let result = cel_iter(qc, p, g, cc, ss, em, kk);

        assert_eq!(result.len(), 2);
        assert!(result[0].is_finite());
        assert!(result[1].is_finite());
        assert!(result[0] > 0.0);
        assert!(result[1] > 0.0);
    }

    #[test]
    fn test_cel_iter_consistency() {
        // Test that vectorized and scalar versions give same results
        let qc_val = 0.5;
        let p_val = 1.0;
        let g_val = 1.0;
        let cc_val = 1.0;
        let ss_val = 1.0;
        let em_val = 1.0;
        let kk_val = 0.25;

        let scalar_result = cel_iter_scalar(qc_val, p_val, g_val, cc_val, ss_val, em_val, kk_val);

        let vec_result = cel_iter(
            array![qc_val],
            array![p_val],
            array![g_val],
            array![cc_val],
            array![ss_val],
            array![em_val],
            array![kk_val],
        );

        assert_relative_eq!(scalar_result, vec_result[0], epsilon = 1e-10);
    }
}
