//! Complete elliptic integrals using Bulirsch's CEL algorithm
//!
//! Direct port from Python magpylib special_cel.py
//! Used for computing magnetic fields from circular current loops.

use ndarray::Array1;
use std::f64::consts::PI;

/// Half pi constant (π/2)
const HALF_PI: f64 = PI / 2.0;

/// Tolerance for convergence
const TOLERANCE: f64 = 1e-8;

/// Maximum iterations to prevent infinite loops
const MAX_ITER: usize = 100;

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
            if (g[i] - qc[i]).abs() >= qc[i] * TOLERANCE {
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

    while (g - qc).abs() >= qc * TOLERANCE && iter_count < MAX_ITER {
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
