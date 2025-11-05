//! Magnetic field from homogeneously magnetized cuboid
//!
//! Direct port from Python magpylib field_BH_cuboid.py
//!
//! References:
//! - Yang, Supercond. Sci. Technol. 3(12):591, 1990
//! - Engel-Herbert, J. Appl. Phys. 97(7):074504-4, 2005
//! - Cichon, IEEE Sensors Journal 19(7):2509, 2019

use crate::utils::FOUR_PI;
use ndarray::{Array1, Array2, Axis};

/// Compute B-field from homogeneously magnetized cuboids
///
/// Cuboid sides are parallel to coordinate axes, geometric center at origin.
/// Uses magnetic surface charge density method.
///
/// # Arguments
/// * `observers` - Observer positions (n, 3) in meters
/// * `dimensions` - Cuboid dimensions (n, 3) as (a, b, c) in meters
/// * `polarizations` - Magnetic polarization (n, 3) as (Jx, Jy, Jz) in Tesla
///
/// # Returns
/// B-field at observer positions (n, 3) in Tesla
///
/// # Algorithm
/// To avoid indeterminate forms at positions along edge extensions, we exploit
/// symmetry by mapping all positions to the bottom-negative-x/positive-y/positive-z
/// octant (bottom Q4), then apply appropriate sign corrections.
pub fn magnet_cuboid_bfield(
    observers: &Array2<f64>,
    dimensions: &Array2<f64>,
    polarizations: &Array2<f64>,
) -> Array2<f64> {
    let n = observers.nrows();

    // Extract polarization components
    let pol_x = polarizations.column(0).to_owned();
    let pol_y = polarizations.column(1).to_owned();
    let pol_z = polarizations.column(2).to_owned();

    // Half-dimensions
    let a = &dimensions.column(0) / 2.0;
    let b = &dimensions.column(1) / 2.0;
    let c = &dimensions.column(2) / 2.0;

    // Observer coordinates (will be modified for symmetry)
    let mut x = observers.column(0).to_owned();
    let mut y = observers.column(1).to_owned();
    let mut z = observers.column(2).to_owned();

    // === Symmetry mapping to bottom Q4 ===
    // Create masks for coordinate signs
    let maskx: Vec<bool> = x.iter().map(|&xi| xi < 0.0).collect();
    let masky: Vec<bool> = y.iter().map(|&yi| yi > 0.0).collect();
    let maskz: Vec<bool> = z.iter().map(|&zi| zi > 0.0).collect();

    // Flip coordinates to map to bottom Q4
    for i in 0..n {
        if maskx[i] {
            x[i] = -x[i];
        }
        if masky[i] {
            y[i] = -y[i];
        }
        if maskz[i] {
            z[i] = -z[i];
        }
    }

    // Create sign flip arrays (3x3 matrices per observer)
    // qsigns[i] is a 3x3 matrix for observer i
    let mut qsigns = vec![[[1.0_f64; 3]; 3]; n];

    let qs_flipx = [[1.0, -1.0, -1.0], [-1.0, 1.0, 1.0], [-1.0, 1.0, 1.0]];
    let qs_flipy = [[1.0, -1.0, 1.0], [-1.0, 1.0, -1.0], [1.0, -1.0, 1.0]];
    let qs_flipz = [[1.0, 1.0, -1.0], [1.0, 1.0, -1.0], [-1.0, -1.0, 1.0]];

    // Apply sign flips
    for i in 0..n {
        if maskx[i] {
            for j in 0..3 {
                for k in 0..3 {
                    qsigns[i][j][k] *= qs_flipx[j][k];
                }
            }
        }
        if masky[i] {
            for j in 0..3 {
                for k in 0..3 {
                    qsigns[i][j][k] *= qs_flipy[j][k];
                }
            }
        }
        if maskz[i] {
            for j in 0..3 {
                for k in 0..3 {
                    qsigns[i][j][k] *= qs_flipz[j][k];
                }
            }
        }
    }

    // === Field computation ===
    // Compute relative positions to corners
    let xma = &x - &a;
    let xpa = &x + &a;
    let ymb = &y - &b;
    let ypb = &y + &b;
    let zmc = &z - &c;
    let zpc = &z + &c;

    // Squares
    let xma2 = &xma * &xma;
    let xpa2 = &xpa * &xpa;
    let ymb2 = &ymb * &ymb;
    let ypb2 = &ypb * &ypb;
    let zmc2 = &zmc * &zmc;
    let zpc2 = &zpc * &zpc;

    // Distances to 8 corners (naming: m=minus, p=plus)
    let mmm = (&xma2 + &ymb2 + &zmc2).mapv(f64::sqrt);
    let pmp = (&xpa2 + &ymb2 + &zpc2).mapv(f64::sqrt);
    let pmm = (&xpa2 + &ymb2 + &zmc2).mapv(f64::sqrt);
    let mmp = (&xma2 + &ymb2 + &zpc2).mapv(f64::sqrt);
    let mpm = (&xma2 + &ypb2 + &zmc2).mapv(f64::sqrt);
    let ppp = (&xpa2 + &ypb2 + &zpc2).mapv(f64::sqrt);
    let ppm = (&xpa2 + &ypb2 + &zmc2).mapv(f64::sqrt);
    let mpp = (&xma2 + &ypb2 + &zpc2).mapv(f64::sqrt);

    // Logarithmic terms (ff2x, ff2y, ff2z)
    let ff2x = ((&xma + &mmm) * (&xpa + &ppm) * (&xpa + &pmp) * (&xma + &mpp)
        / ((&xpa + &pmm) * (&xma + &mpm) * (&xma + &mmp) * (&xpa + &ppp)))
        .mapv(|v| v.ln());

    let ff2y = (((-&ymb + &mmm) * (-&ypb + &ppm) * (-&ymb + &pmp) * (-&ypb + &mpp))
        / ((-&ymb + &pmm) * (-&ypb + &mpm) * (&ymb - &mmp) * (&ypb - &ppp)))
        .mapv(|v| v.ln());

    let ff2z = (((-&zmc + &mmm) * (-&zmc + &ppm) * (-&zpc + &pmp) * (-&zpc + &mpp))
        / ((-&zmc + &pmm) * (&zmc - &mpm) * (-&zpc + &mmp) * (&zpc - &ppp)))
        .mapv(|v| v.ln());

    // Arctangent terms (ff1x, ff1y, ff1z)
    let mut ff1x = Array1::zeros(n);
    for i in 0..n {
        ff1x[i] = f64::atan2(ymb[i] * zmc[i], xma[i] * mmm[i])
            - f64::atan2(ymb[i] * zmc[i], xpa[i] * pmm[i])
            - f64::atan2(ypb[i] * zmc[i], xma[i] * mpm[i])
            + f64::atan2(ypb[i] * zmc[i], xpa[i] * ppm[i])
            - f64::atan2(ymb[i] * zpc[i], xma[i] * mmp[i])
            + f64::atan2(ymb[i] * zpc[i], xpa[i] * pmp[i])
            + f64::atan2(ypb[i] * zpc[i], xma[i] * mpp[i])
            - f64::atan2(ypb[i] * zpc[i], xpa[i] * ppp[i]);
    }

    let mut ff1y = Array1::zeros(n);
    for i in 0..n {
        ff1y[i] = f64::atan2(xma[i] * zmc[i], ymb[i] * mmm[i])
            - f64::atan2(xpa[i] * zmc[i], ymb[i] * pmm[i])
            - f64::atan2(xma[i] * zmc[i], ypb[i] * mpm[i])
            + f64::atan2(xpa[i] * zmc[i], ypb[i] * ppm[i])
            - f64::atan2(xma[i] * zpc[i], ymb[i] * mmp[i])
            + f64::atan2(xpa[i] * zpc[i], ymb[i] * pmp[i])
            + f64::atan2(xma[i] * zpc[i], ypb[i] * mpp[i])
            - f64::atan2(xpa[i] * zpc[i], ypb[i] * ppp[i]);
    }

    let mut ff1z = Array1::zeros(n);
    for i in 0..n {
        ff1z[i] = f64::atan2(xma[i] * ymb[i], zmc[i] * mmm[i])
            - f64::atan2(xpa[i] * ymb[i], zmc[i] * pmm[i])
            - f64::atan2(xma[i] * ypb[i], zmc[i] * mpm[i])
            + f64::atan2(xpa[i] * ypb[i], zmc[i] * ppm[i])
            - f64::atan2(xma[i] * ymb[i], zpc[i] * mmp[i])
            + f64::atan2(xpa[i] * ymb[i], zpc[i] * pmp[i])
            + f64::atan2(xma[i] * ypb[i], zpc[i] * mpp[i])
            - f64::atan2(xpa[i] * ypb[i], zpc[i] * ppp[i]);
    }

    // === Combine contributions from each polarization component ===
    // Extract sign corrections
    let qs00: Array1<f64> = (0..n).map(|i| qsigns[i][0][0]).collect();
    let qs01: Array1<f64> = (0..n).map(|i| qsigns[i][0][1]).collect();
    let qs02: Array1<f64> = (0..n).map(|i| qsigns[i][0][2]).collect();
    let qs10: Array1<f64> = (0..n).map(|i| qsigns[i][1][0]).collect();
    let qs11: Array1<f64> = (0..n).map(|i| qsigns[i][1][1]).collect();
    let qs12: Array1<f64> = (0..n).map(|i| qsigns[i][1][2]).collect();
    let qs20: Array1<f64> = (0..n).map(|i| qsigns[i][2][0]).collect();
    let qs21: Array1<f64> = (0..n).map(|i| qsigns[i][2][1]).collect();
    let qs22: Array1<f64> = (0..n).map(|i| qsigns[i][2][2]).collect();

    // Contributions from x-polarization
    let bx_pol_x = &pol_x * &ff1x * &qs00;
    let by_pol_x = &pol_x * &ff2z * &qs01;
    let bz_pol_x = &pol_x * &ff2y * &qs02;

    // Contributions from y-polarization
    let bx_pol_y = &pol_y * &ff2z * &qs10;
    let by_pol_y = &pol_y * &ff1y * &qs11;
    let bz_pol_y = -&pol_y * &ff2x * &qs12;

    // Contributions from z-polarization
    let bx_pol_z = &pol_z * &ff2y * &qs20;
    let by_pol_z = -&pol_z * &ff2x * &qs21;
    let bz_pol_z = &pol_z * &ff1z * &qs22;

    // Sum all contributions
    let bx_tot = bx_pol_x + bx_pol_y + bx_pol_z;
    let by_tot = by_pol_x + by_pol_y + by_pol_z;
    let bz_tot = bz_pol_x + bz_pol_y + bz_pol_z;

    // Stack into (n, 3) array
    let mut result = Array2::zeros((n, 3));
    result.column_mut(0).assign(&bx_tot);
    result.column_mut(1).assign(&by_tot);
    result.column_mut(2).assign(&bz_tot);

    // Divide by 4π
    result / FOUR_PI
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use ndarray::array;

    #[test]
    fn test_cuboid_on_axis() {
        // Test a simple case: observer on z-axis, z-polarization
        let observers = array![[0.0, 0.0, 0.02]];
        let dimensions = array![[0.01, 0.01, 0.01]];
        let polarizations = array![[0.0, 0.0, 1.0]];

        let b_field = magnet_cuboid_bfield(&observers, &dimensions, &polarizations);

        // Should get field in -z direction (outside magnet, north pole)
        assert_eq!(b_field.shape(), &[1, 3]);
        assert!(b_field[[0, 2]].is_finite());

        // Bx and By should be nearly zero due to symmetry
        assert_relative_eq!(b_field[[0, 0]], 0.0, epsilon = 1e-10);
        assert_relative_eq!(b_field[[0, 1]], 0.0, epsilon = 1e-10);
    }

    #[test]
    fn test_cuboid_symmetry() {
        // Test that field has correct symmetry
        let dimensions = array![[0.01, 0.01, 0.01]];
        let polarizations = array![[0.0, 0.0, 1.0]];

        // Point above and below
        let obs_above = array![[0.0, 0.0, 0.02]];
        let obs_below = array![[0.0, 0.0, -0.02]];

        let b_above = magnet_cuboid_bfield(&obs_above, &dimensions, &polarizations);
        let b_below = magnet_cuboid_bfield(&obs_below, &dimensions, &polarizations);

        // Due to symmetry exploitation in the algorithm, both points
        // give the same result - this is expected behavior from the Python implementation
        assert_relative_eq!(b_above[[0, 2]], b_below[[0, 2]], epsilon = 1e-10);
    }
}
