//! B-field computation for polyline currents (line segments with current)
//!
//! Algorithm from Python magpylib field_BH_polyline.py
//!
//! Implements the Biot-Savart law for straight line segments carrying current.
//! Field computation via: http://www.phys.uri.edu/gerhard/PHY204/tsl216.pdf

use ndarray::{Array1, Array2, Axis};
use rayon::prelude::*;
use std::f64::consts::PI;

/// Permeability of free space (H/m)
const MU0: f64 = 1.25663706212e-6;

/// Tolerance for on-line detection
const ON_LINE_TOL: f64 = 1e-15;

/// Compute H-field for a single line segment carrying current
///
/// # Arguments
/// * `observer` - Observer position [x, y, z]
/// * `segment_start` - Segment start position [x, y, z]
/// * `segment_end` - Segment end position [x, y, z]
/// * `current` - Electrical current (Amperes)
///
/// # Returns
/// * H-field vector [Hx, Hy, Hz] in A/m, or None if observer is on the line
///
/// # Algorithm
/// Uses Biot-Savart law with dimensionless coordinates (normalized by segment length)
/// to avoid numerical issues with large/small inputs.
fn compute_segment_h_field(
    observer: &[f64; 3],
    segment_start: &[f64; 3],
    segment_end: &[f64; 3],
    current: f64,
) -> Option<[f64; 3]> {
    // Segment vector p1 -> p2
    let dx = segment_start[0] - segment_end[0];
    let dy = segment_start[1] - segment_end[1];
    let dz = segment_start[2] - segment_end[2];

    // Segment length (characteristic length scale)
    let norm_12 = (dx * dx + dy * dy + dz * dz).sqrt();

    if norm_12 < 1e-15 {
        // Zero-length segment
        return Some([0.0, 0.0, 0.0]);
    }

    // Make dimensionless by normalizing by segment length
    let p1 = [
        segment_start[0] / norm_12,
        segment_start[1] / norm_12,
        segment_start[2] / norm_12,
    ];
    let p2 = [
        segment_end[0] / norm_12,
        segment_end[1] / norm_12,
        segment_end[2] / norm_12,
    ];
    let po = [
        observer[0] / norm_12,
        observer[1] / norm_12,
        observer[2] / norm_12,
    ];

    // p4 = projection of observer onto line p1-p2
    // t = (po - p1) · (p1 - p2)
    let t = (po[0] - p1[0]) * (p1[0] - p2[0])
        + (po[1] - p1[1]) * (p1[1] - p2[1])
        + (po[2] - p1[2]) * (p1[2] - p2[2]);

    let p4 = [
        p1[0] + t * (p1[0] - p2[0]),
        p1[1] + t * (p1[1] - p2[1]),
        p1[2] + t * (p1[2] - p2[2]),
    ];

    // Distance of observer from line
    let do4_x = po[0] - p4[0];
    let do4_y = po[1] - p4[1];
    let do4_z = po[2] - p4[2];
    let norm_o4 = (do4_x * do4_x + do4_y * do4_y + do4_z * do4_z).sqrt();

    // If observer is on the line, field = 0
    if norm_o4 < ON_LINE_TOL {
        return None;  // Signal that point is on line
    }

    // Field direction: cross product (p2 - p1) × (po - p4)
    let cross_x = (p2[1] - p1[1]) * do4_z - (p2[2] - p1[2]) * do4_y;
    let cross_y = (p2[2] - p1[2]) * do4_x - (p2[0] - p1[0]) * do4_z;
    let cross_z = (p2[0] - p1[0]) * do4_y - (p2[1] - p1[1]) * do4_x;

    let norm_cross = (cross_x * cross_x + cross_y * cross_y + cross_z * cross_z).sqrt();

    if norm_cross < 1e-15 {
        // Degenerate case
        return Some([0.0, 0.0, 0.0]);
    }

    // Unit vector in field direction
    let eb_x = cross_x / norm_cross;
    let eb_y = cross_y / norm_cross;
    let eb_z = cross_z / norm_cross;

    // Compute angles (distances to p1 and p2)
    let do1_x = po[0] - p1[0];
    let do1_y = po[1] - p1[1];
    let do1_z = po[2] - p1[2];
    let norm_o1 = (do1_x * do1_x + do1_y * do1_y + do1_z * do1_z).sqrt();

    let do2_x = po[0] - p2[0];
    let do2_y = po[1] - p2[1];
    let do2_z = po[2] - p2[2];
    let norm_o2 = (do2_x * do2_x + do2_y * do2_y + do2_z * do2_z).sqrt();

    let d41_x = p4[0] - p1[0];
    let d41_y = p4[1] - p1[1];
    let d41_z = p4[2] - p1[2];
    let norm_41 = (d41_x * d41_x + d41_y * d41_y + d41_z * d41_z).sqrt();

    let d42_x = p4[0] - p2[0];
    let d42_y = p4[1] - p2[1];
    let d42_z = p4[2] - p2[2];
    let norm_42 = (d42_x * d42_x + d42_y * d42_y + d42_z * d42_z).sqrt();

    let sin_th1 = norm_41 / norm_o1;
    let sin_th2 = norm_42 / norm_o2;

    // Determine how p1, p2, p4 are sorted on the line (to get correct sign)
    let delta_sin = if norm_41 > 1.0 && norm_41 > norm_42 {
        // Both points below p4
        (sin_th1 - sin_th2).abs()
    } else if norm_42 > 1.0 && norm_42 > norm_41 {
        // Both points above p4
        (sin_th2 - sin_th1).abs()
    } else {
        // One above, one below, or one equals p4
        (sin_th1 + sin_th2).abs()
    };

    // H-field magnitude and direction
    // H = (delta_sin / norm_o4) * eB / norm_12 * current / (4π)
    let h_magnitude = (delta_sin / norm_o4) * current / (norm_12 * 4.0 * PI);

    Some([
        h_magnitude * eb_x,
        h_magnitude * eb_y,
        h_magnitude * eb_z,
    ])
}

/// Compute B-field for polyline current segments
///
/// # Arguments
/// * `observers` - Observer positions (n_points, 3)
/// * `segment_starts` - Segment start positions (n_segments, 3)
/// * `segment_ends` - Segment end positions (n_segments, 3)
/// * `currents` - Electrical currents for each segment (n_segments) in Amperes
///
/// # Returns
/// * B-field vectors (n_points, 3) in Tesla
///
/// # Notes
/// - Field is zero for points on the line segments
/// - Zero-length segments contribute zero field
/// - Uses Biot-Savart law with dimensionless normalization
pub fn compute_polyline_field(
    observers: &Array2<f64>,
    segment_starts: &Array2<f64>,
    segment_ends: &Array2<f64>,
    currents: &Array1<f64>,
) -> Array2<f64> {
    let n_points = observers.nrows();
    let n_segments = segment_starts.nrows();

    // Use parallel iteration over observers for performance
    // Collect results into a Vec first, then convert to Array2
    let field_results: Vec<[f64; 3]> = (0..n_points)
        .into_par_iter()
        .map(|i| {
            let observer = [
                observers[[i, 0]],
                observers[[i, 1]],
                observers[[i, 2]],
            ];

            let mut total_field = [0.0, 0.0, 0.0];

            // Sum contributions from all segments
            for j in 0..n_segments {
                let start = [
                    segment_starts[[j, 0]],
                    segment_starts[[j, 1]],
                    segment_starts[[j, 2]],
                ];
                let end = [
                    segment_ends[[j, 0]],
                    segment_ends[[j, 1]],
                    segment_ends[[j, 2]],
                ];
                let current = currents[j];

                // Compute H-field for this segment
                if let Some(h_field) = compute_segment_h_field(&observer, &start, &end, current) {
                    // Convert H to B (B = μ₀ * H) and accumulate
                    total_field[0] += h_field[0] * MU0;
                    total_field[1] += h_field[1] * MU0;
                    total_field[2] += h_field[2] * MU0;
                }
                // If None (on line), contribution is zero, so nothing to add
            }

            total_field
        })
        .collect();

    // Convert results back to Array2
    let mut b_field = Array2::zeros((n_points, 3));
    for (i, field) in field_results.iter().enumerate() {
        b_field[[i, 0]] = field[0];
        b_field[[i, 1]] = field[1];
        b_field[[i, 2]] = field[2];
    }

    b_field
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use ndarray::array;

    #[test]
    fn test_segment_on_line() {
        // Observer exactly on the line segment
        let observer = [0.5, 0.0, 0.0];
        let start = [0.0, 0.0, 0.0];
        let end = [1.0, 0.0, 0.0];
        let current = 1.0;

        let result = compute_segment_h_field(&observer, &start, &end, current);
        assert!(result.is_none());  // Should return None for on-line case
    }

    #[test]
    fn test_segment_zero_length() {
        // Zero-length segment
        let observer = [1.0, 1.0, 0.0];
        let start = [0.0, 0.0, 0.0];
        let end = [0.0, 0.0, 0.0];  // Same as start
        let current = 1.0;

        let result = compute_segment_h_field(&observer, &start, &end, current);
        assert!(result.is_some());
        let h = result.unwrap();
        assert_relative_eq!(h[0], 0.0, epsilon = 1e-10);
        assert_relative_eq!(h[1], 0.0, epsilon = 1e-10);
        assert_relative_eq!(h[2], 0.0, epsilon = 1e-10);
    }

    #[test]
    fn test_segment_perpendicular() {
        // Observer perpendicular to segment center
        // Segment: x-axis from (0,0,0) to (1,0,0)
        // Observer: (0.5, 1, 0) - perpendicular above center
        let observer = [0.5, 1.0, 0.0];
        let start = [0.0, 0.0, 0.0];
        let end = [1.0, 0.0, 0.0];
        let current = 1.0;

        let result = compute_segment_h_field(&observer, &start, &end, current);
        assert!(result.is_some());
        let h = result.unwrap();

        // Field should be in +z direction (right-hand rule)
        assert_relative_eq!(h[0], 0.0, epsilon = 1e-10);
        assert_relative_eq!(h[1], 0.0, epsilon = 1e-10);
        assert!(h[2] > 0.0);  // Positive z
    }

    #[test]
    fn test_polyline_single_segment() {
        // Test polyline with single segment
        let observers = array![[0.5, 1.0, 0.0]];
        let starts = array![[0.0, 0.0, 0.0]];
        let ends = array![[1.0, 0.0, 0.0]];
        let currents = array![1.0];

        let b_field = compute_polyline_field(&observers, &starts, &ends, &currents);

        // Should match single segment result (converted to B-field)
        let observer = [0.5, 1.0, 0.0];
        let start = [0.0, 0.0, 0.0];
        let end = [1.0, 0.0, 0.0];
        let h = compute_segment_h_field(&observer, &start, &end, 1.0).unwrap();

        assert_relative_eq!(b_field[[0, 0]], h[0] * MU0, epsilon = 1e-15);
        assert_relative_eq!(b_field[[0, 1]], h[1] * MU0, epsilon = 1e-15);
        assert_relative_eq!(b_field[[0, 2]], h[2] * MU0, epsilon = 1e-15);
    }

    #[test]
    fn test_polyline_multiple_segments() {
        // Test polyline with multiple segments (rectangular loop)
        let observers = array![[0.0, 0.0, 1.0]];  // Above center
        let starts = array![
            [1.0, 1.0, 0.0],   // Segment 1
            [-1.0, 1.0, 0.0],  // Segment 2
            [-1.0, -1.0, 0.0], // Segment 3
            [1.0, -1.0, 0.0],  // Segment 4
        ];
        let ends = array![
            [-1.0, 1.0, 0.0],  // Segment 1 end
            [-1.0, -1.0, 0.0], // Segment 2 end
            [1.0, -1.0, 0.0],  // Segment 3 end
            [1.0, 1.0, 0.0],   // Segment 4 end (closes loop)
        ];
        let currents = array![1.0, 1.0, 1.0, 1.0];

        let b_field = compute_polyline_field(&observers, &starts, &ends, &currents);

        // For a square loop with observer on axis, field should be primarily in z
        // and symmetric (x, y components should be small/zero by symmetry)
        assert!(b_field[[0, 0]].abs() < 1e-10);  // Symmetric -> ~0
        assert!(b_field[[0, 1]].abs() < 1e-10);  // Symmetric -> ~0
        assert!(b_field[[0, 2]] > 0.0);  // Positive z (field points up)
    }

    #[test]
    fn test_polyline_current_direction() {
        // Test that reversing current reverses field
        let observers = array![[0.5, 1.0, 0.0]];
        let starts = array![[0.0, 0.0, 0.0]];
        let ends = array![[1.0, 0.0, 0.0]];

        let currents_pos = array![1.0];
        let currents_neg = array![-1.0];

        let b_pos = compute_polyline_field(&observers, &starts, &ends, &currents_pos);
        let b_neg = compute_polyline_field(&observers, &starts, &ends, &currents_neg);

        // Fields should be opposite
        assert_relative_eq!(b_pos[[0, 0]], -b_neg[[0, 0]], epsilon = 1e-15);
        assert_relative_eq!(b_pos[[0, 1]], -b_neg[[0, 1]], epsilon = 1e-15);
        assert_relative_eq!(b_pos[[0, 2]], -b_neg[[0, 2]], epsilon = 1e-15);
    }
}
