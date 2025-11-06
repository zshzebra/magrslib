//! Magnetic field computation
//!
//! This module provides the main API for computing magnetic fields from sources.
//! Field computation happens in source-local coordinates, with coordinate transforms
//! applied automatically.

pub mod currents;
pub mod magnets;
mod transforms;

use crate::sources::{AnySource, FieldFuncId, Source};
use crate::types::Observers;
use ndarray::{Array1, Array2};
use rayon::prelude::*;

pub use transforms::{transform_field_to_global, transform_observers_to_local};

/// Compute B-field from sources at observer positions
///
/// This is the main user-facing API. It handles coordinate transformations
/// and dispatches to appropriate field computation functions.
///
/// # Arguments
/// * `sources` - Collection of magnetic field sources
/// * `observers` - Positions where fields should be evaluated
///
/// # Returns
/// B-field vectors at observer positions (n_observers, 3) in Tesla
///
/// # Algorithm
/// 1. Initialize result array with zeros
/// 2. For each source:
///    a. Transform observers to source-local coordinates
///    b. Call appropriate field function
///    c. Transform field back to global coordinates
///    d. Accumulate into result
/// 3. Return total field
///
/// # Example
/// ```no_run
/// use magrslib::{get_b, Observers};
/// use magrslib::sources::magnets::Cuboid;
///
/// let magnet = Cuboid::builder()
///     .dimension([0.01, 0.01, 0.01])
///     .polarization([0.0, 0.0, 1.0])
///     .position([0.0, 0.0, 0.0])
///     .build()
///     .unwrap();
///
/// let observers = Observers::from_vec3s(vec![[0.0, 0.0, 0.02]]);
/// let b_field = get_b(&[magnet.into()], &observers);
/// ```
pub fn get_b(sources: &[AnySource], observers: &Observers) -> Array2<f64> {
    get_b_at_time(sources, observers, 0.0)
}

/// Compute B-field from sources at observer positions at specific time
///
/// This enables time-varying current sources and dynamic field analysis.
/// For static sources, time parameter is ignored.
///
/// # Arguments
/// * `sources` - Collection of magnetic field sources
/// * `observers` - Positions where fields should be evaluated
/// * `time` - Time in seconds for evaluating time-varying currents
///
/// # Returns
/// B-field vectors at observer positions (n_observers, 3) in Tesla
pub fn get_b_at_time(sources: &[AnySource], observers: &Observers, time: f64) -> Array2<f64> {
    let n_obs = observers.len();

    // Parallel processing across sources - each source is computed independently
    // then the results are summed
    let field_components: Vec<Array2<f64>> = sources
        .par_iter()
        .map(|source| {
            // Get source properties
            let path = source.path();
            let orientation = source.orientation();

            // For now, only use the first position in the path
            // (full path/animation support can be added later)
            let position = path.get(0).expect("Source path must have at least one position");

            // Transform observers to source-local coordinates
            let local_observers = transforms::transform_observers_to_local(
                observers.as_array(),
                position,
                orientation,
            );

            // Compute field in source-local coordinates
            let local_field = compute_field_for_source_at_time(source, &local_observers, time);

            // Transform field back to global coordinates
            transforms::transform_field_to_global(&local_field, orientation)
        })
        .collect();

    // Sum all field contributions
    let mut total_field = Array2::zeros((n_obs, 3));
    for field in field_components {
        total_field = total_field + field;
    }

    total_field
}

/// Compute B-field time series from sources at observer positions
///
/// Enables analysis of field evolution with time-varying currents such as
/// PWM signals, switching transients, and AC waveforms.
///
/// # Arguments
/// * `sources` - Collection of magnetic field sources
/// * `observers` - Positions where fields should be evaluated
/// * `time_start` - Start time in seconds
/// * `time_end` - End time in seconds
/// * `time_step` - Time step in seconds
///
/// # Returns
/// Vector of (time, field) pairs where field is (n_observers, 3) in Tesla
///
/// # Example
/// ```no_run
/// use magrslib::{get_b_time_series, Observers};
/// use magrslib::sources::currents::Polyline;
/// use magrslib::sources::current_functions::PWMCurrent;
///
/// // PWM solenoid simulation
/// let pwm = Box::new(PWMCurrent::new(1000.0, 0.5, 10.0)); // 1kHz, 50%, 10A
/// let solenoid = Polyline::builder()
///     .vertices(helical_vertices)
///     .current_function(pwm)
///     .build()
///     .unwrap();
///
/// let observers = Observers::from_vec3s(vec![[0.0, 0.0, 0.0]]);
/// let time_series = get_b_time_series(
///     &[solenoid.into()],
///     &observers,
///     0.0,    // Start time
///     0.002,  // End time (2ms = 2 PWM cycles)
///     0.00001 // Time step (10μs)
/// );
/// ```
pub fn get_b_time_series(
    sources: &[AnySource],
    observers: &Observers,
    time_start: f64,
    time_end: f64,
    time_step: f64,
) -> Vec<(f64, Array2<f64>)> {
    let num_steps = ((time_end - time_start) / time_step).ceil() as usize + 1;
    let times: Vec<f64> = (0..num_steps)
        .map(|i| time_start + i as f64 * time_step)
        .collect();

    // Parallel computation across time steps
    times
        .par_iter()
        .map(|&time| {
            let field = get_b_at_time(sources, observers, time);
            (time, field)
        })
        .collect()
}

/// Compute field for a single source in source-local coordinates
///
/// Dispatches to the appropriate field computation function based on source type.
/// Backward compatibility wrapper that calls compute_field_for_source_at_time with time=0.
///
/// # Arguments
/// * `source` - The source object
/// * `observers` - Observer positions in source-local coordinates (n, 3)
///
/// # Returns
/// Field vectors in source-local coordinates (n, 3)
fn compute_field_for_source(source: &AnySource, observers: &Array2<f64>) -> Array2<f64> {
    compute_field_for_source_at_time(source, observers, 0.0)
}

/// Compute field for a single source in source-local coordinates at specific time
///
/// Dispatches to the appropriate field computation function based on source type.
/// Supports time-varying current sources.
///
/// # Arguments
/// * `source` - The source object
/// * `observers` - Observer positions in source-local coordinates (n, 3)
/// * `time` - Time in seconds for evaluating time-varying currents
///
/// # Returns
/// Field vectors in source-local coordinates (n, 3)
fn compute_field_for_source_at_time(source: &AnySource, observers: &Array2<f64>, time: f64) -> Array2<f64> {
    let n = observers.nrows();
    let props = source.get_properties();

    match source.field_func_id() {
        FieldFuncId::MagnetCuboid => {
            // Extract cuboid parameters
            let dimension = props.dimension.expect("Cuboid must have dimension");
            let polarization = props.polarization.expect("Cuboid must have polarization");

            // Field function expects (n, 3) arrays for dimensions and polarizations
            // We need to repeat the source parameters for each observer
            let dimensions = tile_row(&dimension, n);
            let polarizations = tile_row(&polarization, n);

            magnets::cuboid::magnet_cuboid_bfield(&observers, &dimensions, &polarizations)
        }

        FieldFuncId::MagnetCylinder => {
            // TODO: Implement cylinder field computation
            unimplemented!("Cylinder field computation not yet implemented")
        }

        FieldFuncId::MagnetSphere => {
            // TODO: Implement sphere field computation
            unimplemented!("Sphere field computation not yet implemented")
        }

        FieldFuncId::CurrentCircle => {
            // Extract Circle source properties
            let props = source.get_properties();
            let diameter = props.diameter.expect("Circle source missing diameter");
            let current = props.evaluate_current(time);

            currents::circle::compute_circle_field(&observers, diameter, current)
        }

        FieldFuncId::CurrentPolyline => {
            let vertices = props.vertices.as_ref().expect("Polyline source missing vertices");
            let current = props.evaluate_current(time);

            // Convert consecutive vertices into line segments
            if vertices.len() < 2 {
                // Not enough vertices for segments
                return Array2::zeros((observers.nrows(), 3));
            }

            let n_segments = vertices.len() - 1;
            let mut segment_starts = Array2::zeros((n_segments, 3));
            let mut segment_ends = Array2::zeros((n_segments, 3));
            let mut currents = Array1::zeros(n_segments);

            for i in 0..n_segments {
                // Start of segment i
                segment_starts[[i, 0]] = vertices[i][0];
                segment_starts[[i, 1]] = vertices[i][1];
                segment_starts[[i, 2]] = vertices[i][2];

                // End of segment i (start of segment i+1)
                segment_ends[[i, 0]] = vertices[i + 1][0];
                segment_ends[[i, 1]] = vertices[i + 1][1];
                segment_ends[[i, 2]] = vertices[i + 1][2];

                // All segments carry the same current
                currents[i] = current;
            }

            currents::polyline::compute_polyline_field(
                &observers,
                &segment_starts,
                &segment_ends,
                &currents,
            )
        }

        FieldFuncId::Dipole => {
            // TODO: Implement dipole field computation
            unimplemented!("Dipole field computation not yet implemented")
        }
    }
}

/// Tile a Vec3 into an (n, 3) array by repeating it n times
///
/// This is equivalent to numpy's `np.tile(vec, (n, 1))`
///
/// # Arguments
/// * `vec` - Input vector [x, y, z]
/// * `n` - Number of times to repeat
///
/// # Returns
/// Array with shape (n, 3) where each row is a copy of `vec`
fn tile_row(vec: &[f64; 3], n: usize) -> Array2<f64> {
    let mut result = Array2::zeros((n, 3));
    for i in 0..n {
        result[[i, 0]] = vec[0];
        result[[i, 1]] = vec[1];
        result[[i, 2]] = vec[2];
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sources::magnets::Cuboid;
    use crate::types::Observers;
    use approx::assert_relative_eq;

    #[test]
    fn test_get_b_single_cuboid() {
        // Simple test: cuboid at origin, observer on z-axis
        let magnet = Cuboid::builder()
            .dimension([0.01, 0.01, 0.01])
            .polarization([0.0, 0.0, 1.0])
            .position([0.0, 0.0, 0.0])
            .build()
            .unwrap();

        let observers = Observers::from_vec3s(vec![[0.0, 0.0, 0.02]]);
        let b_field = get_b(&[magnet.into()], &observers);

        // Field should be finite
        assert!(b_field[[0, 0]].is_finite());
        assert!(b_field[[0, 1]].is_finite());
        assert!(b_field[[0, 2]].is_finite());

        // Due to symmetry, Bx and By should be nearly zero
        assert_relative_eq!(b_field[[0, 0]], 0.0, epsilon = 1e-10);
        assert_relative_eq!(b_field[[0, 1]], 0.0, epsilon = 1e-10);

        // Bz should be positive
        assert!(b_field[[0, 2]] > 0.0);
    }

    #[test]
    fn test_get_b_multiple_observers() {
        // Test with multiple observer positions
        let magnet = Cuboid::builder()
            .dimension([0.01, 0.01, 0.01])
            .polarization([0.0, 0.0, 1.0])
            .position([0.0, 0.0, 0.0])
            .build()
            .unwrap();

        let observers = Observers::from_vec3s(vec![
            [0.0, 0.0, 0.02],
            [0.0, 0.0, -0.02],
        ]);
        let b_field = get_b(&[magnet.into()], &observers);

        assert_eq!(b_field.shape(), &[2, 3]);

        // Both fields should be finite
        assert!(b_field[[0, 2]].is_finite());
        assert!(b_field[[1, 2]].is_finite());

        // Due to the symmetry mapping algorithm in cuboid computation (direct port from Python),
        // both points map to the same computed position. This is the expected behavior.
        assert_relative_eq!(b_field[[0, 2]], b_field[[1, 2]], epsilon = 1e-10);
    }

    #[test]
    fn test_get_b_multiple_sources() {
        // Test field superposition from multiple sources
        let magnet1 = Cuboid::builder()
            .dimension([0.01, 0.01, 0.01])
            .polarization([0.0, 0.0, 1.0])
            .position([0.0, 0.0, 0.01])
            .build()
            .unwrap();

        let magnet2 = Cuboid::builder()
            .dimension([0.01, 0.01, 0.01])
            .polarization([0.0, 0.0, 1.0])
            .position([0.0, 0.0, -0.01])
            .build()
            .unwrap();

        let observers = Observers::from_vec3s(vec![[0.0, 0.0, 0.05]]);

        // Clone magnet1 before using it since .into() consumes the value
        let b_field = get_b(&[magnet1.clone().into(), magnet2.into()], &observers);

        // Field should be finite
        assert!(b_field[[0, 0]].is_finite());
        assert!(b_field[[0, 1]].is_finite());
        assert!(b_field[[0, 2]].is_finite());

        // Due to symmetry, field from two identical magnets should be stronger
        // than from a single magnet
        let b_single = get_b(&[magnet1.into()], &observers);
        assert!(b_field[[0, 2]].abs() > b_single[[0, 2]].abs());
    }

    #[test]
    fn test_tile_row() {
        let vec = [1.0, 2.0, 3.0];
        let tiled = tile_row(&vec, 3);

        assert_eq!(tiled.shape(), &[3, 3]);
        assert_eq!(tiled[[0, 0]], 1.0);
        assert_eq!(tiled[[0, 1]], 2.0);
        assert_eq!(tiled[[0, 2]], 3.0);
        assert_eq!(tiled[[1, 0]], 1.0);
        assert_eq!(tiled[[2, 2]], 3.0);
    }
}
