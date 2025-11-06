//! Unified visualization API for magnetic field sources and field data
//!
//! This module provides a high-level, fluent API for creating comprehensive
//! 3D visualizations of magnetic field systems. It handles:
//! - Multiple source types (magnets, current loops, etc.)
//! - Adaptive grid generation based on source geometry
//! - Export of both field data and source geometry to VTK
//! - Metadata generation for Python visualization scripts
//!
//! # Example
//! ```no_run
//! use magrslib::export::visualization::VisualizationBuilder;
//! use magrslib::sources::magnets::Cuboid;
//! use magrslib::sources::currents::Circle;
//!
//! let magnet = Cuboid::builder()
//!     .dimension([0.01, 0.01, 0.01])
//!     .polarization([0.0, 0.0, 1.0])
//!     .build().unwrap();
//!
//! let current_loop = Circle::builder()
//!     .diameter(0.02)
//!     .current(1.0)
//!     .position([0.0, 0.0, 0.02])
//!     .build().unwrap();
//!
//! VisualizationBuilder::new()
//!     .add_source(magnet)
//!     .add_source(current_loop)
//!     .grid_resolution(64)
//!     .export_to_directory("./output")?;
//! ```

use crate::export::vtk::{export_bfield_to_vtk, export_source_geometries_to_vtk, VtkGrid3D};
use crate::get_b;
use crate::sources::{AnySource, Source};
use crate::types::Observers;
use ndarray::Array2;
use serde::{Deserialize, Serialize};
use std::fs;
use std::io;
use std::path::Path;

/// Metadata describing a visualization scene
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct VisualizationMetadata {
    /// Sources in the scene with their properties
    pub sources: Vec<SourceMetadata>,
    /// Grid configuration
    pub grid: GridMetadata,
    /// Field statistics
    pub field_stats: FieldStatistics,
    /// File paths for generated data
    pub files: FileMetadata,
    /// Recommended visualization settings
    pub visualization_hints: VisualizationHints,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SourceMetadata {
    /// Source type (e.g., "MagnetCuboid", "CurrentCircle")
    pub source_type: String,
    /// Position in global coordinates
    pub position: [f64; 3],
    /// Orientation quaternion (w, x, y, z)
    pub orientation: [f64; 4],
    /// Source-specific parameters
    pub parameters: serde_json::Value,
    /// Characteristic size for visualization
    pub characteristic_size: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GridMetadata {
    /// Grid bounds: [x_min, x_max, y_min, y_max, z_min, z_max]
    pub bounds: [f64; 6],
    /// Grid resolution [nx, ny, nz]
    pub resolution: [usize; 3],
    /// Grid spacing [dx, dy, dz]
    pub spacing: [f64; 3],
    /// Total number of grid points
    pub n_points: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FieldStatistics {
    /// Field magnitude statistics
    pub magnitude: StatisticsData,
    /// Component statistics [Bx, By, Bz]
    pub components: [StatisticsData; 3],
    /// Number of invalid (NaN/Inf) field values
    pub invalid_count: usize,
    /// Total computation time in seconds
    pub computation_time: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct StatisticsData {
    pub min: f64,
    pub max: f64,
    pub mean: f64,
    pub median: f64,
    pub std_dev: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FileMetadata {
    /// Relative path to VTK field data file
    pub field_vtk: String,
    /// Relative path to VTK geometry file (if generated)
    pub geometry_vtk: Option<String>,
    /// Relative path to this metadata file
    pub metadata_json: String,
    /// Absolute path to VTK field data file
    pub field_vtk_absolute: String,
    /// Absolute path to VTK geometry file (if generated)
    pub geometry_vtk_absolute: Option<String>,
    /// Absolute path to this metadata file
    pub metadata_json_absolute: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct VisualizationHints {
    /// Recommended camera position [x, y, z]
    pub camera_position: [f64; 3],
    /// Recommended camera focal point [x, y, z]
    pub camera_focal_point: [f64; 3],
    /// Recommended field magnitude colormap range [min, max]
    pub colormap_range: [f64; 2],
    /// Suggested streamline seed strategy
    pub streamline_strategy: String,
    /// Scene scale factor for visualization
    pub scene_scale: f64,
}

/// Builder for creating comprehensive field visualizations
#[derive(Debug)]
pub struct VisualizationBuilder {
    sources: Vec<AnySource>,
    grid_bounds: Option<([f64; 2], [f64; 2], [f64; 2])>,
    grid_resolution: usize,
    adaptive_resolution: bool,
    export_geometry: bool,
    title: Option<String>,
}

impl VisualizationBuilder {
    /// Create a new visualization builder
    pub fn new() -> Self {
        Self {
            sources: Vec::new(),
            grid_bounds: None,
            grid_resolution: 64,
            adaptive_resolution: true,
            export_geometry: true,
            title: None,
        }
    }

    /// Add a source to the visualization
    pub fn add_source<S: Into<AnySource>>(mut self, source: S) -> Self {
        self.sources.push(source.into());
        self
    }

    /// Add multiple sources to the visualization
    pub fn add_sources<I, S>(mut self, sources: I) -> Self
    where
        I: IntoIterator<Item = S>,
        S: Into<AnySource>,
    {
        self.sources.extend(sources.into_iter().map(|s| s.into()));
        self
    }

    /// Set explicit grid bounds
    ///
    /// # Arguments
    /// * `x_range` - (min, max) bounds in x direction
    /// * `y_range` - (min, max) bounds in y direction
    /// * `z_range` - (min, max) bounds in z direction
    pub fn grid_bounds(
        mut self,
        x_range: (f64, f64),
        y_range: (f64, f64),
        z_range: (f64, f64),
    ) -> Self {
        self.grid_bounds = Some((
            [x_range.0, x_range.1],
            [y_range.0, y_range.1],
            [z_range.0, z_range.1],
        ));
        self
    }

    /// Set grid resolution (number of points per dimension)
    pub fn grid_resolution(mut self, resolution: usize) -> Self {
        self.grid_resolution = resolution;
        self
    }

    /// Enable/disable adaptive resolution based on source sizes
    pub fn adaptive_resolution(mut self, enabled: bool) -> Self {
        self.adaptive_resolution = enabled;
        self
    }

    /// Enable/disable geometry export
    pub fn export_geometry(mut self, enabled: bool) -> Self {
        self.export_geometry = enabled;
        self
    }

    /// Set visualization title
    pub fn title<S: Into<String>>(mut self, title: S) -> Self {
        self.title = Some(title.into());
        self
    }

    /// Export visualization to a directory
    ///
    /// Creates:
    /// - `field.vtk` - B-field data
    /// - `geometry.vtk` - Source geometry (if enabled)
    /// - `metadata.json` - Scene metadata
    pub fn export_to_directory<P: AsRef<Path>>(
        self,
        output_dir: P,
    ) -> io::Result<VisualizationMetadata> {
        let output_dir = output_dir.as_ref();

        // Create output directory if it doesn't exist
        if !output_dir.exists() {
            fs::create_dir_all(output_dir)?;
        }

        // Generate grid bounds if not specified
        let grid_bounds = self
            .grid_bounds
            .unwrap_or_else(|| self.compute_adaptive_bounds());

        // Create grid
        let vtk_grid = VtkGrid3D::from_ranges(
            (grid_bounds.0[0], grid_bounds.0[1], self.grid_resolution),
            (grid_bounds.1[0], grid_bounds.1[1], self.grid_resolution),
            (grid_bounds.2[0], grid_bounds.2[1], self.grid_resolution),
        );

        // Generate observer positions
        let observers = Observers::from_grid_3d(
            (grid_bounds.0[0], grid_bounds.0[1], self.grid_resolution),
            (grid_bounds.1[0], grid_bounds.1[1], self.grid_resolution),
            (grid_bounds.2[0], grid_bounds.2[1], self.grid_resolution),
        );

        // Compute B-field
        println!("Computing B-field at {} points...", observers.len());
        let start_time = std::time::Instant::now();
        let b_field = get_b(&self.sources, &observers);
        let computation_time = start_time.elapsed().as_secs_f64();
        println!("Field computation completed in {:.2}s", computation_time);

        // Clean field data (replace NaN/Inf with zeros for visualization)
        let b_field_clean = self.clean_field_data(&b_field);

        // Compute field statistics
        let field_stats = self.compute_field_statistics(&b_field, computation_time);

        // Export field data
        let field_vtk_path = output_dir.join("field.vtk");
        let title = self
            .title
            .as_deref()
            .unwrap_or("Magnetic field visualization");
        export_bfield_to_vtk(&field_vtk_path, &vtk_grid, &b_field_clean, Some(title))?;

        // Export geometry if requested
        let geometry_vtk_path = if self.export_geometry {
            let geometry_path = output_dir.join("geometry.vtk");
            self.export_source_geometry(&geometry_path)?;
            Some("geometry.vtk".to_string())
        } else {
            None
        };

        // Define file paths
        let metadata_path = output_dir.join("metadata.json");

        // Generate metadata
        let metadata = VisualizationMetadata {
            sources: self.generate_source_metadata(),
            grid: GridMetadata {
                bounds: [
                    grid_bounds.0[0],
                    grid_bounds.0[1],
                    grid_bounds.1[0],
                    grid_bounds.1[1],
                    grid_bounds.2[0],
                    grid_bounds.2[1],
                ],
                resolution: [self.grid_resolution; 3],
                spacing: [vtk_grid.spacing.0, vtk_grid.spacing.1, vtk_grid.spacing.2],
                n_points: vtk_grid.n_points(),
            },
            field_stats,
            files: FileMetadata {
                field_vtk: "field.vtk".to_string(),
                geometry_vtk: geometry_vtk_path.clone(),
                metadata_json: "metadata.json".to_string(),
                field_vtk_absolute: field_vtk_path.to_string_lossy().to_string(),
                geometry_vtk_absolute: geometry_vtk_path.as_ref().map(|_| {
                    output_dir.join("geometry.vtk").to_string_lossy().to_string()
                }),
                metadata_json_absolute: metadata_path.to_string_lossy().to_string(),
            },
            visualization_hints: self.generate_visualization_hints(&grid_bounds, &b_field_clean),
        };

        // Save metadata
        let metadata_json = serde_json::to_string_pretty(&metadata)?;
        fs::write(metadata_path, metadata_json)?;

        println!("Visualization exported to: {}", output_dir.display());
        println!("  - Field data: field.vtk");
        if metadata.files.geometry_vtk.is_some() {
            println!("  - Geometry: geometry.vtk");
        }
        println!("  - Metadata: metadata.json");

        Ok(metadata)
    }

    /// Compute adaptive grid bounds based on source positions and sizes
    fn compute_adaptive_bounds(&self) -> ([f64; 2], [f64; 2], [f64; 2]) {
        if self.sources.is_empty() {
            return ([-0.02, 0.02], [-0.02, 0.02], [-0.02, 0.02]);
        }

        let mut x_min = f64::INFINITY;
        let mut x_max = f64::NEG_INFINITY;
        let mut y_min = f64::INFINITY;
        let mut y_max = f64::NEG_INFINITY;
        let mut z_min = f64::INFINITY;
        let mut z_max = f64::NEG_INFINITY;

        let mut max_source_size: f64 = 0.0;

        for source in &self.sources {
            // Get source position (use first position in path)
            if let Some(position) = source.path().get(0) {
                let pos = position.coords;

                // Estimate source size based on properties
                let source_size = self.estimate_source_size(source);
                max_source_size = max_source_size.max(source_size);

                // Update bounds including source extent
                x_min = x_min.min(pos[0] - source_size);
                x_max = x_max.max(pos[0] + source_size);
                y_min = y_min.min(pos[1] - source_size);
                y_max = y_max.max(pos[1] + source_size);
                z_min = z_min.min(pos[2] - source_size);
                z_max = z_max.max(pos[2] + source_size);
            }
        }

        // Add padding (3x max source size)
        let padding = 3.0 * max_source_size;
        (
            [x_min - padding, x_max + padding],
            [y_min - padding, y_max + padding],
            [z_min - padding, z_max + padding],
        )
    }

    /// Estimate characteristic size of a source for adaptive bounds
    fn estimate_source_size(&self, source: &AnySource) -> f64 {
        let props = source.get_properties();

        // For magnets, use dimension
        if let Some(dim) = props.dimension {
            return dim.iter().fold(0.0_f64, |acc, &x| acc.max(x)) * 0.5;
        }

        // For current loops, use diameter
        if let Some(diameter) = props.diameter {
            return diameter * 0.5;
        }

        // Default fallback
        0.01
    }

    /// Clean field data by replacing NaN/Inf with zeros
    fn clean_field_data(&self, b_field: &Array2<f64>) -> Array2<f64> {
        let mut cleaned = b_field.clone();
        for i in 0..cleaned.nrows() {
            for j in 0..3 {
                if !cleaned[[i, j]].is_finite() {
                    cleaned[[i, j]] = 0.0;
                }
            }
        }
        cleaned
    }

    /// Compute comprehensive field statistics
    fn compute_field_statistics(
        &self,
        b_field: &Array2<f64>,
        computation_time: f64,
    ) -> FieldStatistics {
        let mut invalid_count = 0;
        let mut magnitudes = Vec::new();
        let mut components = [Vec::new(), Vec::new(), Vec::new()];

        for i in 0..b_field.nrows() {
            let bx = b_field[[i, 0]];
            let by = b_field[[i, 1]];
            let bz = b_field[[i, 2]];

            // Check for invalid values
            if !bx.is_finite() || !by.is_finite() || !bz.is_finite() {
                invalid_count += 1;
                continue;
            }

            // Collect valid values
            components[0].push(bx);
            components[1].push(by);
            components[2].push(bz);
            magnitudes.push((bx * bx + by * by + bz * bz).sqrt());
        }

        FieldStatistics {
            magnitude: Self::compute_statistics(&magnitudes),
            components: [
                Self::compute_statistics(&components[0]),
                Self::compute_statistics(&components[1]),
                Self::compute_statistics(&components[2]),
            ],
            invalid_count,
            computation_time,
        }
    }

    /// Compute statistics for a vector of values
    fn compute_statistics(values: &[f64]) -> StatisticsData {
        if values.is_empty() {
            return StatisticsData {
                min: 0.0,
                max: 0.0,
                mean: 0.0,
                median: 0.0,
                std_dev: 0.0,
            };
        }

        let mut sorted_values = values.to_vec();
        sorted_values.sort_by(|a, b| a.partial_cmp(b).unwrap());

        let min = sorted_values[0];
        let max = sorted_values[sorted_values.len() - 1];
        let mean = sorted_values.iter().sum::<f64>() / sorted_values.len() as f64;
        let median = sorted_values[sorted_values.len() / 2];

        let variance = sorted_values
            .iter()
            .map(|x| (x - mean).powi(2))
            .sum::<f64>()
            / sorted_values.len() as f64;
        let std_dev = variance.sqrt();

        StatisticsData {
            min,
            max,
            mean,
            median,
            std_dev,
        }
    }

    /// Generate metadata for all sources in the scene
    fn generate_source_metadata(&self) -> Vec<SourceMetadata> {
        self.sources
            .iter()
            .map(|source| {
                let position = if let Some(pos) = source.path().get(0) {
                    pos.coords
                } else {
                    [0.0, 0.0, 0.0]
                };

                let orientation = {
                    let o = source.orientation();
                    [
                        o.quaternion()[0],
                        o.quaternion()[1],
                        o.quaternion()[2],
                        o.quaternion()[3],
                    ]
                };

                let props = source.get_properties();
                // Create a simplified parameters object for serialization
                let parameters = serde_json::json!({
                    "dimension": props.dimension,
                    "polarization": props.polarization,
                    "diameter": props.diameter,
                    "current": props.current,
                    "vertices": props.vertices
                });

                SourceMetadata {
                    source_type: format!("{:?}", source.field_func_id()),
                    position,
                    orientation,
                    parameters,
                    characteristic_size: self.estimate_source_size(source),
                }
            })
            .collect()
    }

    /// Generate visualization hints for optimal rendering
    fn generate_visualization_hints(
        &self,
        grid_bounds: &([f64; 2], [f64; 2], [f64; 2]),
        b_field: &Array2<f64>,
    ) -> VisualizationHints {
        // Compute scene center and size
        let scene_center = [
            (grid_bounds.0[0] + grid_bounds.0[1]) * 0.5,
            (grid_bounds.1[0] + grid_bounds.1[1]) * 0.5,
            (grid_bounds.2[0] + grid_bounds.2[1]) * 0.5,
        ];

        let scene_size = [
            grid_bounds.0[1] - grid_bounds.0[0],
            grid_bounds.1[1] - grid_bounds.1[0],
            grid_bounds.2[1] - grid_bounds.2[0],
        ];

        let max_scene_size = scene_size.iter().fold(0.0_f64, |acc, &x| acc.max(x));

        // Camera position: offset from center at 45-degree angle
        let camera_distance = max_scene_size * 1.5;
        let camera_position = [
            scene_center[0] + camera_distance * 0.7071,
            scene_center[1] + camera_distance * 0.7071,
            scene_center[2] + camera_distance * 0.5,
        ];

        // Compute field magnitude range for colormap
        let magnitudes: Vec<f64> = (0..b_field.nrows())
            .map(|i| {
                let bx = b_field[[i, 0]];
                let by = b_field[[i, 1]];
                let bz = b_field[[i, 2]];
                (bx * bx + by * by + bz * bz).sqrt()
            })
            .filter(|&m| m.is_finite() && m > 0.0)
            .collect();

        let (colormap_min, colormap_max) = if magnitudes.is_empty() {
            (0.0, 1.0)
        } else {
            let min_val = magnitudes.iter().fold(f64::INFINITY, |acc, &x| acc.min(x));
            let max_val = magnitudes.iter().fold(0.0_f64, |acc, &x| acc.max(x));
            (min_val, max_val)
        };

        VisualizationHints {
            camera_position,
            camera_focal_point: scene_center,
            colormap_range: [colormap_min, colormap_max],
            streamline_strategy: "adaptive".to_string(),
            scene_scale: max_scene_size,
        }
    }

    /// Export source geometry to VTK (placeholder for now)
    fn export_source_geometry<P: AsRef<Path>>(&self, path: P) -> io::Result<()> {
        if self.sources.is_empty() {
            return Ok(()); // No sources to export
        }

        let title = format!("Source geometries for {} sources", self.sources.len());

        export_source_geometries_to_vtk(path, &self.sources, Some(&title))?;

        println!("  ✓ Exported {} source geometries", self.sources.len());
        Ok(())
    }
}

impl Default for VisualizationBuilder {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sources::currents::Circle;
    use crate::sources::magnets::Cuboid;

    #[test]
    fn test_visualization_builder() {
        let magnet = Cuboid::builder()
            .dimension([0.01, 0.01, 0.01])
            .polarization([0.0, 0.0, 1.0])
            .build()
            .unwrap();

        let builder = VisualizationBuilder::new()
            .add_source(magnet)
            .grid_resolution(8) // Small for testing
            .title("Test visualization");

        // Test that builder accepts sources and configurations
        assert_eq!(builder.sources.len(), 1);
        assert_eq!(builder.grid_resolution, 8);
        assert_eq!(builder.title, Some("Test visualization".to_string()));
    }

    #[test]
    fn test_adaptive_bounds() {
        let magnet = Cuboid::builder()
            .dimension([0.01, 0.01, 0.01])
            .polarization([0.0, 0.0, 1.0])
            .position([0.02, 0.03, 0.04])
            .build()
            .unwrap();

        let builder = VisualizationBuilder::new().add_source(magnet);
        let bounds = builder.compute_adaptive_bounds();

        // Bounds should include the magnet position with padding
        assert!(bounds.0[0] < 0.02); // x_min < magnet x
        assert!(bounds.0[1] > 0.02); // x_max > magnet x
        assert!(bounds.1[0] < 0.03); // y_min < magnet y
        assert!(bounds.1[1] > 0.03); // y_max > magnet y
        assert!(bounds.2[0] < 0.04); // z_min < magnet z
        assert!(bounds.2[1] > 0.04); // z_max > magnet z
    }

    #[test]
    fn test_statistics_computation() {
        let values = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let stats = VisualizationBuilder::compute_statistics(&values);

        assert_eq!(stats.min, 1.0);
        assert_eq!(stats.max, 5.0);
        assert_eq!(stats.mean, 3.0);
        assert_eq!(stats.median, 3.0);
        assert!((stats.std_dev - 1.4142135623730951).abs() < 1e-10);
    }
}
