//! Data export functionality for visualization
//!
//! This module provides utilities for exporting field computation results
//! to various formats for visualization with external tools.

pub mod vtk;

use ndarray::Array2;
use serde::{Deserialize, Serialize};
use std::fs::File;
use std::io::Write;
use std::path::Path;

/// Field data structure for export
///
/// Contains grid information and computed field values, suitable for
/// visualization with matplotlib or other plotting libraries.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FieldData {
    /// Description of the data
    pub description: String,

    /// Grid parameters
    pub grid: GridInfo,

    /// Observer positions (n_points, 3) - flattened row-major
    pub positions: Vec<Vec<f64>>,

    /// B-field values (n_points, 3) - flattened row-major
    /// Components: [Bx, By, Bz] in Tesla
    pub b_field: Vec<Vec<f64>>,

    /// Magnet geometry information (for overlay visualization)
    pub magnets: Vec<MagnetInfo>,
}

/// Grid information for structured data
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GridInfo {
    /// Which plane: "xy", "xz", or "yz"
    pub plane: String,

    /// Number of points in first axis
    pub nx: usize,

    /// Number of points in second axis
    pub ny: usize,

    /// Range for first axis (min, max)
    pub x_range: (f64, f64),

    /// Range for second axis (min, max)
    pub y_range: (f64, f64),

    /// Fixed coordinate value
    pub fixed_coord: f64,
}

/// Magnet geometry information for visualization
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MagnetInfo {
    /// Magnet type: "cuboid", "cylinder", "sphere"
    pub magnet_type: String,

    /// Position [x, y, z] in meters
    pub position: Vec<f64>,

    /// Dimensions (interpretation depends on type)
    /// Cuboid: [width, height, depth]
    /// Cylinder: [diameter, height]
    /// Sphere: [diameter]
    pub dimensions: Vec<f64>,

    /// Polarization [Jx, Jy, Jz] in Tesla
    pub polarization: Vec<f64>,
}

impl FieldData {
    /// Create a new FieldData structure
    pub fn new(
        description: String,
        grid: GridInfo,
        positions: &Array2<f64>,
        b_field: &Array2<f64>,
    ) -> Self {
        // Convert ndarray to Vec<Vec<f64>>
        let positions_vec: Vec<Vec<f64>> = (0..positions.nrows())
            .map(|i| vec![positions[[i, 0]], positions[[i, 1]], positions[[i, 2]]])
            .collect();

        let b_field_vec: Vec<Vec<f64>> = (0..b_field.nrows())
            .map(|i| vec![b_field[[i, 0]], b_field[[i, 1]], b_field[[i, 2]]])
            .collect();

        Self {
            description,
            grid,
            positions: positions_vec,
            b_field: b_field_vec,
            magnets: Vec::new(),
        }
    }

    /// Add magnet information for visualization overlay
    pub fn add_magnet(&mut self, magnet: MagnetInfo) {
        self.magnets.push(magnet);
    }

    /// Export to JSON file
    ///
    /// # Arguments
    /// * `path` - Output file path
    ///
    /// # Returns
    /// Result indicating success or error
    ///
    /// # Example
    /// ```no_run
    /// use magrslib::export::{FieldData, GridInfo};
    /// use ndarray::Array2;
    ///
    /// let grid = GridInfo {
    ///     plane: "xy".to_string(),
    ///     nx: 10,
    ///     ny: 10,
    ///     x_range: (-0.05, 0.05),
    ///     y_range: (-0.05, 0.05),
    ///     fixed_coord: 0.0,
    /// };
    ///
    /// let positions = Array2::zeros((100, 3));
    /// let b_field = Array2::zeros((100, 3));
    ///
    /// let data = FieldData::new(
    ///     "Test data".to_string(),
    ///     grid,
    ///     &positions,
    ///     &b_field,
    /// );
    ///
    /// data.export_to_json("output.json").unwrap();
    /// ```
    pub fn export_to_json<P: AsRef<Path>>(&self, path: P) -> std::io::Result<()> {
        let json = serde_json::to_string_pretty(self)
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;

        let mut file = File::create(path)?;
        file.write_all(json.as_bytes())?;

        Ok(())
    }

    /// Load from JSON file
    pub fn load_from_json<P: AsRef<Path>>(path: P) -> std::io::Result<Self> {
        let file = File::open(path)?;
        let data = serde_json::from_reader(file)
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;
        Ok(data)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::array;

    #[test]
    fn test_field_data_creation() {
        let grid = GridInfo {
            plane: "xy".to_string(),
            nx: 2,
            ny: 2,
            x_range: (-1.0, 1.0),
            y_range: (-1.0, 1.0),
            fixed_coord: 0.0,
        };

        let positions = array![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let b_field = array![[0.1, 0.2, 0.3], [0.4, 0.5, 0.6]];

        let data = FieldData::new("Test".to_string(), grid, &positions, &b_field);

        assert_eq!(data.positions.len(), 2);
        assert_eq!(data.b_field.len(), 2);
        assert_eq!(data.positions[0], vec![0.0, 0.0, 0.0]);
        assert_eq!(data.b_field[0], vec![0.1, 0.2, 0.3]);
    }

    #[test]
    fn test_add_magnet() {
        let grid = GridInfo {
            plane: "xy".to_string(),
            nx: 2,
            ny: 2,
            x_range: (-1.0, 1.0),
            y_range: (-1.0, 1.0),
            fixed_coord: 0.0,
        };

        let positions = array![[0.0, 0.0, 0.0]];
        let b_field = array![[0.1, 0.2, 0.3]];

        let mut data = FieldData::new("Test".to_string(), grid, &positions, &b_field);

        data.add_magnet(MagnetInfo {
            magnet_type: "cuboid".to_string(),
            position: vec![0.0, 0.0, 0.0],
            dimensions: vec![0.01, 0.01, 0.01],
            polarization: vec![0.0, 0.0, 1.0],
        });

        assert_eq!(data.magnets.len(), 1);
        assert_eq!(data.magnets[0].magnet_type, "cuboid");
    }
}
