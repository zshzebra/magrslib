//! VTK file format export for 3D visualization
//!
//! Exports magnetic field data using the vtkio crate,
//! compatible with PyVista, Paraview, and other VTK-based tools.

use ndarray::Array2;
use std::io;
use std::path::Path;
use vtkio::model::*;
use vtkio::IOBuffer;

/// VTK structured grid information
#[derive(Debug, Clone)]
pub struct VtkGrid3D {
    /// Number of points in x, y, z directions
    pub dimensions: (usize, usize, usize),

    /// Origin point (x0, y0, z0)
    pub origin: (f64, f64, f64),

    /// Spacing between points (dx, dy, dz)
    pub spacing: (f64, f64, f64),
}

impl VtkGrid3D {
    /// Create VTK grid from ranges
    pub fn from_ranges(
        x_range: (f64, f64, usize),
        y_range: (f64, f64, usize),
        z_range: (f64, f64, usize),
    ) -> Self {
        let (x_start, x_end, nx) = x_range;
        let (y_start, y_end, ny) = y_range;
        let (z_start, z_end, nz) = z_range;

        let dx = if nx > 1 { (x_end - x_start) / ((nx - 1) as f64) } else { 0.0 };
        let dy = if ny > 1 { (y_end - y_start) / ((ny - 1) as f64) } else { 0.0 };
        let dz = if nz > 1 { (z_end - z_start) / ((nz - 1) as f64) } else { 0.0 };

        Self {
            dimensions: (nx, ny, nz),
            origin: (x_start, y_start, z_start),
            spacing: (dx, dy, dz),
        }
    }

    /// Total number of points
    pub fn n_points(&self) -> usize {
        self.dimensions.0 * self.dimensions.1 * self.dimensions.2
    }
}

/// Export B-field data to VTK format for 3D visualization
///
/// Creates a VTK ImageData (structured grid) file with B-field vector data.
/// Uses the vtkio crate for robust VTK file generation.
///
/// # Arguments
/// * `path` - Output file path
/// * `grid` - Grid information (dimensions, origin, spacing)
/// * `b_field` - B-field vectors (n_points, 3) in Tesla
/// * `description` - Optional description for the file header
///
/// # Example
/// ```no_run
/// use magrslib::export::vtk::{VtkGrid3D, export_bfield_to_vtk};
/// use ndarray::Array2;
///
/// let grid = VtkGrid3D::from_ranges(
///     (-0.02, 0.02, 21),
///     (-0.02, 0.02, 21),
///     (-0.02, 0.02, 21),
/// );
///
/// let b_field = Array2::zeros((21*21*21, 3));
///
/// export_bfield_to_vtk(
///     "bfield_3d.vtk",
///     &grid,
///     &b_field,
///     Some("B-field from cuboid magnet")
/// ).unwrap();
/// ```
pub fn export_bfield_to_vtk<P: AsRef<Path>>(
    path: P,
    grid: &VtkGrid3D,
    b_field: &Array2<f64>,
    description: Option<&str>,
) -> io::Result<()> {
    // Validate dimensions
    if b_field.nrows() != grid.n_points() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!(
                "B-field has {} points but grid expects {}",
                b_field.nrows(),
                grid.n_points()
            ),
        ));
    }

    if b_field.ncols() != 3 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "B-field must have 3 components (Bx, By, Bz)",
        ));
    }

    // Convert B-field to flat Vec for vtkio
    let mut b_vectors = Vec::with_capacity(b_field.nrows() * 3);
    for i in 0..b_field.nrows() {
        b_vectors.push(b_field[[i, 0]] as f32); // Bx
        b_vectors.push(b_field[[i, 1]] as f32); // By
        b_vectors.push(b_field[[i, 2]] as f32); // Bz
    }

    // Create VTK file with ImageData dataset
    let vtk_file = Vtk {
        version: Version::new((4, 2)),
        title: description.unwrap_or("B-field data").to_string(),
        byte_order: ByteOrder::BigEndian,
        data: DataSet::ImageData {
            extent: Extent::Dims([
                grid.dimensions.0 as u32,
                grid.dimensions.1 as u32,
                grid.dimensions.2 as u32,
            ]),
            origin: [
                grid.origin.0 as f32,
                grid.origin.1 as f32,
                grid.origin.2 as f32,
            ],
            spacing: [
                grid.spacing.0 as f32,
                grid.spacing.1 as f32,
                grid.spacing.2 as f32,
            ],
            meta: None,
            pieces: vec![Piece::Inline(Box::new(ImageDataPiece {
                extent: Extent::Dims([
                    grid.dimensions.0 as u32,
                    grid.dimensions.1 as u32,
                    grid.dimensions.2 as u32,
                ]),
                data: Attributes {
                    point: vec![Attribute::DataArray(DataArray {
                        name: "B".to_string(),
                        elem: ElementType::Vectors,
                        data: IOBuffer::F32(b_vectors),
                    })],
                    cell: vec![],
                },
            }))],
        },
        file_path: None,
    };

    // Write to file
    vtk_file
        .export(path)
        .map_err(|e| io::Error::new(io::ErrorKind::Other, e))?;

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::array;

    #[test]
    fn test_vtk_grid_creation() {
        let grid = VtkGrid3D::from_ranges(
            (-0.02, 0.02, 5),
            (-0.01, 0.01, 3),
            (0.0, 0.04, 9),
        );

        assert_eq!(grid.dimensions, (5, 3, 9));
        assert_eq!(grid.origin, (-0.02, -0.01, 0.0));
        assert_eq!(grid.n_points(), 5 * 3 * 9);

        // Check spacing
        assert!((grid.spacing.0 - 0.01).abs() < 1e-10);
        assert!((grid.spacing.1 - 0.01).abs() < 1e-10);
        assert!((grid.spacing.2 - 0.005).abs() < 1e-10);
    }

    #[test]
    fn test_vtk_export() {
        let grid = VtkGrid3D::from_ranges(
            (-0.01, 0.01, 3),
            (-0.01, 0.01, 3),
            (-0.01, 0.01, 3),
        );

        // Create dummy field data (3x3x3 = 27 points)
        let mut b_field = Array2::zeros((27, 3));
        for i in 0..27 {
            b_field[[i, 0]] = i as f64 * 0.001;
            b_field[[i, 1]] = i as f64 * 0.002;
            b_field[[i, 2]] = i as f64 * 0.003;
        }

        // Export to temp file
        let temp_path = "/tmp/test_bfield.vtk";
        let result = export_bfield_to_vtk(
            temp_path,
            &grid,
            &b_field,
            Some("Test B-field"),
        );

        assert!(result.is_ok());

        // Verify file exists and has content
        let metadata = std::fs::metadata(temp_path);
        assert!(metadata.is_ok());
        assert!(metadata.unwrap().len() > 0);
    }
}
