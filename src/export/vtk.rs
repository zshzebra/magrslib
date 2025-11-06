//! VTK file format export for 3D visualization
//!
//! Exports magnetic field data and source geometries using the vtkio crate,
//! compatible with PyVista, Paraview, and other VTK-based tools.

use crate::sources::{AnySource, Source};
use crate::types::Vec3;
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

/// Mesh geometry for VTK export
#[derive(Debug, Clone)]
pub struct VtkGeometry {
    /// Points in 3D space [x, y, z, x, y, z, ...]
    pub points: Vec<f64>,
    /// Cell connectivity (flat array)
    pub cell_connectivity: Vec<u32>,
    /// Cell types
    pub cell_types: Vec<CellType>,
    /// Source type labels for each cell
    pub source_types: Vec<String>,
    /// Source IDs for each cell
    pub source_ids: Vec<u32>,
}

impl VtkGeometry {
    pub fn new() -> Self {
        Self {
            points: Vec::new(),
            cell_connectivity: Vec::new(),
            cell_types: Vec::new(),
            source_types: Vec::new(),
            source_ids: Vec::new(),
        }
    }

    /// Add points to the geometry
    pub fn add_points(&mut self, new_points: &[f64]) {
        self.points.extend_from_slice(new_points);
    }

    /// Add a cell to the geometry
    pub fn add_cell(&mut self, cell_type: CellType, connectivity: &[u32], source_type: String, source_id: u32) {
        // Add number of vertices for this cell, then the vertex indices
        self.cell_connectivity.push(connectivity.len() as u32);
        self.cell_connectivity.extend_from_slice(connectivity);

        self.cell_types.push(cell_type);
        self.source_types.push(source_type);
        self.source_ids.push(source_id);
    }
}

/// Create a cuboid mesh for VTK export
fn create_cuboid_mesh(
    position: Vec3,
    dimension: Vec3,
    source_type: String,
    source_id: u32,
) -> VtkGeometry {
    let mut geometry = VtkGeometry::new();

    // Create 8 vertices of a cuboid centered at position
    let (a, b, c) = (dimension[0] / 2.0, dimension[1] / 2.0, dimension[2] / 2.0);
    let (x, y, z) = (position[0], position[1], position[2]);

    #[rustfmt::skip]
    let vertices = vec![
        x - a, y - b, z - c,  // 0: min corner
        x + a, y - b, z - c,  // 1
        x + a, y + b, z - c,  // 2
        x - a, y + b, z - c,  // 3
        x - a, y - b, z + c,  // 4
        x + a, y - b, z + c,  // 5
        x + a, y + b, z + c,  // 6: max corner
        x - a, y + b, z + c,  // 7
    ];

    geometry.add_points(&vertices);

    // Create hexahedron cell (8 vertices)
    let hex_connectivity = vec![0, 1, 2, 3, 4, 5, 6, 7];
    geometry.add_cell(CellType::Hexahedron, &hex_connectivity, source_type, source_id);

    geometry
}

/// Create a cylinder mesh for VTK export
fn create_cylinder_mesh(
    position: Vec3,
    diameter: f64,
    height: f64,
    n_sides: usize,
    source_type: String,
    source_id: u32,
) -> VtkGeometry {
    let mut geometry = VtkGeometry::new();
    let radius = diameter / 2.0;
    let half_height = height / 2.0;
    let (x, y, z) = (position[0], position[1], position[2]);

    // Create vertices for top and bottom circles
    let mut vertices = Vec::new();

    // Bottom circle center
    vertices.extend_from_slice(&[x, y, z - half_height]);
    // Top circle center
    vertices.extend_from_slice(&[x, y, z + half_height]);

    // Bottom circle vertices
    for i in 0..n_sides {
        let theta = 2.0 * std::f64::consts::PI * i as f64 / n_sides as f64;
        vertices.extend_from_slice(&[
            x + radius * theta.cos(),
            y + radius * theta.sin(),
            z - half_height,
        ]);
    }

    // Top circle vertices
    for i in 0..n_sides {
        let theta = 2.0 * std::f64::consts::PI * i as f64 / n_sides as f64;
        vertices.extend_from_slice(&[
            x + radius * theta.cos(),
            y + radius * theta.sin(),
            z + half_height,
        ]);
    }

    geometry.add_points(&vertices);

    // Create cells for cylinder faces
    // Bottom face
    let mut bottom_face = Vec::new();
    for i in 0..n_sides {
        bottom_face.push(2 + i as u32); // Bottom circle vertices
    }
    geometry.add_cell(CellType::Polygon, &bottom_face, source_type.clone(), source_id);

    // Top face
    let mut top_face = Vec::new();
    for i in 0..n_sides {
        top_face.push(2 + n_sides as u32 + i as u32); // Top circle vertices
    }
    geometry.add_cell(CellType::Polygon, &top_face, source_type.clone(), source_id);

    // Side faces (quads)
    for i in 0..n_sides {
        let next_i = (i + 1) % n_sides;
        let quad = vec![
            2 + i as u32,                          // Bottom current
            2 + next_i as u32,                     // Bottom next
            2 + n_sides as u32 + next_i as u32,    // Top next
            2 + n_sides as u32 + i as u32,         // Top current
        ];
        geometry.add_cell(CellType::Quad, &quad, source_type.clone(), source_id);
    }

    geometry
}

/// Create a sphere mesh for VTK export (using icosphere subdivision)
fn create_sphere_mesh(
    position: Vec3,
    diameter: f64,
    source_type: String,
    source_id: u32,
) -> VtkGeometry {
    let mut geometry = VtkGeometry::new();
    let radius = diameter / 2.0;
    let (x, y, z) = (position[0], position[1], position[2]);

    // Simple octahedron-based sphere (8 triangular faces)
    // For production, you'd want to use proper sphere tessellation
    let sqrt3 = 3.0_f64.sqrt();

    #[rustfmt::skip]
    let vertices = vec![
        x + radius, y,         z,          // 0: +X
        x - radius, y,         z,          // 1: -X
        x,          y + radius, z,         // 2: +Y
        x,          y - radius, z,         // 3: -Y
        x,          y,         z + radius, // 4: +Z
        x,          y,         z - radius, // 5: -Z
    ];

    geometry.add_points(&vertices);

    // Create triangular faces of octahedron
    let faces = vec![
        [0, 2, 4], [0, 4, 3], [0, 3, 5], [0, 5, 2], // +X vertex
        [1, 4, 2], [1, 3, 4], [1, 5, 3], [1, 2, 5], // -X vertex
    ];

    for face in faces {
        geometry.add_cell(CellType::Triangle, &face, source_type.clone(), source_id);
    }

    geometry
}

/// Create a circle (current loop) mesh for VTK export
fn create_circle_mesh(
    position: Vec3,
    diameter: f64,
    n_segments: usize,
    source_type: String,
    source_id: u32,
) -> VtkGeometry {
    let mut geometry = VtkGeometry::new();
    let radius = diameter / 2.0;
    let (x, y, z) = (position[0], position[1], position[2]);

    // Create vertices around the circle in xy-plane
    let mut vertices = Vec::new();
    for i in 0..n_segments {
        let theta = 2.0 * std::f64::consts::PI * i as f64 / n_segments as f64;
        vertices.extend_from_slice(&[
            x + radius * theta.cos(),
            y + radius * theta.sin(),
            z,
        ]);
    }

    geometry.add_points(&vertices);

    // Create line segments
    for i in 0..n_segments {
        let next_i = (i + 1) % n_segments;
        let line_connectivity = vec![i as u32, next_i as u32];
        geometry.add_cell(CellType::Line, &line_connectivity, source_type.clone(), source_id);
    }

    geometry
}

/// Create geometry mesh for an electromagnetic source
pub fn create_source_geometry(source: &AnySource, source_id: u32) -> VtkGeometry {
    let position = if source.path().len() > 0 {
        source.path().get(0).unwrap().coords
    } else {
        [0.0, 0.0, 0.0] // Default position
    };
    let properties = source.get_properties();
    let source_type = format!("{:?}", source.field_func_id());

    match source.field_func_id() {
        crate::sources::FieldFuncId::MagnetCuboid => {
            if let Some(dimension) = properties.dimension {
                create_cuboid_mesh(position, dimension, source_type, source_id)
            } else {
                VtkGeometry::new()
            }
        }
        crate::sources::FieldFuncId::MagnetCylinder => {
            if let (Some(diameter), Some(height)) = (properties.diameter, properties.height) {
                create_cylinder_mesh(position, diameter, height, 16, source_type, source_id)
            } else {
                VtkGeometry::new()
            }
        }
        crate::sources::FieldFuncId::MagnetSphere => {
            if let Some(diameter) = properties.diameter {
                create_sphere_mesh(position, diameter, source_type, source_id)
            } else {
                VtkGeometry::new()
            }
        }
        crate::sources::FieldFuncId::CurrentCircle => {
            if let Some(diameter) = properties.diameter {
                create_circle_mesh(position, diameter, 32, source_type, source_id)
            } else {
                VtkGeometry::new()
            }
        }
        crate::sources::FieldFuncId::CurrentPolyline => {
            // For polylines, create line segments between vertices
            if let Some(vertices) = properties.vertices {
                let mut geometry = VtkGeometry::new();
                let mut points = Vec::new();

                // Add all vertices
                for vertex in &vertices {
                    points.extend_from_slice(&[vertex[0], vertex[1], vertex[2]]);
                }
                geometry.add_points(&points);

                // Create line segments
                for i in 0..(vertices.len() - 1) {
                    let line_connectivity = vec![i as u32, (i + 1) as u32];
                    geometry.add_cell(CellType::Line, &line_connectivity, source_type.clone(), source_id);
                }

                geometry
            } else {
                VtkGeometry::new()
            }
        }
        crate::sources::FieldFuncId::Dipole => {
            // Create a small arrow/cone for dipole
            let mut geometry = VtkGeometry::new();
            let scale = 0.001; // 1mm scale

            // Simple arrow as line segment
            let points = vec![
                position[0], position[1], position[2] - scale/2.0,
                position[0], position[1], position[2] + scale/2.0,
            ];
            geometry.add_points(&points);

            let line_connectivity = vec![0, 1];
            geometry.add_cell(CellType::Line, &line_connectivity, source_type, source_id);

            geometry
        }
    }
}

/// Export source geometries to VTK format
pub fn export_source_geometries_to_vtk<P: AsRef<Path>>(
    path: P,
    sources: &[AnySource],
    description: Option<&str>,
) -> io::Result<()> {
    let mut all_points = Vec::new();
    let mut all_connectivity = Vec::new();
    let mut all_cell_types = Vec::new();
    let mut all_source_ids = Vec::new();
    let mut point_offset = 0u32;

    // Combine all source geometries
    for (source_id, source) in sources.iter().enumerate() {
        let geometry = create_source_geometry(source, source_id as u32);

        // Add points
        all_points.extend_from_slice(&geometry.points);

        // Add cells with offset point indices
        let mut offset_connectivity = geometry.cell_connectivity.clone();
        // Offset all vertex indices in the connectivity array
        for i in 0..offset_connectivity.len() {
            // Skip the count entries (every num_vertices + 1 entries)
            if i == 0 || offset_connectivity[i] as usize <= offset_connectivity.len() {
                // This is likely a vertex index, not a count
                if i > 0 && offset_connectivity[i] < 1000 { // Simple heuristic to distinguish indices from counts
                    offset_connectivity[i] += point_offset;
                }
            }
        }
        all_connectivity.extend(offset_connectivity);

        // Add cell types and metadata
        all_cell_types.extend(geometry.cell_types);
        all_source_ids.extend(geometry.source_ids);

        // Update point offset
        point_offset += (geometry.points.len() / 3) as u32;
    }

    if all_points.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "No geometry data to export",
        ));
    }

    // Convert points to Vec<f32> for vtkio
    let points_f32: Vec<f32> = all_points.iter().map(|&x| x as f32).collect();

    // Convert source IDs to cell data
    let source_id_data: Vec<u32> = all_source_ids;

    // Create VTK unstructured grid
    let vtk_file = Vtk {
        version: Version::new((4, 2)),
        title: description.unwrap_or("Source geometries").to_string(),
        byte_order: ByteOrder::BigEndian,
        data: DataSet::UnstructuredGrid {
            meta: None,
            pieces: vec![Piece::Inline(Box::new(UnstructuredGridPiece {
                points: IOBuffer::F32(points_f32),
                cells: Cells {
                    cell_verts: VertexNumbers::Legacy {
                        num_cells: all_cell_types.len() as u32,
                        vertices: all_connectivity,
                    },
                    types: all_cell_types,
                },
                data: Attributes {
                    point: vec![],
                    cell: vec![
                        Attribute::DataArray(DataArray {
                            name: "SourceID".to_string(),
                            elem: ElementType::Scalars { num_comp: 1, lookup_table: None },
                            data: IOBuffer::U32(source_id_data),
                        }),
                    ],
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
