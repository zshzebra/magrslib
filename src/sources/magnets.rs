//! Magnet source implementations (Cuboid, Cylinder, Sphere)

use crate::sources::{FieldFuncId, Source, SourceProperties};
use crate::types::{Orientation, Path, Position, Vec3};
use crate::utils::MU_0;

// =============================================================================
// Cuboid Magnet
// =============================================================================

/// Rectangular cuboid magnet with homogeneous magnetization
///
/// Equivalent to Python's `magpylib.magnet.Cuboid`
#[derive(Debug, Clone)]
pub struct Cuboid {
    dimension: Vec3,         // (a, b, c) in meters
    polarization: Vec3,      // (Jx, Jy, Jz) in Tesla
    path: Path,              // Time-varying positions
    orientation: Orientation, // Rotation
}

impl Cuboid {
    /// Create a builder for constructing a Cuboid magnet
    pub fn builder() -> CuboidBuilder {
        CuboidBuilder::new()
    }
}

impl Source for Cuboid {
    fn field_func_id(&self) -> FieldFuncId {
        FieldFuncId::MagnetCuboid
    }

    fn path(&self) -> &Path {
        &self.path
    }

    fn orientation(&self) -> &Orientation {
        &self.orientation
    }

    fn get_properties(&self) -> SourceProperties {
        SourceProperties {
            dimension: Some(self.dimension),
            polarization: Some(self.polarization),
            ..Default::default()
        }
    }
}

/// Builder for Cuboid magnets
#[derive(Debug, Default)]
pub struct CuboidBuilder {
    dimension: Option<Vec3>,
    polarization: Option<Vec3>,
    magnetization: Option<Vec3>,
    path: Option<Path>,
    orientation: Option<Orientation>,
}

impl CuboidBuilder {
    /// Create a new builder
    pub fn new() -> Self {
        Self::default()
    }

    /// Set the cuboid dimensions (a, b, c) in meters
    pub fn dimension(mut self, dim: Vec3) -> Self {
        self.dimension = Some(dim);
        self
    }

    /// Set the polarization (Jx, Jy, Jz) in Tesla
    pub fn polarization(mut self, pol: Vec3) -> Self {
        self.polarization = Some(pol);
        self
    }

    /// Set the magnetization (Mx, My, Mz) in A/m
    /// Will be converted to polarization: J = μ₀ * M
    pub fn magnetization(mut self, mag: Vec3) -> Self {
        self.magnetization = Some(mag);
        self
    }

    /// Set a single position
    pub fn position(mut self, pos: impl Into<Position>) -> Self {
        self.path = Some(Path::from_position(pos.into()));
        self
    }

    /// Set a path (multiple positions)
    pub fn path(mut self, path: impl Into<Path>) -> Self {
        self.path = Some(path.into());
        self
    }

    /// Set the orientation
    pub fn orientation(mut self, orient: Orientation) -> Self {
        self.orientation = Some(orient);
        self
    }

    /// Build the Cuboid magnet
    pub fn build(self) -> Result<Cuboid, BuildError> {
        let dimension = self.dimension.ok_or(BuildError::MissingDimension)?;

        // Convert magnetization to polarization if provided
        let polarization = match (self.polarization, self.magnetization) {
            (Some(pol), None) => pol,
            (None, Some(mag)) => [mag[0] * MU_0, mag[1] * MU_0, mag[2] * MU_0],
            (Some(_), Some(_)) => return Err(BuildError::BothPolarizationAndMagnetization),
            (None, None) => return Err(BuildError::MissingPolarization),
        };

        let path = self
            .path
            .unwrap_or_else(|| Path::from_position(Position::origin()));
        let orientation = self.orientation.unwrap_or_else(Orientation::identity);

        Ok(Cuboid {
            dimension,
            polarization,
            path,
            orientation,
        })
    }
}

// =============================================================================
// Cylinder Magnet
// =============================================================================

/// Cylindrical magnet with homogeneous magnetization
///
/// Equivalent to Python's `magpylib.magnet.Cylinder`
#[derive(Debug, Clone)]
pub struct Cylinder {
    diameter: f64,        // Diameter in meters
    height: f64,          // Height in meters
    polarization: Vec3,   // (Jx, Jy, Jz) in Tesla
    path: Path,           // Time-varying positions
    orientation: Orientation, // Rotation
}

impl Cylinder {
    /// Create a builder for constructing a Cylinder magnet
    pub fn builder() -> CylinderBuilder {
        CylinderBuilder::new()
    }
}

impl Source for Cylinder {
    fn field_func_id(&self) -> FieldFuncId {
        FieldFuncId::MagnetCylinder
    }

    fn path(&self) -> &Path {
        &self.path
    }

    fn orientation(&self) -> &Orientation {
        &self.orientation
    }

    fn get_properties(&self) -> SourceProperties {
        SourceProperties {
            diameter: Some(self.diameter),
            height: Some(self.height),
            polarization: Some(self.polarization),
            ..Default::default()
        }
    }
}

/// Builder for Cylinder magnets
#[derive(Debug, Default)]
pub struct CylinderBuilder {
    diameter: Option<f64>,
    height: Option<f64>,
    polarization: Option<Vec3>,
    magnetization: Option<Vec3>,
    path: Option<Path>,
    orientation: Option<Orientation>,
}

impl CylinderBuilder {
    /// Create a new builder
    pub fn new() -> Self {
        Self::default()
    }

    /// Set the cylinder diameter in meters
    pub fn diameter(mut self, d: f64) -> Self {
        self.diameter = Some(d);
        self
    }

    /// Set the cylinder height in meters
    pub fn height(mut self, h: f64) -> Self {
        self.height = Some(h);
        self
    }

    /// Set the polarization (Jx, Jy, Jz) in Tesla
    pub fn polarization(mut self, pol: Vec3) -> Self {
        self.polarization = Some(pol);
        self
    }

    /// Set the magnetization (Mx, My, Mz) in A/m
    pub fn magnetization(mut self, mag: Vec3) -> Self {
        self.magnetization = Some(mag);
        self
    }

    /// Set a single position
    pub fn position(mut self, pos: impl Into<Position>) -> Self {
        self.path = Some(Path::from_position(pos.into()));
        self
    }

    /// Set a path (multiple positions)
    pub fn path(mut self, path: impl Into<Path>) -> Self {
        self.path = Some(path.into());
        self
    }

    /// Set the orientation
    pub fn orientation(mut self, orient: Orientation) -> Self {
        self.orientation = Some(orient);
        self
    }

    /// Build the Cylinder magnet
    pub fn build(self) -> Result<Cylinder, BuildError> {
        let diameter = self.diameter.ok_or(BuildError::MissingDiameter)?;
        let height = self.height.ok_or(BuildError::MissingHeight)?;

        let polarization = match (self.polarization, self.magnetization) {
            (Some(pol), None) => pol,
            (None, Some(mag)) => [mag[0] * MU_0, mag[1] * MU_0, mag[2] * MU_0],
            (Some(_), Some(_)) => return Err(BuildError::BothPolarizationAndMagnetization),
            (None, None) => return Err(BuildError::MissingPolarization),
        };

        let path = self
            .path
            .unwrap_or_else(|| Path::from_position(Position::origin()));
        let orientation = self.orientation.unwrap_or_else(Orientation::identity);

        Ok(Cylinder {
            diameter,
            height,
            polarization,
            path,
            orientation,
        })
    }
}

// =============================================================================
// Sphere Magnet
// =============================================================================

/// Spherical magnet with homogeneous magnetization
///
/// Equivalent to Python's `magpylib.magnet.Sphere`
#[derive(Debug, Clone)]
pub struct Sphere {
    diameter: f64,        // Diameter in meters
    polarization: Vec3,   // (Jx, Jy, Jz) in Tesla
    path: Path,           // Time-varying positions
    orientation: Orientation, // Rotation (not physically meaningful for sphere, but kept for API consistency)
}

impl Sphere {
    /// Create a builder for constructing a Sphere magnet
    pub fn builder() -> SphereBuilder {
        SphereBuilder::new()
    }
}

impl Source for Sphere {
    fn field_func_id(&self) -> FieldFuncId {
        FieldFuncId::MagnetSphere
    }

    fn path(&self) -> &Path {
        &self.path
    }

    fn orientation(&self) -> &Orientation {
        &self.orientation
    }

    fn get_properties(&self) -> SourceProperties {
        SourceProperties {
            diameter: Some(self.diameter),
            polarization: Some(self.polarization),
            ..Default::default()
        }
    }
}

/// Builder for Sphere magnets
#[derive(Debug, Default)]
pub struct SphereBuilder {
    diameter: Option<f64>,
    polarization: Option<Vec3>,
    magnetization: Option<Vec3>,
    path: Option<Path>,
    orientation: Option<Orientation>,
}

impl SphereBuilder {
    /// Create a new builder
    pub fn new() -> Self {
        Self::default()
    }

    /// Set the sphere diameter in meters
    pub fn diameter(mut self, d: f64) -> Self {
        self.diameter = Some(d);
        self
    }

    /// Set the polarization (Jx, Jy, Jz) in Tesla
    pub fn polarization(mut self, pol: Vec3) -> Self {
        self.polarization = Some(pol);
        self
    }

    /// Set the magnetization (Mx, My, Mz) in A/m
    pub fn magnetization(mut self, mag: Vec3) -> Self {
        self.magnetization = Some(mag);
        self
    }

    /// Set a single position
    pub fn position(mut self, pos: impl Into<Position>) -> Self {
        self.path = Some(Path::from_position(pos.into()));
        self
    }

    /// Set a path (multiple positions)
    pub fn path(mut self, path: impl Into<Path>) -> Self {
        self.path = Some(path.into());
        self
    }

    /// Set the orientation (not physically meaningful for sphere)
    pub fn orientation(mut self, orient: Orientation) -> Self {
        self.orientation = Some(orient);
        self
    }

    /// Build the Sphere magnet
    pub fn build(self) -> Result<Sphere, BuildError> {
        let diameter = self.diameter.ok_or(BuildError::MissingDiameter)?;

        let polarization = match (self.polarization, self.magnetization) {
            (Some(pol), None) => pol,
            (None, Some(mag)) => [mag[0] * MU_0, mag[1] * MU_0, mag[2] * MU_0],
            (Some(_), Some(_)) => return Err(BuildError::BothPolarizationAndMagnetization),
            (None, None) => return Err(BuildError::MissingPolarization),
        };

        let path = self
            .path
            .unwrap_or_else(|| Path::from_position(Position::origin()));
        let orientation = self.orientation.unwrap_or_else(Orientation::identity);

        Ok(Sphere {
            diameter,
            polarization,
            path,
            orientation,
        })
    }
}

// =============================================================================
// Build Errors
// =============================================================================

/// Errors that can occur when building sources
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum BuildError {
    MissingDimension,
    MissingDiameter,
    MissingHeight,
    MissingPolarization,
    MissingCurrent,
    MissingVertices,
    MissingMoment,
    BothPolarizationAndMagnetization,
}

impl std::fmt::Display for BuildError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            BuildError::MissingDimension => write!(f, "Missing dimension parameter"),
            BuildError::MissingDiameter => write!(f, "Missing diameter parameter"),
            BuildError::MissingHeight => write!(f, "Missing height parameter"),
            BuildError::MissingPolarization => {
                write!(f, "Missing polarization or magnetization parameter")
            }
            BuildError::MissingCurrent => write!(f, "Missing current parameter"),
            BuildError::MissingVertices => write!(f, "Missing vertices parameter"),
            BuildError::MissingMoment => write!(f, "Missing moment parameter"),
            BuildError::BothPolarizationAndMagnetization => {
                write!(f, "Cannot specify both polarization and magnetization")
            }
        }
    }
}

impl std::error::Error for BuildError {}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_cuboid_builder() {
        let cuboid = Cuboid::builder()
            .dimension([0.01, 0.01, 0.01])
            .polarization([0.0, 0.0, 1.0])
            .position([0.0, 0.0, 0.0])
            .build()
            .unwrap();

        assert_eq!(cuboid.dimension, [0.01, 0.01, 0.01]);
        assert_eq!(cuboid.polarization, [0.0, 0.0, 1.0]);
    }

    #[test]
    fn test_cylinder_builder() {
        let cylinder = Cylinder::builder()
            .diameter(0.01)
            .height(0.02)
            .polarization([0.0, 0.0, 1.0])
            .build()
            .unwrap();

        assert_eq!(cylinder.diameter, 0.01);
        assert_eq!(cylinder.height, 0.02);
    }

    #[test]
    fn test_sphere_builder() {
        let sphere = Sphere::builder()
            .diameter(0.01)
            .polarization([0.0, 0.0, 1.0])
            .build()
            .unwrap();

        assert_eq!(sphere.diameter, 0.01);
    }

    #[test]
    fn test_magnetization_conversion() {
        let cuboid = Cuboid::builder()
            .dimension([0.01, 0.01, 0.01])
            .magnetization([0.0, 0.0, 1.0e6]) // 1 MA/m
            .build()
            .unwrap();

        // J = μ₀ * M ≈ 1.257 T
        assert!((cuboid.polarization[2] - MU_0 * 1.0e6).abs() < 1e-6);
    }

    #[test]
    fn test_path_support() {
        let positions = vec![
            Position::new(0.0, 0.0, 0.0),
            Position::new(0.01, 0.0, 0.0),
            Position::new(0.02, 0.0, 0.0),
        ];

        let cuboid = Cuboid::builder()
            .dimension([0.01, 0.01, 0.01])
            .polarization([0.0, 0.0, 1.0])
            .path(positions)
            .build()
            .unwrap();

        assert_eq!(cuboid.path.len(), 3);
    }
}
