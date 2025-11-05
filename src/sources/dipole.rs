//! Dipole source implementation

use crate::sources::{FieldFuncId, Source, SourceProperties};
use crate::types::{Orientation, Path, Position, Vec3};

/// Point magnetic dipole
///
/// A point dipole with magnetic moment m (in A·m²).
/// Equivalent to Python's `magpylib.misc.Dipole`
#[derive(Debug, Clone)]
pub struct Dipole {
    moment: Vec3,         // Magnetic moment (mx, my, mz) in A·m²
    path: Path,           // Time-varying positions
    orientation: Orientation, // Rotation
}

impl Dipole {
    /// Create a builder for constructing a Dipole source
    pub fn builder() -> DipoleBuilder {
        DipoleBuilder::new()
    }
}

impl Source for Dipole {
    fn field_func_id(&self) -> FieldFuncId {
        FieldFuncId::Dipole
    }

    fn path(&self) -> &Path {
        &self.path
    }

    fn orientation(&self) -> &Orientation {
        &self.orientation
    }

    fn get_properties(&self) -> SourceProperties {
        SourceProperties {
            moment: Some(self.moment),
            ..Default::default()
        }
    }
}

/// Builder for Dipole sources
#[derive(Debug, Default)]
pub struct DipoleBuilder {
    moment: Option<Vec3>,
    path: Option<Path>,
    orientation: Option<Orientation>,
}

impl DipoleBuilder {
    /// Create a new builder
    pub fn new() -> Self {
        Self::default()
    }

    /// Set the magnetic moment (mx, my, mz) in A·m²
    pub fn moment(mut self, m: Vec3) -> Self {
        self.moment = Some(m);
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

    /// Build the Dipole source
    pub fn build(self) -> Result<Dipole, super::magnets::BuildError> {
        let moment = self
            .moment
            .ok_or(super::magnets::BuildError::MissingMoment)?;

        let path = self
            .path
            .unwrap_or_else(|| Path::from_position(Position::origin()));
        let orientation = self.orientation.unwrap_or_else(Orientation::identity);

        Ok(Dipole {
            moment,
            path,
            orientation,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_dipole_builder() {
        let dipole = Dipole::builder()
            .moment([0.0, 0.0, 1.0e-6])
            .position([0.0, 0.0, 0.0])
            .build()
            .unwrap();

        assert_eq!(dipole.moment, [0.0, 0.0, 1.0e-6]);
    }

    #[test]
    fn test_dipole_with_path() {
        let positions = vec![
            Position::new(0.0, 0.0, 0.0),
            Position::new(0.01, 0.0, 0.0),
        ];

        let dipole = Dipole::builder()
            .moment([0.0, 0.0, 1.0e-6])
            .path(positions)
            .build()
            .unwrap();

        assert_eq!(dipole.path.len(), 2);
    }
}
