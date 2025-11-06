//! Current source implementations (Circle, Polyline)

use crate::sources::{FieldFuncId, Source, SourceProperties};
use crate::sources::current_functions::CurrentFunction;
use crate::types::{Orientation, Path, Position, Vec3};

// =============================================================================
// Circle Current Source
// =============================================================================

/// Circular current loop
///
/// A current loop lying in the xy-plane (centered at origin, normal in +z).
/// Equivalent to Python's `magpylib.current.Circle`
#[derive(Debug, Clone)]
pub struct Circle {
    diameter: f64,        // Diameter in meters
    current: f64,         // Current in Amperes
    path: Path,           // Time-varying positions
    orientation: Orientation, // Rotation
}

impl Circle {
    /// Create a builder for constructing a Circle current source
    pub fn builder() -> CircleBuilder {
        CircleBuilder::new()
    }
}

impl Source for Circle {
    fn field_func_id(&self) -> FieldFuncId {
        FieldFuncId::CurrentCircle
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
            current: Some(self.current),
            ..Default::default()
        }
    }
}

/// Builder for Circle current sources
#[derive(Debug, Default)]
pub struct CircleBuilder {
    diameter: Option<f64>,
    current: Option<f64>,
    path: Option<Path>,
    orientation: Option<Orientation>,
}

impl CircleBuilder {
    /// Create a new builder
    pub fn new() -> Self {
        Self::default()
    }

    /// Set the circle diameter in meters
    pub fn diameter(mut self, d: f64) -> Self {
        self.diameter = Some(d);
        self
    }

    /// Set the current in Amperes
    pub fn current(mut self, i: f64) -> Self {
        self.current = Some(i);
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

    /// Build the Circle current source
    pub fn build(self) -> Result<Circle, super::magnets::BuildError> {
        let diameter = self
            .diameter
            .ok_or(super::magnets::BuildError::MissingDiameter)?;
        let current = self
            .current
            .ok_or(super::magnets::BuildError::MissingCurrent)?;

        let path = self
            .path
            .unwrap_or_else(|| Path::from_position(Position::origin()));
        let orientation = self.orientation.unwrap_or_else(Orientation::identity);

        Ok(Circle {
            diameter,
            current,
            path,
            orientation,
        })
    }
}

// =============================================================================
// Polyline Current Source
// =============================================================================

/// Piecewise linear current path (line segments)
///
/// Current flows through connected line segments defined by vertices.
/// Supports both static and time-varying current I(t).
/// Equivalent to Python's `magpylib.current.Polyline`
#[derive(Debug, Clone)]
pub struct Polyline {
    vertices: Vec<Vec3>,  // Vertices defining the line segments
    current: f64,         // Static current in Amperes (for backward compatibility)
    current_function: Option<Box<dyn CurrentFunction>>, // Time-varying current I(t)
    path: Path,           // Time-varying positions (applies to entire polyline)
    orientation: Orientation, // Rotation
}

impl Polyline {
    /// Create a builder for constructing a Polyline current source
    pub fn builder() -> PolylineBuilder {
        PolylineBuilder::new()
    }
}

impl Source for Polyline {
    fn field_func_id(&self) -> FieldFuncId {
        FieldFuncId::CurrentPolyline
    }

    fn path(&self) -> &Path {
        &self.path
    }

    fn orientation(&self) -> &Orientation {
        &self.orientation
    }

    fn get_properties(&self) -> SourceProperties {
        SourceProperties {
            vertices: Some(self.vertices.clone()),
            current: Some(self.current),
            current_function: self.current_function.clone(),
            ..Default::default()
        }
    }
}

/// Builder for Polyline current sources
#[derive(Debug, Default)]
pub struct PolylineBuilder {
    vertices: Option<Vec<Vec3>>,
    current: Option<f64>,
    current_function: Option<Box<dyn CurrentFunction>>,
    path: Option<Path>,
    orientation: Option<Orientation>,
}

impl PolylineBuilder {
    /// Create a new builder
    pub fn new() -> Self {
        Self::default()
    }

    /// Set the vertices defining the line segments
    pub fn vertices(mut self, verts: Vec<Vec3>) -> Self {
        self.vertices = Some(verts);
        self
    }

    /// Set the current in Amperes (static/DC current)
    pub fn current(mut self, i: f64) -> Self {
        self.current = Some(i);
        self.current_function = None; // Clear any current function
        self
    }

    /// Set a time-varying current function I(t)
    pub fn current_function(mut self, func: Box<dyn CurrentFunction>) -> Self {
        self.current_function = Some(func);
        self.current = Some(0.0); // Set fallback value for backward compatibility
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

    /// Build the Polyline current source
    pub fn build(self) -> Result<Polyline, super::magnets::BuildError> {
        let vertices = self
            .vertices
            .ok_or(super::magnets::BuildError::MissingVertices)?;

        // Require either static current or current function
        let current = self.current.unwrap_or(0.0);
        if current == 0.0 && self.current_function.is_none() {
            return Err(super::magnets::BuildError::MissingCurrent);
        }

        if vertices.len() < 2 {
            return Err(super::magnets::BuildError::MissingVertices);
        }

        let path = self
            .path
            .unwrap_or_else(|| Path::from_position(Position::origin()));
        let orientation = self.orientation.unwrap_or_else(Orientation::identity);

        Ok(Polyline {
            vertices,
            current,
            current_function: self.current_function,
            path,
            orientation,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_circle_builder() {
        let circle = Circle::builder()
            .diameter(0.01)
            .current(1.0)
            .position([0.0, 0.0, 0.0])
            .build()
            .unwrap();

        assert_eq!(circle.diameter, 0.01);
        assert_eq!(circle.current, 1.0);
    }

    #[test]
    fn test_polyline_builder() {
        let vertices = vec![[0.0, 0.0, 0.0], [0.01, 0.0, 0.0], [0.01, 0.01, 0.0]];

        let polyline = Polyline::builder()
            .vertices(vertices.clone())
            .current(1.0)
            .build()
            .unwrap();

        assert_eq!(polyline.vertices.len(), 3);
        assert_eq!(polyline.current, 1.0);
    }

    #[test]
    fn test_polyline_minimum_vertices() {
        // Should fail with fewer than 2 vertices
        let result = Polyline::builder()
            .vertices(vec![[0.0, 0.0, 0.0]])
            .current(1.0)
            .build();

        assert!(result.is_err());
    }
}
