//! Source objects (magnets, current sources, dipoles)
//!
//! This module defines all source types that generate magnetic fields.
//! Sources implement the `Source` trait and can be collected in `AnySource` enum
//! for heterogeneous collections.

pub mod currents;
pub mod dipole;
pub mod magnets;

use crate::types::{Orientation, Path, Vec3};

/// Identifier for field computation functions
///
/// Used to group sources by type for vectorized computation
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum FieldFuncId {
    MagnetCuboid,
    MagnetCylinder,
    MagnetSphere,
    CurrentCircle,
    CurrentPolyline,
    Dipole,
}

/// Properties extracted from sources for vectorized field computation
///
/// This struct contains all the geometry-specific parameters needed
/// for field computation. Not all fields are used by all source types.
#[derive(Debug, Clone)]
pub struct SourceProperties {
    // Magnet properties
    pub dimension: Option<Vec3>,     // (a, b, c) for cuboid
    pub polarization: Option<Vec3>,  // (Jx, Jy, Jz) in Tesla
    pub magnetization: Option<Vec3>, // (Mx, My, Mz) in A/m
    pub diameter: Option<f64>,       // For cylinders, circles
    pub height: Option<f64>,         // For cylinders

    // Current properties
    pub current: Option<f64>,     // Current in Amperes
    pub vertices: Option<Vec<Vec3>>, // For polylines

    // Dipole properties
    pub moment: Option<Vec3>, // Magnetic moment in A·m²
}

impl SourceProperties {
    /// Create empty properties
    pub fn new() -> Self {
        Self {
            dimension: None,
            polarization: None,
            magnetization: None,
            diameter: None,
            height: None,
            current: None,
            vertices: None,
            moment: None,
        }
    }
}

impl Default for SourceProperties {
    fn default() -> Self {
        Self::new()
    }
}

/// Common interface for all magnetic field sources
pub trait Source {
    /// Get the field function identifier for this source type
    fn field_func_id(&self) -> FieldFuncId;

    /// Get the path (time-varying positions)
    fn path(&self) -> &Path;

    /// Get the orientation
    fn orientation(&self) -> &Orientation;

    /// Extract properties needed for field computation
    fn get_properties(&self) -> SourceProperties;
}

/// Type-erased container for heterogeneous source collections
///
/// This enum allows storing different source types in a single Vec
/// while maintaining type safety and enabling dynamic dispatch.
#[derive(Debug, Clone)]
pub enum AnySource {
    Cuboid(magnets::Cuboid),
    Cylinder(magnets::Cylinder),
    Sphere(magnets::Sphere),
    Circle(currents::Circle),
    Polyline(currents::Polyline),
    Dipole(dipole::Dipole),
}

impl Source for AnySource {
    fn field_func_id(&self) -> FieldFuncId {
        match self {
            AnySource::Cuboid(s) => s.field_func_id(),
            AnySource::Cylinder(s) => s.field_func_id(),
            AnySource::Sphere(s) => s.field_func_id(),
            AnySource::Circle(s) => s.field_func_id(),
            AnySource::Polyline(s) => s.field_func_id(),
            AnySource::Dipole(s) => s.field_func_id(),
        }
    }

    fn path(&self) -> &Path {
        match self {
            AnySource::Cuboid(s) => s.path(),
            AnySource::Cylinder(s) => s.path(),
            AnySource::Sphere(s) => s.path(),
            AnySource::Circle(s) => s.path(),
            AnySource::Polyline(s) => s.path(),
            AnySource::Dipole(s) => s.path(),
        }
    }

    fn orientation(&self) -> &Orientation {
        match self {
            AnySource::Cuboid(s) => s.orientation(),
            AnySource::Cylinder(s) => s.orientation(),
            AnySource::Sphere(s) => s.orientation(),
            AnySource::Circle(s) => s.orientation(),
            AnySource::Polyline(s) => s.orientation(),
            AnySource::Dipole(s) => s.orientation(),
        }
    }

    fn get_properties(&self) -> SourceProperties {
        match self {
            AnySource::Cuboid(s) => s.get_properties(),
            AnySource::Cylinder(s) => s.get_properties(),
            AnySource::Sphere(s) => s.get_properties(),
            AnySource::Circle(s) => s.get_properties(),
            AnySource::Polyline(s) => s.get_properties(),
            AnySource::Dipole(s) => s.get_properties(),
        }
    }
}

// Implement From conversions for easy construction
impl From<magnets::Cuboid> for AnySource {
    fn from(s: magnets::Cuboid) -> Self {
        AnySource::Cuboid(s)
    }
}

impl From<magnets::Cylinder> for AnySource {
    fn from(s: magnets::Cylinder) -> Self {
        AnySource::Cylinder(s)
    }
}

impl From<magnets::Sphere> for AnySource {
    fn from(s: magnets::Sphere) -> Self {
        AnySource::Sphere(s)
    }
}

impl From<currents::Circle> for AnySource {
    fn from(s: currents::Circle) -> Self {
        AnySource::Circle(s)
    }
}

impl From<currents::Polyline> for AnySource {
    fn from(s: currents::Polyline) -> Self {
        AnySource::Polyline(s)
    }
}

impl From<dipole::Dipole> for AnySource {
    fn from(s: dipole::Dipole) -> Self {
        AnySource::Dipole(s)
    }
}
