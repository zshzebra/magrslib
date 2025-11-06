//! Export functionality for visualization and data interchange
//!
//! This module provides tools for exporting magnetic field data and source
//! geometry in various formats for visualization and analysis.

pub mod vtk;
pub mod visualization;

// Re-export main APIs for convenience
pub use vtk::{export_bfield_to_vtk, VtkGrid3D};
pub use visualization::{VisualizationBuilder, VisualizationMetadata};