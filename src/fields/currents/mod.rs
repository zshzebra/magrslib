//! Magnetic field computation for current sources
//!
//! This module contains field computation functions for various current geometries.

pub mod circle;
pub mod polyline;

// Re-export for convenience
pub use circle::compute_circle_field;
pub use polyline::compute_polyline_field;
