//! Magnetic field computation for magnet sources
//!
//! This module contains field computation functions for various magnet geometries.

pub mod cuboid;
pub mod sphere;

// Re-export for convenience
pub use cuboid::magnet_cuboid_bfield;
pub use sphere::compute_sphere_field;
