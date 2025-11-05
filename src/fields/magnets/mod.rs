//! Magnetic field computation for magnet sources
//!
//! This module contains field computation functions for various magnet geometries.

pub mod cuboid;

// Re-export for convenience
pub use cuboid::magnet_cuboid_bfield;
