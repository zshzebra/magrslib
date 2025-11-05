//! magrslib - Magnetic field computation library
//!
//! A high-performance Rust port of magpylib's computational core.
//! Computes magnetic fields from various source geometries (magnets, currents, dipoles)
//! at observer positions.
//!
//! # Example
//! ```no_run
//! use magrslib::sources::magnets::Cuboid;
//! use magrslib::types::{Observers, Position};
//!
//! // Create a cuboid magnet
//! let magnet = Cuboid::builder()
//!     .dimension([0.01, 0.01, 0.01])  // 1cm cube
//!     .polarization([0.0, 0.0, 1.0])  // 1 Tesla in z-direction
//!     .position([0.0, 0.0, 0.0])
//!     .build()
//!     .unwrap();
//!
//! // Define observer positions
//! let observers = Observers::from_vec3s(vec![
//!     [0.0, 0.0, 0.02],  // 2cm above magnet
//! ]);
//!
//! // Compute B-field
//! use magrslib::get_b;
//! let b_field = get_b(&[magnet.into()], &observers);
//! ```

pub mod export;
pub mod fields;
pub mod sources;
pub mod special;
pub mod types;
pub mod utils;

// Re-export commonly used types
pub use fields::get_b;
pub use types::{Observers, Orientation, Path, Position, Vec3};