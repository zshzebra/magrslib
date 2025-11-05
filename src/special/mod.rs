//! Special mathematical functions
//!
//! This module contains special functions needed for field computations,
//! particularly complete elliptic integrals used in circular current loop calculations.

pub mod elliptic;

pub use elliptic::{cel_iter, cel_iter_scalar};
