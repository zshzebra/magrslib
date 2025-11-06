//! Special mathematical functions
//!
//! This module contains special functions needed for field computations,
//! particularly complete elliptic integrals used in circular current loop calculations.

pub mod elliptic;
pub mod standard_elliptic;

pub use elliptic::{cel, cel0, celv, cel_iter, cel_iter_scalar};
pub use standard_elliptic::{ellipk, ellipe, ellippi, ellipk_vec, ellipe_vec, ellippi_vec};
