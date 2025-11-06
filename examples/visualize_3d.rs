//! Example: Generate 3D B-field data for PyVista/Paraview visualization
//!
//! This example creates a cuboid magnet, computes the B-field on a 3D grid,
//! and exports to VTK format for advanced 3D visualization.
//!
//! The output can be visualized with:
//! - PyVista (Python): Interactive 3D streamlines, volume rendering
//! - Paraview: Professional scientific visualization
//! - Any VTK-compatible tool

use magrslib::export::vtk::{export_bfield_to_vtk, VtkGrid3D};
use magrslib::get_b;
use magrslib::sources::magnets::Cuboid;
use magrslib::types::Observers;

fn main() {
    println!("==============================================");
    println!("  3D B-Field Visualization Example");
    println!("==============================================\n");

    // Create a cuboid magnet
    let magnet = Cuboid::builder()
        .dimension([0.010, 0.010, 0.01]) // 1cm x 1cm x 4mm (similar to example)
        .polarization([0.0, 0.0, 1.0]) // 1 Tesla in +z direction
        .position([0.0, 0.0, 0.0]) // Centered at origin
        .build()
        .unwrap();

    println!("Magnet configuration:");
    println!("  Type: Cuboid");
    println!("  Dimensions: 10mm x 10mm x 4mm");
    println!("  Polarization: 1.0 T in +z direction");
    println!("  Position: origin\n");

    // Create 3D grid
    // Grid spacing of 1mm, covering -20mm to +20mm in all directions
    let grid_min = -0.050; // -20mm
    let grid_max = 0.050; // +20mm
    let grid_points = 128; // 41^3 = 68,921 points

    println!("Generating 3D observation grid:");
    println!("  Range: {}m to {}m", grid_min, grid_max);
    println!(
        "  Resolution: {}x{}x{} = {} points",
        grid_points,
        grid_points,
        grid_points,
        grid_points * grid_points * grid_points
    );
    println!("  Grid spacing: 1mm\n");

    let observers = Observers::from_grid_3d(
        (grid_min, grid_max, grid_points),
        (grid_min, grid_max, grid_points),
        (grid_min, grid_max, grid_points),
    );

    // Compute B-field
    println!("Computing B-field at {} points...", observers.len());
    let start = std::time::Instant::now();
    let b_field = get_b(&[magnet.into()], &observers);
    let elapsed = start.elapsed();
    println!("Done! Computation took {:.2}s\n", elapsed.as_secs_f64());

    // === DIAGNOSTICS ===
    println!("=== FIELD DIAGNOSTICS ===\n");

    // Check for NaN/Inf
    let mut nan_count = 0;
    let mut inf_count = 0;
    for i in 0..b_field.nrows() {
        for j in 0..3 {
            let val = b_field[[i, j]];
            if val.is_nan() {
                nan_count += 1;
            }
            if val.is_infinite() {
                inf_count += 1;
            }
        }
    }

    if nan_count > 0 || inf_count > 0 {
        println!("⚠️  WARNING: Invalid values detected!");
        println!("  NaN values: {}", nan_count);
        println!("  Inf values: {}\n", inf_count);
    } else {
        println!("✓ All field values are finite\n");
    }

    // Compute statistics (filter out NaN values)
    let mut magnitudes: Vec<f64> = (0..b_field.nrows())
        .map(|i| {
            let bx = b_field[[i, 0]];
            let by = b_field[[i, 1]];
            let bz = b_field[[i, 2]];
            (bx * bx + by * by + bz * bz).sqrt()
        })
        .filter(|&m| m.is_finite()) // Filter out NaN and Inf
        .collect();

    if magnitudes.is_empty() {
        println!("⚠️  ERROR: All field values are NaN/Inf!\n");
    } else {
        magnitudes.sort_by(|a, b| a.partial_cmp(b).unwrap());

        println!("Field magnitude statistics (excluding NaN/Inf):");
        println!("  Valid points: {} / {}", magnitudes.len(), b_field.nrows());
        println!("  Min:    {:.3e} T", magnitudes[0]);
        println!("  Max:    {:.3e} T", magnitudes[magnitudes.len() - 1]);
        println!(
            "  Mean:   {:.3e} T",
            magnitudes.iter().sum::<f64>() / magnitudes.len() as f64
        );
        println!("  Median: {:.3e} T\n", magnitudes[magnitudes.len() / 2]);
    }

    // Memory usage estimate
    let memory_mb = (b_field.nrows() * b_field.ncols() * 8) as f64 / 1024.0 / 1024.0;
    println!("Memory usage: {:.2} MB\n", memory_mb);

    println!("=========================\n");

    // Replace NaN/Inf values with zeros for visualization
    // (Points inside magnet produce NaN)
    let mut b_field_clean = b_field.clone();
    for i in 0..b_field_clean.nrows() {
        for j in 0..3 {
            if !b_field_clean[[i, j]].is_finite() {
                b_field_clean[[i, j]] = 0.0;
            }
        }
    }

    // Create VTK grid structure
    let vtk_grid = VtkGrid3D::from_ranges(
        (grid_min, grid_max, grid_points),
        (grid_min, grid_max, grid_points),
        (grid_min, grid_max, grid_points),
    );

    // Export to VTK
    let output_file = "bfield_3d.vtk";
    println!("Exporting to VTK format: {}...", output_file);
    println!("  (NaN values replaced with zeros for visualization)");

    let export_start = std::time::Instant::now();
    export_bfield_to_vtk(
        output_file,
        &vtk_grid,
        &b_field_clean, // Use cleaned data
        Some("B-field from cuboid magnet (10x10x4mm, 1T z-polarization)"),
    )
    .expect("Failed to export to VTK");

    let export_elapsed = export_start.elapsed();
    println!("Export took {:.2}s\n", export_elapsed.as_secs_f64());

    // Check file size
    let metadata = std::fs::metadata(output_file).unwrap();
    let file_size_mb = metadata.len() as f64 / 1024.0 / 1024.0;
    println!("Output file size: {:.2} MB\n", file_size_mb);

    println!("==============================================");
    println!("  Visualization Instructions");
    println!("==============================================\n");
    println!("To visualize the 3D field, run:");
    println!("  python3 visualize_3d.py\n");
    println!("This will create an interactive 3D visualization with:");
    println!("  - Streamlines following the magnetic field");
    println!("  - Color-coded field magnitude");
    println!("  - Magnet geometry overlay");
    println!("  - Interactive rotation/zoom\n");
    println!("Requirements:");
    println!("  pip install pyvista numpy\n");
}
