//! Example: Generate B-field data for visualization with matplotlib
//!
//! This example creates a cuboid magnet, computes the B-field on a 2D grid,
//! and exports the data to JSON for visualization with Python.

use magrslib::export::{FieldData, GridInfo, MagnetInfo};
use magrslib::sources::magnets::Cuboid;
use magrslib::types::Observers;
use magrslib::get_b;

fn main() {
    println!("Generating B-field data for visualization...\n");

    // Create a cuboid magnet at the origin
    // Dimensions in meters, polarization in Tesla
    let magnet = Cuboid::builder()
        .dimension([0.02, 0.02, 0.01])      // 2cm x 2cm x 1cm
        .polarization([0.0, 0.0, 1.0])       // 1 Tesla in +z direction
        .position([0.0, 0.0, 0.0])           // Centered at origin
        .build()
        .unwrap();

    println!("Magnet configuration:");
    println!("  Type: Cuboid");
    println!("  Dimensions: 2cm x 2cm x 1cm");
    println!("  Polarization: 1.0 T in +z direction");
    println!("  Position: origin\n");

    // Create a 2D grid of observer positions in the xy-plane
    // Grid from -5cm to +5cm in both x and y, at z=0
    let grid_min = -0.05;  // -5cm
    let grid_max = 0.05;   // +5cm
    let grid_points = 50;  // 50x50 = 2500 points

    println!("Generating observation grid:");
    println!("  Plane: xy (at z=0)");
    println!("  Range: {}m to {}m", grid_min, grid_max);
    println!("  Resolution: {}x{} = {} points\n",
             grid_points, grid_points, grid_points * grid_points);

    let observers = Observers::from_grid_2d(
        (grid_min, grid_max, grid_points),
        (grid_min, grid_max, grid_points),
        "xy",
        0.0  // z-coordinate
    );

    // Compute B-field at all observer positions
    println!("Computing B-field...");
    let b_field = get_b(&[magnet.into()], &observers);
    println!("Done! Computed field at {} points\n", observers.len());

    // === DIAGNOSTIC OUTPUT ===
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

    // Compute field magnitudes
    let mut magnitudes: Vec<f64> = (0..b_field.nrows())
        .map(|i| {
            let bx = b_field[[i, 0]];
            let by = b_field[[i, 1]];
            let bz = b_field[[i, 2]];
            (bx * bx + by * by + bz * bz).sqrt()
        })
        .collect();

    // Statistics
    magnitudes.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let min_b = magnitudes[0];
    let max_b = magnitudes[magnitudes.len() - 1];
    let mean_b = magnitudes.iter().sum::<f64>() / magnitudes.len() as f64;
    let median_b = magnitudes[magnitudes.len() / 2];

    println!("Field magnitude statistics:");
    println!("  Min:    {:.3e} T", min_b);
    println!("  Max:    {:.3e} T", max_b);
    println!("  Mean:   {:.3e} T", mean_b);
    println!("  Median: {:.3e} T\n", median_b);

    // Sample specific points for verification
    println!("Sample field values at specific points:");

    // Find point closest to (0, 0, 0) - should be at magnet center
    let center_idx = (grid_points / 2) * grid_points + (grid_points / 2);
    println!("  At origin (0,0,0):");
    println!("    Bx = {:.6e} T", b_field[[center_idx, 0]]);
    println!("    By = {:.6e} T", b_field[[center_idx, 1]]);
    println!("    Bz = {:.6e} T", b_field[[center_idx, 2]]);
    println!("    |B| = {:.6e} T", magnitudes[center_idx]);

    // Point far from magnet (corner of grid)
    let corner_idx = 0;  // (-0.05, -0.05, 0)
    println!("  At corner (-0.05, -0.05, 0):");
    println!("    Bx = {:.6e} T", b_field[[corner_idx, 0]]);
    println!("    By = {:.6e} T", b_field[[corner_idx, 1]]);
    println!("    Bz = {:.6e} T", b_field[[corner_idx, 2]]);
    println!("    |B| = {:.6e} T", magnitudes[corner_idx]);

    // Verify grid ordering - print first 3 positions
    println!("\nGrid ordering (first 3 points):");
    let obs_array = observers.as_array();
    for i in 0..3 {
        println!("  Point {}: ({:.4}, {:.4}, {:.4})",
                 i, obs_array[[i, 0]], obs_array[[i, 1]], obs_array[[i, 2]]);
    }

    println!("\n=========================\n");

    // Create field data structure for export
    let grid_info = GridInfo {
        plane: "xy".to_string(),
        nx: grid_points,
        ny: grid_points,
        x_range: (grid_min, grid_max),
        y_range: (grid_min, grid_max),
        fixed_coord: 0.0,
    };

    let mut field_data = FieldData::new(
        "B-field from cuboid magnet (2cm x 2cm x 1cm, 1T z-polarization)".to_string(),
        grid_info,
        observers.as_array(),
        &b_field,
    );

    // Add magnet information for overlay visualization
    field_data.add_magnet(MagnetInfo {
        magnet_type: "cuboid".to_string(),
        position: vec![0.0, 0.0, 0.0],
        dimensions: vec![0.02, 0.02, 0.01],
        polarization: vec![0.0, 0.0, 1.0],
    });

    // Export to JSON
    let output_file = "bfield_data.json";
    println!("Exporting data to {}...", output_file);
    field_data.export_to_json(output_file)
        .expect("Failed to export to JSON");

    println!("Done!\n");
    println!("To visualize the data, run:");
    println!("  python3 visualize_bfield.py");
    println!("\nThe Python script will load {} and create a streamplot.", output_file);
}
