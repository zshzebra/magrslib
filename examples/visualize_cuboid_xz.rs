//! Example: Generate B-field data for XZ-plane visualization
//!
//! This creates the same magnet but visualizes in the XZ plane (side view)
//! where we can actually see the field looping from north to south pole.

use magrslib::export::{FieldData, GridInfo, MagnetInfo};
use magrslib::sources::magnets::Cuboid;
use magrslib::types::Observers;
use magrslib::get_b;

fn main() {
    println!("Generating B-field data for XZ-plane visualization...\n");

    // Create a cuboid magnet at the origin
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

    // Create a 2D grid in the XZ-plane (side view)
    // This is where we'll see the field looping from north to south pole!
    let grid_min = -0.05;  // -5cm
    let grid_max = 0.05;   // +5cm
    let grid_points = 50;  // 50x50 = 2500 points

    println!("Generating observation grid:");
    println!("  Plane: xz (side view at y=0)");
    println!("  Range: {}m to {}m", grid_min, grid_max);
    println!("  Resolution: {}x{} = {} points\n",
             grid_points, grid_points, grid_points * grid_points);

    let observers = Observers::from_grid_2d(
        (grid_min, grid_max, grid_points),  // x-range
        (grid_min, grid_max, grid_points),  // z-range
        "xz",  // XZ-plane!
        0.0    // y=0 (cut through the middle)
    );

    // Compute B-field
    println!("Computing B-field...");
    let b_field = get_b(&[magnet.into()], &observers);
    println!("Done! Computed field at {} points\n", observers.len());

    // === DIAGNOSTIC OUTPUT ===
    println!("=== FIELD DIAGNOSTICS ===\n");

    // Check for in-plane components
    let mut bx_nonzero = 0;
    let mut bz_nonzero = 0;
    for i in 0..b_field.nrows() {
        if b_field[[i, 0]].abs() > 1e-10 {
            bx_nonzero += 1;
        }
        if b_field[[i, 2]].abs() > 1e-10 {
            bz_nonzero += 1;
        }
    }

    println!("In-plane field components (for streamplot):");
    println!("  Bx non-zero: {} points ({:.1}%)",
             bx_nonzero, 100.0 * bx_nonzero as f64 / b_field.nrows() as f64);
    println!("  Bz non-zero: {} points ({:.1}%)",
             bz_nonzero, 100.0 * bz_nonzero as f64 / b_field.nrows() as f64);

    // Field magnitude stats
    let mut magnitudes: Vec<f64> = (0..b_field.nrows())
        .map(|i| {
            let bx = b_field[[i, 0]];
            let by = b_field[[i, 1]];
            let bz = b_field[[i, 2]];
            (bx * bx + by * by + bz * bz).sqrt()
        })
        .collect();

    magnitudes.sort_by(|a, b| a.partial_cmp(b).unwrap());
    println!("\nField magnitude:");
    println!("  Min: {:.3e} T", magnitudes[0]);
    println!("  Max: {:.3e} T", magnitudes[magnitudes.len() - 1]);

    println!("\n=========================\n");

    // Create field data structure for export
    let grid_info = GridInfo {
        plane: "xz".to_string(),  // XZ plane!
        nx: grid_points,
        ny: grid_points,
        x_range: (grid_min, grid_max),
        y_range: (grid_min, grid_max),  // Actually z-range for xz plane
        fixed_coord: 0.0,
    };

    let mut field_data = FieldData::new(
        "B-field from cuboid magnet - XZ plane (side view)".to_string(),
        grid_info,
        observers.as_array(),
        &b_field,
    );

    // Add magnet information
    field_data.add_magnet(MagnetInfo {
        magnet_type: "cuboid".to_string(),
        position: vec![0.0, 0.0, 0.0],
        dimensions: vec![0.02, 0.02, 0.01],
        polarization: vec![0.0, 0.0, 1.0],
    });

    // Export to JSON
    let output_file = "bfield_data_xz.json";
    println!("Exporting data to {}...", output_file);
    field_data.export_to_json(output_file)
        .expect("Failed to export to JSON");

    println!("Done!\n");
    println!("To visualize the data, run:");
    println!("  python3 visualize_bfield.py");
    println!("  (it will automatically detect it's an XZ plane)\n");
    println!("Or load {} instead of bfield_data.json", output_file);
}
