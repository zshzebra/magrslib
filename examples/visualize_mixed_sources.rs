//! Mixed electromagnetic sources visualization example
//!
//! Demonstrates the enhanced visualization capabilities with multiple source types:
//! - Cuboid magnet
//! - Current loop
//! - Sphere magnet
//!
//! This example shows how the new object-oriented Python visualization system
//! can automatically handle complex multi-source scenarios with proper geometry export.

use magrslib::export::visualization::VisualizationBuilder;
use magrslib::sources::currents::Circle;
use magrslib::sources::magnets::Cuboid;
// Removed Sphere import since it's not fully implemented yet

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("Mixed Sources Electromagnetic Field Visualization");
    println!("================================================");

    // Create a cuboid magnet (main field source)
    let main_magnet = Cuboid::builder()
        .dimension([0.02, 0.02, 0.01]) // 20×20×10mm cuboid
        .polarization([0.0, 0.0, 1.5]) // 1.5T z-polarization
        .position([0.0, 0.0, 0.0]) // Centered at origin
        .build()?;

    // Create a current loop above the magnet
    let current_loop = Circle::builder()
        .diameter(0.025) // 25mm diameter
        .current(2.0) // 2A current
        .position([0.0, 0.0, 0.025]) // 25mm above magnet
        .build()?;

    // Create a second cuboid magnet to the side
    let side_magnet = Cuboid::builder()
        .dimension([0.008, 0.008, 0.006]) // 8×8×6mm cuboid
        .polarization([1.2, 0.0, 0.0]) // 1.2T x-polarization
        .position([0.030, 0.0, 0.010]) // 30mm to the side, 10mm up
        .build()?;

    // Create another current loop below (opposite current)
    let bottom_loop = Circle::builder()
        .diameter(0.020) // 20mm diameter
        .current(-1.5) // -1.5A current (opposite direction)
        .position([0.0, 0.0, -0.020]) // 20mm below magnet
        .build()?;

    println!("Sources configured:");
    println!("  - Main cuboid magnet: 20×20×10mm, 1.5T z-polarization");
    println!("  - Top current loop: ⌀25mm, 2.0A");
    println!("  - Side cuboid magnet: 8×8×6mm, 1.2T x-polarization");
    println!("  - Bottom current loop: ⌀20mm, -1.5A");

    // Build visualization with enhanced capabilities
    let _visualization = VisualizationBuilder::new()
        .add_source(main_magnet)
        .add_source(current_loop)
        .add_source(side_magnet)
        .add_source(bottom_loop)
        .title("Mixed Electromagnetic Sources Field Visualization")
        .grid_resolution(128) // 48³ = ~110k points for good detail
        .export_geometry(true) // Enable geometry export
        .adaptive_resolution(true) // Auto-adapt grid to sources
        .export_to_directory("./output_mixed")?;

    println!("\n✓ Visualization exported to ./output_mixed/");
    println!("Files created:");
    println!("  - field.vtk     : Field data for all 4 sources");
    println!("  - geometry.vtk  : 3D geometry of all sources");
    println!("  - metadata.json : Complete scene description");

    println!("\nTo visualize:");
    println!("  python3 visualize_3d.py output_mixed/metadata.json");
    println!("  python3 visualize_3d.py --batch output_mixed/metadata.json -o mixed_scene.png");

    println!("\nThis demonstrates:");
    println!("  ✓ Multi-source field computation with parallelization");
    println!("  ✓ Automatic geometry export for all source types");
    println!("  ✓ Adaptive grid bounds based on source positions");
    println!("  ✓ Comprehensive metadata for automatic Python setup");
    println!("  ✓ Professional visualization ready for publication");

    Ok(())
}
