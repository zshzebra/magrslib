//! Solenoid coil electromagnetic field visualization example
//!
//! Demonstrates creating a realistic solenoid using CurrentPolyline with a helical path:
//! - Continuous helical current path (physically accurate)
//! - Uniform magnetic field inside the coil
//! - Dipole-like field pattern outside the coil
//! - Proper solenoid field characteristics
//!
//! This example shows how to model inductors, electromagnets, and solenoid actuators
//! using the continuous current path approach for maximum accuracy.

use magrslib::export::visualization::VisualizationBuilder;
use magrslib::sources::currents::Polyline;
use magrslib::types::Vec3;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("Solenoid Coil Electromagnetic Field Visualization");
    println!("================================================");

    // Solenoid parameters
    let radius = 0.010;          // 10mm radius (20mm diameter)
    let length = 0.040;          // 40mm total length
    let num_turns = 10;          // 10 complete turns
    let current = 5.0;           // 5A current
    let points_per_turn = 20;    // Points per turn for smooth helix

    println!("Solenoid configuration:");
    println!("  Diameter: {:.1}mm", radius * 2000.0);
    println!("  Length: {:.1}mm", length * 1000.0);
    println!("  Number of turns: {}", num_turns);
    println!("  Current: {:.1}A", current);
    println!("  Pitch: {:.1}mm/turn", (length * 1000.0) / num_turns as f64);

    // Generate helical path vertices
    let total_points = num_turns * points_per_turn;
    let mut vertices = Vec::with_capacity(total_points);

    // Calculate pitch (distance between turns)
    let pitch = length / num_turns as f64;
    let angular_step = 2.0 * std::f64::consts::PI / points_per_turn as f64;

    for i in 0..total_points {
        let angle = (i as f64) * angular_step;
        let z_pos = -length / 2.0 + (i as f64 / points_per_turn as f64) * pitch;

        let vertex: Vec3 = [
            radius * angle.cos(),  // x = r * cos(θ)
            radius * angle.sin(),  // y = r * sin(θ)
            z_pos                  // z = linear progression along axis
        ];
        vertices.push(vertex);
    }

    println!("Generated {} vertices for helical path", vertices.len());
    println!("  Angular resolution: {:.1}° per segment", 360.0 / points_per_turn as f64);

    // Create the solenoid as a continuous current polyline
    let solenoid = Polyline::builder()
        .vertices(vertices)
        .current(current)
        .position([0.0, 0.0, 0.0])  // Centered at origin
        .build()?;

    println!("\n✓ Solenoid created as continuous helical current path");

    // Build visualization with enhanced capabilities
    let _visualization = VisualizationBuilder::new()
        .add_source(solenoid)
        .title("Solenoid Coil Electromagnetic Field Visualization")
        .grid_resolution(128)       // 128³ ≈ 2M points for detailed field
        .export_geometry(true)      // Export solenoid geometry
        .adaptive_resolution(true)  // Auto-adapt grid to solenoid bounds
        .export_to_directory("./output_solenoid")?;

    println!("\n✓ Visualization exported to ./output_solenoid/");
    println!("Files created:");
    println!("  - field.vtk     : Magnetic field data (B-vector field)");
    println!("  - geometry.vtk  : Solenoid coil geometry (helical curve)");
    println!("  - metadata.json : Complete scene description and field stats");

    println!("\nTo visualize:");
    println!("  python3 visualize_3d.py output_solenoid/metadata.json");
    println!("  python3 visualize_3d.py --batch output_solenoid/metadata.json -o solenoid.png");

    println!("\n=== EXPECTED FIELD CHARACTERISTICS ===");
    println!("Inside the solenoid:");
    println!("  ✓ Nearly uniform magnetic field parallel to axis");
    println!("  ✓ Field strength ≈ μ₀nI (n = turns per unit length)");
    println!("  ✓ Minimal radial field components");

    let turns_per_meter = num_turns as f64 / length;
    let mu_0 = 4.0 * std::f64::consts::PI * 1e-7;
    let theoretical_field = mu_0 * turns_per_meter * current;

    println!("  ✓ Theoretical internal field: B ≈ {:.2e} T", theoretical_field);

    println!("\nOutside the solenoid:");
    println!("  ✓ Dipole-like field pattern");
    println!("  ✓ Field lines close back through exterior space");
    println!("  ✓ Rapid field decay with distance");

    println!("\nPhysical insights:");
    println!("  • Solenoids are fundamental to electromagnets, inductors, relays");
    println!("  • Internal field uniformity makes them ideal for magnetic sensors");
    println!("  • Field energy storage: U = ½LI² where L ∝ μ₀n²V");
    println!("  • Self-inductance scales with turn density squared");

    println!("\n=== ENGINEERING APPLICATIONS ===");
    println!("This solenoid model applies to:");
    println!("  🔧 Electromagnetic actuators and solenoid valves");
    println!("  🔧 Inductor design for power electronics");
    println!("  🔧 MRI gradient coils and magnetic confinement");
    println!("  🔧 Magnetic field uniformity analysis");
    println!("  🔧 EMC and magnetic shielding studies");

    println!("\nNext steps for advanced modeling:");
    println!("  1. Multi-layer solenoids: Stack multiple helical polylines");
    println!("  2. Finite wire thickness: Add CurrentCircle per turn");
    println!("  3. Core materials: Modify permeability in field computation");
    println!("  4. Time-varying fields: I(t) → B(t) for transient analysis");

    Ok(())
}

/// Helper function to calculate theoretical solenoid field
#[allow(dead_code)]
fn theoretical_solenoid_field(turns_per_length: f64, current: f64) -> f64 {
    /// Internal field of an ideal infinite solenoid: B = μ₀nI
    let mu_0 = 4.0 * std::f64::consts::PI * 1e-7;
    mu_0 * turns_per_length * current
}

/// Helper function for finite solenoid correction
#[allow(dead_code)]
fn finite_solenoid_correction(length: f64, radius: f64, z_position: f64) -> f64 {
    /// Correction factor for finite solenoid length
    /// Returns multiplier for ideal field to account for end effects
    let z1 = (length / 2.0 - z_position) / radius;
    let z2 = (length / 2.0 + z_position) / radius;

    let factor1 = z1 / (1.0 + z1 * z1).sqrt();
    let factor2 = z2 / (1.0 + z2 * z2).sqrt();

    0.5 * (factor1 + factor2)
}