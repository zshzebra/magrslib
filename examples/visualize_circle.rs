//! Example: Generate B-field from circular current loop for visualization
//!
//! This example demonstrates:
//! - Creating a circular current loop
//! - Computing the B-field on a 3D grid using Biot-Savart law and elliptic integrals
//! - Exporting to VTK format with metadata for advanced 3D visualization
//! - Using the unified visualization API
//!
//! The output can be visualized with:
//! - PyVista (Python): Interactive 3D streamlines, current loop geometry
//! - Paraview: Professional scientific visualization
//! - Any VTK-compatible tool
//!
//! This example is ideal for understanding electromagnetic fields around current-carrying
//! conductors and serves as foundation for modeling inductors and transformers.

use magrslib::export::visualization::VisualizationBuilder;
use magrslib::sources::currents::Circle;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("===============================================");
    println!("  Circle Current Loop B-Field Visualization");
    println!("===============================================\n");

    // Create a circular current loop
    let current_loop = Circle::builder()
        .diameter(0.020)  // 20mm diameter (10mm radius)
        .current(5.0)     // 5 Amperes current
        .position([0.0, 0.0, 0.0])  // Centered at origin
        .build()?;

    println!("Current loop configuration:");
    println!("  Type: Circular current loop");
    println!("  Diameter: 20mm (10mm radius)");
    println!("  Current: 5.0 A (counter-clockwise when viewed from +z)");
    println!("  Position: origin");
    println!("  Orientation: xy-plane, normal in +z direction\n");

    // Create visualization with adaptive grid
    println!("Setting up visualization with unified API...");
    let visualization = VisualizationBuilder::new()
        .add_source(current_loop)
        .grid_resolution(48)  // 48³ = 110,592 points - good balance of detail vs speed
        .adaptive_resolution(true)  // Auto-adjust grid around the current loop
        .export_geometry(true)      // Export current loop geometry for visualization
        .title("B-field from 20mm diameter current loop carrying 5A")
        .export_to_directory("./output/circle_field")?;

    println!();
    println!("=== FIELD ANALYSIS ===\n");

    // Display field statistics
    let stats = &visualization.field_stats;
    println!("Field computation completed!");
    println!("  Computation time: {:.2}s", stats.computation_time);
    println!("  Total grid points: {}", visualization.grid.n_points);
    println!("  Grid resolution: {}³", visualization.grid.resolution[0]);
    println!("  Grid bounds: [{:.1}, {:.1}] mm",
             visualization.grid.bounds[0] * 1000.0,
             visualization.grid.bounds[1] * 1000.0);

    if stats.invalid_count > 0 {
        println!("  ⚠️  Invalid values: {} (replaced with zeros for visualization)", stats.invalid_count);
    } else {
        println!("  ✓ All field values are finite");
    }

    println!("\nField magnitude statistics:");
    println!("  Min:    {:.3e} T", stats.magnitude.min);
    println!("  Max:    {:.3e} T", stats.magnitude.max);
    println!("  Mean:   {:.3e} T", stats.magnitude.mean);
    println!("  Median: {:.3e} T", stats.magnitude.median);
    println!("  Std:    {:.3e} T", stats.magnitude.std_dev);

    // Component analysis
    println!("\nField components statistics:");
    let components = ["Bx", "By", "Bz"];
    for (i, comp) in components.iter().enumerate() {
        let comp_stats = &stats.components[i];
        println!("  {} range: [{:.3e}, {:.3e}] T", comp, comp_stats.min, comp_stats.max);
    }

    println!("\n=== THEORETICAL VALIDATION ===\n");

    // Theoretical field at center: B = μ₀I/(2R)
    let radius = 0.010;  // 10mm radius
    let current = 5.0;   // 5A
    let mu_0 = 4.0 * std::f64::consts::PI * 1e-7;  // μ₀ in H/m
    let theoretical_center_field = mu_0 * current / (2.0 * radius);

    println!("Theoretical field at loop center:");
    println!("  B_center = μ₀I/(2R) = {:.3e} T", theoretical_center_field);
    println!("  (Using B = 4π×10⁻⁷ × 5.0 / (2 × 0.01) = {:.3e} T)", theoretical_center_field);

    // Far-field dipole approximation check
    let magnetic_moment = current * std::f64::consts::PI * radius * radius;
    println!("\nMagnetic dipole moment:");
    println!("  m = IA = {:.3e} A⋅m²", magnetic_moment);

    println!("\nFar-field characteristics:");
    println!("  Field follows dipole pattern at distances >> {:.1}mm", radius * 1000.0);
    println!("  Expected axial field decay: ∝ 1/r³");
    println!("  Expected radial field pattern with positive/negative lobes");

    println!("\n=== ELECTROMAGNETIC INSIGHTS ===\n");

    println!("Physical interpretation:");
    println!("  • Current flows counter-clockwise (right-hand rule: B points +z inside loop)");
    println!("  • Field lines form closed loops around the current path");
    println!("  • Maximum field occurs at loop center: {:.1} mT", theoretical_center_field * 1000.0);
    println!("  • Field strength scales linearly with current: B ∝ I");
    println!("  • This geometry is fundamental to inductors, solenoids, and transformers");

    println!("\nField pattern expectations:");
    println!("  • Axial symmetry around z-axis");
    println!("  • Field strength decreases with distance from loop");
    println!("  • Characteristic 'toroidal' field line pattern");
    println!("  • Zero radial field on the z-axis (axial symmetry)");

    println!("\n===============================================");
    println!("  Visualization Instructions");
    println!("===============================================\n");

    println!("To visualize the 3D electromagnetic field:");
    println!("  cd output/circle_field");
    println!("  python3 ../../visualize_3d.py\n");

    println!("The visualization will show:");
    println!("  • Magnetic field streamlines following Biot-Savart law");
    println!("  • Color-coded field magnitude (Tesla)");
    println!("  • Current loop geometry (20mm diameter ring)");
    println!("  • Interactive 3D exploration with rotation/zoom");

    println!("\nAdvanced analysis options:");
    println!("  • Compare with analytical solutions");
    println!("  • Study field patterns at different current values");
    println!("  • Investigate near-field vs far-field behavior");
    println!("  • Model multi-turn coils by scaling current");

    println!("\nFiles generated:");
    println!("  • field.vtk:     B-field data (vector field)");
    println!("  • geometry.vtk:  Current loop geometry (planned)");
    println!("  • metadata.json: Scene description and statistics\n");

    println!("Next steps for inductor modeling:");
    println!("  1. Multi-turn coils: Stack multiple Circle sources");
    println!("  2. Time-varying currents: I(t) → B(t) field evolution");
    println!("  3. Electromagnetic induction: ∇×E = -∂B/∂t");
    println!("  4. Energy storage: U = ½LI² in magnetic field");

    Ok(())
}

// Helper functions for theoretical analysis (could be moved to utils)
#[allow(dead_code)]
fn theoretical_axial_field(current: f64, radius: f64, z: f64) -> f64 {
    /// Axial field from circular current loop: B_z = (μ₀I a²) / (2(z² + a²)^(3/2))
    let mu_0 = 4.0 * std::f64::consts::PI * 1e-7;
    let denominator = (z * z + radius * radius).powf(1.5);
    mu_0 * current * radius * radius / (2.0 * denominator)
}

#[allow(dead_code)]
fn theoretical_dipole_field(magnetic_moment: f64, r: f64, theta: f64) -> (f64, f64) {
    /// Far-field dipole approximation: B_r ∝ cos(θ)/r³, B_θ ∝ sin(θ)/r³
    let mu_0 = 4.0 * std::f64::consts::PI * 1e-7;
    let prefactor = mu_0 * magnetic_moment / (4.0 * std::f64::consts::PI * r.powi(3));

    let b_r = prefactor * 2.0 * theta.cos();      // Radial component
    let b_theta = prefactor * theta.sin();        // Polar component

    (b_r, b_theta)
}