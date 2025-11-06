//! PWM-driven inductor electromagnetic field visualization
//!
//! Demonstrates the complete PWM → Physical Inductor → Current Visualization pipeline:
//! - PWM signal generation with configurable frequency, duty cycle, and amplitude
//! - Physical solenoid model using helical CurrentPolyline
//! - Time-varying field computation with realistic switching behavior
//! - Current-to-color visualization mapping showing field evolution
//!
//! This example shows how to model inductors in switching power supplies,
//! motor drives, and other power electronics applications where current
//! varies dynamically with time.

use magrslib::export::visualization::VisualizationBuilder;
use magrslib::sources::currents::Polyline;
use magrslib::sources::current_functions::PWMCurrent;
use magrslib::types::{Observers, Vec3};
use magrslib::{get_b_time_series, get_b_at_time};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("PWM-Driven Inductor Electromagnetic Field Visualization");
    println!("======================================================");

    // PWM signal parameters
    let pwm_frequency = 50_000.0;     // 50kHz switching frequency (typical for power supplies)
    let duty_cycle = 0.6;             // 60% duty cycle
    let peak_current = 15.0;          // 15A peak current

    // Solenoid parameters (power inductor specification)
    let radius = 0.008;               // 8mm radius (16mm diameter)
    let length = 0.020;               // 20mm total length
    let num_turns = 15;               // 15 turns for higher inductance
    let points_per_turn = 16;         // 16 points per turn (smooth enough)

    println!("PWM Signal Configuration:");
    println!("  Frequency: {:.1}kHz", pwm_frequency / 1000.0);
    println!("  Duty Cycle: {:.1}%", duty_cycle * 100.0);
    println!("  Peak Current: {:.1}A", peak_current);
    println!("  Period: {:.1}μs", 1e6 / pwm_frequency);

    println!("\nInductor Configuration:");
    println!("  Diameter: {:.1}mm", radius * 2000.0);
    println!("  Length: {:.1}mm", length * 1000.0);
    println!("  Number of turns: {}", num_turns);
    println!("  Inductance per turn: L ∝ μ₀n²A/l");

    // Generate helical solenoid geometry
    let total_points = num_turns * points_per_turn;
    let mut vertices = Vec::with_capacity(total_points);

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

    println!("\nGenerated {} vertices for helical inductor geometry", vertices.len());

    // Create PWM current function
    let pwm_current = Box::new(PWMCurrent::new(pwm_frequency, duty_cycle, peak_current));

    println!("✓ PWM current function created with ideal switching (instantaneous rise/fall)");

    // Create the PWM-driven solenoid inductor
    let pwm_inductor = Polyline::builder()
        .vertices(vertices)
        .current_function(pwm_current)
        .position([0.0, 0.0, 0.0])  // Centered at origin
        .build()?;

    println!("✓ PWM-driven solenoid inductor created");

    // Test single time point first (at peak current)
    let observers = Observers::from_vec3s(vec![[0.0, 0.0, 0.0]]);  // Center of inductor

    let field_at_peak = get_b_at_time(&[pwm_inductor.clone().into()], &observers, 0.0);
    let field_at_zero = get_b_at_time(&[pwm_inductor.clone().into()], &observers, duty_cycle / pwm_frequency);

    println!("\nField validation:");
    println!("  B-field at PWM peak (t=0): [{:.2e}, {:.2e}, {:.2e}] T",
             field_at_peak[[0,0]], field_at_peak[[0,1]], field_at_peak[[0,2]]);
    println!("  B-field at PWM zero (t={}μs): [{:.2e}, {:.2e}, {:.2e}] T",
             (duty_cycle / pwm_frequency) * 1e6,
             field_at_zero[[0,0]], field_at_zero[[0,1]], field_at_zero[[0,2]]);

    // Compute time series over one PWM cycle
    let period = 1.0 / pwm_frequency;           // PWM period
    let time_step = period / 100.0;             // 100 points per cycle
    let num_cycles = 2.0;                       // Show 2 full cycles

    println!("\nComputing field time series:");
    println!("  Time span: {:.1} PWM cycles ({:.1}μs)", num_cycles, num_cycles * period * 1e6);
    println!("  Time resolution: {:.1}ns per step", time_step * 1e9);
    println!("  Total time steps: {}", (num_cycles * period / time_step) as usize);

    let time_series = get_b_time_series(
        &[pwm_inductor.clone().into()],
        &observers,
        0.0,                    // Start time
        num_cycles * period,    // End time
        time_step               // Time step
    );

    println!("✓ Field time series computed ({} time points)", time_series.len());

    // Analyze field variation
    let field_magnitudes: Vec<f64> = time_series.iter()
        .map(|(_, field)| {
            let bx = field[[0, 0]];
            let by = field[[0, 1]];
            let bz = field[[0, 2]];
            (bx*bx + by*by + bz*bz).sqrt()
        })
        .collect();

    let max_field = field_magnitudes.iter().fold(0.0_f64, |a, &b| a.max(b));
    let min_field = field_magnitudes.iter().fold(f64::INFINITY, |a, &b| a.min(b));
    let avg_field = field_magnitudes.iter().sum::<f64>() / field_magnitudes.len() as f64;

    println!("\nField magnitude statistics:");
    println!("  Maximum: {:.2e} T (at peak current)", max_field);
    println!("  Minimum: {:.2e} T (at zero current)", min_field);
    println!("  Average: {:.2e} T", avg_field);
    println!("  Dynamic range: {:.1}:1", max_field / min_field);

    // Build enhanced visualization with time-series data
    let _visualization = VisualizationBuilder::new()
        .add_source(pwm_inductor)
        .title("PWM-Driven Inductor Electromagnetic Field Evolution")
        .grid_resolution(64)         // Smaller grid for faster animation
        .export_geometry(true)       // Export inductor geometry
        .adaptive_resolution(true)   // Auto-adapt to inductor bounds
        .export_to_directory("./output_pwm_inductor")?;

    println!("\n✓ Static field visualization exported to ./output_pwm_inductor/");
    println!("Files created:");
    println!("  - field.vtk     : Instantaneous B-field at t=0 (peak current)");
    println!("  - geometry.vtk  : PWM inductor geometry (helical coil)");
    println!("  - metadata.json : Complete scene description and PWM parameters");

    // Export time series data for custom animation
    write_time_series_data(&time_series, "./output_pwm_inductor/time_series.csv")?;
    write_pwm_metadata(pwm_frequency, duty_cycle, peak_current, num_cycles, "./output_pwm_inductor/pwm_parameters.json")?;

    println!("  - time_series.csv    : Field evolution data for animation");
    println!("  - pwm_parameters.json : PWM signal parameters");

    println!("\nTo visualize:");
    println!("  python3 visualize_3d.py output_pwm_inductor/metadata.json");
    println!("  python3 visualize_3d.py --batch output_pwm_inductor/metadata.json -o pwm_inductor.png");

    println!("\n=== PWM INDUCTOR FIELD CHARACTERISTICS ===");
    println!("During PWM ON phase (current = {:.1}A):", peak_current);
    println!("  ✓ Strong uniform magnetic field parallel to inductor axis");
    println!("  ✓ Field magnitude proportional to instantaneous current");
    println!("  ✓ Energy stored: U = ½LI² (maximum)");

    println!("\nDuring PWM OFF phase (current = 0A):");
    println!("  ✓ Magnetic field collapses to zero");
    println!("  ✓ Stored energy released (back-EMF if load present)");
    println!("  ✓ Field transitions create dB/dt for EMI analysis");

    println!("\nSwitching behavior:");
    println!("  ✓ Field rise/fall tracks current waveform exactly");
    println!("  ✓ No hysteresis (ideal linear inductor model)");
    println!("  ✓ Field energy oscillates at PWM frequency");

    println!("\n=== POWER ELECTRONICS INSIGHTS ===");
    println!("This PWM inductor model applies to:");
    println!("  🔧 Buck/boost converter output inductors");
    println!("  🔧 Motor drive current control inductors");
    println!("  🔧 LED driver current regulation");
    println!("  🔧 Battery charger current limiting");
    println!("  🔧 EMI filter design and analysis");

    let theoretical_inductance = estimate_solenoid_inductance(num_turns, radius, length);
    let energy_stored = 0.5 * theoretical_inductance * peak_current * peak_current;

    println!("\nQuantitative analysis:");
    println!("  Estimated inductance: L ≈ {:.1}μH", theoretical_inductance * 1e6);
    println!("  Peak energy stored: U = {:.1}μJ", energy_stored * 1e6);
    println!("  Energy oscillation frequency: {:.1}kHz", pwm_frequency / 1000.0);
    println!("  Power handling: P ∝ LI²f = {:.1}mW", energy_stored * pwm_frequency * 1000.0);

    println!("\nNext steps for advanced modeling:");
    println!("  1. Add finite rise/fall times: PWMCurrent::with_transitions()");
    println!("  2. Include core saturation: B(H) nonlinearity");
    println!("  3. Model parasitic capacitance: LC resonance effects");
    println!("  4. Temperature dependence: R(T) and core losses");
    println!("  5. Couple with circuit simulator: spiceng integration");

    Ok(())
}

/// Estimate solenoid inductance using standard formula
fn estimate_solenoid_inductance(turns: usize, radius: f64, length: f64) -> f64 {
    // L = μ₀ * n² * A / l
    // where n = number of turns, A = cross-sectional area, l = length
    let mu_0 = 4.0 * std::f64::consts::PI * 1e-7;  // H/m
    let area = std::f64::consts::PI * radius * radius;  // m²
    let turns_per_length = turns as f64 / length;  // turns/m

    mu_0 * (turns_per_length * turns_per_length) * area * length
}

/// Write time series data to CSV for external analysis/animation
fn write_time_series_data(
    time_series: &[(f64, ndarray::Array2<f64>)],
    filename: &str
) -> Result<(), Box<dyn std::error::Error>> {
    use std::fs::File;
    use std::io::Write;

    let mut file = File::create(filename)?;
    writeln!(file, "time_s,Bx_T,By_T,Bz_T,magnitude_T")?;

    for (time, field) in time_series {
        let bx = field[[0, 0]];
        let by = field[[0, 1]];
        let bz = field[[0, 2]];
        let magnitude = (bx*bx + by*by + bz*bz).sqrt();

        writeln!(file, "{:.6e},{:.6e},{:.6e},{:.6e},{:.6e}",
                 time, bx, by, bz, magnitude)?;
    }

    Ok(())
}

/// Write PWM parameters to JSON for reference
fn write_pwm_metadata(
    frequency: f64,
    duty_cycle: f64,
    peak_current: f64,
    num_cycles: f64,
    filename: &str
) -> Result<(), Box<dyn std::error::Error>> {
    use std::fs::File;
    use std::io::Write;

    let mut file = File::create(filename)?;
    let json_content = format!(r#"{{
  "pwm_parameters": {{
    "frequency_hz": {:.1},
    "frequency_khz": {:.1},
    "duty_cycle": {:.2},
    "duty_cycle_percent": {:.1},
    "peak_current_a": {:.1},
    "period_us": {:.2},
    "on_time_us": {:.2},
    "off_time_us": {:.2}
  }},
  "simulation_parameters": {{
    "num_cycles": {:.1},
    "total_time_us": {:.2},
    "time_resolution_ns": {:.1}
  }},
  "generated_by": "magrslib PWM inductor example",
  "description": "PWM-driven solenoid inductor electromagnetic field simulation"
}}"#,
        frequency,
        frequency / 1000.0,
        duty_cycle,
        duty_cycle * 100.0,
        peak_current,
        1e6 / frequency,
        duty_cycle * 1e6 / frequency,
        (1.0 - duty_cycle) * 1e6 / frequency,
        num_cycles,
        num_cycles * 1e6 / frequency,
        (1e6 / frequency) / 100.0  // 100 points per cycle
    );

    file.write_all(json_content.as_bytes())?;
    Ok(())
}