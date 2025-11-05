# B-Field Visualization Guide

## Overview

This guide explains how to visualize magnetic field data computed by magrslib using Python and matplotlib.

## What is the B-Field?

### Physical Meaning

The **B-field** (magnetic flux density) represents the actual magnetic field present in space. It's measured in **Tesla (T)**.

**Key Concepts:**
- **B = μ₀(H + M)** - Relationship between B, H, and magnetization M
- **B-field** includes contributions from both external sources and material magnetization
- **H-field** (magnetic field intensity, A/m) is the "external" applied field
- **Polarization J = μ₀M** (Tesla) is what our API uses for magnets

**Physical Interpretation:**
- Field lines show the direction magnetic forces act on moving charges
- Field magnitude indicates strength (1 Tesla is very strong!)
- For reference: Earth's field ≈ 50 μT, neodymium magnet surface ≈ 1-1.4 T

## Quick Start

### 1. Generate Field Data (Rust)

Run the provided example to compute B-field on a 2D grid:

```bash
cargo run --example visualize_cuboid
```

This creates `bfield_data.json` containing:
- 50×50 grid of observer positions in the xy-plane
- B-field vectors at each position
- Magnet geometry information

### 2. Visualize with Python

Make sure you have the required Python packages:

```bash
pip install matplotlib numpy
```

Run the visualization script:

```bash
python3 visualize_bfield.py
```

This creates:
- Interactive plot window showing the field
- Saved image: `bfield_visualization.png`

## Understanding the Visualization

### Streamplot Explanation

The visualization shows:
- **Streamlines**: Follow the direction of the B-field
- **Color**: log₁₀(|B|) - logarithmic scale of field magnitude
  - Bright colors = stronger field (near magnet)
  - Dark colors = weaker field (far from magnet)
- **Dashed box**: Magnet outline
- **Arrows**: Field direction

### Reading the Field

- **Near the magnet**: Strong, complex field patterns
- **Far from magnet**: Weaker, dipole-like pattern
- **Along z-axis** (for z-polarized magnet): Strongest field

## Creating Custom Visualizations

### Example: Different Magnet Configuration

```rust
use magrslib::export::{FieldData, GridInfo, MagnetInfo};
use magrslib::sources::magnets::Cuboid;
use magrslib::types::Observers;
use magrslib::get_b;

fn main() {
    // Create your magnet
    let magnet = Cuboid::builder()
        .dimension([0.03, 0.03, 0.01])      // 3cm × 3cm × 1cm
        .polarization([1.0, 0.0, 0.0])       // 1T in +x direction
        .position([0.0, 0.0, 0.0])
        .build()
        .unwrap();

    // Create observation grid
    let observers = Observers::from_grid_2d(
        (-0.08, 0.08, 80),  // x: -8cm to +8cm, 80 points
        (-0.08, 0.08, 80),  // y: -8cm to +8cm, 80 points
        "xy",               // xy-plane
        0.005               // at z = 5mm above magnet
    );

    // Compute field
    let b_field = get_b(&[magnet.into()], &observers);

    // Export
    let grid_info = GridInfo {
        plane: "xy".to_string(),
        nx: 80,
        ny: 80,
        x_range: (-0.08, 0.08),
        y_range: (-0.08, 0.08),
        fixed_coord: 0.005,
    };

    let mut data = FieldData::new(
        "Custom magnet field".to_string(),
        grid_info,
        observers.as_array(),
        &b_field,
    );

    data.add_magnet(MagnetInfo {
        magnet_type: "cuboid".to_string(),
        position: vec![0.0, 0.0, 0.0],
        dimensions: vec![0.03, 0.03, 0.01],
        polarization: vec![1.0, 0.0, 0.0],
    });

    data.export_to_json("custom_field.json").unwrap();
}
```

### Different Planes

You can visualize in any plane:

**XY-Plane** (top view):
```rust
Observers::from_grid_2d((-0.05, 0.05, 50), (-0.05, 0.05, 50), "xy", 0.0)
```

**XZ-Plane** (side view):
```rust
Observers::from_grid_2d((-0.05, 0.05, 50), (-0.05, 0.05, 50), "xz", 0.0)
```

**YZ-Plane** (front view):
```rust
Observers::from_grid_2d((-0.05, 0.05, 50), (-0.05, 0.05, 50), "yz", 0.0)
```

## Advanced: Multiple Magnets

```rust
// Create multiple magnets
let magnet1 = Cuboid::builder()
    .dimension([0.01, 0.01, 0.01])
    .polarization([0.0, 0.0, 1.0])
    .position([-0.015, 0.0, 0.0])  // Left magnet
    .build()
    .unwrap();

let magnet2 = Cuboid::builder()
    .dimension([0.01, 0.01, 0.01])
    .polarization([0.0, 0.0, -1.0])  // Opposite polarization!
    .position([0.015, 0.0, 0.0])     // Right magnet
    .build()
    .unwrap();

// Compute combined field
let b_field = get_b(&[magnet1.clone().into(), magnet2.clone().into()], &observers);

// Add both magnets to visualization
field_data.add_magnet(/* magnet1 info */);
field_data.add_magnet(/* magnet2 info */);
```

## Field Statistics

The visualization displays:
- **Min |B|**: Minimum field magnitude in the grid
- **Max |B|**: Maximum field magnitude (usually near magnet surface)
- **Mean |B|**: Average field magnitude

## Tips for Good Visualizations

1. **Grid Resolution**:
   - 50×50 = 2,500 points (fast, good for quick checks)
   - 100×100 = 10,000 points (detailed, slower)
   - 200×200 = 40,000 points (high-quality, slow)

2. **Grid Size**:
   - Include at least 2-3× magnet size on each side
   - Too small: Miss far-field behavior
   - Too large: Lose near-field detail

3. **Plane Selection**:
   - Use symmetry plane for clearest view
   - For z-polarized magnet, xy-plane at z=0 shows classic pattern

4. **Multiple Configurations**:
   - Visualize different heights (z-values) to see 3D structure
   - Compare different polarization directions

## Integration with Your Circuit Simulation

For your **circuits → currents → fields → parasitic currents** workflow:

1. Compute current distribution from circuit solver
2. Create Current source objects (Circle, Polyline)
3. Compute B-field at observation points
4. Use B-field to compute induced currents/voltages
5. Feed back into circuit solver

The visualization helps verify:
- Field patterns are physically reasonable
- Proper positioning of sources and observers
- Field magnitude in expected range

## Troubleshooting

**Issue**: "FileNotFoundError: bfield_data.json"
- **Solution**: Run `cargo run --example visualize_cuboid` first

**Issue**: Colors too bright/dark
- **Solution**: Adjust colormap in Python script (`cmap='plasma'` → `cmap='viridis'`)

**Issue**: Field looks wrong
- **Solution**: Check magnet polarization direction and position

**Issue**: Plot is slow
- **Solution**: Reduce grid resolution (nx, ny in Rust code)

## Next Steps

- Implement cylinder/sphere field functions for more geometries
- Add current sources for your circuit integration
- Create 3D visualizations using mayavi or plotly
- Animate time-varying fields for AC circuits
