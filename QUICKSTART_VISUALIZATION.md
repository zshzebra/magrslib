# Quick Start: B-Field Visualization

## TL;DR

```bash
# 1. Generate field data (Rust)
cargo run --example visualize_cuboid

# 2. Visualize (Python)
python3 visualize_bfield.py
```

## What You Get

A streamplot showing magnetic field lines from a cuboid magnet, colored by field strength.

## Understanding the B-Field

**B-field = Magnetic Flux Density (Tesla)**

Think of it as the "magnetic force field" in space:
- **Direction**: Arrows show which way the field points
- **Strength**: Color shows magnitude (bright = strong)
- **North pole**: Field lines emerge
- **South pole**: Field lines enter

For a magnet with +z polarization:
- Above magnet: Field points upward (away from north pole)
- Below magnet: Field points downward (toward south pole)
- Sides: Field loops around

## Customize It

### Change Magnet Size
Edit `examples/visualize_cuboid.rs`:
```rust
.dimension([0.04, 0.04, 0.02])  // 4cm x 4cm x 2cm instead of 2x2x1
```

### Change Polarization Direction
```rust
.polarization([1.0, 0.0, 0.0])  // +x direction instead of +z
```

### Change Grid Resolution
```rust
let grid_points = 100;  // 100x100 = 10,000 points (more detailed)
```

### Different Viewing Plane
```rust
// XZ-plane (side view) instead of XY
let observers = Observers::from_grid_2d(
    (grid_min, grid_max, grid_points),
    (grid_min, grid_max, grid_points),
    "xz",  // ← change this
    0.0
);

// Update grid_info too:
plane: "xz".to_string(),
```

## For Your Circuit Simulation

Your workflow: **circuits → currents → fields → parasitic currents**

Current status:
- ✅ Can compute B-field from magnets
- ✅ Can visualize field patterns
- 🔄 Next: Add current sources (for circuit integration)
- 🔄 Then: Compute induced fields/voltages

Use visualization to:
1. Verify magnet placement is correct
2. Check field magnitude is reasonable
3. Identify regions of strong/weak field
4. Plan where to place sensors/measurements

## Python Script Customization

Edit `visualize_bfield.py` to change:

**Colormap**:
```python
stream = ax.streamplot(X, Y, Bx, By, color=color, cmap='viridis')  # or 'plasma', 'inferno', 'coolwarm'
```

**Streamline Density**:
```python
stream = ax.streamplot(X, Y, Bx, By, density=2.5)  # more lines (default: 1.5)
```

**Figure Size**:
```python
fig, ax = plt.subplots(figsize=(12, 10))  # larger plot
```

## Typical Field Values

To calibrate your intuition:
- **1 T** = Very strong (neodymium magnet surface, MRI scanner)
- **0.1 T** = Strong (small permanent magnet at surface)
- **0.01 T** = Moderate (1 cm away from 1cm magnet)
- **0.001 T (1 mT)** = Weak but measurable
- **0.00005 T (50 μT)** = Earth's magnetic field
- **< 1 μT** = Very weak (far from sources)

## Quick Checks

**Sanity check your results**:
1. Field should be strongest near magnet surface
2. Field should decrease with distance (roughly 1/r³ for dipole)
3. Field direction should follow intuitive north→south pattern
4. Symmetry: Field pattern should match magnet symmetry

**If something looks wrong**:
- Check polarization vector is normalized (or use magnetization with μ₀)
- Verify position units are in meters
- Ensure grid includes enough area around magnet
- Check that plane selection makes sense for your geometry
