# B-Field Visualization Issue - ROOT CAUSE & FIX

## Problem Reported

User reported that B-field streamplot visualization showed:
1. Field directions appearing almost random
2. Field lines getting LONGER with distance (should get shorter!)

## Root Cause Analysis

### What We Found

After adding comprehensive diagnostics, we discovered:

1. **✓ Rust computation is CORRECT** - No NaN/Inf, field decreases with distance
2. **✓ Array ordering is CORRECT** - Grid generation matches Python reshape
3. **✓ Field values are PHYSICALLY CORRECT** - All validations passed

### The ACTUAL Problem

**The issue was the choice of visualization plane, not a bug!**

In the **XY-plane at z=0** (original example):
```
Bx range: [-2.65e-17, 2.65e-17] T   ← Essentially ZERO!
By range: [-3.53e-17, 3.53e-17] T   ← Essentially ZERO!
Bz range: [-0.753, 0.753] T         ← All the field is here
```

**Why this happens (Physics)**:

For a **cuboid magnet with +z polarization**:
- Polarization vector points in +z direction (north pole on top)
- In the **equatorial plane** (xy-plane at z=0), by symmetry:
  - Field is almost entirely in the **z-direction**
  - In-plane components (Bx, By) are nearly zero

**Why streamplot looks broken**:
- `matplotlib.streamplot()` needs **in-plane** field components to draw streamlines
- With Bx ≈ 0 and By ≈ 0, it's trying to draw with ~10^-17 T vectors (numerical noise)
- This causes apparently random directions and strange scaling

## The Fix

### Use XZ or YZ Plane Instead!

In the **XZ-plane** (side view) at y=0:
```
Bx range: [-0.361, 0.361] T   ← Real field vectors!
Bz range: [-0.266, 0.618] T   ← Real field vectors!
```

**100% of points have non-zero in-plane components** → Beautiful streamlines!

## How to Visualize Properly

### Rule of Thumb

**Choose a plane that cuts through the field lines**, not parallel to them:

| Magnet Polarization | Good Visualization Planes | Bad Plane |
|---------------------|---------------------------|-----------|
| +z (vertical) | **XZ** (side), **YZ** (front) | XY (top) |
| +x (horizontal) | **XY** (top), **XZ** (side) | YZ (front) |
| +y (horizontal) | **XY** (top), **YZ** (front) | XZ (side) |

### Quick Fix for Existing Code

**Option 1**: Change plane in Rust example:
```rust
let observers = Observers::from_grid_2d(
    (grid_min, grid_max, grid_points),
    (grid_min, grid_max, grid_points),
    "xz",  // ← Change from "xy" to "xz"
    0.0
);

// Update GridInfo:
let grid_info = GridInfo {
    plane: "xz".to_string(),  // ← Change here too
    // ... rest stays the same
};
```

**Option 2**: Use the new example:
```bash
cargo run --example visualize_cuboid_xz
python3 visualize_xz.py
```

## Verification

Run diagnostics to check if your plane choice is good:

```rust
// Add after computing b_field:
let mut bx_nonzero = 0;
let mut by_nonzero = 0;
for i in 0..b_field.nrows() {
    if b_field[[i, 0]].abs() > 1e-10 { bx_nonzero += 1; }
    if b_field[[i, 1]].abs() > 1e-10 { by_nonzero += 1; }
}
println!("Plane components: Bx: {}, By: {}", bx_nonzero, by_nonzero);
```

**Goal**: Both components should be >80% non-zero for good visualization.

## Understanding the Physics

### Why is the XY-plane field mostly in Z-direction?

Imagine the magnetic field lines:
1. **Exit** from the north pole (top of magnet, +z)
2. **Loop around** through space
3. **Enter** back at the south pole (bottom of magnet, -z)

In the **equatorial plane** (XY at z=0):
- Field is transitioning from "going up" (near north) to "going down" (near south)
- In the middle region, by symmetry, there's no preferred x or y direction
- So field is predominantly z-direction

In a **side view** (XZ or YZ plane):
- You see the field lines **looping**
- Clear x-component (horizontal part of loop)
- Clear z-component (vertical part of loop)
- Perfect for visualization!

## Summary

| Aspect | Status |
|--------|--------|
| Rust field computation | ✓ **CORRECT** |
| Grid generation | ✓ **CORRECT** |
| Array ordering | ✓ **CORRECT** |
| Field physics | ✓ **CORRECT** |
| Original visualization plane choice | ✗ **POOR CHOICE** |
| **Solution** | **Use XZ or YZ plane** |

## Examples

### Bad Visualization (XY plane)
- Bx, By ≈ 10^-17 T (noise)
- Streamlines look random
- No physical insight

### Good Visualization (XZ plane)
- Bx, Bz ≈ 0.1-0.6 T (real vectors)
- Streamlines show field looping from N to S pole
- Clear physical structure

## For Your Circuit Simulation

When computing fields for circuit simulation:
- **Computation plane doesn't matter** - field values are correct in any plane
- **Visualization plane DOES matter** - choose wisely to see field structure
- Use XZ/YZ planes to understand how field couples to your circuit geometry
