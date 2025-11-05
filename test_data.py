#!/usr/bin/env python3
"""
Test script to verify the JSON data structure and field values
"""

import json
import numpy as np

def main():
    # Load data
    with open('bfield_data.json', 'r') as f:
        data = json.load(f)

    print("="*60)
    print("DATA STRUCTURE VERIFICATION")
    print("="*60)

    # Extract data
    grid = data['grid']
    positions = np.array(data['positions'])
    b_field = np.array(data['b_field'])

    nx, ny = grid['nx'], grid['ny']

    print(f"\nGrid: {nx} x {ny} = {nx*ny} points")
    print(f"Data shapes:")
    print(f"  positions: {positions.shape}")
    print(f"  b_field: {b_field.shape}")

    # Check grid ordering
    print(f"\nFirst 5 positions:")
    for i in range(min(5, len(positions))):
        print(f"  {i}: ({positions[i,0]:.4f}, {positions[i,1]:.4f}, {positions[i,2]:.4f})")

    # Reshape
    X = positions[:, 0].reshape(ny, nx)
    Y = positions[:, 1].reshape(ny, nx)
    Bx = b_field[:, 0].reshape(ny, nx)
    By = b_field[:, 1].reshape(ny, nx)
    Bz = b_field[:, 2].reshape(ny, nx)

    print(f"\nReshaped arrays:")
    print(f"  X shape: {X.shape}")
    print(f"  Y shape: {Y.shape}")

    # Verify grid structure
    print(f"\nGrid structure check:")
    print(f"  X[0,:5] (first row, first 5 x-values): {X[0,:5]}")
    print(f"  Y[:5,0] (first col, first 5 y-values): {Y[:5,0]}")
    print(f"  X varies along axis 1 (columns): {np.allclose(X[:, 0], X[:, 1])=}")  # Should be False
    print(f"  Y varies along axis 0 (rows): {np.allclose(Y[0, :], Y[1, :])=}")     # Should be False

    # Field magnitude
    B_mag = np.sqrt(Bx**2 + By**2 + Bz**2)

    print(f"\nField magnitude statistics:")
    print(f"  Min:  {B_mag.min():.3e} T")
    print(f"  Max:  {B_mag.max():.3e} T")
    print(f"  Mean: {B_mag.mean():.3e} T")

    # Check field at center
    center_i, center_j = ny // 2, nx // 2
    print(f"\nField at grid center [{center_i},{center_j}]:")
    print(f"  Position: ({X[center_i, center_j]:.4f}, {Y[center_i, center_j]:.4f}, 0)")
    print(f"  Bx = {Bx[center_i, center_j]:.6e} T")
    print(f"  By = {By[center_i, center_j]:.6e} T")
    print(f"  Bz = {Bz[center_i, center_j]:.6e} T")
    print(f"  |B| = {B_mag[center_i, center_j]:.6e} T")

    # Check field at corners
    print(f"\nField at corners:")
    corners = [(0, 0), (0, -1), (-1, 0), (-1, -1)]
    for i, j in corners:
        print(f"  [{i:2},{j:2}] at ({X[i,j]:6.3f}, {Y[i,j]:6.3f}): |B| = {B_mag[i,j]:.3e} T")

    # Check that field decreases with distance
    center_pos = np.array([X[center_i, center_j], Y[center_i, center_j], 0])
    distances = np.sqrt((X - center_pos[0])**2 + (Y - center_pos[1])**2)

    print(f"\nField vs distance correlation:")
    # Sample a few distance ranges
    for r_min, r_max in [(0, 0.01), (0.01, 0.02), (0.02, 0.03), (0.03, 0.05)]:
        mask = (distances >= r_min) & (distances < r_max)
        if np.any(mask):
            mean_b = B_mag[mask].mean()
            print(f"  Distance {r_min:.3f} to {r_max:.3f} m: mean |B| = {mean_b:.3e} T")

    # Check field direction consistency
    print(f"\nField direction check (Bz component):")
    print(f"  Bz > 0 (upward): {np.sum(Bz > 0)} points ({100*np.sum(Bz > 0)/Bz.size:.1f}%)")
    print(f"  Bz < 0 (downward): {np.sum(Bz < 0)} points ({100*np.sum(Bz < 0)/Bz.size:.1f}%)")
    print(f"  Bz ≈ 0: {np.sum(np.abs(Bz) < 1e-10)} points")

    # Simple visual check - print a small region around center
    print(f"\nBz values in 5x5 region around center:")
    i0, i1 = center_i - 2, center_i + 3
    j0, j1 = center_j - 2, center_j + 3
    print(f"  (showing Bz in units of mT)")
    for i in range(i0, i1):
        row = " ".join([f"{Bz[i,j]*1000:6.1f}" for j in range(j0, j1)])
        marker = " <--" if i == center_i else ""
        print(f"  {row}{marker}")

    print("\n" + "="*60)

if __name__ == '__main__':
    main()
