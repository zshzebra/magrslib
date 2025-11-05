#!/usr/bin/env python3
"""
Visualize B-field data exported from magrslib

This script loads JSON data exported by the Rust library and creates
a streamplot visualization similar to the magpylib example.
"""

import json
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Circle
from numpy.linalg import norm

def load_field_data(filename):
    """Load field data from JSON file"""
    with open(filename, 'r') as f:
        return json.load(f)

def reshape_grid_data(data):
    """Reshape flattened data back to 2D grid format"""
    grid = data['grid']
    nx, ny = grid['nx'], grid['ny']

    # Extract positions
    positions = np.array(data['positions'])  # shape: (n_points, 3)

    # Extract B-field components
    b_field = np.array(data['b_field'])      # shape: (n_points, 3)

    # Determine which components to use based on plane
    plane = grid['plane']
    if plane == 'xy':
        X = positions[:, 0].reshape(ny, nx)
        Y = positions[:, 1].reshape(ny, nx)
        Bx = b_field[:, 0].reshape(ny, nx)
        By = b_field[:, 1].reshape(ny, nx)
        Bz = b_field[:, 2].reshape(ny, nx)
        xlabel, ylabel = 'x-position (m)', 'y-position (m)'
    elif plane == 'xz':
        X = positions[:, 0].reshape(ny, nx)
        Y = positions[:, 2].reshape(ny, nx)
        Bx = b_field[:, 0].reshape(ny, nx)
        By = b_field[:, 2].reshape(ny, nx)
        Bz = b_field[:, 1].reshape(ny, nx)
        xlabel, ylabel = 'x-position (m)', 'z-position (m)'
    elif plane == 'yz':
        X = positions[:, 1].reshape(ny, nx)
        Y = positions[:, 2].reshape(ny, nx)
        Bx = b_field[:, 1].reshape(ny, nx)
        By = b_field[:, 2].reshape(ny, nx)
        Bz = b_field[:, 0].reshape(ny, nx)
        xlabel, ylabel = 'y-position (m)', 'z-position (m)'
    else:
        raise ValueError(f"Unknown plane: {plane}")

    return X, Y, Bx, By, Bz, xlabel, ylabel

def plot_magnet_outline(ax, magnet, plane):
    """Add magnet outline to the plot"""
    mag_type = magnet['magnet_type']
    pos = magnet['position']
    dims = magnet['dimensions']

    if plane == 'xy':
        x, y = pos[0], pos[1]
        if mag_type == 'cuboid':
            width, height = dims[0], dims[1]
            rect = Rectangle((x - width/2, y - height/2), width, height,
                           linewidth=2, edgecolor='black', facecolor='none',
                           linestyle='--', label='Magnet outline')
            ax.add_patch(rect)
        elif mag_type in ['cylinder', 'sphere']:
            radius = dims[0] / 2
            circle = Circle((x, y), radius, linewidth=2, edgecolor='black',
                          facecolor='none', linestyle='--', label='Magnet outline')
            ax.add_patch(circle)
    elif plane == 'xz':
        x, z = pos[0], pos[2]
        if mag_type == 'cuboid':
            width, depth = dims[0], dims[2]
            rect = Rectangle((x - width/2, z - depth/2), width, depth,
                           linewidth=2, edgecolor='black', facecolor='none',
                           linestyle='--', label='Magnet outline')
            ax.add_patch(rect)
        elif mag_type == 'cylinder':
            diameter, height = dims[0], dims[1]
            rect = Rectangle((x - diameter/2, z - height/2), diameter, height,
                           linewidth=2, edgecolor='black', facecolor='none',
                           linestyle='--', label='Magnet outline')
            ax.add_patch(rect)
        elif mag_type == 'sphere':
            radius = dims[0] / 2
            circle = Circle((x, z), radius, linewidth=2, edgecolor='black',
                          facecolor='none', linestyle='--', label='Magnet outline')
            ax.add_patch(circle)
    elif plane == 'yz':
        y, z = pos[1], pos[2]
        if mag_type == 'cuboid':
            height, depth = dims[1], dims[2]
            rect = Rectangle((y - height/2, z - depth/2), height, depth,
                           linewidth=2, edgecolor='black', facecolor='none',
                           linestyle='--', label='Magnet outline')
            ax.add_patch(rect)
        elif mag_type == 'cylinder':
            diameter, height = dims[0], dims[1]
            rect = Rectangle((y - diameter/2, z - height/2), diameter, height,
                           linewidth=2, edgecolor='black', facecolor='none',
                           linestyle='--', label='Magnet outline')
            ax.add_patch(rect)
        elif mag_type == 'sphere':
            radius = dims[0] / 2
            circle = Circle((y, z), radius, linewidth=2, edgecolor='black',
                          facecolor='none', linestyle='--', label='Magnet outline')
            ax.add_patch(circle)

def main():
    # Load data
    print("Loading field data from bfield_data.json...")
    data = load_field_data('bfield_data.json')

    print(f"Description: {data['description']}")
    print(f"Grid: {data['grid']['nx']}x{data['grid']['ny']} points in {data['grid']['plane']}-plane")
    print(f"Number of magnets: {len(data['magnets'])}")

    # Reshape data for plotting
    X, Y, Bx, By, Bz, xlabel, ylabel = reshape_grid_data(data)

    # Compute field magnitude
    B_magnitude = np.sqrt(Bx**2 + By**2 + Bz**2)

    # Create figure
    fig, ax = plt.subplots(figsize=(10, 9))

    # Create streamplot
    # Use log scale for color to see both near and far field
    # Add small epsilon to avoid log(0)
    color = np.log10(B_magnitude + 1e-10)

    stream = ax.streamplot(X, Y, Bx, By, color=color, cmap='plasma',
                          density=1.5, linewidth=1, arrowsize=1.5)

    # Add colorbar
    cbar = plt.colorbar(stream.lines, ax=ax, label='log₁₀(|B|) [log₁₀(T)]')

    # Plot magnet outlines
    for magnet in data['magnets']:
        plot_magnet_outline(ax, magnet, data['grid']['plane'])

    # Labels and title
    ax.set_xlabel(xlabel, fontsize=12)
    ax.set_ylabel(ylabel, fontsize=12)
    ax.set_title(f"B-Field Visualization\n{data['description']}", fontsize=14)
    ax.set_aspect('equal')

    # Set axis limits from grid range
    grid = data['grid']
    ax.set_xlim(grid['x_range'])
    ax.set_ylim(grid['y_range'])

    # Add grid
    ax.grid(True, alpha=0.3, linestyle=':')

    # Add statistics text
    stats_text = f"Min |B|: {B_magnitude.min():.2e} T\n"
    stats_text += f"Max |B|: {B_magnitude.max():.2e} T\n"
    stats_text += f"Mean |B|: {B_magnitude.mean():.2e} T"
    ax.text(0.02, 0.98, stats_text, transform=ax.transAxes,
            fontsize=10, verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

    plt.tight_layout()

    # Save figure
    output_filename = 'bfield_visualization.png'
    plt.savefig(output_filename, dpi=150, bbox_inches='tight')
    print(f"\nVisualization saved to {output_filename}")

    # Show interactive plot
    print("Displaying interactive plot...")
    plt.show()

if __name__ == '__main__':
    main()
