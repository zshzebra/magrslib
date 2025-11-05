import matplotlib.pyplot as plt
import numpy as np
import json
from matplotlib.patches import Rectangle

# Load data
with open('bfield_data_xz.json', 'r') as f:
    data = json.load(f)

grid = data['grid']
positions = np.array(data['positions'])
b_field = np.array(data['b_field'])

nx, ny = grid['nx'], grid['ny']

# Reshape for XZ plane
X = positions[:, 0].reshape(ny, nx)  # x-coordinates
Z = positions[:, 2].reshape(ny, nx)  # z-coordinates  
Bx = b_field[:, 0].reshape(ny, nx)
Bz = b_field[:, 2].reshape(ny, nx)

B_mag = np.sqrt(Bx**2 + Bz**2 + b_field[:, 1].reshape(ny, nx)**2)

print(f"Bx range: [{Bx.min():.3e}, {Bx.max():.3e}]")
print(f"Bz range: [{Bz.min():.3e}, {Bz.max():.3e}]")
print(f"|B| range: [{B_mag.min():.3e}, {B_mag.max():.3e}]")

# Create plot
fig, ax = plt.subplots(figsize=(10, 9))

# Streamplot with log color scale
color = np.log10(B_mag + 1e-10)
stream = ax.streamplot(X, Z, Bx, Bz, color=color, cmap='plasma',
                      density=2.0, linewidth=1, arrowsize=1.5)

plt.colorbar(stream.lines, ax=ax, label='log₁₀(|B|) [log₁₀(T)]')

# Add magnet outline (XZ plane view)
magnet_width = 0.02  # x-dimension
magnet_height = 0.01  # z-dimension
rect = Rectangle((-magnet_width/2, -magnet_height/2), magnet_width, magnet_height,
                linewidth=2, edgecolor='black', facecolor='lightgray',
                linestyle='--', alpha=0.5, label='Magnet')
ax.add_patch(rect)

# Add labels
ax.set_xlabel('x-position (m)', fontsize=12)
ax.set_ylabel('z-position (m)', fontsize=12)
ax.set_title('B-Field in XZ-Plane (Side View)\nCuboid: 2cm×2cm×1cm, Polarization: +z (1T)', fontsize=14)
ax.set_aspect('equal')
ax.grid(True, alpha=0.3, linestyle=':')
ax.legend(loc='upper right')

# Add statistics
stats_text = f"Min |B|: {B_mag.min():.2e} T\n"
stats_text += f"Max |B|: {B_mag.max():.2e} T"
ax.text(0.02, 0.98, stats_text, transform=ax.transAxes,
        fontsize=10, verticalalignment='top',
        bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

plt.tight_layout()
plt.savefig('bfield_xz_visualization.png', dpi=150, bbox_inches='tight')
print("\nSaved bfield_xz_visualization.png")
plt.show()
