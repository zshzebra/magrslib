import matplotlib.pyplot as plt
import numpy as np
import json

# Load data
with open('bfield_data.json', 'r') as f:
    data = json.load(f)

grid = data['grid']
positions = np.array(data['positions'])
b_field = np.array(data['b_field'])

nx, ny = grid['nx'], grid['ny']

# Reshape
X = positions[:, 0].reshape(ny, nx)
Y = positions[:, 1].reshape(ny, nx)
Bx = b_field[:, 0].reshape(ny, nx)
By = b_field[:, 1].reshape(ny, nx)
Bz = b_field[:, 2].reshape(ny, nx)

B_mag = np.sqrt(Bx**2 + By**2 + Bz**2)

print(f"X shape: {X.shape}, range: [{X.min()}, {X.max()}]")
print(f"Y shape: {Y.shape}, range: [{Y.min()}, {Y.max()}]")
print(f"Bx shape: {Bx.shape}, range: [{Bx.min()}, {Bx.max()}]")
print(f"By shape: {By.shape}, range: [{By.min()}, {By.max()}]")

# Check if X and Y are meshgrids
print(f"\nX[0,0] = {X[0,0]}, X[0,1] = {X[0,1]}, X[1,0] = {X[1,0]}")
print(f"Y[0,0] = {Y[0,0]}, Y[1,0] = {Y[1,0]}, Y[0,1] = {Y[0,1]}")
print(f"\nDoes X vary along rows (axis 0)? {not np.allclose(X[0,:], X[1,:])}")
print(f"Does X vary along cols (axis 1)? {not np.allclose(X[:,0], X[:,1])}")
print(f"Does Y vary along rows (axis 0)? {not np.allclose(Y[0,:], Y[1,:])}")
print(f"Does Y vary along cols (axis 1)? {not np.allclose(Y[:,0], Y[:,1])}")

#  For streamplot, we need:
# - X: meshgrid of x-coordinates, varying along axis 1 (columns)
# - Y: meshgrid of y-coordinates, varying along axis 0 (rows)

# Create a simple test plot
fig, ax = plt.subplots(figsize=(8, 8))

color = np.log10(B_mag + 1e-10)
stream = ax.streamplot(X, Y, Bx, By, color=color, cmap='plasma', density=1.5)

ax.set_xlabel('x (m)')
ax.set_ylabel('y (m)')
ax.set_title('B-field streamlines')
ax.set_aspect('equal')
ax.grid(True, alpha=0.3)

# Add magnet outline
from matplotlib.patches import Rectangle
magnet_width = 0.02
magnet_height = 0.02
rect = Rectangle((-magnet_width/2, -magnet_height/2), magnet_width, magnet_height,
                linewidth=2, edgecolor='black', facecolor='none', linestyle='--')
ax.add_patch(rect)

plt.colorbar(stream.lines, label='log₁₀(|B|) [T]')
plt.tight_layout()
plt.savefig('test_streamplot.png', dpi=150)
print("\nSaved test_streamplot.png")
print("Check if streamlines look reasonable!")
