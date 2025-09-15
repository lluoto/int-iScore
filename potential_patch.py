from chimera import openModels
from chimera import runCommand as rc

import argparse
parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('--input','-i', type=str, help='Input file name')
parser.add_argument('--dx','-d', type=str, help='Input file name')

args = parser.parse_args()
filename=args.input
dx_file=args.dx
print(filename,dx_file)
rc("open " + filename)
rc("open " + dx_file)
rc('surface')
# Get the surface model (assuming one exists)
surface_model = [m for m in openModels.list() if m.__class__.__name__ == 'MSMSModel'][0]

# Get the volume (DX file) model
volume_model = [m for m in openModels.list() if hasattr(m, 'interpolated_values')][0]

# Access the surface vertices (positions in 3D)
vertex_coords = []
for surface_piece in surface_model.surfacePieces:
    vertices = surface_piece.geometry[0]  # Surface vertices (positions in 3D)
    for vertex in vertices:
        vertex_coords.append([vertex[0], vertex[1], vertex[2]])

# Convert the list to a 2D numpy array (Nx3)
# Interpolate the potential values from the DX file at these vertex positions
# The correct usage of interpolate_volume_data is to pass the volume model, coordinates, and an optional parameter
interpolated_values = volume_model.interpolated_values(vertex_coords)

# Print the interpolated potential values at each vertex

# Print the interpolated potential values at each vertex
for vertex, potential in zip(vertex_coords, interpolated_values):
        vertex.append(potential)
        print(vertex)
