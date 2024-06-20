import gmsh
import sys
import math
import random
import numpy as np

# Initialize gmsh:
gmsh.initialize()

# Function to generate the initial mesh for a cell
def Create_Cell_Mesh(N, lc, offset_x, offset_y, Cell_Coordinates, cell_number, radius):
    cell_points = {}
    Theta = 0
    dTheta = 2 * math.pi / N
    Cell_Coordinates[0, cell_number-1, 0] = offset_x
    Cell_Coordinates[1, cell_number-1, 0] = offset_y
    for i in range(N):
        x = radius * math.cos(Theta) + offset_x
        y = radius * math.sin(Theta) + offset_y
        cell_points[str(i+1)] = gmsh.model.geo.add_point(x, y, 0, lc)
        Cell_Coordinates[0, cell_number-1, i+1] = x
        Cell_Coordinates[1, cell_number-1, i+1] = y
        Theta += dTheta
    return cell_points

# Function to connect the points along a cell's mesh
def Connect_Mesh_Points(N, cell_points):
    line_dict = {}
    for i in range(N):
        if (i + 1 != N):
            line_dict[str(i+1)] = gmsh.model.geo.add_line(cell_points[str(i+1)], cell_points[str(i+2)])
        else:
            line_dict[str(i+1)] = gmsh.model.geo.add_line(cell_points[str(i+1)], cell_points[str(1)])
    return line_dict

# Function to generate random offsets ensuring no overlap for initial cell placement
def generate_random_offset(used_offsets, spacing, max_attempts=100):
    Size_Of_Grid = 25
    for _ in range(max_attempts):
        offset_x = random.uniform(-Size_Of_Grid, Size_Of_Grid)
        offset_y = random.uniform(-Size_Of_Grid, Size_Of_Grid)
        if all(math.sqrt((offset_x - x)**2 + (offset_y - y)**2) >= spacing for x, y in used_offsets):
            used_offsets.append((offset_x, offset_y))
            return offset_x, offset_y
    raise ValueError("Could not find a non-overlapping position after multiple attempts")

def run_gmesh():
    # Create the relevant Gmsh data structures from Gmsh model
    gmsh.model.geo.synchronize()

    # Generate mesh
    gmsh.model.mesh.generate(2)

    # Write mesh data
    gmsh.write("GFG.msh")

    # Creates graphical user interface
    if 'close' not in sys.argv:
        gmsh.fltk.run()

    # Finalize the Gmsh API
    gmsh.finalize()

# Circle points:
lc = 1e-2
Node_Amount = 30
Cell_Amount = 5

# Create a dictionary to store points of all cells
Outer_Cell_Points = {}
Inner_Cell_Points = {}
Outer_Cell_Lines = {}
Inner_Cell_Lines = {}
Outer_Cell_Loops = {}
Inner_Cell_Loops = {}

# Parameters for cell spacing to avoid overlap
cell_spacing = 3  # Minimum distance between cell centers to avoid overlap
used_offsets = []
Cell_Membrane_Coordinates = np.zeros((2, Cell_Amount + 1, Node_Amount + 1)) #x dimension is coordinates y dimension is cell # and z dimension is node #
Cell_Nucleus_Coordinates = np.zeros((2, Cell_Amount + 1, Node_Amount + 1))

# Generate mesh for each cell with a random offset
for i in range(1, Cell_Amount + 1):
    offset_x, offset_y = generate_random_offset(used_offsets, cell_spacing)
    Outer_Cell_Points[str(i)] = Create_Cell_Mesh(Node_Amount, lc, offset_x, offset_y, Cell_Membrane_Coordinates, i, 1)
    Inner_Cell_Points[str(i)] = Create_Cell_Mesh(Node_Amount, lc, offset_x, offset_y, Cell_Nucleus_Coordinates, i, 0.25)

# Connect mesh points for each cell
for i in range(1, Cell_Amount + 1):
    Outer_Cell_Lines[str(i)] = Connect_Mesh_Points(Node_Amount, Outer_Cell_Points[str(i)])
    Inner_Cell_Lines[str(i)] = Connect_Mesh_Points(Node_Amount, Inner_Cell_Points[str(i)])
    
    # Create a curve loop for the outer and inner circles
    Outer_Cell_Loops[str(i)] = gmsh.model.geo.add_curve_loop(list(Outer_Cell_Lines[str(i)].values()))
    Inner_Cell_Loops[str(i)] = gmsh.model.geo.add_curve_loop(list(Inner_Cell_Lines[str(i)].values()))
    
    # Define the annular region between the inner and outer circles
    annular_surface = gmsh.model.geo.add_plane_surface([Outer_Cell_Loops[str(i)], Inner_Cell_Loops[str(i)]])
    
    # Define the inner circle surface
    inner_surface = gmsh.model.geo.add_plane_surface([Inner_Cell_Loops[str(i)]])
    
run_gmesh()
