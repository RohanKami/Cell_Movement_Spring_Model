import gmsh
import sys
import math
import random
import numpy

# Initialize gmsh:
gmsh.initialize()

# Function used to generate the initial mesh for a cell
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

# Function used to connect the points along a cell's mesh
def Connect_Mesh_Points(N, cell_points):
    line_dict = {}
    for i in range(N):
        if (i + 1 != N):
            line_dict[str(i+1)] = gmsh.model.geo.add_line(cell_points[str(i+1)], cell_points[str(i+2)])
        else:
            line_dict[str(i+1)] = gmsh.model.geo.add_line(cell_points[str(i+1)], cell_points[str(1)])
    return line_dict

def create_center_to_inner_lines(N, center_point_id, inner_points):
    lines = {}
    for i in range(1, N + 1):
        lines[str(i)] = gmsh.model.geo.add_line(center_point_id, inner_points[str(i)])
    return lines

def create_inner_to_outer_radial_lines(N, inner_points, outer_points):
    lines = {}
    for i in range(1, N + 1):
        lines[str(i)] = gmsh.model.geo.add_line(inner_points[str(i)], outer_points[str(i)])
    return lines

# Function to generate random offsets ensuring no overlap for initial cell placement 
def generate_random_offset(used_offsets, spacing, max_attempts=100):
    Size_Of_Grid = 50
    for _ in range(max_attempts):
        offset_x = random.uniform(-Size_Of_Grid, Size_Of_Grid)
        offset_y = random.uniform(-Size_Of_Grid, Size_Of_Grid)
        if all(math.sqrt((offset_x - x)**2 + (offset_y - y)**2) >= spacing for x, y in used_offsets):
            used_offsets.append((offset_x, offset_y))
            return offset_x, offset_y
    raise ValueError("Could not find a non-overlapping position after multiple attempts")

#This function will take in a specified cell number and node number alongside the nodes updated x and y positions and then redraw the cell
def update_inner_cell_shape(cell_number, node_number, Cell_Points, Cell_Lines, Center_ID, Cell_Radial_Lines, new_x, new_y):

    for line_id in Cell_Lines[str(cell_number)].values():
        gmsh.model.geo.remove([(1, line_id)])  
    for line_id in Cell_Radial_Lines[str(cell_number)].values():
        gmsh.model.geo.remove([(1, line_id)])  

    old_point_id = Cell_Points[str(cell_number)][str(node_number)]
    gmsh.model.geo.remove([(0, old_point_id)])
    
    # Add a new point with updated coordinates
    new_point_id = gmsh.model.geo.add_point(new_x, new_y, 0, lc)
    # Update the Cell_Points dictionary with the new point ID
    Cell_Points[str(cell_number)][str(node_number)] = new_point_id
    
    Cell_Lines[str(cell_number)] = Connect_Mesh_Points(Node_Amount, Cell_Points[str(cell_number)])
    Cell_Radial_Lines[str(cell_number)] = create_center_to_inner_lines(Node_Amount, Center_ID, Cell_Points[str(cell_number)])
    

# This function will take in a specified cell number and node number alongside the nodes updated x and y positions and then redraw the outer cell
def update_outer_cell_shape(cell_number, node_number, Cell_Points, Cell_Lines, Inner_Cell_Points, Inner_Cell_Radial_Lines, new_x, new_y):
    for line_id in Cell_Lines[str(cell_number)].values():
        gmsh.model.geo.remove([(1, line_id)])  
    for line_id in Inner_Cell_Radial_Lines[str(cell_number)].values():
        gmsh.model.geo.remove([(1, line_id)])  
        
    old_point_id = Cell_Points[str(cell_number)][str(node_number)]
    gmsh.model.geo.remove([(0, old_point_id)])
    
    # Add a new point with updated coordinates
    new_point_id = gmsh.model.geo.add_point(new_x, new_y, 0, lc)
    # Update the Cell_Points dictionary with the new point ID
    Cell_Points[str(cell_number)][str(node_number)] = new_point_id
    

    Cell_Lines[str(cell_number)] = Connect_Mesh_Points(Node_Amount, Cell_Points[str(cell_number)])
    Inner_Cell_Radial_Lines[str(cell_number)] = create_inner_to_outer_radial_lines(Node_Amount, Inner_Cell_Points[str(cell_number)], Cell_Points[str(cell_number)])
    
def solve_Mx(node_amount, membrane_cell_coordinates, nucleus_cell_coordinates, cell_number, rho_1, rho_2):
    sum1 = 0
    sum2 = 0
    for i in range(1, node_amount):
        sum1+= (membrane_cell_coordinates[0][cell_number][i] - membrane_cell_coordinates[0][cell_number][i + 1]) * ((membrane_cell_coordinates[1][cell_number][i] + 2*membrane_cell_coordinates[1][cell_number][i + 1]) * membrane_cell_coordinates[0][cell_number][i+1] + (2*membrane_cell_coordinates[1][cell_number][i] + membrane_cell_coordinates[1][cell_number][i+1]) * membrane_cell_coordinates[0][cell_number][i])
        sum2+= (nucleus_cell_coordinates[0][cell_number][i] - nucleus_cell_coordinates[0][cell_number][i + 1]) * ((nucleus_cell_coordinates[1][cell_number][i] + 2*nucleus_cell_coordinates[1][cell_number][i + 1]) * nucleus_cell_coordinates[0][cell_number][i+1] + (2*nucleus_cell_coordinates[1][cell_number][i] + nucleus_cell_coordinates[1][cell_number][i+1]) * nucleus_cell_coordinates[0][cell_number][i])
    
    sum1+= (membrane_cell_coordinates[0][cell_number][node_amount] - membrane_cell_coordinates[0][cell_number][1]) * ((membrane_cell_coordinates[1][cell_number][node_amount] + 2*membrane_cell_coordinates[1][cell_number][1]) * membrane_cell_coordinates[0][cell_number][1] + (2*membrane_cell_coordinates[1][cell_number][node_amount] + membrane_cell_coordinates[1][cell_number][1]) * membrane_cell_coordinates[0][cell_number][node_amount])
    sum2+= (nucleus_cell_coordinates[0][cell_number][node_amount] - nucleus_cell_coordinates[0][cell_number][1]) * ((nucleus_cell_coordinates[1][cell_number][node_amount] + 2*nucleus_cell_coordinates[1][cell_number][1]) * nucleus_cell_coordinates[0][cell_number][1] + (2*nucleus_cell_coordinates[1][cell_number][node_amount] + nucleus_cell_coordinates[1][cell_number][1]) * nucleus_cell_coordinates[0][cell_number][node_amount])
    Mx = rho_1/6 * sum1 + (rho_2-rho_1)/6 *  sum2
    return Mx

def solve_My(node_amount, membrane_cell_coordinates, nucleus_cell_coordinates, cell_number, rho_1, rho_2):
    sum1 = 0
    sum2 = 0
    for i in range(1, node_amount):
        sum1+= (membrane_cell_coordinates[0][cell_number][i] - membrane_cell_coordinates[0][cell_number][i+1]) * (membrane_cell_coordinates[1][cell_number][i]**2 + membrane_cell_coordinates[1][cell_number][i] * membrane_cell_coordinates[1][cell_number][i+1] + membrane_cell_coordinates[1][cell_number][i+1]**2)
        sum2+= (nucleus_cell_coordinates[0][cell_number][i] - nucleus_cell_coordinates[0][cell_number][i+1]) * (nucleus_cell_coordinates[1][cell_number][i]**2 + nucleus_cell_coordinates[1][cell_number][i] * nucleus_cell_coordinates[1][cell_number][i+1] + nucleus_cell_coordinates[1][cell_number][i+1]**2)

    sum1+= (membrane_cell_coordinates[0][cell_number][node_amount] - membrane_cell_coordinates[0][cell_number][1]) * (membrane_cell_coordinates[1][cell_number][node_amount]**2 + membrane_cell_coordinates[1][cell_number][node_amount] * membrane_cell_coordinates[1][cell_number][1] + membrane_cell_coordinates[1][cell_number][1]**2)
    sum2+= (nucleus_cell_coordinates[0][cell_number][node_amount] - nucleus_cell_coordinates[0][cell_number][1]) * (nucleus_cell_coordinates[1][cell_number][node_amount]**2 + nucleus_cell_coordinates[1][cell_number][node_amount] * nucleus_cell_coordinates[1][cell_number][1] + nucleus_cell_coordinates[1][cell_number][1]**2)
    My = rho_1/6 * sum1 + (rho_2-rho_1)/6 *  sum2
    return My

def solve_A(node_amount, membrane_cell_coordinates, nucleus_cell_coordinates, cell_number, rho_1, rho_2):
    sum1 = 0
    sum2 = 0
    for i in range(1, node_amount):
        sum1+=(membrane_cell_coordinates[1][cell_number][i]+membrane_cell_coordinates[1][cell_number][i+1]) * (membrane_cell_coordinates[0][cell_number][i]-membrane_cell_coordinates[0][cell_number][i+1])
        sum2+=(nucleus_cell_coordinates[1][cell_number][i]+nucleus_cell_coordinates[1][cell_number][i+1]) * (nucleus_cell_coordinates[0][cell_number][i]-nucleus_cell_coordinates[0][cell_number][i+1])
   
    sum1+=(membrane_cell_coordinates[1][cell_number][node_amount]+membrane_cell_coordinates[1][cell_number][1]) * (membrane_cell_coordinates[0][cell_number][node_amount]-membrane_cell_coordinates[0][cell_number][1])    
    sum2+=(nucleus_cell_coordinates[1][cell_number][node_amount]+nucleus_cell_coordinates[1][cell_number][1]) * (nucleus_cell_coordinates[0][cell_number][node_amount]-nucleus_cell_coordinates[0][cell_number][1])
    A = rho_1/2 * sum1 + (rho_2-rho_1)/2 * sum2
    return A

def solve_center_of_mass(Mx,My,A):
    x_c = Mx/A
    y_c = My/A
    return x_c, y_c

def run_gmesh():
    # Create the relevant Gmsh data structures from Gmsh model
    gmsh.model.geo.synchronize()

    # Generate mesh
    gmsh.model.mesh.generate()

    # Write mesh data
    gmsh.write("GFG.msh")

    # Creates graphical user interface
    if 'close' not in sys.argv:
        gmsh.fltk.run()

    # Finalize the Gmsh API
    gmsh.finalize()
    
def test_center_of_mass(node_amount, membrane_cell_coordinates, nucleus_cell_coordinates, cell_number, rho_1, rho_2):
    Mx = solve_Mx(Node_Amount, Cell_Membrane_Coordinates, Cell_Nucleus_Coordinates, 1, rho_1, rho_2)
    My = solve_My(Node_Amount, Cell_Membrane_Coordinates, Cell_Nucleus_Coordinates, 1, rho_1, rho_2)
    A1 = solve_A(Node_Amount, Cell_Membrane_Coordinates, Cell_Nucleus_Coordinates, 1, rho_1, rho_2)
    xc, yc = solve_center_of_mass(Mx,My, A1)
    print(Cell_Nucleus_Coordinates[0][1][0])
    print(Cell_Nucleus_Coordinates[1][1][0])
    print(xc)
    print(yc)

# Circle points:
lc = 1e-2
Node_Amount = 20
Cell_Amount = 2
rho_1 = 2
rho_2 = 4

# Create a dictionary to store points of all cells
Center_Cell_Points = {}
Outer_Cell_Points = {}
Inner_Cell_Points = {}
Outer_Cell_Lines = {}
Inner_Cell_Lines = {}
Center_to_Inner_Lines = {}
Inner_to_Outer_Lines = {}



# Parameters for cell spacing to avoid overlap
cell_spacing = 3  # Minimum distance between cell centers to avoid overlap
used_offsets = []
Cell_Membrane_Coordinates = numpy.zeros((2, Cell_Amount + 1, Node_Amount + 1)) #x dimension is coordinates y dimension is cell # and z dimension is node #
Cell_Nucleus_Coordinates = numpy.zeros((2, Cell_Amount + 1, Node_Amount + 1))

# Generate mesh for each cell with a random offset
for i in range(1, Cell_Amount + 1):
    offset_x, offset_y = generate_random_offset(used_offsets, cell_spacing)
    Center_Cell_Points[str(i)] = gmsh.model.geo.add_point(offset_x, offset_y, 0, lc)
    Outer_Cell_Points[str(i)] = Create_Cell_Mesh(Node_Amount, lc, offset_x, offset_y, Cell_Membrane_Coordinates, i, 1)
    Inner_Cell_Points[str(i)] = Create_Cell_Mesh(Node_Amount, lc, offset_x, offset_y, Cell_Nucleus_Coordinates, i, 0.25)
    Center_to_Inner_Lines[str(i)] = create_center_to_inner_lines(Node_Amount, Center_Cell_Points[str(i)], Inner_Cell_Points[str(i)])
    Inner_to_Outer_Lines[str(i)] = create_inner_to_outer_radial_lines(Node_Amount, Inner_Cell_Points[str(i)], Outer_Cell_Points[str(i)])
    
    
test_center_of_mass(Node_Amount, Cell_Membrane_Coordinates, Cell_Nucleus_Coordinates, 0, rho_1, rho_2)
    

# Connect mesh points for each cell
for i in range(1, Cell_Amount + 1):
    Outer_Cell_Lines[str(i)] = Connect_Mesh_Points(Node_Amount, Outer_Cell_Points[str(i)])
    Inner_Cell_Lines[str(i)] = Connect_Mesh_Points(Node_Amount, Inner_Cell_Points[str(i)])
    
    

gmsh.model.geo.synchronize() 

# Define new coordinates for the point
new_x = 1.0
new_y = 2.0
# update_inner_cell_shape(1,1, Inner_Cell_Points, Inner_Cell_Lines, Center_Cell_Points[str(1)], Center_to_Inner_Lines, new_x, new_y)
#update_outer_cell_shape(1, 1, Outer_Cell_Points, Outer_Cell_Lines, Inner_Cell_Points, Inner_to_Outer_Lines, 1.5, 1.5)



#run_gmesh()
