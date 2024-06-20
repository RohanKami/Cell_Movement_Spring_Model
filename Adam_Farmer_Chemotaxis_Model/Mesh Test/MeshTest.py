
# Import modules:
import gmsh
import sys
 
# Initialize gmsh:
gmsh.initialize()
 
# cube points:
lc = 1e-2
point1 = gmsh.model.geo.add_point(0, 0, 0, lc)
point2 = gmsh.model.geo.add_point(1, 0, 0, lc)
point3 = gmsh.model.geo.add_point(1, 1, 0, lc)
point4 = gmsh.model.geo.add_point(0, 1, 0, lc)
point5 = gmsh.model.geo.add_point(0, 1, 1, lc)
point6 = gmsh.model.geo.add_point(0, 0, 1, lc)
point7 = gmsh.model.geo.add_point(1, 0, 1, lc)
point8 = gmsh.model.geo.add_point(1, 1, 1, lc)
 
# Create the relevant Gmsh data structures 
# from Gmsh model.
gmsh.model.geo.synchronize()
 
# Generate mesh:
gmsh.model.mesh.generate()
 
# Write mesh data:
gmsh.write("GFG.msh")
 
# Creates  graphical user interface
if 'close' not in sys.argv:
    gmsh.fltk.run()
 
# It finalize the Gmsh API
gmsh.finalize()