import pyvista as pv
import numpy as np
import sys

# Load STL file
mesh = pv.read(sys.argv[1])

# Create a plotter
plotter = pv.Plotter()
plotter.add_mesh(mesh, color="lightgray", opacity=0.5, show_edges=True, label="Original Mesh")

# Highlight holes (open edges)
edges = mesh.extract_feature_edges(
    boundary_edges=True,
    non_manifold_edges=False,
    feature_edges=False,
    manifold_edges=False
)

if edges.n_points > 0:
    plotter.add_mesh(edges, color="red", line_width=3, label="Open Edges")
    print(f"Found {edges.n_cells} open edge segments.")
else:
    print("No open edges found — mesh might be watertight.")

# Show the plot with a legend
plotter.add_legend()
plotter.show()

