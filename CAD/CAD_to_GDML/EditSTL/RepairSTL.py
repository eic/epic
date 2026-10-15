import argparse
import os
import trimesh
import trimesh.repair
import sys

parser = argparse.ArgumentParser(description="Check and clean an STL mesh.")
parser.add_argument("input", help="Input mesh file")
parser.add_argument(
    "-o", "--output",
    help="Output filename (default: input filename with _CLEANED suffix)"
)
args = parser.parse_args()

file = args.input

# Default output: input filename with "_CLEANED" before the extension
if args.output:
    output_file = args.output
else:
    base, ext = os.path.splitext(file)
    output_file = base + "_CLEANED" + ext

print("Checking file " + file)
mesh = trimesh.load_mesh(file)
print("Watertight? " + str(mesh.is_watertight))
components = mesh.split(only_watertight=False)

if not mesh.is_watertight or len(components) > 1:
    print("  Not watertight - cleaning file")
    print("Merge vertices")
    mesh.merge_vertices()
    print("Watertight? " + str(mesh.is_watertight))
    print("Fill holes")
    mesh.fill_holes()
    print("Watertight? " + str(mesh.is_watertight))
    
    broken = trimesh.repair.broken_faces(mesh)
    print("Broken faces:", len(broken))
    
    print("Watertight:", mesh.is_watertight)
    print("Winding consistent:", mesh.is_winding_consistent)
    print("Euler number:", mesh.euler_number) 
    print("Connected components:", len(mesh.split(only_watertight=False)))
    
    trimesh.repair.fix_winding(mesh)
    trimesh.repair.fix_normals(mesh) 

    print("Watertight? " + str(mesh.is_watertight))


    print("Writing " + output_file)
    mesh.export(output_file)
    sys.exit(1)
else:
    sys.exit(2)
