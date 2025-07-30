import trimesh
import sys

file = sys.argv[1]
print("Checking file "+file)
mesh = trimesh.load_mesh(file)
components = mesh.split(only_watertight=False)

if not mesh.is_watertight or len(components) > 1:
    print("  Not watertight - cleaning file")
    mesh.merge_vertices()
    mesh.fill_holes()
    mesh.export(file[:-4]+"_CLEANED.stl")
    
