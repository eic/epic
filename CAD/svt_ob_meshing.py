# this is a python script

#=======================================================================
#   Copyright (C) 2024 Univ. of Bham  All rights reserved.
#   
#   		FileName：		svt_ob_meshing.py
#   	 	Author：		LongLI <long.l@cern.ch>
#   		Time：			2024.12.02
#   		Description：
#
#======================================================================

import Mesh
import FreeCAD as App
import Part
import os
import sys

'''
    This python scripts is for FreeCAD python API(freecad.cmd/freecadcmd)
    usage:
        freecad.cmd sub_ob_meshing.py Freecadfile.FCStd to convert the FCStd file 
        in to Binary stl file.
'''

global freecad_file


def get_global_placement(obj):
    
    global_placement = obj.Placement

    # Traverse up the hierarchy
    while obj.InList: # parent object
        # Assume single parent scenario
        parent = obj.InList[0]

        # Multiply the current global placement by the parent's placement
        global_placement = parent.Placement.multiply(global_placement)

        # Move to the next parent
        obj = parent

    return global_placement


def mesh_object(component, outfile, tessellation_level = 0.1):
    if not outfile.endswith('.stl'): outfile += '.stl'

    shape = component.Shape
    
    #transform the component from the local to the global
    placement = get_global_placement(component)
    trans_shape = shape.copy()
    trans_shape.Placement = placement

    #check the shape is valid
    if shape.isNull():
        raise ValueError("The shape is not valid!")

    #Tessellate the shape
    try:
        mesh_data = trans_shape.tessellate(tessellation_level)
    except Exception as e:
        raise RuntimeError(f"fail to tessellate the shape {shape}")
    
    # check the mesh is valid
    if not mesh_data or len(mesh_data) == 0:
        raise RuntimeError("Tessellation results in empty mesh")
    
    #create Mesh:Feature
    mesh_obj = App.ActiveDocument.addObject("Mesh::Feature", "Mesh")
    mesh_obj.Mesh = Mesh.Mesh(mesh_data)

    # save the mesh in the binary stl file
    save_dir = freecad_file.split('.')[0]
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)

    try:
        mesh_obj.Mesh.write(os.path.join(save_dir, outfile))
        print(f'Mesh saved in {outfile}!')

    except Exception as e:
        raise RuntimeError(f"Fail to write STL file {outfile}")


def convert_all_components(doc):
    components = [obj for obj in doc.Objects if obj.TypeId in ["Part::Feature"]]
    print(f"{len(components)} objects will be processed!")

    default_name = 'Unknown_obj'
    for component in components:
        if component.Label:
            name = component.Label
        else:
            name = default_name

        name = name.replace(' ', '_')
        name = name.replace('(', '')
        name = name.replace(')', '')
        
        print("Meshing", component.Name, component.Label, '...')

        mesh_object(component, name)



if len(sys.argv) < 2:
        print("Usage: python refine_step_for_ePIC.py <arg1> where <arg1> is the name of the FreeCAD(.FCStd) file.")
        sys.exit(1)
freecad_file = str(sys.argv[2])
print(freecad_file)
doc = App.openDocument(freecad_file)

if doc is None:
    print("File: " + freecad_file + " not found!!")
convert_all_components(doc)
# Optionally, close the document
App.closeDocument(doc.Name)

