import Mesh
import FreeCAD as App
import Part
import os
import sys

freecad_file = ""

def get_global_placement(obj): #this method I created with chatgpt
    """
    Calculate the global placement of a given object by considering its parent hierarchy.

    Parameters:
    obj (FreeCAD object): The object for which the global placement is calculated.

    Returns:
    FreeCAD.Placement: The global placement of the object.
    """
    # Initialize global placement with the object's own placement
    global_placement = obj.Placement

    # Traverse up the hierarchy
    while obj.InList:
        # Assume single parent scenario
        parent = obj.InList[0]

        # Multiply the current global placement by the parent's placement
        global_placement = parent.Placement.multiply(global_placement)

        # Move to the next parent
        obj = parent

    return global_placement

def export_component_as_mesh(component, filename, tessellation_level = 0.1): #filename must contain .stl extension
    # Ensure the filename ends with .stl
    if not filename.lower().endswith(".stl"):
        raise ValueError("Filename must have a .stl extension")
    
    # Get the shape from the component
    shape = component.Shape

    # Apply the component's placement to transform it to the global coordinate system (THIS PART IS VERY IMPORTANT)
    #placement = component.Placement
    placement = get_global_placement(component)
    transformed_shape = shape.copy()
    transformed_shape.Placement = placement
    
    # Check if the shape is valid
    if shape.isNull():
        raise ValueError("Shape is null or invalid")
    
    # Recursion Case
    
    # Base Case
    # Tessellate the shape
    try:
        mesh_data = transformed_shape.tessellate(tessellation_level)
    except Exception as e:
        raise RuntimeError(f"Failed to tessellate shape: {e}")

    # Check if tessellation produced a valid mesh
    if not mesh_data or len(mesh_data) == 0:
        raise RuntimeError("Tessellation resulted in an empty mesh")

    # Create a Mesh::Feature object and set its mesh
    mesh_obj = App.ActiveDocument.addObject("Mesh::Feature", "Mesh")
    mesh_obj.Mesh = Mesh.Mesh(mesh_data)
    
    # Save the mesh as an STL file
    global freecad_file
    dir_path = os.path.splitext(freecad_file)[0] #remove the freecad extension from the freecad filename
    dir_path += "/"
    # Check if the directory exists
    if not os.path.exists(dir_path):
        try:
            # Create the directory
            os.makedirs(dir_path, exist_ok=True)
            print(f"Directory '{dir_path}' created successfully.")
        except OSError as e:
            print(f"Error creating directory: {e}")
    try:
        mesh_obj.Mesh.write(dir_path + filename)
        print(f"Mesh saved as {filename}")
    except Exception as e:
        raise RuntimeError(f"Failed to write STL file: {e}")


def process_component(doc, component, parent_placement=App.Placement(), default_name = "unnamed_mesh"):
    if component.Label:
        name = component.Label
    else:
        name = default_name
    #note, given name is the name passed above from recursion
    print(f"Processing Component: {component.Name} (Label: {component.Label})")
    #If the component is a basic component (e.g., a "Part::Feature" or "PartDesign::Body"), transform it
    if component.TypeId in ["Part::Feature"]:
        export_component_as_mesh(component, name  + ".stl")

# Main function to process all components in the document
def transform_all_components_to_global(doc):
    components = [obj for obj in doc.Objects if obj.TypeId in ["Part::Feature"]]
    print(f"{len(components)} bodies found.")
    for component in components:
        process_component(doc, component)
        
# Run the script

# Load the FreeCAD file, the directory is the first argument of the script
if len(sys.argv) < 2:
        print("Usage: python refine_step_for_ePIC.py <arg1> where <arg1> is the name of the FreeCAD(.FCStd) file.")
        sys.exit(1)
freecad_file = str(sys.argv[2])
print(freecad_file)
doc = App.openDocument(freecad_file)
print("dbg")
if doc is None:
    print("File: " + freecad_file + " not found!!")
transform_all_components_to_global(doc)
# Optionally, close the document
App.closeDocument(doc.Name)