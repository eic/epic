import os

import shutil

from mesh_simplifier import MeshSimplifier



# Config paths

input_folder = "simplifyIP"

output_stl_folder = "simplifiedOP"




# os.makedirs(output_stl_folder, exist_ok=True)

# os.makedirs(output_gdml_folder, exist_ok=True)



# Instantiate simplifiers

simplifiers = [

    MeshSimplifier(keyword="Readout", target_faces=12),

    MeshSimplifier(keyword="Biasing", target_faces=12),

    MeshSimplifier(keyword="kapton_sensor", target_faces=12),

    MeshSimplifier(keyword="Pads", target_faces=12),

    MeshSimplifier(keyword="Bridge", target_percentage=0.1),

    MeshSimplifier(keyword="Ultem", target_percentage = 0.2),
    
    MeshSimplifier(keyword="K9", target_percentage = 0.8)
]



default_simplifier = MeshSimplifier(keyword="", target_percentage=0.5)



for filename in os.listdir(input_folder):

    if not filename.endswith('.stl'):

        continue



    filepath = os.path.join(input_folder, filename)

    matched = False



    for simplifier in simplifiers:

        if simplifier.matches(filename):

            simplifier.simplify(filepath, output_stl_folder)

            matched = True

            break  # only apply one simplifier per file



    if not matched:

        file_size_kb = os.path.getsize(filepath) / 1024

        threshold_kb = 30

        if file_size_kb > threshold_kb:

            default_simplifier.simplify(filepath, output_stl_folder)

        # No simplifier matched: copy original file as is

        dest_stl = os.path.join(output_stl_folder, filename)

        shutil.copy2(filepath, dest_stl)

        print(f"Copied without changes: {filename}")
