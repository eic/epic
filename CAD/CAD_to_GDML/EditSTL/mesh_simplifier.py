import os
import shutil
import pymeshlab

class MeshSimplifier:

    def __init__(self, keyword: str, target_faces: int = None, target_percentage: float = None,
                 quality_threshold: float = 0.5, preserve_normal: bool = True, preserve_topology: bool = True,
                 preserve_boundary: bool = True):
        self.keyword = keyword
        self.target_faces = target_faces
        self.target_percentage = target_percentage
        self.quality_threshold = quality_threshold
        self.preserve_normal = preserve_normal
        self.preserve_topology = preserve_topology
        self.preserve_boundary = preserve_boundary

    def matches(self, filename: str) -> bool:
        return filename.endswith('.stl') and self.keyword in filename

    def simplify(self, input_file: str, output_stl_folder: str):

        filename = os.path.basename(input_file)
        name_root, _ = os.path.splitext(filename)
        output_stl_path = os.path.join(output_stl_folder, name_root + ".stl")
        
        if not self.matches(filename):
            print("Error: Probably not STL")
            return

        print(f"  Simplifying: {filename}")

        ms = pymeshlab.MeshSet()
        ms.load_new_mesh(input_file)

        original_faces = ms.current_mesh().face_number()

        if self.target_faces:
            target_faces = self.target_faces
        elif self.target_percentage:
            target_faces = int(original_faces * self.target_percentage)
        else:
            print(f"!!  No simplification target set for: {filename}")
            return

        # Decimate the mesh
        ms.meshing_decimation_quadric_edge_collapse(
            targetfacenum=target_faces,
            preservenormal=self.preserve_normal,
            preservetopology=self.preserve_topology,
            preserveboundary=self.preserve_boundary,
            optimalplacement=True,
            qualitythr=self.quality_threshold,
            planarquadric=True,
        )

        # Save the simplified STL
        ms.save_current_mesh(output_stl_path)
        print(f" Saved STL to: {output_stl_path}")
