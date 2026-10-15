#!/usr/bin/env python3
"""
intersect_stl.py

Compute the portion of the first STL mesh that lies inside
the second STL mesh, and write it to a new STL file.

Usage:
    python intersect_stl.py first.stl second.stl intersection.stl

Optional:
    python intersect_stl.py first.stl second.stl intersection.stl --engine manifold
"""

import argparse
import sys
from pathlib import Path

import trimesh


def load_mesh(path):
    mesh = trimesh.load_mesh(path, force="mesh", process=True)

    if not isinstance(mesh, trimesh.Trimesh):
        raise ValueError("{} did not load as a single mesh.".format(path))

    return mesh


def combine_result(result):
    """
    Boolean operations can return a Trimesh, Scene, or collection
    depending on the trimesh version/backend. Convert everything
    to one Trimesh.
    """
    if result is None:
        return None

    if isinstance(result, trimesh.Trimesh):
        return result

    if isinstance(result, trimesh.Scene):
        meshes = [
            geom for geom in result.geometry.values()
            if isinstance(geom, trimesh.Trimesh)
        ]
        if not meshes:
            return None
        return trimesh.util.concatenate(meshes)

    if isinstance(result, (list, tuple)):
        meshes = [
            geom for geom in result
            if isinstance(geom, trimesh.Trimesh)
        ]
        if not meshes:
            return None
        return trimesh.util.concatenate(meshes)

    raise TypeError(
        "Unsupported boolean result type: {}".format(type(result))
    )


def meshes_intersection(mesh_a, mesh_b, engine=None):
    """
    Calculate the geometric intersection of mesh_a and mesh_b.
    """
    result = trimesh.boolean.intersection(
        [mesh_a, mesh_b],
        engine=engine
    )

    return combine_result(result)


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Write the portion of the first STL mesh that "
            "intersects the second STL mesh."
        )
    )

    parser.add_argument(
        "input_a",
        help="First STL mesh"
    )

    parser.add_argument(
        "input_b",
        help="Second STL mesh"
    )

    parser.add_argument(
        "output",
        help="Output STL containing the intersection"
    )

    parser.add_argument(
        "--engine",
        default=None,
        help="Boolean engine: manifold, blender, or scad. Default: auto."
    )

    args = parser.parse_args()

    mesh_a = load_mesh(args.input_a)
    mesh_b = load_mesh(args.input_b)

    if not mesh_a.is_watertight:
        print(
            "Warning: first mesh is not watertight.",
            file=sys.stderr
        )

    if not mesh_b.is_watertight:
        print(
            "Warning: second mesh is not watertight.",
            file=sys.stderr
        )

    try:
        result = meshes_intersection(
            mesh_a,
            mesh_b,
            engine=args.engine
        )
    except Exception as exc:
        print(
            "Boolean intersection failed: {}".format(exc),
            file=sys.stderr
        )
        print(
            "Make sure both meshes are closed/watertight and that "
            "a boolean backend such as manifold3d is installed.",
            file=sys.stderr
        )
        sys.exit(2)

    # No intersection
    if result is None or len(result.faces) == 0:
        print("The meshes do not have a volumetric intersection.")
        print("No output mesh was written.")
        sys.exit(1)

    # Remove any empty/invalid geometry
    result.remove_duplicate_faces()
    result.remove_degenerate_faces()
    result.remove_unreferenced_vertices()
    result.merge_vertices()

    output_path = Path(args.output)
    result.export(str(output_path))

    print("Intersection mesh written to:")
    print("  {}".format(output_path))
    print("Vertices: {}".format(len(result.vertices)))
    print("Faces:    {}".format(len(result.faces)))
    print("Volume:   {}".format(result.volume))


if __name__ == "__main__":
    main()
