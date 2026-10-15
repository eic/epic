#!/usr/bin/env python3

"""
scale_stl.py

Scale an STL mesh about its centre and write the result to a new STL.

Usage:
    python scale_stl.py input.stl output.stl --scale 1.01

Examples:
    --scale 1.01    Increase dimensions by 1%
    --scale 0.99    Decrease dimensions by 1%
    --scale 2.0     Double the dimensions

The centre of the mesh remains unchanged.
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import trimesh


def load_mesh(path):
    mesh = trimesh.load_mesh(
        path,
        force="mesh",
        process=True
    )

    if not isinstance(mesh, trimesh.Trimesh):
        raise ValueError(
            "{} did not load as a single mesh.".format(path)
        )

    return mesh


def scale_mesh_about_centre(mesh, scale):
    """
    Scale mesh about its centroid.

    x_new = centre + scale * (x - centre)
    """

    if scale <= 0:
        raise ValueError("Scale factor must be greater than zero.")

    scaled = mesh.copy()

    centre = np.asarray(mesh.centroid, dtype=float)

    scaled.vertices = (
        centre
        + scale * (mesh.vertices - centre)
    )

    return scaled, centre


def main():

    parser = argparse.ArgumentParser(
        description=(
            "Scale an STL mesh about its centre while "
            "keeping the centre coordinate unchanged."
        )
    )

    parser.add_argument(
        "input",
        help="Input STL file"
    )

    parser.add_argument(
        "output",
        help="Output STL file"
    )

    parser.add_argument(
        "--scale",
        type=float,
        required=True,
        help=(
            "Scale factor. For example, 1.01 increases the "
            "dimensions by 1 percent."
        )
    )

    args = parser.parse_args()

    # Load input mesh
    try:
        mesh = load_mesh(args.input)
    except Exception as exc:
        print(
            "Error loading mesh: {}".format(exc),
            file=sys.stderr
        )
        sys.exit(1)

    print("Input mesh:")
    print("  vertices:   {}".format(len(mesh.vertices)))
    print("  faces:      {}".format(len(mesh.faces)))
    print("  watertight: {}".format(mesh.is_watertight))
    print("  volume:     {}".format(abs(mesh.volume)))

    # Scale
    try:
        scaled, centre = scale_mesh_about_centre(
            mesh,
            args.scale
        )
    except Exception as exc:
        print(
            "Error scaling mesh: {}".format(exc),
            file=sys.stderr
        )
        sys.exit(1)

    print()
    print("Scaling:")
    print("  scale factor: {}".format(args.scale))
    print(
        "  centre: ({:.6g}, {:.6g}, {:.6g})".format(
            centre[0],
            centre[1],
            centre[2]
        )
    )

    print()
    print("Output mesh:")
    print("  vertices:   {}".format(len(scaled.vertices)))
    print("  faces:      {}".format(len(scaled.faces)))
    print("  watertight: {}".format(scaled.is_watertight))
    print("  volume:     {}".format(abs(scaled.volume)))

    # Write output
    output_path = Path(args.output)

    try:
        scaled.export(str(output_path))
    except Exception as exc:
        print(
            "Error writing output: {}".format(exc),
            file=sys.stderr
        )
        sys.exit(1)

    print()
    print("Wrote:")
    print("  {}".format(output_path))


if __name__ == "__main__":
    main()
