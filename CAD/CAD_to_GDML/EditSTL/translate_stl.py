#!/usr/bin/env python3

"""
translate_stl.py

Translate all points in an STL mesh by a specified displacement
(dx, dy, dz) and write the result to a new STL file.

Usage:
    python translate_stl.py input.stl output.stl --dx 10 --dy -5 --dz 2.5

The displacement values use the same coordinate units as the STL file.
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import trimesh


def load_mesh(path):
    """Load an STL file as a Trimesh."""

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


def translate_mesh(mesh, dx, dy, dz):
    """
    Translate every vertex by (dx, dy, dz).

    Returns a copy of the original mesh.
    """

    translated = mesh.copy()

    displacement = np.array(
        [dx, dy, dz],
        dtype=float
    )

    translated.vertices += displacement

    return translated


def main():

    parser = argparse.ArgumentParser(
        description=(
            "Translate an STL mesh by a specified "
            "X, Y and Z displacement."
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
        "--dx",
        type=float,
        default=0.0,
        help="Displacement in X (default: 0)"
    )

    parser.add_argument(
        "--dy",
        type=float,
        default=0.0,
        help="Displacement in Y (default: 0)"
    )

    parser.add_argument(
        "--dz",
        type=float,
        default=0.0,
        help="Displacement in Z (default: 0)"
    )

    args = parser.parse_args()

    # Load mesh
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
    print(
        "  centroid:   ({:.6g}, {:.6g}, {:.6g})".format(
            mesh.centroid[0],
            mesh.centroid[1],
            mesh.centroid[2]
        )
    )

    print()
    print("Displacement:")
    print("  dx: {}".format(args.dx))
    print("  dy: {}".format(args.dy))
    print("  dz: {}".format(args.dz))

    # Translate mesh
    translated = translate_mesh(
        mesh,
        args.dx,
        args.dy,
        args.dz
    )

    print()
    print(
        "New centroid: ({:.6g}, {:.6g}, {:.6g})".format(
            translated.centroid[0],
            translated.centroid[1],
            translated.centroid[2]
        )
    )

    # Write output
    output_path = Path(args.output)

    try:
        translated.export(str(output_path))

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
