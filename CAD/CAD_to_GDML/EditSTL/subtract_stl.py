#!/usr/bin/env python3

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


def scale_mesh_from_clearance(mesh, clearance):
    """
    Radially expand a mesh about its centre so that the vertex
    furthest from the centre moves outward by `clearance`.

    The transformation is:

        x_new = centre + scale * (x - centre)

    where:

        scale = 1 + clearance / maximum_radius

    Thus the largest radial displacement is exactly `clearance`.

    Note that the actual displacement of other vertices is smaller
    than `clearance` because it is proportional to their distance
    from the centre.
    """

    if clearance < 0:
        raise ValueError("Clearance must be >= 0.")

    expanded = mesh.copy()

    if clearance == 0:
        return expanded, 1.0, mesh.centroid, 0.0

    # Use the volume centroid for a watertight solid.
    # Fall back to the bounding-box centre if necessary.
    try:
        centre = np.asarray(mesh.centroid, dtype=float)
    except Exception:
        centre = None

    if centre is None or not np.all(np.isfinite(centre)):
        centre = mesh.bounds.mean(axis=0)

    # Distance of every vertex from the centre.
    distances = np.linalg.norm(
        mesh.vertices - centre,
        axis=1
    )

    maximum_radius = float(np.max(distances))

    if maximum_radius <= 0:
        raise ValueError(
            "Mesh has zero extent from its centre."
        )

    # Calculate radial scale required for the outermost vertex
    # to move by exactly `clearance`.
    scale = 1.0 + clearance / maximum_radius

    expanded.vertices = (
        centre
        + scale * (mesh.vertices - centre)
    )

    return expanded, scale, centre, maximum_radius


def combine_result(result):
    if result is None:
        return None

    if isinstance(result, trimesh.Trimesh):
        return result

    if isinstance(result, trimesh.Scene):
        meshes = [
            geom
            for geom in result.geometry.values()
            if isinstance(geom, trimesh.Trimesh)
        ]

        if not meshes:
            return None

        return trimesh.util.concatenate(meshes)

    if isinstance(result, (list, tuple)):
        meshes = [
            geom
            for geom in result
            if isinstance(geom, trimesh.Trimesh)
        ]

        if not meshes:
            return None

        return trimesh.util.concatenate(meshes)

    raise TypeError(
        "Unsupported boolean result type: {}".format(type(result))
    )


def subtract_mesh(mesh_a, mesh_b, engine=None):
    result = trimesh.boolean.difference(
        [mesh_a, mesh_b],
        engine=engine
    )

    return combine_result(result)


def main():

    parser = argparse.ArgumentParser(
        description=(
            "Subtract mesh B from mesh A, optionally expanding B "
            "radially by a specified absolute clearance."
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
        help="Output STL mesh"
    )

    parser.add_argument(
        "--clearance",
        type=float,
        default=0.0,
        help=(
            "Radial clearance to add to the furthest point of mesh B "
            "before subtraction. Uses STL coordinate units. "
            "Default: 0"
        )
    )

    parser.add_argument(
        "--engine",
        default=None,
        help=(
            "Boolean engine: manifold, blender, or scad. "
            "Default: auto"
        )
    )

    args = parser.parse_args()

    # Load meshes.
    mesh_a = load_mesh(args.input_a)
    mesh_b = load_mesh(args.input_b)

    print("Mesh A:")
    print("  vertices:   {}".format(len(mesh_a.vertices)))
    print("  faces:      {}".format(len(mesh_a.faces)))
    print("  watertight: {}".format(mesh_a.is_watertight))
    print("  volume:     {}".format(abs(mesh_a.volume)))

    print()
    print("Mesh B:")
    print("  vertices:   {}".format(len(mesh_b.vertices)))
    print("  faces:      {}".format(len(mesh_b.faces)))
    print("  watertight: {}".format(mesh_b.is_watertight))
    print("  volume:     {}".format(abs(mesh_b.volume)))

    # Expand B if requested.
    if args.clearance > 0:

        mesh_b, scale, centre, maximum_radius = (
            scale_mesh_from_clearance(
                mesh_b,
                args.clearance
            )
        )

        print()
        print("Expanding mesh B:")
        print("  centre: ({:.6g}, {:.6g}, {:.6g})".format(
            centre[0],
            centre[1],
            centre[2]
        ))
        print("  maximum radius: {:.6g}".format(maximum_radius))
        print("  requested clearance: {:.6g}".format(args.clearance))
        print("  calculated scale: {:.12g}".format(scale))
        print("  new volume: {}".format(abs(mesh_b.volume)))

    # Perform boolean subtraction.
    try:
        result = subtract_mesh(
            mesh_a,
            mesh_b,
            engine=args.engine
        )

    except Exception as exc:
        print(
            "Boolean subtraction failed: {}".format(exc),
            file=sys.stderr
        )
        sys.exit(2)

    if result is None or len(result.faces) == 0:
        print("Subtraction produced an empty mesh.")
        sys.exit(1)

    # Cleanup result.
    result.remove_duplicate_faces()
    result.remove_degenerate_faces()
    result.remove_unreferenced_vertices()

    output_path = Path(args.output)
    result.export(str(output_path))

    print()
    print("Result:")
    print("  output:     {}".format(output_path))
    print("  vertices:   {}".format(len(result.vertices)))
    print("  faces:      {}".format(len(result.faces)))
    print("  watertight: {}".format(result.is_watertight))
    print("  volume:     {}".format(abs(result.volume)))


if __name__ == "__main__":
    main()
