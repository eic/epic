#!/usr/bin/env python3

"""
plot_stl_intersection.py

Plot two STL meshes and highlight their intersection in red.

Usage:
    python plot_stl_intersection.py A.stl B.stl

Optional:
    --output figure.png
        Save the plot to a PNG file instead of only displaying it.

    --engine manifold
        Boolean engine to use for calculating the intersection.
"""

import argparse
import sys

import matplotlib.pyplot as plt
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


def combine_result(result):
    """Convert a boolean result into a single Trimesh."""
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


def calculate_intersection(mesh_a, mesh_b, engine=None):
    """Calculate A ∩ B."""
    result = trimesh.boolean.intersection(
        [mesh_a, mesh_b],
        engine=engine
    )

    return combine_result(result)


def add_mesh_to_plot(ax, mesh, color, alpha, label):
    """Plot a triangular mesh using Matplotlib."""
    if mesh is None or len(mesh.faces) == 0:
        return

    vertices = mesh.vertices
    faces = mesh.faces

    ax.plot_trisurf(
        vertices[:, 0],
        vertices[:, 1],
        vertices[:, 2],
        triangles=faces,
        color=color,
        alpha=alpha,
        linewidth=0.1,
        edgecolor="none",
        label=label
    )


def set_equal_aspect(ax, meshes):
    """Set equal XYZ scaling based on all supplied meshes."""
    mins = []
    maxs = []

    for mesh in meshes:
        if mesh is not None and len(mesh.vertices) > 0:
            mins.append(mesh.bounds[0])
            maxs.append(mesh.bounds[1])

    if not mins:
        return

    overall_min = np.min(np.array(mins), axis=0)
    overall_max = np.max(np.array(maxs), axis=0)

    centre = (overall_min + overall_max) / 2.0
    half_range = np.max(overall_max - overall_min) / 2.0

    if half_range == 0:
        half_range = 1.0

    ax.set_xlim(
        centre[0] - half_range,
        centre[0] + half_range
    )

    ax.set_ylim(
        centre[1] - half_range,
        centre[1] + half_range
    )

    ax.set_zlim(
        centre[2] - half_range,
        centre[2] + half_range
    )


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Plot two STL volumes and highlight their intersection in red."
        )
    )

    parser.add_argument(
        "mesh_a",
        help="First STL file"
    )

    parser.add_argument(
        "mesh_b",
        help="Second STL file"
    )

    parser.add_argument(
        "--output",
        help="Save the figure to this file, e.g. plot.png"
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

    print("Loading meshes...")

    try:
        mesh_a = load_mesh(args.mesh_a)
        mesh_b = load_mesh(args.mesh_b)
    except Exception as exc:
        print(
            "Error loading meshes: {}".format(exc),
            file=sys.stderr
        )
        sys.exit(1)

    print()
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

    # Calculate intersection.
    print()
    print("Calculating intersection...")

    try:
        intersection = calculate_intersection(
            mesh_a,
            mesh_b,
            engine=args.engine
        )
    except Exception as exc:
        print(
            "Error calculating intersection: {}".format(exc),
            file=sys.stderr
        )
        sys.exit(2)

    if intersection is None or len(intersection.faces) == 0:
        intersection_volume = 0.0
        print("No volumetric intersection found.")
    else:
        intersection_volume = abs(float(intersection.volume))

        print("Intersection:")
        print("  vertices:   {}".format(len(intersection.vertices)))
        print("  faces:      {}".format(len(intersection.faces)))
        print("  volume:     {}".format(intersection_volume))

    # Create plot.
    fig = plt.figure(figsize=(11, 9))
    ax = fig.add_subplot(111, projection="3d")

    # Plot the two original meshes.
    add_mesh_to_plot(
        ax,
        mesh_a,
        color="blue",
        alpha=0.25,
        label="Mesh A"
    )

    add_mesh_to_plot(
        ax,
        mesh_b,
        color="green",
        alpha=0.25,
        label="Mesh B"
    )

    # Plot intersection on top in red.
    if intersection is not None and len(intersection.faces) > 0:
        add_mesh_to_plot(
            ax,
            intersection,
            color="red",
            alpha=0.9,
            label="Intersection"
        )

    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")

    ax.set_title(
        "STL Mesh Intersection\n"
        "A volume = {:.6g}, "
        "B volume = {:.6g}, "
        "Intersection = {:.6g}".format(
            abs(mesh_a.volume),
            abs(mesh_b.volume),
            intersection_volume
        )
    )

    set_equal_aspect(
        ax,
        [mesh_a, mesh_b, intersection]
    )

    # Matplotlib's 3D axes don't always handle legend labels
    # from plot_trisurf automatically, so create explicit handles.
    from matplotlib.patches import Patch

    legend_handles = [
        Patch(
            facecolor="blue",
            alpha=0.25,
            label="Mesh A"
        ),
        Patch(
            facecolor="green",
            alpha=0.25,
            label="Mesh B"
        ),
    ]

    if intersection is not None and len(intersection.faces) > 0:
        legend_handles.append(
            Patch(
                facecolor="red",
                alpha=0.9,
                label="Intersection"
            )
        )

    ax.legend(handles=legend_handles)

    plt.tight_layout()

    if args.output:
        plt.savefig(
            args.output,
            dpi=300,
            bbox_inches="tight"
        )
        print()
        print("Figure saved to: {}".format(args.output))

    plt.show()


if __name__ == "__main__":
    main()
