#!/usr/bin/env python3

"""
plot_all_stl_xy_sections_overlay.py

Plot XY cross-sections from all STL files in a folder on one
interactive Matplotlib plot.

Each STL file has its own colour. All requested Z cross-sections
are super-imposed.

No intersection or boolean calculations are performed.

Usage:
    python plot_all_stl_xy_sections_overlay.py ./meshes \
        --zmin 0 \
        --zmax 100 \
        --dz 5

Specify an initial X/Y plotting range:

    python plot_all_stl_xy_sections_overlay.py ./meshes \
        --zmin 0 \
        --zmax 100 \
        --dz 5 \
        --xmin 10 \
        --xmax 50 \
        --ymin -20 \
        --ymax 20

Save to an image:

    python plot_all_stl_xy_sections_overlay.py ./meshes \
        --zmin 0 \
        --zmax 100 \
        --dz 5 \
        --output sections.png

Matplotlib provides zoom and pan through its normal GUI toolbar.
No WebGL is required.
"""

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import trimesh


def load_mesh(path):
    """Load an STL file as a single Trimesh."""

    mesh = trimesh.load_mesh(
        str(path),
        force="mesh",
        process=True
    )

    if not isinstance(mesh, trimesh.Trimesh):
        raise ValueError(
            "{} did not load as a single mesh.".format(path)
        )

    return mesh


def get_xy_section(mesh, z):
    """
    Calculate an XY cross-section at absolute Z coordinate z.

    Returned curves retain the original STL X/Y coordinates.
    """

    if z < mesh.bounds[0, 2] or z > mesh.bounds[1, 2]:
        return []

    section = mesh.section(
        plane_origin=[0.0, 0.0, z],
        plane_normal=[0.0, 0.0, 1.0]
    )

    if section is None:
        return []

    return section.discrete


def make_z_values(zmin, zmax, dz):
    """Generate Z values from zmin to zmax."""

    if dz <= 0:
        raise ValueError("--dz must be greater than zero.")

    if zmax < zmin:
        raise ValueError(
            "--zmax must be greater than or equal to --zmin."
        )

    values = []

    z = zmin
    tolerance = abs(dz) * 1.0e-9

    while z <= zmax + tolerance:
        values.append(z)
        z += dz

    return values


def get_plot_ranges(meshes, xmin, xmax, ymin, ymax):
    """
    Determine initial X/Y ranges.

    Explicit user values take precedence. Missing values are
    obtained from the STL bounding boxes.
    """

    if xmin is None:
        xmin = min(mesh.bounds[0, 0] for mesh in meshes)

    if xmax is None:
        xmax = max(mesh.bounds[1, 0] for mesh in meshes)

    if ymin is None:
        ymin = min(mesh.bounds[0, 1] for mesh in meshes)

    if ymax is None:
        ymax = max(mesh.bounds[1, 1] for mesh in meshes)

    if xmax <= xmin:
        raise ValueError("xmax must be greater than xmin.")

    if ymax <= ymin:
        raise ValueError("ymax must be greater than ymin.")

    return xmin, xmax, ymin, ymax


def main():

    parser = argparse.ArgumentParser(
        description=(
            "Plot super-imposed XY cross-sections of all STL files "
            "using Matplotlib."
        )
    )

    parser.add_argument(
        "directory",
        help="Directory containing STL files"
    )

    parser.add_argument(
        "--zmin",
        type=float,
        required=True,
        help="Minimum absolute Z coordinate"
    )

    parser.add_argument(
        "--zmax",
        type=float,
        required=True,
        help="Maximum absolute Z coordinate"
    )

    parser.add_argument(
        "--dz",
        type=float,
        required=True,
        help="Z spacing between cross-sections"
    )

    parser.add_argument(
        "--xmin",
        type=float,
        default=None,
        help="Initial minimum X coordinate"
    )

    parser.add_argument(
        "--xmax",
        type=float,
        default=None,
        help="Initial maximum X coordinate"
    )

    parser.add_argument(
        "--ymin",
        type=float,
        default=None,
        help="Initial minimum Y coordinate"
    )

    parser.add_argument(
        "--ymax",
        type=float,
        default=None,
        help="Initial maximum Y coordinate"
    )

    parser.add_argument(
        "--output",
        default=None,
        help="Optional PNG output filename"
    )

    parser.add_argument(
        "--linewidth",
        type=float,
        default=1.0,
        help="Line width (default: 1.0)"
    )

    parser.add_argument(
        "--alpha",
        type=float,
        default=0.5,
        help="Line transparency (default: 0.5)"
    )

    args = parser.parse_args()

    directory = Path(args.directory)

    if not directory.is_dir():
        print(
            "Error: directory does not exist: {}".format(
                directory
            ),
            file=sys.stderr
        )
        sys.exit(1)

    # ------------------------------------------------------------
    # Find STL files
    # ------------------------------------------------------------

    stl_files = sorted(
        path
        for path in directory.iterdir()
        if path.is_file()
        and path.suffix.lower() == ".stl"
    )

    if not stl_files:
        print(
            "No STL files found in: {}".format(directory),
            file=sys.stderr
        )
        sys.exit(1)

    print(
        "Found {} STL files.".format(
            len(stl_files)
        )
    )

    # ------------------------------------------------------------
    # Load meshes
    # ------------------------------------------------------------

    meshes = []
    names = []

    for path in stl_files:

        print("Loading {}".format(path.name))

        try:
            mesh = load_mesh(path)

        except Exception as exc:
            print(
                "  WARNING: {}".format(exc),
                file=sys.stderr
            )
            continue

        meshes.append(mesh)
        names.append(path.stem)

    if not meshes:
        print(
            "No meshes could be loaded.",
            file=sys.stderr
        )
        sys.exit(1)

    # ------------------------------------------------------------
    # Generate Z values
    # ------------------------------------------------------------

    try:
        z_values = make_z_values(
            args.zmin,
            args.zmax,
            args.dz
        )

        xmin, xmax, ymin, ymax = get_plot_ranges(
            meshes,
            args.xmin,
            args.xmax,
            args.ymin,
            args.ymax
        )

    except ValueError as exc:
        print(
            "Error: {}".format(exc),
            file=sys.stderr
        )
        sys.exit(1)

    print()
    print(
        "Z range: {:.6g} to {:.6g}".format(
            args.zmin,
            args.zmax
        )
    )

    print(
        "Z step: {:.6g}".format(
            args.dz
        )
    )

    print()
    print("Initial X/Y range:")
    print(
        "  X: {:.6g} to {:.6g}".format(
            xmin,
            xmax
        )
    )
    print(
        "  Y: {:.6g} to {:.6g}".format(
            ymin,
            ymax
        )
    )

    # ------------------------------------------------------------
    # Create figure
    # ------------------------------------------------------------

    fig, ax = plt.subplots(
        figsize=(11, 11)
    )

    # tab20 gives up to 20 easily distinguishable colours.
    cmap = plt.get_cmap("tab20")

    total_curves = 0

    # ------------------------------------------------------------
    # Plot all sections
    # ------------------------------------------------------------

    for mesh_index, (mesh, name) in enumerate(
        zip(meshes, names)
    ):

        colour = cmap(
            mesh_index % cmap.N
        )

        first_curve = True

        print()
        print(
            "Processing {}".format(name)
        )

        for z in z_values:

            section = get_xy_section(
                mesh,
                z
            )

            if not section:
                continue

            for curve in section:

                curve = np.asarray(curve)

                if (
                    curve.ndim != 2
                    or curve.shape[1] < 2
                ):
                    continue

                ax.plot(
                    curve[:, 0],
                    curve[:, 1],
                    color=colour,
                    linewidth=args.linewidth,
                    alpha=args.alpha,
                    label=(
                        name
                        if first_curve
                        else None
                    )
                )

                first_curve = False
                total_curves += 1

    # ------------------------------------------------------------
    # Plot formatting
    # ------------------------------------------------------------

    ax.set_xlim(
        xmin,
        xmax
    )

    ax.set_ylim(
        ymin,
        ymax
    )

    ax.set_aspect(
        "equal",
        adjustable="box"
    )

    ax.set_xlabel("X")
    ax.set_ylabel("Y")

    ax.set_title(
        "Super-imposed XY Cross-Sections\n"
        "Z = {:.6g} to {:.6g}, step = {:.6g}".format(
            args.zmin,
            args.zmax,
            args.dz
        )
    )

    ax.grid(
        True,
        alpha=0.25
    )

    ax.legend(
        loc="best",
        fontsize="small"
    )

    plt.tight_layout()

    # ------------------------------------------------------------
    # Save image if requested
    # ------------------------------------------------------------

    if args.output:

        plt.savefig(
            args.output,
            dpi=300,
            bbox_inches="tight"
        )

        print()
        print(
            "Saved image to: {}".format(
                args.output
            )
        )

    print()
    print(
        "Total cross-section curves plotted: {}".format(
            total_curves
        )
    )

    print()
    print(
        "Use the Matplotlib toolbar to zoom and pan."
    )

    plt.show()


if __name__ == "__main__":
    main()
