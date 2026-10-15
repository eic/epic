#!/usr/bin/env python3

import argparse
from pathlib import Path


def write_cuboid_stl(filename, length, width, height, center):
    """
    Write an ASCII STL file describing a cuboid.

    Parameters
    ----------
    filename : str or Path
        Output STL filename.
    length : float
        Size along the X axis.
    width : float
        Size along the Y axis.
    height : float
        Size along the Z axis.
    center : tuple[float, float, float]
        Centre coordinate (x, y, z).
    """

    cx, cy, cz = center

    # Half-dimensions
    hx = length / 2.0
    hy = width / 2.0
    hz = height / 2.0

    # 8 vertices
    vertices = {
        "000": (cx - hx, cy - hy, cz - hz),
        "100": (cx + hx, cy - hy, cz - hz),
        "110": (cx + hx, cy + hy, cz - hz),
        "010": (cx - hx, cy + hy, cz - hz),
        "001": (cx - hx, cy - hy, cz + hz),
        "101": (cx + hx, cy - hy, cz + hz),
        "111": (cx + hx, cy + hy, cz + hz),
        "011": (cx - hx, cy + hy, cz + hz),
    }

    # Each face is made from two triangles.
    # The vertex order is chosen so the normals point outward.
    triangles = [
        # Bottom (-Z)
        (vertices["000"], vertices["110"], vertices["100"]),
        (vertices["000"], vertices["010"], vertices["110"]),

        # Top (+Z)
        (vertices["001"], vertices["101"], vertices["111"]),
        (vertices["001"], vertices["111"], vertices["011"]),

        # Front (-Y)
        (vertices["000"], vertices["100"], vertices["101"]),
        (vertices["000"], vertices["101"], vertices["001"]),

        # Back (+Y)
        (vertices["010"], vertices["011"], vertices["111"]),
        (vertices["010"], vertices["111"], vertices["110"]),

        # Left (-X)
        (vertices["000"], vertices["001"], vertices["011"]),
        (vertices["000"], vertices["011"], vertices["010"]),

        # Right (+X)
        (vertices["100"], vertices["110"], vertices["111"]),
        (vertices["100"], vertices["111"], vertices["101"]),
    ]

    # Corresponding outward-pointing normals
    normals = [
        (0, 0, -1),
        (0, 0, -1),

        (0, 0, 1),
        (0, 0, 1),

        (0, -1, 0),
        (0, -1, 0),

        (0, 1, 0),
        (0, 1, 0),

        (-1, 0, 0),
        (-1, 0, 0),

        (1, 0, 0),
        (1, 0, 0),
    ]

    filename = Path(filename)

    with filename.open("w", encoding="ascii") as f:
        f.write("solid cuboid\n")

        for normal, triangle in zip(normals, triangles):
            nx, ny, nz = normal

            f.write(f"  facet normal {nx} {ny} {nz}\n")
            f.write("    outer loop\n")

            for x, y, z in triangle:
                f.write(f"      vertex {x:.9g} {y:.9g} {z:.9g}\n")

            f.write("    endloop\n")
            f.write("  endfacet\n")

        f.write("endsolid cuboid\n")


def main():
    parser = argparse.ArgumentParser(
        description="Generate an STL cuboid centred on a specified coordinate."
    )

    parser.add_argument("output", help="Output STL filename")
    parser.add_argument("length", type=float, help="Cuboid length along X")
    parser.add_argument("width", type=float, help="Cuboid width along Y")
    parser.add_argument("height", type=float, help="Cuboid height along Z")
    parser.add_argument(
        "--center",
        nargs=3,
        type=float,
        metavar=("X", "Y", "Z"),
        default=(0.0, 0.0, 0.0),
        help="Centre coordinate (default: 0 0 0)",
    )

    args = parser.parse_args()

    if args.length <= 0 or args.width <= 0 or args.height <= 0:
        parser.error("Length, width, and height must all be positive.")

    write_cuboid_stl(
        filename=args.output,
        length=args.length,
        width=args.width,
        height=args.height,
        center=tuple(args.center),
    )

    print(f"Wrote cuboid STL to {args.output}")


if __name__ == "__main__":
    main()
