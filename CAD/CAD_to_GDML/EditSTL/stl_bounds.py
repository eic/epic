#!/usr/bin/env python3

import argparse
import struct
from pathlib import Path


def update_bounds(vertex, mins, maxs):
    """Update coordinate bounds with a single (x, y, z) vertex."""
    for i in range(3):
        mins[i] = min(mins[i], vertex[i])
        maxs[i] = max(maxs[i], vertex[i])


def bounds_ascii(filename):
    """Find bounds of an ASCII STL file."""
    mins = [float("inf")] * 3
    maxs = [float("-inf")] * 3
    vertex_count = 0

    with open(filename, "r", encoding="ascii", errors="strict") as f:
        for line in f:
            parts = line.split()

            if len(parts) == 4 and parts[0].lower() == "vertex":
                vertex = tuple(map(float, parts[1:4]))
                update_bounds(vertex, mins, maxs)
                vertex_count += 1

    if vertex_count == 0:
        raise ValueError("No vertices found in ASCII STL file.")

    return mins, maxs


def bounds_binary(filename):
    """Find bounds of a binary STL file."""
    mins = [float("inf")] * 3
    maxs = [float("-inf")] * 3

    with open(filename, "rb") as f:
        # Binary STL:
        # 80-byte header
        # 4-byte unsigned integer containing triangle count
        header = f.read(80)

        if len(header) != 80:
            raise ValueError("Invalid binary STL header.")

        count_data = f.read(4)

        if len(count_data) != 4:
            raise ValueError("Invalid binary STL triangle count.")

        triangle_count = struct.unpack("<I", count_data)[0]

        for _ in range(triangle_count):
            # Each triangle is 50 bytes:
            # normal:    3 x float32 = 12 bytes
            # vertices:  9 x float32 = 36 bytes
            # attribute: uint16      =  2 bytes
            data = f.read(50)

            if len(data) != 50:
                raise ValueError("Unexpected end of binary STL file.")

            values = struct.unpack("<12fH", data)

            # values[0:3] = normal
            # values[3:12] = three vertices
            for i in range(3):
                vertex = values[3 + i * 3 : 6 + i * 3]
                update_bounds(vertex, mins, maxs)

    return mins, maxs


def is_binary_stl(filename):
    """
    Determine whether an STL file is binary by comparing its size
    with the size implied by the binary STL triangle count.
    """
    filename = Path(filename)
    file_size = filename.stat().st_size

    if file_size < 84:
        return False

    with filename.open("rb") as f:
        f.seek(80)
        count_data = f.read(4)

    triangle_count = struct.unpack("<I", count_data)[0]

    expected_size = 84 + triangle_count * 50

    return file_size == expected_size


def stl_bounds(filename):
    """Return minimum and maximum XYZ coordinates from an STL file."""
    if is_binary_stl(filename):
        return bounds_binary(filename)
    else:
        return bounds_ascii(filename)


def main():
    parser = argparse.ArgumentParser(
        description="Print the minimum and maximum XYZ coordinates of an STL file."
    )

    parser.add_argument("filename", help="STL file to inspect")

    args = parser.parse_args()

    mins, maxs = stl_bounds(args.filename)

    print(f"X: min = {mins[0]:.9g}, max = {maxs[0]:.9g}")
    print(f"Y: min = {mins[1]:.9g}, max = {maxs[1]:.9g}")
    print(f"Z: min = {mins[2]:.9g}, max = {maxs[2]:.9g}")


if __name__ == "__main__":
    main()
