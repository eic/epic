# -*- encoding: utf-8 -*-
# This script comes from https://github.com/IsseiMori/binary-stl-toASCII/blob/master/BinaryToASCII.py
# Converted to Python 3

import struct
import sys

file = sys.argv[1]
file_with_no_extension = file[:-4]

with open(sys.argv[1], 'rb') as infile, open(file_with_no_extension + '_ASCII.stl', 'w') as out:
    data = infile.read()

    out.write("solid ")

    for x in range(80):
        if data[x] != 0:
            out.write(data[x:x+1].decode('utf-8'))
    out.write("\n")

    number = data[80:84]
    faces = struct.unpack('I', number)[0]

    for x in range(faces):
        out.write("facet normal ")
        base = 84 + x * 50

        nx = struct.unpack('f', data[base:base+4])[0]
        ny = struct.unpack('f', data[base+4:base+8])[0]
        nz = struct.unpack('f', data[base+8:base+12])[0]

        out.write(f"{nx} {ny} {nz}\n")
        out.write("outer loop\n")

        for y in range(3):
            vertex_base = base + 12 + y * 12
            vx = struct.unpack('f', data[vertex_base:vertex_base+4])[0]
            vy = struct.unpack('f', data[vertex_base+4:vertex_base+8])[0]
            vz = struct.unpack('f', data[vertex_base+8:vertex_base+12])[0]
            out.write(f"vertex {vx} {vy} {vz}\n")

        out.write("endloop\n")
        out.write("endfacet\n")

    out.write("endsolid\n")

print("end")
