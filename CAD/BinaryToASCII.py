# -*- encoding: utf-8 -*-
import struct
import sys

infile = open(sys.argv[1], 'rb') #import file

file = sys.argv[1]
file_with_no_extension = file[:-4]


out = open(file_with_no_extension + '_ASCII.stl', 'w') #export file

data = infile.read()


out.write("solid ")
'''
for x in range(0,80):
    if not ord(str(data[x])) == 0:
        out.write(struct.unpack('c', data[x])[0])
    else:
        pass'''
out.write("\n")

number = data[80] + data[81] + data[82] + data[83] 
faces = struct.unpack('I',number)[0]

for x in range(0,faces):
    out.write("facet normal ")

    xc = data[84+x*50] + data[85+x*50] + data[86+x*50] + data[87+x*50]
    yc = data[88+x*50] + data[89+x*50] + data[90+x*50] + data[91+x*50]
    zc = data[92+x*50] + data[93+x*50] + data[94+x*50] + data[95+x*50]

    out.write(str(struct.unpack('f',xc)[0]) + " ")
    out.write(str(struct.unpack('f',yc)[0]) + " ")
    out.write(str(struct.unpack('f',zc)[0]) + "\n")

    out.write("outer loop\n")

    for y in range(1,4):
        out.write("vertex ")

        xc = data[84+y*12+x*50] + data[85+y*12+x*50] + data[86+y*12+x*50] + data[87+y*12+x*50]
        yc = data[88+y*12+x*50] + data[89+y*12+x*50] + data[90+y*12+x*50] + data[91+y*12+x*50]
        zc = data[92+y*12+x*50] + data[93+y*12+x*50] + data[94+y*12+x*50] + data[95+y*12+x*50]

        out.write(str(struct.unpack('f',xc)[0]) + " ")
        out.write(str(struct.unpack('f',yc)[0]) + " ")
        out.write(str(struct.unpack('f',zc)[0]) + "\n")

    out.write("endloop\n")
    out.write("endfacet\n")

out.close()
print("end")
