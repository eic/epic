import os
import glob
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from lxml import etree
# written by chatgpt - doesn't work
def parse_gdml_vertices(gdml_file):
    vertices = []
    tree = etree.parse(gdml_file)
    root = tree.getroot()

 # Find all vertices (assumes vertices are stored in <position> elements)
    for position in root.findall(".//position"):
        x = position.attrib.get('x')
        y = position.attrib.get('y')
        z = position.attrib.get('z')
#        print(str(x)+", "+str(y))
        vertices.append((x, y, z))

    return vertices

def read_all_gdml_vertices(directory):
    all_vertices = []
    gdml_files = glob.glob(os.path.join(directory, "*.gdml"))
    for gdml_file in gdml_files:
        print(f"Parsing: {gdml_file}")
        vertices = parse_gdml_vertices(gdml_file)
        all_vertices.extend(vertices)
    return all_vertices


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Plot X-Y vertex positions from GDML files.")
    parser.add_argument("directory", help="Directory containing GDML files")
    args = parser.parse_args()

    all_vertices = read_all_gdml_vertices(args.directory)
    print("Write to "+args.directory+".txt")
    f  = open(args.directory+".txt","w")
    for line in all_vertices:
    	f.write(str(line[0]+"\t"+str(line[1])+"\t"+str(line[2]))+"\n")
    f.close()    
    

