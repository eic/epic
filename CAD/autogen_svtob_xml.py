# this is a python script

#=======================================================================
#   Copyright (C) 2024 Univ. of Bham  All rights reserved.
#   
#   		FileName：		autogen_svtob_xml.py
#   	 	Author：		LongLI <l.li.8@bham.ac.uk> <long.l@cern.ch>
#   		Time：			2024.11.12
#   		Description：
#
#======================================================================

'''
    This work is an integration of Tuna's work on auto-generation
    of xml description file of SVTOB from FCStd file
'''

import os
import sys
import argparse
import  stl
from stl import mesh
import stl_gdml
import create_xml 
from pathlib import Path



def generate_xml(args):
    
    # convert the FCStd file in to stl file
    args.fcstd = sorted(args.fcstd)
    for fcfile in args.fcstd:
        print(f'Processing {fcfile}')
        
        os.system(f'/snap/bin/freecad.cmd svt_ob_meshing.py {fcfile}') # for Long's laptop setting
        for rootdir, _, files in os.walk(fcfile.replace('.FCStd', '/')):
            for file in files:
                file = os.path.join(rootdir, file)
                if not file.endswith('ASCII.stl') and file.endswith('.stl'):
                    print(f'Converting {file} to ADCII file...')
                    binarystl = mesh.Mesh.from_file(file)
                    outfile = file.replace('.stl', 'ASCII.stl')
                    binarystl.save(outfile, mode=stl.Mode.ASCII)

                    # convert the ASCII stl file to gdml file
                    
                    print(f'Converting {outfile} to gdml file...')
                    stl_gdml.stl_to_gdml(outfile)
                

    if len(args.fcstd) > 1:
        save_path = os.path.join(os.path.dirname(__file__), 'xml')
        if not os.path.exists(save_path): os.makedirs(save_path)

        create_xml.generate_xml(args.fcstd[0].replace('.FCStd', '/'), args.fcstd[1].replace('.FCStd', '/'), os.path.join(save_path, args.output))

                
                
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='This is a parameter description for the auto-generation of SVTOB xml file')
    parser.add_argument('--fcstd', '-f', type=str, nargs='+', default=['L3_stave_Long.FCStd', 'L4_stave_Long.FCStd'], help='input FCStd file from FreeCAD')
    parser.add_argument('--output', '-o', type=str, default='silicon_barrel.xml', help='the output file name eg. silicon_barrel.xml')
    
    args = parser.parse_args()
    
    generate_xml(args)