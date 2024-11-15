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
import struct
import stl_gdml 
import create_xml 

def convert_ASCIISTL(stlfile):
    
    if 'ASCII' in stlfile: pass
    
    else:
        print(f'convert {stlfile} in to ASCII!') 
        # can not implement in python3 (byte buffer issues, need to solve)
        os.system(f'/usr/bin/python2 BinaryToASCII.py {stlfile}')   

        # convert to GDML file
        stlfile = stlfile.replace('.stl', '_ASCII.stl')
        stl_gdml.creat_gdml_bundle('./out.gdml', [stlfile])


def generate_xml(args):
    
    stave_paths = []
    # convert the FCStd file in to stl file
    for fcfile in args.fcstd:
        print(f'Processing {fcfile}')
        
        os.system(f'/snap/bin/freecad.cmd refine_step_for_ePIC.py {fcfile}') # for Long's laptop setting
        #rename the stl file for linux OS
        tarpath = fcfile.split('.')[0]
        for currentPath, _, files in os.walk(tarpath):
            for stlfile in files:
                newstlfile = stlfile.replace(' ', '').replace('(', '').replace(')', '')
                oldstlfile = stlfile.replace(' ', '\ ').replace('(', '\(').replace(')', '\)')
                
                os.system(f'mv {tarpath}/{oldstlfile} {tarpath}/{newstlfile}')
                              
                # convert stl file to ASCII STL
                convert_ASCIISTL(os.path.join(currentPath, newstlfile))
            
            stave_paths.append(tarpath)
    
    # create xml file
    save_path = os.path.join(os.path.dirname(__file__), 'xml')
    if not os.path.exists(save_path):
        os.makedirs(save_path)
    

    save_file = os.path.join(save_path, args.output)
    print(save_file)
    print(stave_paths)
    create_xml.generate_xml(stave_paths[0], stave_paths[1], save_file)

                
                
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='This is a parameter description for the auto-generation of SVTOB xml file')
    parser.add_argument('--fcstd', '-f', type=str, nargs='+', default=['L3_stave_Long.FCStd', 'L4_stave_Long.FCStd'], help='input FCStd file from FreeCAD')
    parser.add_argument('--output', '-o', type=str, default='silicon_barrel.xml', help='the output file name eg. silicon_barrel.xml')
    
    args = parser.parse_args()
    
    generate_xml(args)