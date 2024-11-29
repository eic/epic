# this is a python script

#=======================================================================
#   Copyright (C) 2024 Univ. of Bham  All rights reserved.
#   
#   		FileName：		GeoTest.py
#   	 	Author：		LongLI <long.l@cern.ch>
#   		Time：			2024.11.26
#   		Description：
#
#======================================================================

import FreeCAD as fc
import Part   
import argparse


def geo_check():
    doc = fc.openDocument('data/L3_stave.FCStd')
    first = True
    for obj in doc.Objects:
        # print(obj.Label, obj.TypeId)
        if obj.TypeId not in "Part::Feature":
            continue

        # for Part::Feature
        hierarchy = 0
        while(obj.InList):
            parent = obj.InList[0]
            print(hierarchy)
            placement = parent.Placement
            print("before multiply:", placement)

            placement = parent.Placement.multiply(placement)

            print("after multiply:", placement)
            obj = parent


        
            
            

# if __name__ == '__main__':
#     parser = argparse.ArgumentParser(description='Geo check for OB stave')
#     parser.add_argument('--fcstd', '-f', type=str, default='./data/L4_LAS_only.FCStd', help='FreeCAD file to check')

#     args = parser.parse_args()


geo_check()
    