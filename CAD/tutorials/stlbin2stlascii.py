# this is a python script

#=======================================================================
#   Copyright (C) 2024 Univ. of Bham  All rights reserved.
#   
#   		FileName：		stlbin2stlascii.py
#   	 	Author：		LongLI <long.l@cern.ch>
#   		Time：			2024.11.27
#   		Description：
#
#======================================================================


import stl
from stl import mesh

binarystl = mesh.Mesh.from_file('data/L3_SolidCarbonFoam.stl')
binarystl.save('data/L3_SolidCarbonFoam_ASCII.stl', mode=stl.Mode.ASCII)