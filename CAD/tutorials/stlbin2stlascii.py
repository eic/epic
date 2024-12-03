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

binarystl = mesh.Mesh.from_file('data/L3_stave/Pads_033Mirror_1001.stl')
binarystl.save('data/Pads_033Mirror_1001_ASCII.stl', mode=stl.Mode.ASCII)
