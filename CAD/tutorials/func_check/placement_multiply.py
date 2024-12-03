# this is a python script

#=======================================================================
#   Copyright (C) 2024 Univ. of Bham  All rights reserved.
#   
#   		FileName：		placement_multiply.py
#   	 	Author：		LongLI <long.l@cern.ch>
#   		Time：			2024.12.02
#   		Description：
#
#======================================================================


import FreeCAD

# Create two placement objects
placement1 = FreeCAD.Placement(FreeCAD.Vector(10, 0, 0), FreeCAD.Rotation(FreeCAD.Vector(0, 0, 1), 45))  # Translation + Rotation
placement2 = FreeCAD.Placement(FreeCAD.Vector(0, 10, 0), FreeCAD.Rotation(FreeCAD.Vector(1, 0, 0), 30))  # Another transformation

# Combine the placements

print("Placement1:", placement1)
print("Placement2:", placement2)
result_placement = placement1.multiply(placement2)

# Access the combined base and rotation
print("P1 multiply P2:", result_placement)
print("Resulting Base:", result_placement.Base)
print("Resulting Rotation:", result_placement.Rotation)

