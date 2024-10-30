Author: Tuna Tasali, epicUK
Date: 19/08/24

Description:
This CAD/ directory contains essential utilities to convert the components of a .STEP file into .gdml files
which can then be used by dd4hep to define geometry. The scripts here will produce two outputs under the directory CAD/
    1- the silicon_barrel.xml file
    2- a Folder with a user defined name containing the .gdml files which are referenced by the silicon_barrel.xml file
The conversion relies on the FreeCAD API so the utility needs some configuration before it can be used.

************************************************************************************************************************

To configure the CAD utility:

1-  Download a FreeCAD AppImage and place it in a directory of your choice
2-  Extract the AppImage: run the appimage (as if it is an executable) with the option '--appimage-extract'
3-  A directory named 'squashfs-root' should have appeared next to the AppImage
4-  Navigate to the directory of the freecadcmd executable, it is normally under 'squashfs-root/usr/bin/'
5-  Get that directory with 'pwd' and paste it as the first line of the file 'freecadcmd_dir'
6-  You are now ready to produce some step files. Yay.

************************************************************************************************************************

To start using the CAD utility:

1-  Open your STEP file in FreeCAD and save it as a FreeCAD(.FCStd) file with a name of your choice.
2-  Source the script with arguments: source FCStd_to_gdml.sh <FreeCAD_file_name>
    NOTE: this file contains a few parameters, adjust them to your own needs: 
        MAX_COMPONENT_NUMBER: maximum number of sensitive <module_component/>s per <module/>
        dict_list: a map allowing file names to be mapped to materials and sensitive/insensitive etc. 
3-  A folder containing the .gdml files named <FreeCAD_file_name> will be generated.

************************************************************************************************************************
You can also produce the silicon_barrel.xml automatically. For this process you need to have produced two folders contining gdml files 
for components of the L3 stave and L4 stave. 

1-  Run the create_xml.py script with arguments <L3_stave_gdml_folder_name> <L4_stave_gdml_folder_name>
2-  Move the silicon_barrel.xml file to the ../compact/tracking folder (make sure to backup the older version with a git commit!!)
3-  Source the script: source ../recompile.sh to get the .xml file copied over to $DETECTOR_PATH and then you are done!
