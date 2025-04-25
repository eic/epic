[![CI status](https://github.com/eic/epic/actions/workflows/linux-eic-shell.yml/badge.svg)](https://github.com/eic/epic/actions/workflows/linux-eic-shell.yml)

Overview
--------

This branch contains a version of the epic detector geometry code, in which a detailed model of the SVT outer barrels is implemented from .gdml files created from the CAD drawings. This was originally developed by Tuna Tasali at the University of Oxford in August 2024. His version is available here: https://github.com/Tunat66/epicUK/tree/UK_Contribution

This version has been edited by Long Li (University of Birmingham) and Sam Henry (University of Oxford) to fix bugs, simplify, and add further tools.

This is intended for design studies of the outer barrels to investigate the impact of the design of detector performances and verify the simple geometry used for routine simulations is accurate. It may be of interested to other groups wishing to do similar studies on other detector components.



Getting Started
---------------

Get a copy of this branch from this repository:
```bash
git clone -b SVTOB_UK https://github.com/eic/epic.git
```

Configure, build, and install in eic-shell in the usual way:
```bash
cd epic
cmake -B build -S . -DCMAKE_INSTALL_PREFIX=install
cmake --build build
cmake --install build
```
Load this geometry:
```bash
source install/bin/thisepic.sh
```

Run simulations using ddsim as described in the ePIC software tutorials.


Description of new SVTOB geometry code
--------------------------------------------

The assembly of the L3 (Sagitta) and L4 (Outer) barrels is done in the file *BarrelTrackerOuter_standardized_geo.cpp*. This reads parameters from *silicon_barrel.xml*. The barrels are assembled from 46 (L3) or 70 (L4) staves. Each stave is divided into a large number of module components with the vertex points in a .gdml file. The .xml files contains links to the .gdml for each component and other information including the material and whether the component is a sensitive detector component. The class to read GDML files is defined in TGDMLParseBiggerFiles.cpp.

*BarrelTrackerOuter_standardized_geo.cpp* assembles the stave, and then places copies to build the cylinder at the radii specified in *silicon_barrel.xml* and *definitions_craterlake.xml*.

The GDML files are inside the CAD directory.

Converting CAD files to GDML
------------------------------

This section outlines how to convert a CAD file to a folder of GDML files and produce the *silicon_barrel.xml file*

To configure the CAD utility:
<ol>
<li>Download a FreeCAD AppImage and place it in a directory of your choice</li>
<li>Extract the AppImage: run the appimage (as if it is an executable) with the option '--appimage-extract'</li>
<li>A directory named 'squashfs-root' should have appeared next to the AppImage</li>
<li>Navigate to the directory of the freecadcmd executable, it is normally under 'squashfs-root/usr/bin/'</li>
<li>Get that directory with 'pwd' and paste it as the first line of the file 'freecadcmd_dir'</li>
</ol>
Copy files in CAD/CAD_to_GDML to the directory containing the .STEP, as well as freecadcmd_dir

To start using the CAD utility:
<ol>
<li>Open your STEP file in FreeCAD and save it as a FreeCAD(.FCStd) file with a name of your choice.</li>
<li>Source the script with arguments: source FCStd_to_gdml.sh <FreeCAD_file_name>
    NOTE: this file contains a few parameters, adjust them to your own needs: 
        MAX_COMPONENT_NUMBER: maximum number of sensitive <module_component/>s per <module/>
        dict_list: a map allowing file names to be mapped to materials and sensitive/insensitive etc. </li>
<li>A folder containing the .gdml files named <FreeCAD_file_name> will be generated.</li>
</ol>

*FCStd_to_gdml.sh* runs the following scripts:
<ul>
<li>refine_step_for_ePIC.py - to convert the .FCStd file to .stl files</li>
<li>BinaryToASCII.py - to convert the binary .stl files to ASCII</li>
<li>stl_gdml.py - to convert the ASCII .stl files to .gdml</li>
</ul>

On offset on the coordinate system can be set on lines 51-53 of *stl_gdml.py*

You can also produce the *silicon_barrel.xml* automatically. For this process you need to have produced two folders contining gdml files 
for components of the L3 stave and L4 stave. 
<ol>
<li>Run the create_xml.py script with arguments <L3_stave_gdml_folder_name> <L4_stave_gdml_folder_name>. Not it may be useful to use the full path length for the directories so it can be used in any location.</li>
<li>Move the silicon_barrel.xml file to the ../compact/tracking folder (make sure to backup the older version with a git commit!!)</li>
<li>Source the script: source ../recompile.sh to get the .xml file copied over to $DETECTOR_PATH and then you are done!</li>
</ol>



Visualization and debugging tools
-----------------------------------

To help debug code and verify the correct placement of the staves, I wrote some visualization tools which are in the *tools* directory.

*ExtractVertexCoordinates.C* is a ROOT C++ macro to plot the x,y coordinates of all vertices and draw circles at the radii of the barrels. This reads the *detector_geometry.root* file, which can be created by:
```
dd_web_display --export $DETECTOR_PATH/epic_silicon_barrel_only.xml
```
*epic_silicon_barrel_only.xml* is created from *epic_inner_detector_only.xml* by removing all components except the barrels. It can be copied from *tools* to the $DETECTOR_PATH directory.

Due to the very large number of vertices, it is not possible to view the detector_geometry.root file using the ROOT geoviewer unless the number of the staves is reduced to 1.

*read_edm4hep_SiBarrel.C* plots the x,y coordinartes of the silicon barrel hits in an edm4hep.root file. This can be produced by running ddsim. For example, to simulate 100000 muons:
```
ddsim --compactFile $DETECTOR_PATH/epic_craterlake_tracking_only.xml --outputFile TestHitMap.edm4hep.root --numberOfEvents 100000 --enableGun --gun
.thetaMin 63*deg --gun.thetaMax 117*deg --gun.particle mu- --gun.momentumMin 10*GeV --gun.momentumMax 11*GeV --gun.distribution eta --gun.multiplicity 1 
--random.seed 100000
```

*gdml_xy_vertices.py* will read all the .gdml files in a directory and write the x,y,z coordinates to a text file. These can be plotted using *plotXY.C*

*BarrelTrackerOuter_standardized_geo.cpp* contains a function *ExportVolumePoints(TGeoVolume* vol, string file_suffix)* which will extract the x,y,z coordinates of all vertices in a TGeoVolume object and write these to a text file.


Material Scans
---------------

Tracking performance benchmark tests
-----------------------------------------
