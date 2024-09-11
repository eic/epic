#!/bin/bash
set -e
# script for material map validation with ACTS python bindings
# run as : ./run_material_map_validation.sh --nevents 1000
# Shujie Li, 03. 2024 (https://github.com/eic/snippets/pull/3)

# Check if DETECTOR_PATH and DETECTOR_CONFIG are set
if [[ -z ${DETECTOR_PATH} || -z ${DETECTOR_CONFIG} ]] ; then
  echo "You must set \$DETECTOR_PATH and \$DETECTOR_CONFIG before running this script."
  exit -1
fi

# Download required Acts files
ACTS_VERSION="00591a593a648430820e980b031301d25c18f1c7" # v33.1.0
ACTS_URL="https://github.com/acts-project/acts/raw/"
ACTS_FILES=(
  "Examples/Scripts/Python/geometry.py"
  "Examples/Scripts/Python/material_mapping.py"
  "Examples/Scripts/Python/material_recording.py"
  "Examples/Scripts/Python/material_validation.py"
  "Examples/Scripts/MaterialMapping/writeMapConfig.py"
  "Examples/Scripts/MaterialMapping/configureMap.py"
  "Examples/Scripts/MaterialMapping/Mat_map.C"
  "Examples/Scripts/MaterialMapping/Mat_map_surface_plot.C"
  "Examples/Scripts/MaterialMapping/Mat_map_surface_plot_ratio.C"
  "Examples/Scripts/MaterialMapping/Mat_map_surface_plot_dist.C"
  "Examples/Scripts/MaterialMapping/Mat_map_surface_plot_1D.C"
  "Examples/Scripts/MaterialMapping/materialPlotHelper.cpp"
  "Examples/Scripts/MaterialMapping/materialPlotHelper.hpp"
)
for file in ${ACTS_FILES[@]} ; do
  if [ ! -f ${file} ] ; then
    curl --silent --location --create-dirs --output ${file} ${ACTS_URL}/${ACTS_VERSION}/${file}
  fi
done
export PYTHONPATH=$PWD/Examples/Scripts/Python:$PYTHONPATH

# Default arguments
nevents=1000
nparticles=1000
while [[ $# -gt 1 ]]
do
  key="$1"
  case $key in
    --nevents)
      nevents=$2
      shift # past value
      shift
      ;;
    --nparticles)
      nparticles=$2
      shift # past value
      shift
      ;;
    *)    # unknown option
      #POSITIONAL+=("$1") # save it in an array for later
      echo "unknown option $1"
      print_the_help
      shift # past argument
      ;;
  esac
done
set -- "${POSITIONAL[@]}" # restore positional parameters

recordingFile=geant4_material_tracks.root
geoFile=geometry-map.json
matFile=material-map.json
trackFile=material-map_tracks.root
propFile=propagation_material

echo "::group::----GEANTINO SCAN------"
# output geant4_material_tracks.root
# The result of the geantino scan will be a root file containing material tracks. Those contain the direction and production vertex of the geantino, the total material accumulated and all the interaction points in the detector.
python material_recording_epic.py -i ${DETECTOR_PATH}/${DETECTOR_CONFIG}.xml -n ${nevents} -t ${nparticles} -o ${recordingFile}
echo "::endgroup::"

echo "::group::-----MAPPING Configuration-----"
# map geometry to geometry-map.json
python geometry_epic.py -i ${DETECTOR_PATH}/${DETECTOR_CONFIG}.xml

# take geometry-map.json and read out config-map.json
python Examples/Scripts/MaterialMapping/writeMapConfig.py ${geoFile} config-map.json

# turn on approaches and beampipe surfaces for material mapping
# you can always manually adjust the mapmaterial flag and binnings in config-map.json
python materialmap_config.py -i config-map.json -o config-map_new.json

# turn config-map.json into modified geometry-map.json
python Examples/Scripts/MaterialMapping/configureMap.py ${geoFile} config-map_new.json
echo "::endgroup::"

echo "::group::----MAPPING------------"
# input: geant4_material_tracks.root, geometry-map.json
# output: material-maps.json or cbor. This is the material map that you want to provide to EICrecon, i.e.  -Pacts:MaterialMap=XXX  .Please --matFile to specify the name and type
#         material-maps_tracks.root(recorded steps from geantino, for validation purpose)
python material_mapping_epic.py --xmlFile ${DETECTOR_PATH}/${DETECTOR_CONFIG}.xml --geoFile ${geoFile} --matFile ${matFile}
echo "::endgroup::"

echo "::group::----Prepare validation rootfile--------"
# output propagation-material.root
python material_validation_epic.py --xmlFile ${DETECTOR_PATH}/${DETECTOR_CONFIG}.xml --outputName ${propFile}_new --matFile ${matFile} -n ${nevents}
python material_validation_epic.py --xmlFile ${DETECTOR_PATH}/${DETECTOR_CONFIG}.xml --outputName ${propFile}_old --matFile "calibrations/materials-map.cbor" -n ${nevents}
echo "::endgroup::"

echo "::group::-------Comparison plots---------"
rm -rf Validation/new
mkdir -p Validation/new
root -l -b -q Examples/Scripts/MaterialMapping/Mat_map.C'("'${propFile}_new'.root","'${trackFile}'","Validation/new")'
rm -rf Validation/old
mkdir -p Validation/old
root -l -b -q Examples/Scripts/MaterialMapping/Mat_map.C'("'${propFile}_old'.root","'${trackFile}'","Validation/old")'

rm -rf Surfaces
mkdir -p Surfaces/new/ratio_plot
mkdir -p Surfaces/new/prop_plot
mkdir -p Surfaces/new/map_plot
mkdir -p Surfaces/old/ratio_plot
mkdir -p Surfaces/old/prop_plot
mkdir -p Surfaces/old/map_plot
mkdir -p Surfaces/dist_plot
mkdir -p Surfaces/1D_plot

root -l -b -q Examples/Scripts/MaterialMapping/Mat_map_surface_plot_ratio.C'("'${propFile}_new'.root","'${trackFile}'",-1,"Surfaces/new/ratio_plot","Surfaces/new/prop_plot","Surfaces/new/map_plot")'
root -l -b -q Examples/Scripts/MaterialMapping/Mat_map_surface_plot_ratio.C'("'${propFile}_old'.root","'${trackFile}'",-1,"Surfaces/old/ratio_plot","Surfaces/old/prop_plot","Surfaces/old/map_plot")'
root -l -b -q Examples/Scripts/MaterialMapping/Mat_map_surface_plot_dist.C'("'${trackFile}'",-1,"Surfaces/dist_plot")'
root -l -b -q Examples/Scripts/MaterialMapping/Mat_map_surface_plot_1D.C'("'${trackFile}'",-1,"Surfaces/1D_plot")'
echo "::endgroup::"
