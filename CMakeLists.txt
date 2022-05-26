cmake_minimum_required(VERSION 3.12 FATAL_ERROR)

# CMP0074: find_package() uses <PackageName>_ROOT variables
cmake_policy(SET CMP0074 NEW)

PROJECT(ecce
  VERSION 0.0.1
  LANGUAGES CXX
  DESCRIPTION "A template dd4hep+acts detector"
  )

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror")

find_package(DD4hep REQUIRED COMPONENTS DDCore DDG4)
find_package(ActsDD4hep REQUIRED)
find_package(fmt REQUIRED)

#-----------------------------------------------------------------------------------
set(a_lib_name ecce)

dd4hep_configure_output()

dd4hep_add_plugin(${a_lib_name} SOURCES 
  src/BarrelBarDetectorWithSideFrame_geo.cpp
  src/BarrelCalorimeter_geo.cpp
  src/BarrelCalorimeterInterlayers_geo.cpp
  src/BarrelTrackerWithFrame_geo.cpp
  src/CompositeTracker_geo.cpp
  src/CylinderTrackerBarrel_geo.cpp
  src/DIRC_geo.cpp
  src/DRICH_geo.cpp
  src/PFRICH_geo.cpp
  src/FileLoader.cpp
  src/FieldMapBrBz.cpp
  src/GaseousRICH_geo.cpp
  src/GeometryHelpers.cpp
  src/HomogeneousCalorimeter_geo.cpp
  src/HybridCalorimeter_geo.cpp
  src/MRich_geo.cpp
  src/PolyhedraEndcapCalorimeter2_geo.cpp
  src/ScFiCalorimeter_geo.cpp
  src/SciGlassCalorimeter_geo.cpp
  src/ShashlikCalorimeter_geo.cpp
  src/SimpleDiskDetector_geo.cpp
  src/Solenoid_geo.cpp
  src/TrapEndcapTracker_geo.cpp
  src/TrackerSupport_geo.cpp
  USES ActsDD4hep::ActsDD4hep
  )
target_link_libraries(${a_lib_name}
  PUBLIC DD4hep::DDCore DD4hep::DDRec fmt::fmt
  )

#-----------------------------------------------------------------------------------
# Install the detector description files.
install(DIRECTORY compact/
    DESTINATION share/${PROJECT_NAME}/${a_lib_name}
    FILES_MATCHING PATTERN "*.xml"
    )
install(DIRECTORY compact/
    DESTINATION share/${PROJECT_NAME}/compact
    FILES_MATCHING PATTERN "*.xml"
    )
install(FILES ecce.xml
    DESTINATION share/${PROJECT_NAME}/
    )
#-----------------------------------------------------------------------------------
# Install the detector calibration files.
install(DIRECTORY calibrations/
    DESTINATION share/${PROJECT_NAME}/calibrations
    )
