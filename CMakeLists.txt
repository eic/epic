cmake_minimum_required(VERSION 3.12 FATAL_ERROR)

# CMP0074: find_package() uses <PackageName>_ROOT variables
cmake_policy(SET CMP0074 NEW)

PROJECT(ecce
  LANGUAGES CXX
  DESCRIPTION "A template dd4hep+acts detector"
  )

# C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD 17 CACHE STRING "Set the C++ standard to be used")
if(NOT CMAKE_CXX_STANDARD MATCHES "17|20|23")
  message(FATAL_ERROR "Unsupported C++ standard: ${CMAKE_CXX_STANDARD}")
endif()
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Build type
if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "Build type configuration" FORCE)
  message(STATUS "Setting default build type: ${CMAKE_BUILD_TYPE}")
endif()

# Error on all warnings
if(NOT CMAKE_BUILD_TYPE STREQUAL "Release")
  add_compile_options(-Wall -Wextra -Werror -pedantic)
endif()

# Export compile commands as json for run-clang-tidy
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

find_package(DD4hep REQUIRED COMPONENTS DDCore DDG4)
find_package(ActsDD4hep)
if(ActsDD4hep_FOUND)
  add_compile_definitions(USE_ACTSDD4HEP)
  set(ActsDD4hep ActsDD4hep::ActsDD4hep)
else()
  find_package(Acts REQUIRED COMPONENTS Core PluginIdentification PluginTGeo PluginDD4hep)
  set(ActsDD4hep ActsCore ActsPluginDD4hep)
endif()
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
  USES ${ActsDD4hep}
  )
target_link_libraries(${a_lib_name}
  PUBLIC DD4hep::DDCore DD4hep::DDRec fmt::fmt
  )

#-----------------------------------------------------------------------------------
# Parse jinja templates: once by default, and once for all json files

set(TEMPLATE ${CMAKE_CURRENT_SOURCE_DIR}/templates/${PROJECT_NAME}.xml.jinja2)

add_custom_target(${PROJECT_NAME}.xml ALL)
add_custom_command(TARGET ${PROJECT_NAME}.xml
    COMMAND jinja2 ${TEMPLATE} -o ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}.xml
    COMMENT "Creating default ${PROJECT_NAME}.xml geometry"
    )
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}.xml
    DESTINATION share/${PROJECT_NAME}/
    )

file(GLOB CONFIG_JSONS ${CMAKE_CURRENT_SOURCE_DIR}/configurations/*.json)
foreach(config_json ${CONFIG_JSONS})
    get_filename_component(config ${config_json} NAME_WE)
    add_custom_target(${config}.xml ALL)
    add_custom_command(TARGET ${config}.xml
        COMMAND jinja2 ${TEMPLATE} ${config_json} -o ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}_${config}.xml
        COMMENT "Creating configuration ${config} for ${PROJECT_NAME}"
        )
    install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}_${config}.xml
        DESTINATION share/${PROJECT_NAME}/
        )
endforeach()

#-----------------------------------------------------------------------------------
# Install the detector description files.
install(DIRECTORY compact/
    DESTINATION share/${PROJECT_NAME}/compact
    FILES_MATCHING PATTERN "*.xml"
    )

#-----------------------------------------------------------------------------------
# Install the detector calibration files.
install(DIRECTORY calibrations/
    DESTINATION share/${PROJECT_NAME}/calibrations
    )
