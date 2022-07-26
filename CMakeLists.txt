cmake_minimum_required(VERSION 3.12 FATAL_ERROR)

# CMP0074: find_package() uses <PackageName>_ROOT variables
cmake_policy(SET CMP0074 NEW)

PROJECT(epic
  LANGUAGES CXX
  DESCRIPTION "DD4hep Geometry Description of the EPIC Experiment"
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

# Set required DD4hep components based on user input
set(USE_DDG4 "TRUE" CACHE BOOL "Require DD4hep DDG4 component")
set(DD4hep_required_components DDCore DDRec)
if(${USE_DDG4})
  list(APPEND DD4hep_required_components DDG4)
endif()

# The items passed to target_link_libraries below need DD4hep:: prepended
set(DD4hep_required_libraries ${DD4hep_required_components})
list(TRANSFORM DD4hep_required_libraries PREPEND DD4hep::)

# Dependencies
find_package(DD4hep REQUIRED COMPONENTS ${DD4hep_required_components})
find_package(ActsDD4hep)
if(ActsDD4hep_FOUND)
  add_compile_definitions(USE_ACTSDD4HEP)
  set(ActsDD4hep ActsDD4hep::ActsDD4hep)
else()
  find_package(Acts REQUIRED COMPONENTS Core PluginIdentification PluginTGeo PluginDD4hep)
  set(ActsDD4hep ActsCore ActsPluginDD4hep)
endif()
find_package(fmt REQUIRED)

if(DEFINED IRT_AUXFILE)
  add_compile_definitions(IRT_AUXFILE)
  find_package(IRT REQUIRED)
endif()

#-----------------------------------------------------------------------------------
set(a_lib_name ${PROJECT_NAME})

dd4hep_configure_output()

file(GLOB sources CONFIGURE_DEPENDS src/*.cpp)
dd4hep_add_plugin(${a_lib_name}
  SOURCES ${sources}
  USES ${ActsDD4hep}
  )
target_link_libraries(${a_lib_name}
  PUBLIC ${DD4hep_required_libraries} fmt::fmt
  )
if(DEFINED IRT_AUXFILE)
  target_link_libraries(${a_lib_name} PUBLIC IRT)
endif()

#-----------------------------------------------------------------------------------
# Parse jinja templates: once by default, and once for all yml files
#
# FIXME: should not write rendered templates to ${CMAKE_CURRENT_SOURCE_DIR}

set(TEMPLATE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/templates)
set(TEMPLATE_XML ${PROJECT_NAME}.xml.jinja2)

add_custom_target(${PROJECT_NAME}.xml ALL)
add_custom_command(TARGET ${PROJECT_NAME}.xml
    COMMAND
        ${CMAKE_CURRENT_SOURCE_DIR}/bin/make_detector_configuration
            --dir ${TEMPLATE_DIR}
            --template ${TEMPLATE_XML}
            --output ${CMAKE_CURRENT_SOURCE_DIR}/${PROJECT_NAME}.xml
    COMMENT "Creating default ${PROJECT_NAME}.xml geometry"
    )
install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/${PROJECT_NAME}.xml
    DESTINATION share/${PROJECT_NAME}/
    )

file(GLOB CONFIG_YMLS ${CMAKE_CURRENT_SOURCE_DIR}/configurations/*.yml)
foreach(config_yml ${CONFIG_YMLS})
    get_filename_component(config ${config_yml} NAME_WE)
    add_custom_target(${config}.xml ALL)
    add_custom_command(TARGET ${config}.xml
        COMMAND
            ${CMAKE_CURRENT_SOURCE_DIR}/bin/make_detector_configuration
                --dir ${TEMPLATE_DIR}
                --template ${TEMPLATE_XML}
                --config ${config_yml}
                --output ${CMAKE_CURRENT_SOURCE_DIR}/${PROJECT_NAME}_${config}.xml
        COMMENT "Creating configuration ${config} for ${PROJECT_NAME}"
        )
    install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/${PROJECT_NAME}_${config}.xml
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

#-----------------------------------------------------------------------------------
# Configure and install detector setup script
#
execute_process(
  COMMAND git rev-parse --abbrev-ref HEAD
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE GIT_BRANCH
  OUTPUT_STRIP_TRAILING_WHITESPACE
)
configure_file(templates/setup.sh.in ${CMAKE_CURRENT_BINARY_DIR}/setup.sh @ONLY)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/setup.sh
    DESTINATION ${CMAKE_INSTALL_PREFIX}
)
