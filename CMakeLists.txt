cmake_minimum_required(VERSION 3.12 FATAL_ERROR)

# CMP0074: find_package() uses <PackageName>_ROOT variables
cmake_policy(SET CMP0074 NEW)

PROJECT(ip6
  LANGUAGES CXX
  DESCRIPTION "A template dd4hep+acts beamline"
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

# Install into GNU standard directories
include(GNUInstallDirs)

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
if(${DD4hep_VERSION} VERSION_LESS 1.21)
  message(WARNING "DD4hep before 1.21 does not write collections correctly \n"
                  "More info at https://github.com/AIDASoft/DD4hep/pull/922")
endif()
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
set(a_lib_name ${PROJECT_NAME})

dd4hep_configure_output()

file(GLOB sources CONFIGURE_DEPENDS src/*.cpp)
dd4hep_add_plugin(${a_lib_name}
  SOURCES ${sources}
  USES ${ActsDD4hep} ROOT::Core ROOT::Gdml
  )
target_link_libraries(${a_lib_name}
  PUBLIC ${DD4hep_required_libraries}
  )

#-----------------------------------------------------------------------------------
install(DIRECTORY ip6/
    DESTINATION ${CMAKE_INSTALL_DATADIR}/${PROJECT_NAME}/${a_lib_name}
    FILES_MATCHING PATTERN "*.xml"
    )
install(DIRECTORY ip6/
    DESTINATION ${CMAKE_INSTALL_DATADIR}/${PROJECT_NAME}/ip6
    FILES_MATCHING PATTERN "*.xml"
    )
install(FILES ip6.xml
    DESTINATION ${CMAKE_INSTALL_DATADIR}/${PROJECT_NAME}/
    )

#-----------------------------------------------------------------------------------
# Configure and install beamline setup script
#
execute_process(
  COMMAND git rev-parse --abbrev-ref HEAD
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE GIT_BRANCH
  OUTPUT_STRIP_TRAILING_WHITESPACE
)
configure_file(templates/setup.sh.in ${CMAKE_CURRENT_BINARY_DIR}/setup.sh @ONLY)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/setup.sh
    DESTINATION ${CMAKE_INSTALL_DATADIR}/ip6/ip6
)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/setup.sh
    DESTINATION ${CMAKE_INSTALL_PREFIX}
)
