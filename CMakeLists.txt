cmake_minimum_required(VERSION 3.12 FATAL_ERROR)

# CMP0074: find_package() uses <PackageName>_ROOT variables
cmake_policy(SET CMP0074 NEW)

PROJECT(ip6
  LANGUAGES CXX
  DESCRIPTION "A template dd4hep+acts beamline"
  )

SET(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Werror")

find_package(DD4hep REQUIRED COMPONENTS DDCore DDG4)
find_package(Acts REQUIRED COMPONENTS Core PluginTGeo PluginDD4hep )

#-----------------------------------------------------------------------------------
set(a_lib_name IP6)
dd4hep_configure_output()

dd4hep_add_plugin(${a_lib_name} 
  SOURCES 
  src/B0Tracker_geo.cpp
  src/B0Preshower_geo.cpp
  src/BackwardsBeamPipe_geo.cpp
  src/BackwardsCollimator.cpp
  src/BackwardsLumiVac_geo.cpp
  src/BackwardsTagWin_geo.cpp
  src/BackwardsVacBox_geo.cpp
  src/BackwardsVacDrift_geo.cpp
  src/CylindricalDipoleMagnet_geo.cpp
  src/DD4hep_GdmlDetector.cpp
  src/ForwardRomanPot_geo.cpp
  src/hadronDownstreamBeamPipe.cpp
  src/IP6BeamPipe.cpp
  src/magnetVacuumFF.cpp
  src/OffMomentumTracker_geo.cpp
  src/TaggerCalWSi_geo.cpp
  src/TaggerTracker_geo.cpp
  src/Timepix_geo.cpp
  src/ZeroDegreeCalorimeterEcalWSciFi_geo.cpp
  src/ZeroDegreeCalorimeterEcal_geo.cpp
  src/ZeroDegreeCalorimeterSampling_geo.cpp
  USES ActsCore ActsPluginDD4hep ROOT::Core ROOT::Gdml
  )
target_link_libraries(${a_lib_name}
  PUBLIC DD4hep::DDCore  DD4hep::DDRec DD4hep::DDG4
  )

#-----------------------------------------------------------------------------------
install(DIRECTORY ip6/
    DESTINATION share/${PROJECT_NAME}/${a_lib_name}
    FILES_MATCHING PATTERN "*.xml"
    )
install(DIRECTORY ip6/
    DESTINATION share/${PROJECT_NAME}/ip6
    FILES_MATCHING PATTERN "*.xml"
    )
install(FILES eic_ip6.xml
    DESTINATION share/${PROJECT_NAME}/
    )
