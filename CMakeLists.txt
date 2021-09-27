cmake_minimum_required(VERSION 3.3 FATAL_ERROR)
PROJECT(ip6
  LANGUAGES CXX
  DESCRIPTION "A template dd4hep+acts beamline"
  )

SET(CMAKE_CXX_STANDARD 17)
find_package(Acts REQUIRED COMPONENTS Core PluginTGeo PluginDD4hep )

#-----------------------------------------------------------------------------------
set(a_lib_name IP6)
dd4hep_configure_output()

dd4hep_add_plugin(${a_lib_name} 
  SOURCES 
  src/B0Tracker_geo.cpp
  src/B0Preshower_geo.cpp
  src/DD4hep_GdmlDetector.cpp
  src/ForwardRomanPot_geo.cpp
  src/hadronDownstreamBeamPipe.cpp
  src/IP6BeamPipe.cpp
  src/magnetVacuumFF.cpp
  src/OffMomentumTracker_geo.cpp
  src/SimpleDiskTracker_geo.cpp
  src/SimpleRectangularTracker_geo.cpp
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
