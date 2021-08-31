cmake_minimum_required(VERSION 3.3 FATAL_ERROR)

SET(CMAKE_CXX_STANDARD 17)

find_package(Acts REQUIRED COMPONENTS Core PluginIdentification PluginTGeo PluginDD4hep )

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
install(DIRECTORY compact/
    DESTINATION share/${PROJECT_NAME}/${a_lib_name}
    FILES_MATCHING PATTERN "*.xml"
    )

