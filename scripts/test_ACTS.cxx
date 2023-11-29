// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Whitney Armstrong

#include "DD4hep/DetElement.h"
#include "DD4hep/Detector.h"
#include "DD4hep/Objects.h"
#include "DDG4/Geant4Data.h"
#include "DDRec/CellIDPositionConverter.h"
#include "DDRec/Surface.h"
#include "DDRec/SurfaceManager.h"

#include "TCanvas.h"
#include "TChain.h"

#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Plugins/DD4hep/ConvertDD4hepDetector.hpp"

/** Example loading ACTs.
 *
 *
 */
void test_ACTS(const char* compact = "epic.xml")
{

  using namespace ROOT::Math;
  // -------------------------
  // Get the DD4hep instance
  // Load the compact XML file
  // Initialize the position converter tool
  dd4hep::Detector& detector = dd4hep::Detector::getInstance();
  detector.fromCompact(compact);
  dd4hep::rec::CellIDPositionConverter cellid_converter(detector);

  auto logger = Acts::getDefaultLogger("Acts", Acts::Logging::Level::VERBOSE);
  auto acts_tracking_geometry = Acts::convertDD4hepDetector(detector.world(), *logger);

  if (acts_tracking_geometry) {
    std::cout << "success?\n";
  }
  //  if(acts_tracking_geometry->highestTrackingVolume()) {
  //    std::cout << " volume name \n ";
  //    std::cout << acts_tracking_geometry->highestTrackingVolume()->volumeName() << std::endl;
  //  } else {
  //    std::cout << "derp\n";
  //  }
  //}
}
