// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 - 2024, Whitney Armstrong, Wouter Deconinck

#include "DD4hep/Detector.h"

#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Plugins/DD4hep/ConvertDD4hepDetector.hpp"

/** Example loading ACTs.
 *
 *
 */
void test_ACTS(const char* compact = "epic.xml") {
  // -------------------------
  // Get the DD4hep instance
  // Load the compact XML file
  // Initialize the position converter tool
  auto detector = dd4hep::Detector::make_unique("");
  detector->fromCompact(compact);

  auto logger                 = Acts::getDefaultLogger("Acts", Acts::Logging::Level::VERBOSE);
  auto acts_tracking_geometry = Acts::convertDD4hepDetector(detector->world(), *logger);

  // Visit all surfaces
  acts_tracking_geometry->visitSurfaces([](const Acts::Surface* surface) {});
}
