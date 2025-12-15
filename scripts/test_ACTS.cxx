// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 - 2024, Whitney Armstrong, Wouter Deconinck

#include "DD4hep/Detector.h"

#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include <Acts/Visualization/GeometryView3D.hpp>
#include <Acts/Visualization/ObjVisualization3D.hpp>
#include <Acts/Visualization/PlyVisualization3D.hpp>
#include <Acts/Visualization/ViewConfig.hpp>

#if __has_include(<ActsPlugins/DD4hep/ConvertDD4hepDetector.hpp>)
// Acts_MAJOR_VERSION >= 44
#include <ActsPlugins/DD4hep/ConvertDD4hepDetector.hpp>
using ActsPlugins::convertDD4hepDetector;
#else
// Acts_MAJOR_VERSION < 44
#include <Acts/Plugins/DD4hep/ConvertDD4hepDetector.hpp>
using Acts::convertDD4hepDetector;
#endif

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
  auto acts_tracking_geometry = convertDD4hepDetector(detector->world(), *logger);

  // Visit all surfaces
  acts_tracking_geometry->visitSurfaces([](const Acts::Surface* surface) {});

  // Export to obj+mtl and ply collections
  Acts::ViewConfig containerView{.color = {220, 220, 220}}; // alto
  Acts::ViewConfig volumeView{.color = {220, 220, 0}};      // barberry yellow
  Acts::ViewConfig sensitiveView{.color = {0, 180, 240}};   // picton blue
  Acts::ViewConfig passiveView{.color = {240, 180, 0}};     // lightning yellow
  Acts::ViewConfig gridView{.color = {220, 0, 0}};          // scarlet red
  Acts::GeometryContext trackingGeoCtx;
  const Acts::TrackingVolume* world = acts_tracking_geometry->highestTrackingVolume();
  // Export to obj+mtl
  Acts::ObjVisualization3D objVis;
  Acts::GeometryView3D::drawTrackingVolume(objVis, *world, trackingGeoCtx, containerView,
                                           volumeView, passiveView, sensitiveView, gridView, true,
                                           "", "");
  // Export to ply
  Acts::PlyVisualization3D plyVis;
  Acts::GeometryView3D::drawTrackingVolume(plyVis, *world, trackingGeoCtx, containerView,
                                           volumeView, passiveView, sensitiveView, gridView, true,
                                           "", "");
}
