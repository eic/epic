// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 - 2024, Whitney Armstrong, Wouter Deconinck

#include <DD4hep/Detector.h>

#include <Acts/Detector/CylindricalContainerBuilder.hpp>
#include <Acts/Detector/Detector.hpp>
#include <Acts/Detector/DetectorBuilder.hpp>
#include <Acts/Detector/detail/BlueprintDrawer.hpp>
#include <Acts/Detector/detail/BlueprintHelper.hpp>
#include <Acts/Geometry/TrackingGeometry.hpp>
#include <Acts/Geometry/TrackingVolume.hpp>
#include <Acts/Surfaces/Surface.hpp>
#include <Acts/Surfaces/SurfaceBounds.hpp>
#include <Acts/Utilities/Enumerate.hpp>
#include <Acts/Utilities/Logger.hpp>
#include <Acts/Visualization/GeometryView3D.hpp>
#include <Acts/Visualization/ObjVisualization3D.hpp>
#include <Acts/Visualization/PlyVisualization3D.hpp>
#include <Acts/Visualization/ViewConfig.hpp>

#include <ActsPlugins/DD4hep/DD4hepBlueprintFactory.hpp>
#include <ActsPlugins/DD4hep/DD4hepDetectorStructure.hpp>
#include <ActsPlugins/DD4hep/DD4hepDetectorSurfaceFactory.hpp>
#include <ActsPlugins/DD4hep/DD4hepLayerStructure.hpp>

/** Example loading ACTs.
 *
 *
 */
void test_ACTS_gen3(const char* compact = "epic.xml") {
  // -------------------------
  // Get the DD4hep instance
  // Load the compact XML file
  // Initialize the position converter tool
  auto dd4hep_detector = dd4hep::Detector::make_unique("");
  dd4hep_detector->fromCompact(compact);
  auto dd4hep_world = dd4hep_detector->world();

  Acts::GeometryContext trackingGeoCtx;

  ActsPlugins::DD4hepDetectorSurfaceFactory::Config surfaceFactoryConfig;
  auto surfaceFactory = std::make_shared<ActsPlugins::DD4hepDetectorSurfaceFactory>(
      surfaceFactoryConfig,
      Acts::getDefaultLogger("DD4hepDetectorSurfaceFactory", Acts::Logging::VERBOSE));
  auto layerStructure = std::make_shared<ActsPlugins::DD4hepLayerStructure>(
      std::move(surfaceFactory),
      Acts::getDefaultLogger("DD4hepLayerStructure", Acts::Logging::VERBOSE));
  ActsPlugins::DD4hepBlueprintFactory::Config blueprintConfig{layerStructure};
  ActsPlugins::DD4hepBlueprintFactory::Cache blueprintCache;
  ActsPlugins::DD4hepBlueprintFactory blueprint(
      blueprintConfig, Acts::getDefaultLogger("DD4hepBlueprintFactory", Acts::Logging::VERBOSE));
  auto dd4hepBlueprint = blueprint.create(blueprintCache, trackingGeoCtx, dd4hep_world);

  // Now fill the gaps
  Acts::Experimental::detail::BlueprintHelper::fillGaps(*dd4hepBlueprint);

  // dot -P -Tsvg -o plugins.svg
  std::ofstream cbp("cylindrical_detector_dd4hep.dot");
  Acts::Experimental::detail::BlueprintDrawer::dotStream(cbp, *dd4hepBlueprint);
  cbp.close();

  // Create a Cylindrical detector builder from this blueprint
  auto detectorBuilder = std::make_shared<Acts::Experimental::CylindricalContainerBuilder>(
      *dd4hepBlueprint, Acts::Logging::VERBOSE);

  // Detector builder
  Acts::Experimental::DetectorBuilder::Config detectorBuilderConfig;
  detectorBuilderConfig.auxiliary = "*** Test : auto generated cylindrical detector builder  ***";
  detectorBuilderConfig.name      = "Cylindrical detector from blueprint";
  detectorBuilderConfig.builder   = detectorBuilder;
  detectorBuilderConfig.geoIdGenerator = dd4hepBlueprint->geoIdGenerator;
  auto detector =
      Acts::Experimental::DetectorBuilder(detectorBuilderConfig).construct(trackingGeoCtx);

  // Export to obj+mtl and ply collections
  Acts::ViewConfig containerView{.color = {220, 220, 220}}; // alto
  Acts::ViewConfig volumeView{.color = {220, 220, 0}};      // barberry yellow
  Acts::ViewConfig sensitiveView{.color = {0, 180, 240}};   // picton blue
  Acts::ViewConfig passiveView{.color = {240, 180, 0}};     // lightning yellow
  Acts::ViewConfig gridView{.color = {220, 0, 0}};          // scarlet red

  // Get the root detector volumes
  for (const auto* volume : detector->rootVolumes()) {

    // Export to obj+mtl
    Acts::ObjVisualization3D objVis;
    Acts::GeometryView3D::drawDetectorVolume(objVis, *volume, trackingGeoCtx);
    // Export to ply
    Acts::PlyVisualization3D plyVis;
    Acts::GeometryView3D::drawDetectorVolume(plyVis, *volume, trackingGeoCtx);
  }
}
