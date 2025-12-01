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

#include <fstream>

/**
 * Example: Loading and testing ACTS Gen3 geometry from a DD4hep compact file.
 *
 * Gen3 geometry in ACTS refers to the third-generation detector geometry description,
 * which introduces improved modularity, flexibility, and support for advanced features
 * compared to Gen1. Gen3 geometry typically uses blueprints and container builders to
 * construct detector volumes, allowing for more complex and realistic detector layouts.
 *
 * This test function demonstrates how to:
 *   - Load a DD4hep compact XML file describing the detector geometry.
 *   - Convert the DD4hep geometry into ACTS Gen3 blueprints.
 *   - Fill gaps in the blueprint to ensure completeness.
 *   - Export the geometry to DOT, OBJ, and PLY formats for visualization.
 *   - Build a cylindrical detector using the Gen3 blueprint and visualize its volumes.
 *
 * Compared to Gen1, Gen3 geometry provides a more robust and extensible framework for
 * detector description, making it suitable for modern experiments and simulation workflows.
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

  // Generate dot file (can be converted to SVG with: dot -Tsvg cylindrical_detector_dd4hep.dot -o output.svg)
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
