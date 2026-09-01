// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026
//
// Forward MPGD endcap with four independently placed quadrant assemblies.
// The overlap panels are declared in compact XML; no boundary coordinates are
// encoded in this driver.

#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/Shapes.h"
#include "DDRec/DetectorData.h"
#include "DDRec/Surface.h"
#include "DD4hepDetectorHelper.h"
#include "XML/Utilities.h"

#include <cmath>
#include <map>
#include <string>
#include <utility>
#include <vector>
#include <stdexcept>

using namespace dd4hep;
using namespace dd4hep::detail;
using namespace dd4hep::rec;

namespace {

using Placements = std::vector<PlacedVolume>;

struct ModulePrototype {
  Volume volume;
  Placements sensitivePlacements;
  std::vector<VolPlane> surfaces;
};

struct ModulePair {
  ModulePrototype normal;
  ModulePrototype overlap;
};

ModulePrototype buildPrototype(Detector& description, SensitiveDetector& sensitiveDetector,
                               xml_comp_t xModule, const std::string& suffix, bool isOverlap) {
  const std::string moduleName = xModule.nameStr() + suffix;
  const xml_comp_t trd         = xModule.trd();
  const double x1              = trd.x1();
  const double x2              = trd.x2();
  const double halfLength      = trd.z();

  double totalThickness = 0.0;
  for (xml_coll_t component(xModule, _U(module_component)); component; ++component) {
    totalThickness += xml_comp_t(component).thickness();
  }

  Solid moduleSolid;
  double overlapHalfWidth  = 0.0;
  double overlapHalfLength = 0.0;
  if (isOverlap) {
    xml_comp_t xOverlap = xModule.child(_U(overlap));
    overlapHalfWidth    = xOverlap.x();
    overlapHalfLength   = xOverlap.z();
    moduleSolid         = Box(overlapHalfWidth, totalThickness / 2.0, overlapHalfLength);
  } else {
    moduleSolid = Trapezoid(x1, x2, totalThickness / 2.0, totalThickness / 2.0, halfLength);
  }

  ModulePrototype prototype;
  prototype.volume = Volume(moduleName, moduleSolid, description.vacuum());
  prototype.volume.setVisAttributes(description.visAttributes(xModule.visStr()));

  double positionY       = -totalThickness / 2.0;
  double thicknessBefore = 0.0;
  int componentId        = 0;
  int sensorId           = 1;
  for (xml_coll_t component(xModule, _U(module_component)); component; ++component, ++componentId) {
    xml_comp_t xComponent  = component;
    const double thickness = xComponent.thickness();
    const std::string componentName =
        getAttrOrDefault(xComponent, _Unicode(name), _toString(componentId, "component%d"));

    Solid componentSolid;
    if (isOverlap) {
      componentSolid = Box(overlapHalfWidth, thickness / 2.0, overlapHalfLength);
    } else {
      const double componentX1         = getAttrOrDefault(xComponent, _Unicode(x1), x1);
      const double componentX2         = getAttrOrDefault(xComponent, _Unicode(x2), x2);
      const double componentHalfLength = getAttrOrDefault(xComponent, _Unicode(height), halfLength);
      componentSolid = Trapezoid(componentX1, componentX2, thickness / 2.0, thickness / 2.0,
                                 componentHalfLength);
    }

    Volume componentVolume(componentName + suffix, componentSolid,
                           description.material(xComponent.materialStr()));
    componentVolume.setVisAttributes(description.visAttributes(xComponent.visStr()));
    PlacedVolume componentPlacement = prototype.volume.placeVolume(
        componentVolume, Position(0.0, positionY + thickness / 2.0, 0.0));

    if (xComponent.isSensitive()) {
      if (sensorId > 2) {
        throw std::runtime_error(
            "MPGDQuadrantEndcapTracker: at most two sensitive components are supported");
      }
      componentPlacement.addPhysVolID("sensor", sensorId++);
      componentVolume.setSensitiveDetector(sensitiveDetector);
      prototype.sensitivePlacements.push_back(componentPlacement);

      const double innerThickness = thicknessBefore + thickness / 2.0;
      const double outerThickness = totalThickness - innerThickness;
      const SurfaceType type(SurfaceType::Sensitive);
      const Vector3D u(0.0, 0.0, -1.0);
      const Vector3D v(-1.0, 0.0, 0.0);
      const Vector3D n(0.0, 1.0, 0.0);
      prototype.surfaces.emplace_back(componentVolume, type, innerThickness, outerThickness, u, v,
                                      n);
    }
    positionY += thickness;
    thicknessBefore += thickness;
  }
  return prototype;
}

void addSensitiveDetElements(DetElement& moduleElement, const ModulePrototype& prototype,
                             int moduleId) {
  for (std::size_t index = 0; index < prototype.sensitivePlacements.size(); ++index) {
    const PlacedVolume& placement = prototype.sensitivePlacements.at(index);
    DetElement sensorElement(moduleElement, placement.volume().name(), moduleId);
    auto& parameters =
        DD4hepDetectorHelper::ensureExtension<dd4hep::rec::VariantParameters>(sensorElement);
    parameters.set<std::string>("axis_definitions", "XZY");
    sensorElement.setPlacement(placement);

    volSurfaceList(sensorElement)->push_back(prototype.surfaces.at(index));
  }
}

} // namespace

static Ref_t create_detector(Detector& description, xml_h e, SensitiveDetector sensitiveDetector) {

  xml_det_t xDetector            = e;
  const int detectorId           = xDetector.id();
  const std::string detectorName = xDetector.nameStr();
  const bool reflect             = xDetector.reflect(false);
  DetElement detector(detectorName, detectorId);
  Assembly detectorAssembly(detectorName);
  detectorAssembly.setVisAttributes(description.invisible());
  sensitiveDetector.setType("tracker");
  dd4hep::xml::setDetectorTypeFlag(xDetector, detector);

  auto& detectorParameters =
      DD4hepDetectorHelper::ensureExtension<dd4hep::rec::VariantParameters>(detector);
  for (xml_coll_t boundaryMaterial(xDetector, _Unicode(boundary_material)); boundaryMaterial;
       ++boundaryMaterial) {
    DD4hepDetectorHelper::xmlToProtoSurfaceMaterial(xml_comp_t(boundaryMaterial),
                                                    detectorParameters, "boundary_material");
  }

  std::map<std::string, ModulePair> modules;
  for (xml_coll_t module(xDetector, _U(module)); module; ++module) {
    xml_comp_t xModule = module;
    ModulePair pair;
    pair.normal  = buildPrototype(description, sensitiveDetector, xModule, "", false);
    pair.overlap = buildPrototype(description, sensitiveDetector, xModule, "_overlap", true);
    modules.emplace(xModule.nameStr(), std::move(pair));
  }

  for (xml_coll_t layer(xDetector, _U(layer)); layer; ++layer) {
    xml_comp_t xLayer               = layer;
    const int layerId               = xLayer.id();
    xml_comp_t xEnvelope            = xLayer.child(_U(envelope));
    const std::string baseLayerName = detectorName + "_layer" + std::to_string(layerId);
    const double layerLength        = xEnvelope.length();
    const double layerStart         = xEnvelope.zstart();
    const double layerCenter        = layerStart + layerLength / 2.0;

    Volume layerVolume(baseLayerName, Tube(xEnvelope.rmin(), xEnvelope.rmax(), layerLength / 2.0),
                       description.material("Air"));
    layerVolume.setVisAttributes(description.visAttributes(xEnvelope.visStr()));
    const Transform3D layerTransform =
        reflect ? Transform3D(RotationZYX(0.0, -M_PI, 0.0), Position(0.0, 0.0, -layerCenter))
                : Transform3D(Rotation3D(), Position(0.0, 0.0, layerCenter));
    PlacedVolume layerPlacement = detectorAssembly.placeVolume(layerVolume, layerTransform);
    layerPlacement.addPhysVolID("layer", layerId);
    DetElement layerElement(detector, baseLayerName + (reflect ? "_N" : "_P"), layerId);
    layerElement.setPlacement(layerPlacement);

    auto& layerParameters =
        DD4hepDetectorHelper::ensureExtension<dd4hep::rec::VariantParameters>(layerElement);
    for (xml_coll_t material(xLayer, _Unicode(layer_material)); material; ++material) {
      DD4hepDetectorHelper::xmlToProtoSurfaceMaterial(xml_comp_t(material), layerParameters,
                                                      "layer_material");
    }

    int nextModuleId = 1;
    for (xml_coll_t ring(xLayer, _U(ring)); ring; ++ring) {
      xml_comp_t xRing             = ring;
      const ModulePair& modulePair = modules.at(xRing.moduleStr());
      const int modulesPerQuadrant = xRing.attr<int>(_Unicode(modules_per_quadrant));
      const double radius          = xRing.r();
      const double phi0            = xRing.phi0(0.0);
      const double localZ          = xRing.zstart();
      const double dz              = xRing.dz(0.0);
      const double quadrantStagger = getAttrOrDefault(xRing, _Unicode(dz_offset), 0.0);
      const double phiStep         = (M_PI / 2.0) / modulesPerQuadrant;

      for (int quadrantId = 0; quadrantId < 4; ++quadrantId) {
        const std::string quadrantName = "quadrant" + std::to_string(quadrantId);
        Assembly quadrantAssembly(baseLayerName + "_" + quadrantName);
        PlacedVolume quadrantPlacement = layerVolume.placeVolume(quadrantAssembly);
        DetElement quadrantElement(layerElement, quadrantName, quadrantId);
        quadrantElement.setPlacement(quadrantPlacement);

        const double zOffset = (quadrantId % 2 == 0) ? quadrantStagger : 0.0;
        for (int inQuadrant = 0; inQuadrant < modulesPerQuadrant; ++inQuadrant) {
          const double phi = phi0 + (quadrantId * modulesPerQuadrant + inQuadrant) * phiStep;
          const double x   = -radius * std::cos(phi);
          const double y   = -radius * std::sin(phi);
          const double z   = reflect ? -localZ - dz - zOffset : localZ + dz + zOffset;
          PlacedVolume modulePlacement = quadrantAssembly.placeVolume(
              modulePair.normal.volume,
              Transform3D(RotationZYX(0.0, -M_PI / 2.0 - phi, -M_PI / 2.0), Position(x, y, z)));
          modulePlacement.addPhysVolID("module", nextModuleId);
          DetElement moduleElement(quadrantElement, "module" + std::to_string(nextModuleId),
                                   nextModuleId);
          moduleElement.setPlacement(modulePlacement);
          addSensitiveDetElements(moduleElement, modulePair.normal, nextModuleId++);
        }

        // Exactly one explicitly declared overlap belongs to each quadrant boundary.
        for (xml_coll_t overlap(xRing, _U(overlap)); overlap; ++overlap) {
          xml_comp_t xOverlap = overlap;
          if (xOverlap.attr<int>(_Unicode(quadrant)) != quadrantId) {
            continue;
          }
          const double z                = reflect ? -localZ - dz - zOffset : localZ + dz + zOffset;
          PlacedVolume overlapPlacement = quadrantAssembly.placeVolume(
              modulePair.overlap.volume,
              Transform3D(RotationZYX(0.0, xOverlap.attr<double>(_Unicode(rotation)), -M_PI / 2.0),
                          Position(xOverlap.x(), xOverlap.y(), z)));
          overlapPlacement.addPhysVolID("module", nextModuleId);
          DetElement overlapElement(quadrantElement, "overlap" + std::to_string(nextModuleId),
                                    nextModuleId);
          overlapElement.setPlacement(overlapPlacement);
          addSensitiveDetElements(overlapElement, modulePair.overlap, nextModuleId++);
        }
      }
    }
  }

  PlacedVolume detectorPlacement = description.pickMotherVolume(detector).placeVolume(
      detectorAssembly, Position(0.0, 0.0, reflect ? -1.0e-9 : 1.0e-9));
  detectorPlacement.addPhysVolID("system", detectorId);
  detector.setPlacement(detectorPlacement);

  return detector;
}

DECLARE_DETELEMENT(epic_MPGDQuadrantEndcapTracker, create_detector)
