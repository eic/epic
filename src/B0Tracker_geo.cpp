// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 Tom Bleher, Igor Korover

/*! B0 Tracker.
 *
 * @author Tom Bleher, Igor Korover
 * Compact elements read under <detector>:
 *   - <module name="TrackingUnit">: <module_component> boxes (<box>, <position>,
 *     material, sensitive) forming one double-sided sensor stave
 *   - <support_stack>: <slice name material thickness> layers, listed upstream to
 *     downstream, extruded along each <support_plate name> outline of <point x y>
 *   - <envelope rmin/rmax/zmin/zmax_tolerance vis> and <layer_material>: ACTS layer
 *     envelopes and material binning, shared by every face
 *   - <module_layout name>: a back and a front <face side>, each listing its
 *     <module x y [rotZ]> placements
 *   - <station id>: <position>, <support ref>, and <layout ref>
 *   - optional <acts_guard gap rmin rmax>: empty ACTS layers outside the first
 *     and last faces
 *
 * Hierarchy station -> face -> module -> sensor, with these invariants:
 *   - Each <face> of a station's layout is one ACTS disc layer, the back or
 *     front sensor stack of that <station>
 *   - Layer id = 2*(station-1) + (1=back | 2=front), so the cellID layer field
 *     is monotonic in z and separates front from back
 *   - Module ids restart at 1 per face, so cellIDs do not depend on the
 *     order of <station> blocks in the compact file
 *   - Layer, module and sensor ids are checked against the readout fields
 *   - TrackingUnit Assembly is built once and reused via placeVolume for
 *     every (face, module position)
 *
 * @{
 */

#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/IDDescriptor.h"
#include "DD4hep/Printout.h"
#include "DD4hep/Readout.h"
#include "DD4hep/Shapes.h"
#include "DD4hepDetectorHelper.h"
#include "DDRec/DetectorData.h"
#include "DDRec/Surface.h"
#include "XML/Utilities.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

using namespace dd4hep;
using namespace dd4hep::rec;

namespace {

// DetElement ids of the guard layers, above the face layer ids
constexpr int kGuardIdOffset = 100;

struct SupportSlice {
  Volume volume;
  Position position;
};

struct SupportPlates {
  std::map<std::string, std::vector<SupportSlice>> slices; // by <support_plate name>
  double thickness = 0.0;                                  // of the whole stack
};

// The shared TrackingUnit, placed at every module position
struct TrackingUnit {
  Assembly assembly;
  std::vector<PlacedVolume> sensors; // sensitive placements, sensor id = index + 1
  std::vector<VolPlane> surfaces;    // ACTS measurement plane of each sensor
  double zMin = +std::numeric_limits<double>::infinity(); // z extent of the stack
  double zMax = -std::numeric_limits<double>::infinity();
};

// Throw if a volume id does not fit its readout field; DD4hep would otherwise
// wrap it silently into the cellID. The limits come from the field width, as
// BitFieldElement::maxValue() is one too high for unsigned fields (DD4hep 1.38)
void checkReadoutField(const std::string& det_name, SensitiveDetector sens,
                       const std::string& field, int value) {
  const auto* element = sens.readout().idSpec().field(field);
  const int maxValue  = element->isSigned() ? element->maxValue() : (1 << element->width()) - 1;
  if (value < element->minValue() || value > maxValue) {
    throw std::runtime_error(det_name + ": " + field + " id " + std::to_string(value) +
                             " does not fit the readout field (" +
                             std::to_string(element->minValue()) + ".." + std::to_string(maxValue) +
                             ")");
  }
}

// The <tag name="..."> element under <detector>, or an invalid handle
xml_h findNamed(xml_det_t x_det, const xml::Strng_t& tag, const std::string& name) {
  for (xml_coll_t it(x_det, tag); it; ++it) {
    xml_comp_t xm = it;
    if (xm.nameStr() == name) {
      return xm;
    }
  }
  return xml_h();
}

TrackingUnit buildTrackingUnit(Detector& description, SensitiveDetector sens,
                               const std::string& det_name, xml_comp_t x_unit) {
  TrackingUnit unit;
  unit.assembly = Assembly(det_name + "_TrackingUnit");
  unit.assembly.setVisAttributes(description,
                                 getAttrOrDefault<std::string>(x_unit, _Unicode(vis), ""));

  // Extent of the stack in z, which sets the thickness of the ACTS measurement
  // surfaces and the mounting on the support plates
  for (xml_coll_t comp(x_unit, _U(module_component)); comp; ++comp) {
    xml_comp_t xc      = comp;
    const double pz    = xc.position().z();
    const double halfZ = xml_dim_t(xc.child(_U(box))).z() / 2.0;
    unit.zMin          = std::min(unit.zMin, pz - halfZ);
    unit.zMax          = std::max(unit.zMax, pz + halfZ);
  }

  for (xml_coll_t comp(x_unit, _U(module_component)); comp; ++comp) {
    xml_comp_t xc   = comp;
    xml_dim_t x_box = xc.child(_U(box));
    xml_dim_t x_pos = xc.position();

    Volume c_vol(det_name + "_" + xc.nameStr(),
                 Box(x_box.x() / 2.0, x_box.y() / 2.0, x_box.z() / 2.0),
                 description.material(xc.materialStr()));
    c_vol.setVisAttributes(description, getAttrOrDefault<std::string>(xc, _Unicode(vis), ""));
    PlacedVolume comp_pv =
        unit.assembly.placeVolume(c_vol, Position(x_pos.x(), x_pos.y(), x_pos.z()));

    if (xc.isSensitive()) {
      c_vol.setSensitiveDetector(sens);
      comp_pv.addPhysVolID("sensor", static_cast<int>(unit.sensors.size()) + 1);
      unit.sensors.push_back(comp_pv);
      // Measurement plane; its inner and outer thicknesses span the whole stack
      unit.surfaces.emplace_back(c_vol, SurfaceType(SurfaceType::Sensitive), x_pos.z() - unit.zMin,
                                 unit.zMax - x_pos.z(), Vector3D(-1.0, 0.0, 0.0),
                                 Vector3D(0.0, -1.0, 0.0), Vector3D(0.0, 0.0, 1.0));
    }
  }
  checkReadoutField(det_name, sens, "sensor", static_cast<int>(unit.sensors.size()));
  return unit;
}

// The <support_stack> slices, stacked upstream to downstream around the plate
// mid-plane and extruded along each <support_plate> outline
SupportPlates buildSupportPlates(Detector& description, const std::string& det_name,
                                 xml_det_t x_det) {
  SupportPlates plates;
  xml_comp_t x_stack         = x_det.child(_Unicode(support_stack));
  const std::string stackVis = getAttrOrDefault<std::string>(x_stack, _Unicode(vis), "");
  for (xml_coll_t sl(x_stack, _U(slice)); sl; ++sl) {
    plates.thickness += xml_comp_t(sl).thickness();
  }

  for (xml_coll_t pl(x_det, _Unicode(support_plate)); pl; ++pl) {
    xml_comp_t x_plate          = pl;
    const std::string plateName = x_plate.nameStr();
    std::vector<double> xVertices;
    std::vector<double> yVertices;
    for (xml_coll_t point(x_plate, _U(point)); point; ++point) {
      xml_comp_t x_point = point;
      xVertices.push_back(x_point.x());
      yVertices.push_back(x_point.y());
    }
    if (xVertices.size() < 3) {
      throw std::runtime_error(det_name + ": " + plateName + " has fewer than three points");
    }

    auto& slices = plates.slices[plateName];
    double zLow  = -plates.thickness / 2.0;
    for (xml_coll_t sl(x_stack, _U(slice)); sl; ++sl) {
      xml_comp_t x_slice     = sl;
      const double thickness = x_slice.thickness();
      Volume sliceVol(plateName + "_" + x_slice.nameStr(),
                      ExtrudedPolygon(xVertices, yVertices, {-thickness / 2., thickness / 2.},
                                      {0., 0.}, {0., 0.}, {1., 1.}),
                      description.material(x_slice.materialStr()));
      sliceVol.setVisAttributes(description,
                                getAttrOrDefault<std::string>(x_slice, _Unicode(vis), stackVis));
      slices.push_back({sliceVol, Position(0.0, 0.0, zLow + thickness / 2.0)});
      zLow += thickness;
    }
  }
  return plates;
}

// Empty ACTS layers `gap` outside the first and last faces, which widen the B0
// tracking volume to the outer sensors of the faces tilted by the crossing angle.
// ACTS takes only rmin and rmax of an empty TGeoTubeSeg layer and builds a full
// disc; the half disc only shapes the Geant4 volume and the event display
void addActsGuardLayers(Detector& description, DetElement sdet, Assembly assembly,
                        xml_comp_t x_guard, const Position& firstFacePos,
                        const Position& lastFacePos) {
  const std::string det_name = sdet.name();
  if (!std::isfinite(firstFacePos.z())) {
    throw std::runtime_error(det_name + ": <acts_guard> needs at least one <face>");
  }
  const Position gap(0.0, 0.0, x_guard.attr<double>(_Unicode(gap)));
  const double rmin = x_guard.rmin();
  const double rmax = x_guard.rmax();
  int guardID       = kGuardIdOffset;
  for (const auto& [name, guardPos] :
       {std::pair{"upstream", firstFacePos - gap}, std::pair{"downstream", lastFacePos + gap}}) {
    const std::string guardName = det_name + "_guard_" + name;
    // Half disc on the side away from the electron beam pipe
    Volume guardVol(guardName, Tube(rmin, rmax, 0.5 * dd4hep::um, 0.5 * M_PI, 1.5 * M_PI),
                    description.vacuum());
    guardVol.setVisAttributes(description.invisible());
    PlacedVolume guardPV = assembly.placeVolume(guardVol, guardPos);
    DetElement guardDE(sdet, guardName + "_P", guardID++);
    guardDE.setPlacement(guardPV);
    DD4hepDetectorHelper::ensureExtension<VariantParameters>(guardDE);
  }
}

} // namespace

static Ref_t create_B0Tracker(Detector& description, xml_h e, SensitiveDetector sens) {
  xml_det_t x_det            = e;
  const int det_id           = x_det.id();
  const std::string det_name = x_det.nameStr();

  DetElement sdet(det_name, det_id);
  Assembly assembly(det_name);

  Volume motherVol   = description.pickMotherVolume(sdet);
  xml::Component pos = x_det.position();
  xml::Component rot = x_det.rotation();
  Transform3D posAndRot(RotationZYX(rot.z(), rot.y(), rot.x()),
                        Position(pos.x(), pos.y(), pos.z()));

  // Set detector type flag
  xml::setDetectorTypeFlag(x_det, sdet);
  auto& detParams = DD4hepDetectorHelper::ensureExtension<VariantParameters>(sdet);

  // Add the volume boundary material if configured
  for (xml_coll_t bmat(x_det, _Unicode(boundary_material)); bmat; ++bmat) {
    xml_comp_t x_boundary_material = bmat;
    DD4hepDetectorHelper::xmlToProtoSurfaceMaterial(x_boundary_material, detParams,
                                                    "boundary_material");
  }

  assembly.setAttributes(description, x_det.regionStr(), x_det.limitsStr(), x_det.visStr());
  sens.setType("tracker");

  xml_comp_t x_unit = findNamed(x_det, _U(module), "TrackingUnit");
  if (!x_unit.ptr()) {
    throw std::runtime_error(det_name +
                             ": <module name=\"TrackingUnit\"> not found under <detector>");
  }
  const TrackingUnit unit    = buildTrackingUnit(description, sens, det_name, x_unit);
  const SupportPlates plates = buildSupportPlates(description, det_name, x_det);

  // Distance from the support mid-plane to the tracking unit center, with the
  // unit's lowest-z face flush on the support
  const double moduleOffset = plates.thickness / 2.0 - unit.zMin;

  // ACTS layer settings, shared by every face
  xml_comp_t x_env = x_det.child(_U(envelope), false);
  if (!x_env.ptr()) {
    throw std::runtime_error(det_name + ": <envelope> not found under <detector>");
  }
  const double env_rmin_tol = getAttrOrDefault<double>(x_env, _Unicode(rmin_tolerance), 0.0);
  const double env_rmax_tol = getAttrOrDefault<double>(x_env, _Unicode(rmax_tolerance), 0.0);
  const double env_zmin_tol = getAttrOrDefault<double>(x_env, _Unicode(zmin_tolerance), 0.0);
  const double env_zmax_tol = getAttrOrDefault<double>(x_env, _Unicode(zmax_tolerance), 0.0);
  if (env_zmin_tol <= 0.0 || env_zmax_tol <= 0.0) {
    printout(WARNING, det_name,
             "Non-positive envelope z tolerance; the ACTS approach surfaces will collapse onto "
             "the layer surfaces");
  }
  const std::string env_vis = getAttrOrDefault<std::string>(x_env, _Unicode(vis), "");

  // First and last tracking faces, for the ACTS guard layers
  Position firstFacePos(0.0, 0.0, +std::numeric_limits<double>::infinity());
  Position lastFacePos(0.0, 0.0, -std::numeric_limits<double>::infinity());

  for (xml_coll_t st(x_det, _Unicode(station)); st; ++st) {
    xml_comp_t x_station = st;
    const int station    = x_station.id();
    if (station < 1) {
      throw std::runtime_error(det_name + ": station id " + std::to_string(station) +
                               " must be positive");
    }

    // Station origin, shared by its support plate and both faces
    xml_dim_t x_station_pos = x_station.position();
    const Position stationPos(x_station_pos.x(), x_station_pos.y(), x_station_pos.z());

    // The ACTS measurement layers are the faces, so material maps project the
    // support material onto the adjacent layer surfaces
    for (xml_coll_t sup(x_station, _Unicode(support)); sup; ++sup) {
      const std::string ref = xml_comp_t(sup).attr<std::string>(_Unicode(ref));
      const auto plateIt    = plates.slices.find(ref);
      if (plateIt == plates.slices.end()) {
        throw std::runtime_error(det_name + ": station " + std::to_string(station) +
                                 " has unknown <support ref=\"" + ref + "\">");
      }
      for (const auto& slice : plateIt->second) {
        assembly.placeVolume(slice.volume, stationPos + slice.position);
      }
    }

    const std::string layoutRef =
        xml_comp_t(x_station.child(_Unicode(layout))).attr<std::string>(_Unicode(ref));
    const xml_h x_layout = findNamed(x_det, _Unicode(module_layout), layoutRef);
    if (!x_layout.ptr()) {
      throw std::runtime_error(det_name + ": station " + std::to_string(station) +
                               " has unknown <layout ref=\"" + layoutRef + "\">");
    }

    for (xml_coll_t fc(x_layout, _Unicode(face)); fc; ++fc) {
      xml_comp_t x_face      = fc;
      const std::string side = x_face.attr<std::string>(_Unicode(side));
      if (side != "front" && side != "back") {
        throw std::runtime_error(det_name + ": station " + std::to_string(station) +
                                 " has side=\"" + side + "\"; expected \"front\" or \"back\"");
      }
      const bool isFront = side == "front";

      // cellID layer field, monotonic in z and separating front from back
      const int layerID = 2 * (station - 1) + (isFront ? 2 : 1);
      checkReadoutField(det_name, sens, "layer", layerID);
      checkReadoutField(det_name, sens, "module",
                        static_cast<int>(xml_coll_t(x_face, _U(module)).size()));

      // Face origin at the TrackingUnit center, midway between its sensor planes:
      // ACTS centres the disc layer on this origin with the full thickness of its
      // surfaces, so an off-centre origin would shift the layer off its sensors
      const Position facePos(stationPos.x(), stationPos.y(),
                             stationPos.z() + (isFront ? moduleOffset : -moduleOffset));

      const std::string faceName = det_name + "_station" + std::to_string(station) + "_" + side;
      Assembly faceVol(faceName);
      if (!env_vis.empty()) {
        faceVol.setVisAttributes(description.visAttributes(env_vis));
      }
      PlacedVolume facePV = assembly.placeVolume(faceVol, facePos);
      facePV.addPhysVolID("layer", layerID);
      if (facePos.z() < firstFacePos.z()) {
        firstFacePos = facePos;
      }
      if (facePos.z() > lastFacePos.z()) {
        lastFacePos = facePos;
      }

      DetElement faceDE(sdet, faceName + "_P", layerID);
      faceDE.setPlacement(facePV);

      // Place the shared TrackingUnit at each <module>; back modules are
      // flipped about x so their lowest-z side also sits on the support
      int moduleID = 1;
      for (xml_coll_t mp(x_face, _U(module)); mp; ++mp, ++moduleID) {
        xml_comp_t xm = mp;
        RotationZYX modRot(getAttrOrDefault<double>(xm, _Unicode(rotZ), 0.0), 0.0,
                           isFront ? 0.0 : M_PI);
        PlacedVolume mod_pv =
            faceVol.placeVolume(unit.assembly, Transform3D(modRot, Position(xm.x(), xm.y(), 0.0)));
        mod_pv.addPhysVolID("module", moduleID);

        DetElement modDE(faceDE, _toString(moduleID, "module%d"), moduleID);
        modDE.setPlacement(mod_pv);

        for (size_t ic = 0; ic < unit.sensors.size(); ++ic) {
          const int sensorID = static_cast<int>(ic) + 1;
          DetElement comp_de(modDE, _toString(sensorID, "sensor%d"), sensorID);
          comp_de.setPlacement(unit.sensors[ic]);

          auto& comp_de_params = DD4hepDetectorHelper::ensureExtension<VariantParameters>(comp_de);
          comp_de_params.set<std::string>("axis_definitions", "XYZ");

          volSurfaceList(comp_de)->push_back(unit.surfaces[ic]);
        }
      }

      auto& faceParams = DD4hepDetectorHelper::ensureExtension<VariantParameters>(faceDE);
      faceParams.set<double>("envelope_r_min", env_rmin_tol / dd4hep::mm);
      faceParams.set<double>("envelope_r_max", env_rmax_tol / dd4hep::mm);
      faceParams.set<double>("envelope_z_min", env_zmin_tol / dd4hep::mm);
      faceParams.set<double>("envelope_z_max", env_zmax_tol / dd4hep::mm);

      for (xml_coll_t lmat(x_det, _Unicode(layer_material)); lmat; ++lmat) {
        xml_comp_t x_layer_material = lmat;
        DD4hepDetectorHelper::xmlToProtoSurfaceMaterial(x_layer_material, faceParams,
                                                        "layer_material");
      }

      printout(DEBUG, det_name, "Layer %d (station %d %s) z=%8.3f mm", layerID, station,
               side.c_str(), facePos.z() / dd4hep::mm);
    }
  }

  if (xml_comp_t x_guard = x_det.child(_Unicode(acts_guard), false); x_guard.ptr()) {
    addActsGuardLayers(description, sdet, assembly, x_guard, firstFacePos, lastFacePos);
  }

  PlacedVolume pv = motherVol.placeVolume(assembly, posAndRot);
  pv.addPhysVolID("system", det_id);
  sdet.setPlacement(pv);

  return sdet;
}

//@}
// clang-format off
DECLARE_DETELEMENT(ip6_B0Tracker, create_B0Tracker)
