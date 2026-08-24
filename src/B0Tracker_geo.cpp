// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 - 2026 Whitney Armstrong, Igor Korover, Tom Bleher

/** \addtogroup Trackers Trackers
 * \brief Type: **B0Tracker**.
 * \author W. Armstrong, I. Korover, T. Bleher
 *
 * \ingroup trackers
 *
 * Hierarchy layer -> module -> sensor, with these invariants:
 *   - One compact <layer> is one ACTS disc layer, front or back
 *   - Layer id = 2*(station-1) + (1=back | 2=front), so the cellID layer field
 *     is monotonic in z and separates front from back
 *   - Module ids restart at 1 per layer, so cellIDs do not depend on the
 *     order of <layer> blocks in the compact file
 *   - TrackingUnit Assembly is built once and reused via placeVolume for
 *     every (layer, module position)
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

#include <cmath>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

using namespace dd4hep;
using namespace dd4hep::rec;

namespace {
struct ModuleComponentDef {
  std::string name;
  std::string material;
  std::string vis;
  double dx{0.0};
  double dy{0.0};
  double dz{0.0};
  double px{0.0};
  double py{0.0};
  double pz{0.0};
  bool sensitive{false};
  double inner{0.0};
  double outer{0.0};
};

struct SupportComponentDef {
  Volume volume;
  Position position;
};

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
  dd4hep::xml::setDetectorTypeFlag(x_det, sdet);
  auto& detParams = DD4hepDetectorHelper::ensureExtension<dd4hep::rec::VariantParameters>(sdet);

  // Add the volume boundary material if configured
  for (xml_coll_t bmat(x_det, _Unicode(boundary_material)); bmat; ++bmat) {
    xml_comp_t x_boundary_material = bmat;
    DD4hepDetectorHelper::xmlToProtoSurfaceMaterial(x_boundary_material, detParams,
                                                    "boundary_material");
  }

  assembly.setVisAttributes("InvisibleWithDaughters");
  sens.setType("tracker");

  // Read the shared <module> definitions declared under <detector>
  auto findModule = [&x_det](const std::string& name) {
    xml_h found;
    for (xml_coll_t it(x_det, _U(module)); it; ++it) {
      xml_comp_t xm = it;
      if (xm.nameStr() == name) {
        found = xm;
        break;
      }
    }
    return found;
  };

  xml_h trackingUnit = findModule("TrackingUnit");
  if (!trackingUnit.ptr()) {
    throw std::runtime_error(
        "B0Tracker: <module name=\"TrackingUnit\"> not found under <detector>");
  }

  std::vector<ModuleComponentDef> moduleComponents;

  double zMin          = +std::numeric_limits<double>::infinity();
  double zMax          = -std::numeric_limits<double>::infinity();
  double sensitiveZMin = +std::numeric_limits<double>::infinity();
  double sensitiveZMax = -std::numeric_limits<double>::infinity();

  for (xml_coll_t comp(trackingUnit, _U(module_component)); comp; ++comp) {
    xml_comp_t xc   = comp;
    xml_dim_t x_box = xc.child(_U(box));

    const double dx = x_box.x();
    const double dy = x_box.y();
    const double dz = x_box.z();

    xml::Component cpos = xc.position();
    const double px     = cpos.x();
    const double py     = cpos.y();
    const double pz     = cpos.z();

    zMin = std::min(zMin, pz - dz / 2.0);
    zMax = std::max(zMax, pz + dz / 2.0);

    ModuleComponentDef cdef;
    cdef.name      = xc.nameStr();
    cdef.material  = xc.attr<std::string>(_Unicode(material));
    cdef.vis       = getAttrOrDefault<std::string>(xc, _Unicode(vis), "");
    cdef.dx        = dx;
    cdef.dy        = dy;
    cdef.dz        = dz;
    cdef.px        = px;
    cdef.py        = py;
    cdef.pz        = pz;
    cdef.sensitive = xc.isSensitive();

    if (cdef.sensitive) {
      sensitiveZMin = std::min(sensitiveZMin, pz - dz / 2.0);
      sensitiveZMax = std::max(sensitiveZMax, pz + dz / 2.0);
    }

    moduleComponents.push_back(cdef);
  }

  const double sensitiveCenterZ = 0.5 * (sensitiveZMin + sensitiveZMax);

  // compute the inner and outer thicknesses that need to be assigned to the tracking
  // surface; for B0 they span the full TrackingUnit stack for every sensor
  for (auto& cdef : moduleComponents) {
    if (!cdef.sensitive) {
      continue;
    }
    cdef.inner = cdef.pz - zMin;
    cdef.outer = zMax - cdef.pz;
  }

  // Collect the module's sensitive PlacedVolumes and ACTS VolPlane surfaces,
  // reused at every module position
  Assembly moduleAsm("TrackingUnit");
  std::vector<PlacedVolume> moduleSensVols;
  std::vector<VolPlane> moduleSensSurfs;

  moduleAsm.setVisAttributes(
      description, getAttrOrDefault<std::string>(xml_comp_t(trackingUnit), _Unicode(vis), ""));

  {
    int sensorIndex = 1;
    for (const auto& cdef : moduleComponents) {
      Material mat = description.material(cdef.material);
      Box shape(cdef.dx / 2.0, cdef.dy / 2.0, cdef.dz / 2.0);
      Volume c_vol(cdef.name, shape, mat);

      c_vol.setVisAttributes(description, cdef.vis);
      if (cdef.sensitive) {
        c_vol.setSensitiveDetector(sens);
      }

      PlacedVolume comp_pv = moduleAsm.placeVolume(c_vol, Position(cdef.px, cdef.py, cdef.pz));

      if (cdef.sensitive) {
        comp_pv.addPhysVolID("sensor", sensorIndex);
        moduleSensVols.push_back(comp_pv);

        // Measurement plane attached to the sensitive volume
        Vector3D u(-1.0, 0.0, 0.0);
        Vector3D v(0.0, -1.0, 0.0);
        Vector3D n(0.0, 0.0, 1.0);
        SurfaceType type(SurfaceType::Sensitive);
        VolPlane surf(c_vol, type, cdef.inner, cdef.outer, u, v, n);
        moduleSensSurfs.push_back(surf);

        ++sensorIndex;
      }
    }
  }

  // Guard cellID field capacities against the readout definition itself,
  // so a readout or module change cannot silently overflow a bit field
  const dd4hep::IDDescriptor idSpec = sens.readout().idSpec();
  const auto maxSensorID            = idSpec.field("sensor")->maxValue();
  if (static_cast<long long>(moduleSensVols.size()) > static_cast<long long>(maxSensorID)) {
    throw std::runtime_error("B0Tracker: TrackingUnit has " +
                             std::to_string(moduleSensVols.size()) +
                             " sensitive components; the 'sensor' readout field holds at most " +
                             std::to_string(maxSensorID));
  }

  // Support plates: every <module> carrying an <extruded_polygon> outline, built
  // once here and placed per station by name from the <layer> <component ref=...>
  std::map<std::string, std::vector<SupportComponentDef>> supportComponents;
  for (xml_coll_t it(x_det, _U(module)); it; ++it) {
    xml_comp_t x_support_module = it;

    // One CAD-traced outline per plate, shared by its skins; each
    // <module_component> contributes its own thickness, material and offset
    xml_comp_t x_polygon = x_support_module.child(dd4hep::xml::Strng_t("extruded_polygon"), false);
    if (!x_polygon.ptr()) {
      continue;
    }
    const std::string supportName = x_support_module.nameStr();
    std::vector<double> xVertices;
    std::vector<double> yVertices;
    for (xml_coll_t point(x_polygon, _U(point)); point; ++point) {
      xml_comp_t x_point = point;
      xVertices.push_back(x_point.x());
      yVertices.push_back(x_point.y());
    }
    if (xVertices.size() < 3) {
      throw std::runtime_error("B0Tracker: " + supportName +
                               " has fewer than three polygon points");
    }

    auto& components = supportComponents[supportName];
    for (xml_coll_t comp(x_support_module, _U(module_component)); comp; ++comp) {
      xml_comp_t x_support_component = comp;
      const std::string supportVis   = getAttrOrDefault<std::string>(
          x_support_component, _Unicode(vis),
          getAttrOrDefault<std::string>(x_support_module, _Unicode(vis), ""));
      Material supportMat    = description.material(x_support_component.materialStr());
      const double thickness = x_support_component.attr<double>(_Unicode(thickness));
      ExtrudedPolygon supportSolid(xVertices, yVertices, {-thickness / 2., thickness / 2.},
                                   {0., 0.}, {0., 0.}, {1., 1.});
      Volume supportVol(supportName + "_" + x_support_component.nameStr(), supportSolid,
                        supportMat);
      if (!supportVis.empty()) {
        supportVol.setVisAttributes(description.visAttributes(supportVis));
      }
      Position supportPosition;
      if (x_support_component.hasChild(_U(position))) {
        xml_dim_t x_support_position = x_support_component.child(_U(position));
        supportPosition =
            Position(x_support_position.x(), x_support_position.y(), x_support_position.z());
      }
      components.push_back({supportVol, supportPosition});
    }
  }

  // Distance from the support mid-plane to the tracking unit center
  const double moduleOffset = xml_comp_t(trackingUnit).attr<double>(_Unicode(offset_from_support));

  // now build the layers
  for (xml_coll_t layer(x_det, _U(layer)); layer; ++layer) {
    xml_comp_t x_layer = layer;
    const int station  = x_layer.attr<int>(_Unicode(station));

    const std::string side = x_layer.attr<std::string>(_Unicode(side));
    if (side != "front" && side != "back") {
      throw std::runtime_error("B0Tracker: station " + std::to_string(station) + " has side=\"" +
                               side + "\"; expected \"front\" or \"back\"");
    }
    const bool isFront = side == "front";

    // cellID layer field, monotonic in z and separating front from back
    const int layerID = 2 * (station - 1) + (isFront ? 2 : 1);

    xml_comp_t x_env = x_layer.child(_U(envelope), false);
    if (!x_env.ptr()) {
      throw std::runtime_error("B0Tracker: layer " + std::to_string(layerID) +
                               " is missing its <envelope> element");
    }

    const double env_rmin_tol = getAttrOrDefault<double>(x_env, _Unicode(rmin_tolerance), 0.0);
    const double env_rmax_tol = getAttrOrDefault<double>(x_env, _Unicode(rmax_tolerance), 0.0);
    const double env_zmin_tol = getAttrOrDefault<double>(x_env, _Unicode(zmin_tolerance), 0.0);
    const double env_zmax_tol = getAttrOrDefault<double>(x_env, _Unicode(zmax_tolerance), 0.0);
    if (env_zmin_tol <= 0.0 || env_zmax_tol <= 0.0) {
      printout(WARNING, det_name,
               "Layer %d has a non-positive envelope z tolerance; the ACTS approach surfaces "
               "will collapse onto the layer surface",
               layerID);
    }

    std::string env_vis;
    if (x_env.hasAttr(_Unicode(vis))) {
      env_vis = x_env.attr<std::string>(_Unicode(vis));
    }

    // Layer origin from XML: the station origin, shared by front and
    // back layers of a station. Mechanical support stays at this origin
    xml_dim_t lp        = x_layer.child(_U(position));
    const double layerX = lp.x();
    const double layerY = lp.y();
    const double layerZ = lp.z();

    // The ACTS measurement layers are the sensor stacks, so material maps
    // project this support material onto the adjacent layer surfaces
    for (xml_coll_t comp(x_layer, _U(component)); comp; ++comp) {
      xml_comp_t xc         = comp;
      const std::string ref = xc.attr<std::string>(_Unicode(ref));
      const auto supportIt  = supportComponents.find(ref);
      if (supportIt == supportComponents.end()) {
        throw std::runtime_error("B0Tracker: layer " + std::to_string(layerID) +
                                 " has unsupported <component ref=\"" + ref + "\">");
      }
      xml_dim_t sp = xc.child(_U(position));
      for (const auto& support : supportIt->second) {
        assembly.placeVolume(support.volume, Position(layerX + sp.x() + support.position.x(),
                                                      layerY + sp.y() + support.position.y(),
                                                      layerZ + sp.z() + support.position.z()));
      }
    }

    xml_comp_t mpos = x_layer.child(_Unicode(module_positions), false);
    if (!mpos.ptr()) {
      throw std::runtime_error("B0Tracker: layer " + std::to_string(layerID) +
                               " is missing its <module_positions> element");
    }

    // Layer assembly, centered on the sensor stack so the binned ACTS
    // measurement surfaces sit at the physical sensor planes
    const double sideSign   = isFront ? 1.0 : -1.0;
    const double sideLayerZ = sideSign * (moduleOffset + sensitiveCenterZ);

    const std::string sideLayerName = det_name + "_layer" + std::to_string(station) + "_" + side;
    Assembly sideVol(sideLayerName);
    if (!env_vis.empty()) {
      sideVol.setVisAttributes(description.visAttributes(env_vis));
    }

    PlacedVolume sidePV =
        assembly.placeVolume(sideVol, Position(layerX, layerY, layerZ + sideLayerZ));
    sidePV.addPhysVolID("layer", layerID);

    DetElement sideDE(sdet, sideLayerName + "_P", layerID);
    sideDE.setPlacement(sidePV);

    // Place the shared TrackingUnit Assembly at each <module_positions> entry
    int moduleID = 1;

    for (xml_coll_t mp(mpos, _U(module)); mp; ++mp, ++moduleID) {
      xml_comp_t xm = mp;

      const double modX    = xm.attr<double>(_Unicode(posX));
      const double modY    = xm.attr<double>(_Unicode(posY));
      const double modRotZ = xm.attr<double>(_Unicode(rotZ));
      // Offset each module so its sensitive center lands on the layer plane;
      // back modules are flipped about x so both sides face the support
      const double modZ = -sideSign * sensitiveCenterZ;
      RotationZYX rotLocal(modRotZ, 0.0, isFront ? 0.0 : M_PI);
      Transform3D modTr(rotLocal, Position(modX, modY, modZ));

      PlacedVolume mod_pv = sideVol.placeVolume(moduleAsm, modTr);
      mod_pv.addPhysVolID("module", moduleID);

      std::string m_base = _toString(layerID, "layer%d") + _toString(moduleID, "_module%d");
      DetElement modDE(sideDE, m_base, moduleID);
      modDE.setPlacement(mod_pv);

      for (size_t ic = 0; ic < moduleSensVols.size(); ++ic) {
        PlacedVolume sens_pv = moduleSensVols[ic];
        DetElement comp_de(modDE, std::string("de_") + sens_pv.volume().name(), moduleID);
        comp_de.setPlacement(sens_pv);

        auto& comp_de_params =
            DD4hepDetectorHelper::ensureExtension<dd4hep::rec::VariantParameters>(comp_de);
        comp_de_params.set<std::string>("axis_definitions", "XYZ");

        volSurfaceList(comp_de)->push_back(moduleSensSurfs[ic]);
      }
    }
    sideVol->GetShape()->ComputeBBox();

    auto& sideParams =
        DD4hepDetectorHelper::ensureExtension<dd4hep::rec::VariantParameters>(sideDE);

    sideParams.set<double>("envelope_r_min", env_rmin_tol / dd4hep::mm);
    sideParams.set<double>("envelope_r_max", env_rmax_tol / dd4hep::mm);
    sideParams.set<double>("envelope_z_min", env_zmin_tol / dd4hep::mm);
    sideParams.set<double>("envelope_z_max", env_zmax_tol / dd4hep::mm);

    for (xml_coll_t lmat(x_layer, _Unicode(layer_material)); lmat; ++lmat) {
      xml_comp_t x_layer_material = lmat;
      DD4hepDetectorHelper::xmlToProtoSurfaceMaterial(x_layer_material, sideParams,
                                                      "layer_material");
    }

    printout(DEBUG, det_name,
             "Layer %d (station %d %s) z=%8.3f mm "
             "tol(rmin,rmax,zmin,zmax)=(%6.3f,%6.3f,%6.3f,%6.3f) mm",
             layerID, station, side.c_str(), (layerZ + sideLayerZ) / mm, env_rmin_tol / mm,
             env_rmax_tol / mm, env_zmin_tol / mm, env_zmax_tol / mm);
  }

  PlacedVolume pv = motherVol.placeVolume(assembly, posAndRot);
  pv.addPhysVolID("system", det_id);
  sdet.setPlacement(pv);

  return sdet;
}

//@}
// clang-format off
DECLARE_DETELEMENT(ip6_B0Tracker, create_B0Tracker)
