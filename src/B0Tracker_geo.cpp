// SPDX-License-Identifier: LGPL-3.0-or-later
// B0 Tracker - ACTS/DDRec-friendly builder
//
// Flattened per-layer placement version:
//   - Keep top detector Assembly and per-layer Assembly
//   - Do NOT build a TrackingUnit Assembly
//   - For each layer and module position, place each module component
//     directly into the layer assembly as a basic solid
//   - Sensitive components get direct sensor DetElements under the layer
//   - Read layer <envelope> and propagate envelope_* parameters like the
//     original working B0 implementation

#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/Printout.h"
#include "DD4hep/Shapes.h"
#include "DDRec/DetectorData.h"
#include "DDRec/Surface.h"
#include "DD4hepDetectorHelper.h"
#include "XML/Utilities.h"

#include <map>
#include <vector>
#include <string>
#include <cmath>
#include <limits>
#include <stdexcept>

using namespace std;
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
}

static Ref_t create_B0Tracker(Detector& description, xml_h e, SensitiveDetector sens) {
  xml_det_t x_det = e;
  const int det_id = x_det.id();
  const std::string det_name = x_det.nameStr();

  DetElement sdet(det_name, det_id);
  Assembly assembly(det_name);

  // Mother volume and top transform
  Volume motherVol = description.pickMotherVolume(sdet);
  xml::Component pos = x_det.position();
  xml::Component rot = x_det.rotation();
  Transform3D posAndRot(RotationZYX(rot.z(), rot.y(), rot.x()),
                        Position(pos.x(), pos.y(), pos.z()));

  PlacedVolume pv;

  // Set detector type flag + VariantParameters extension
  dd4hep::xml::setDetectorTypeFlag(x_det, sdet);
  auto& detParams = DD4hepDetectorHelper::ensureExtension<dd4hep::rec::VariantParameters>(sdet);

  // Optional boundary material configuration
  for (xml_coll_t bmat(x_det, _Unicode(boundary_material)); bmat; ++bmat) {
    xml_comp_t x_boundary_material = bmat;
    DD4hepDetectorHelper::xmlToProtoSurfaceMaterial(x_boundary_material, detParams,
                                                    "boundary_material");
  }

  assembly.setVisAttributes("InvisibleWithDaughters");
  sens.setType("tracker");

  // ------------------------------------------------------------------
  // Read TrackingUnit component definitions ONCE, but do not build an assembly
  // ------------------------------------------------------------------
  xml_h root = x_det.document().root();
  xml_h trackingUnit;
  for (xml_coll_t it(root, _U(module)); it; ++it) {
    xml_comp_t xm = it;
    if (xm.nameStr() == "TrackingUnit") {
      trackingUnit = xm;
      break;
    }
  }
  if (!trackingUnit.ptr()) {
    throw std::runtime_error("FATAL: <module name=\"TrackingUnit\"> not found in compact file");
  }

  std::vector<ModuleComponentDef> moduleComponents;

  double zMin = +std::numeric_limits<double>::infinity();
  double zMax = -std::numeric_limits<double>::infinity();

  for (xml_coll_t comp(trackingUnit, _U(module_component)); comp; ++comp) {
    xml_comp_t xc = comp;
    xml_h x_box = xc.child(_U(box));

    const double dx = x_box.attr<double>(_Unicode(x));
    const double dy = x_box.attr<double>(_Unicode(y));
    const double dz = x_box.attr<double>(_Unicode(z));

    xml::Component cpos = xc.position();
    const double px = cpos.x();
    const double py = cpos.y();
    const double pz = cpos.z();

    zMin = std::min(zMin, pz - dz / 2.0);
    zMax = std::max(zMax, pz + dz / 2.0);

    ModuleComponentDef cdef;
    cdef.name = xc.nameStr();
    cdef.material = xc.attr<std::string>(_Unicode(material));
    cdef.vis = xc.hasAttr(_Unicode(vis)) ? xc.attr<std::string>(_Unicode(vis)) : "";
    cdef.dx = dx;
    cdef.dy = dy;
    cdef.dz = dz;
    cdef.px = px;
    cdef.py = py;
    cdef.pz = pz;
    cdef.sensitive = xc.hasAttr(_Unicode(sensitive)) && xc.attr<bool>(_Unicode(sensitive));

    moduleComponents.push_back(cdef);
  }

  if (!std::isfinite(zMin) || !std::isfinite(zMax) || zMax <= zMin) {
    throw std::runtime_error("FATAL: failed to compute TrackingUnit z-extent (zMin/zMax invalid)");
  }

  // Assign per-sensitive surface thicknesses
  for (auto& cdef : moduleComponents) {
    if (!cdef.sensitive) {
      continue;
    }
    cdef.inner = 0.15 * mm;
    cdef.outer = 0.15 * mm;
  }

  // ------------------------------------------------------------------
  // Support disk volume
  // ------------------------------------------------------------------
  const double rmin  = description.constant<double>("SupportRMin");
  const double rmax  = description.constant<double>("SupportRMax");
  const double thick = description.constant<double>("SupportThickness");
  const double phi0  = description.constant<double>("SupportPhiStart");
  const double dphi  = description.constant<double>("SupportPhiDelta");

  Tube supportSolid(rmin, rmax, thick / 2.0, phi0, phi0 + dphi);
  Material supportMat = description.material("Copper");
  Volume supportVol("B0SupportDiskVol", supportSolid, supportMat);

  const double moduleOffset = description.constant<double>("ModuleOffsetFromSupport");
  const double frontZ = +moduleOffset;
  const double backZ  = -moduleOffset;

  // ------------------------------------------------------------------
  // Layers
  // ------------------------------------------------------------------
  int globalModuleID = 1;

  for (xml_coll_t layer(x_det, _U(layer)); layer; ++layer) {
    xml_comp_t x_layer = layer;
    const int layerID = x_layer.id();

    // --------------------------------------------------------------
    // Read layer envelope like in the working original
    // --------------------------------------------------------------
    xml_comp_t x_env = x_layer.child(_U(envelope), false);

    double env_rmin_tol = 0.0;
    double env_rmax_tol = 0.0;
    double env_zmin_tol = 0.0;
    double env_zmax_tol = 0.0;
    double env_length   = 0.0;
    double env_zstart   = 0.0;
    std::string env_vis;

    if (x_env.ptr()) {
      env_rmin_tol = getAttrOrDefault<double>(x_env, _Unicode(rmin_tolerance), 0.0);
      env_rmax_tol = getAttrOrDefault<double>(x_env, _Unicode(rmax_tolerance), 0.0);
      env_zmin_tol = getAttrOrDefault<double>(x_env, _Unicode(zmin_tolerance), 0.0);
      env_zmax_tol = getAttrOrDefault<double>(x_env, _Unicode(zmax_tolerance), 0.0);
      env_length   = getAttrOrDefault<double>(x_env, _Unicode(length), 0.0);
      env_zstart   = getAttrOrDefault<double>(x_env, _Unicode(zstart), 0.0);

      if (x_env.hasAttr(_Unicode(vis))) {
        env_vis = x_env.attr<std::string>(_Unicode(vis));
      }
    }

    std::string layer_name = det_name + std::string("_layer") + std::to_string(layerID);
    Assembly layer_vol(layer_name);

    if (!env_vis.empty()) {
      layer_vol.setVisAttributes(description.visAttributes(env_vis));
    }

    // Place the layer in the detector assembly using its <position>
    xml_comp_t lp = x_layer.child(_U(position));
    Transform3D layerTr(Rotation3D(),
                        Position(lp.attr<double>(_Unicode(x)),
                                 lp.attr<double>(_Unicode(y)),
                                 lp.attr<double>(_Unicode(z))));

    PlacedVolume layer_pv = assembly.placeVolume(layer_vol, layerTr);
    layer_pv.addPhysVolID("layer", layerID);

    DetElement layerDE(sdet, layer_name + "_P", layerID);
    layerDE.setPlacement(layer_pv);

    auto& layerParams =
        DD4hepDetectorHelper::ensureExtension<dd4hep::rec::VariantParameters>(layerDE);

    for (xml_coll_t lmat(x_layer, _Unicode(layer_material)); lmat; ++lmat) {
      xml_comp_t x_layer_material = lmat;
      DD4hepDetectorHelper::xmlToProtoSurfaceMaterial(x_layer_material, layerParams,
                                                      "layer_material");
    }

    // --------------------------------------------------------------
    // Support disks go directly into the layer
    // --------------------------------------------------------------
    for (xml_coll_t comp(x_layer, _U(component)); comp; ++comp) {
      xml_comp_t xc = comp;
      if (xc.hasAttr(_Unicode(ref)) && xc.attr<std::string>(_Unicode(ref)) == "B0SupportDisk") {
        xml_comp_t sp = xc.child(_U(position));
        Transform3D tr(Rotation3D(),
                       Position(sp.attr<double>(_Unicode(x)),
                                sp.attr<double>(_Unicode(y)),
                                sp.attr<double>(_Unicode(z))));
        PlacedVolume spv = layer_vol.placeVolume(supportVol, tr);
        spv.addPhysVolID("layer", layerID);
      }
    }

    // --------------------------------------------------------------
    // Place every module component directly while looping over modules
    // --------------------------------------------------------------
    xml_comp_t mpos = x_layer.child("module_positions");
    if (!mpos.ptr()) {
      printout(WARNING, det_name, "Layer %d has no <module_positions> - skipping modules", layerID);

      // Even if no modules, still propagate envelope metadata
      layer_vol->GetShape()->ComputeBBox();
      layerParams.set<double>("envelope_r_min", env_rmin_tol / dd4hep::mm);
      layerParams.set<double>("envelope_r_max", env_rmax_tol / dd4hep::mm);
      layerParams.set<double>("envelope_z_min", env_zmin_tol / dd4hep::mm);
      layerParams.set<double>("envelope_z_max", env_zmax_tol / dd4hep::mm);
      continue;
    }

    int moduleIndexInLayer = 1;
    for (xml_coll_t mp(mpos, _U(module)); mp; ++mp, ++moduleIndexInLayer, ++globalModuleID) {
      xml_comp_t xm = mp;

      const double modX = xm.attr<double>(_Unicode(posX));
      const double modY = xm.attr<double>(_Unicode(posY));
      const double modRotZ = xm.attr<double>(_Unicode(rotZ));
      const std::string side = xm.attr<std::string>(_Unicode(side));
      const double modZ = (side == "front" ? frontZ : backZ);

      // Keep your required front/back rotation convention
      RotationZYX rotLocal(modRotZ, 0.0, (side == "back" ? M_PI : 0.0));
      Transform3D modTr(rotLocal, Position(modX, modY, modZ));

      int sensorIndexInModule = 1;

      for (const auto& cdef : moduleComponents) {
        Material mat = description.material(cdef.material);
        Box shape(cdef.dx / 2.0, cdef.dy / 2.0, cdef.dz / 2.0);
        Volume c_vol(cdef.name, shape, mat);

        if (!cdef.vis.empty()) {
          c_vol.setVisAttributes(description.visAttributes(cdef.vis));
        }
        if (cdef.sensitive) {
          c_vol.setSensitiveDetector(sens);
        }

        Transform3D compLocalTr(Rotation3D(), Position(cdef.px, cdef.py, cdef.pz));
        Transform3D compTr = modTr * compLocalTr;

        PlacedVolume comp_pv = layer_vol.placeVolume(c_vol, compTr);
        comp_pv.addPhysVolID("layer", layerID)
               .addPhysVolID("module", globalModuleID);

        if (cdef.sensitive) {
          comp_pv.addPhysVolID("sensor", sensorIndexInModule);

          std::string sensorName =
              _toString(layerID, "layer%d") +
              _toString(globalModuleID, "_module%d") +
              _toString(sensorIndexInModule, "_sensor%d");

          DetElement sensorDE(layerDE, sensorName, globalModuleID * 10 + sensorIndexInModule);
          sensorDE.setPlacement(comp_pv);

          auto& sensorParams =
              DD4hepDetectorHelper::ensureExtension<dd4hep::rec::VariantParameters>(sensorDE);
          sensorParams.set<std::string>("axis_definitions", "XZY");

          Vector3D u(-1., 0., 0.);
          Vector3D v(0., -1., 0.);
          Vector3D n(0., 0., 1.);
          SurfaceType type(SurfaceType::Sensitive);
          VolPlane surf(c_vol, type, cdef.inner, cdef.outer, u, v, n);
          volSurfaceList(sensorDE)->push_back(surf);

          ++sensorIndexInModule;
        }
      }
    }

    // --------------------------------------------------------------
    // Envelope metadata, like the working original
    // --------------------------------------------------------------
    layer_vol->GetShape()->ComputeBBox();

    layerParams.set<double>("envelope_r_min", env_rmin_tol / dd4hep::mm);
    layerParams.set<double>("envelope_r_max", env_rmax_tol / dd4hep::mm);
    layerParams.set<double>("envelope_z_min", env_zmin_tol / dd4hep::mm);
    layerParams.set<double>("envelope_z_max", env_zmax_tol / dd4hep::mm);

    if (x_env.ptr()) {
      printout(INFO, det_name,
               "Layer %d envelope: length=%8.3f mm zstart=%8.3f mm "
               "tol(rmin,rmax,zmin,zmax)=(%6.3f,%6.3f,%6.3f,%6.3f) mm",
               layerID,
               env_length / mm, env_zstart / mm,
               env_rmin_tol / mm, env_rmax_tol / mm,
               env_zmin_tol / mm, env_zmax_tol / mm);
    }
  }

  // ------------------------------------------------------------------
  // Place the full detector assembly into the mother volume LAST
  // ------------------------------------------------------------------
  pv = motherVol.placeVolume(assembly, posAndRot);
  pv.addPhysVolID("system", det_id);
  sdet.setPlacement(pv);

  return sdet;
}

DECLARE_DETELEMENT(ip6_B0Tracker, create_B0Tracker)
