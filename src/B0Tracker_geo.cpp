// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 - 2026 Whitney Armstrong, Tom Bleher, Igor Korover

/*! B0 Tracker: stations, each a support plate with sensor staves on both faces.
 *
 * @author Whitney Armstrong, Tom Bleher, Igor Korover
 *
 * Each <station> places a <support_plate> (the <support_stack> slices extruded
 * along its outline) and a <module_layout> of front and back <face>. Every face
 * is one ACTS disc layer holding the shared TrackingUnit at its <module>
 * positions. The cellID layer is 2*(station-1) + (1 front | 2 back), so it is
 * monotonic in z, and module ids restart per face.
 *
 * @{
 */

#include "DD4hep/DetFactoryHelper.h"
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
#include <vector>

using namespace dd4hep;
using namespace dd4hep::rec;

namespace {

struct SupportSlice {
  Volume volume;
  Position position;
};

struct SupportPlates {
  std::map<std::string, std::vector<SupportSlice>> slices; // by plate name
  double thickness = 0.0;
};

struct TrackingUnit {
  Assembly assembly;
  std::vector<PlacedVolume> sensors; // sensor id = index + 1
  std::vector<VolPlane> surfaces;
  double z_min = +std::numeric_limits<double>::infinity(); // z extent of the whole stack
  double z_max = -std::numeric_limits<double>::infinity();
};

xml_h find_named(xml_det_t x_det, const xml::Strng_t& tag, const std::string& name) {
  for (xml_coll_t it(x_det, tag); it; ++it) {
    xml_comp_t x = it;
    if (x.nameStr() == name) {
      return x;
    }
  }
  return xml_h();
}

TrackingUnit build_tracking_unit(Detector& description, SensitiveDetector sens,
                                 const std::string& det_name, xml_comp_t x_unit) {
  TrackingUnit unit;
  unit.assembly = Assembly(det_name + "_TrackingUnit");
  unit.assembly.setVisAttributes(description,
                                 getAttrOrDefault<std::string>(x_unit, _Unicode(vis), ""));

  // The stack extent is needed first: it sets the sensor surface thicknesses below
  for (xml_coll_t comp(x_unit, _U(module_component)); comp; ++comp) {
    xml_comp_t x_comp   = comp;
    const double pz     = x_comp.position().z();
    const double half_z = xml_dim_t(x_comp.child(_U(box))).z() / 2.0;
    unit.z_min          = std::min(unit.z_min, pz - half_z);
    unit.z_max          = std::max(unit.z_max, pz + half_z);
  }

  const Vector3D u(-1.0, 0.0, 0.0);
  const Vector3D v(0.0, -1.0, 0.0);
  const Vector3D n(0.0, 0.0, 1.0);
  for (xml_coll_t comp(x_unit, _U(module_component)); comp; ++comp) {
    xml_comp_t x_comp = comp;
    xml_dim_t x_box   = x_comp.child(_U(box));
    xml_dim_t x_pos   = x_comp.position();

    Volume vol(det_name + "_" + x_comp.nameStr(),
               Box(x_box.x() / 2.0, x_box.y() / 2.0, x_box.z() / 2.0),
               description.material(x_comp.materialStr()));
    vol.setVisAttributes(description, getAttrOrDefault<std::string>(x_comp, _Unicode(vis), ""));
    PlacedVolume pv = unit.assembly.placeVolume(vol, Position(x_pos.x(), x_pos.y(), x_pos.z()));

    if (x_comp.isSensitive()) {
      vol.setSensitiveDetector(sens);
      pv.addPhysVolID("sensor", static_cast<int>(unit.sensors.size()) + 1);
      unit.sensors.push_back(pv);
      // DDRec inner/outer thickness: the extent of the stave stack on either side of the
      // sensor
      unit.surfaces.emplace_back(vol, SurfaceType(SurfaceType::Sensitive), x_pos.z() - unit.z_min,
                                 unit.z_max - x_pos.z(), u, v, n);
    }
  }
  unit.assembly->GetShape()->ComputeBBox();
  return unit;
}

// Slices are stacked upstream to downstream, centered on the plate mid-plane
SupportPlates build_support_plates(Detector& description, xml_det_t x_det) {
  SupportPlates plates;
  xml_comp_t x_stack          = x_det.child(_Unicode(support_stack));
  const std::string stack_vis = getAttrOrDefault<std::string>(x_stack, _Unicode(vis), "");
  for (xml_coll_t sl(x_stack, _U(slice)); sl; ++sl) {
    plates.thickness += xml_comp_t(sl).thickness();
  }

  for (xml_coll_t pl(x_det, _Unicode(support_plate)); pl; ++pl) {
    xml_comp_t x_plate           = pl;
    const std::string plate_name = x_plate.nameStr();
    std::vector<double> xs;
    std::vector<double> ys;
    for (xml_coll_t point(x_plate, _U(point)); point; ++point) {
      xs.push_back(xml_comp_t(point).x());
      ys.push_back(xml_comp_t(point).y());
    }

    auto& slices = plates.slices[plate_name];
    double z_low = -plates.thickness / 2.0;
    for (xml_coll_t sl(x_stack, _U(slice)); sl; ++sl) {
      xml_comp_t x_slice     = sl;
      const double thickness = x_slice.thickness();
      Volume vol(plate_name + "_" + x_slice.nameStr(),
                 ExtrudedPolygon(xs, ys, {-thickness / 2.0, thickness / 2.0}, {0.0, 0.0},
                                 {0.0, 0.0}, {1.0, 1.0}),
                 description.material(x_slice.materialStr()));
      vol.setVisAttributes(description,
                           getAttrOrDefault<std::string>(x_slice, _Unicode(vis), stack_vis));
      slices.push_back({vol, Position(0.0, 0.0, z_low + thickness / 2.0)});
      z_low += thickness;
    }
  }
  return plates;
}

} // namespace

static Ref_t create_B0Tracker(Detector& description, xml_h e, SensitiveDetector sens) {
  xml_det_t x_det            = e;
  const int det_id           = x_det.id();
  const std::string det_name = x_det.nameStr();

  DetElement sdet(det_name, det_id);
  Assembly assembly(det_name);

  Volume mother_vol  = description.pickMotherVolume(sdet);
  xml::Component pos = x_det.position();
  xml::Component rot = x_det.rotation();
  Transform3D pos_and_rot(RotationZYX(rot.z(), rot.y(), rot.x()),
                          Position(pos.x(), pos.y(), pos.z()));

  xml::setDetectorTypeFlag(x_det, sdet);
  auto& det_params = DD4hepDetectorHelper::ensureExtension<VariantParameters>(sdet);
  for (xml_coll_t bmat(x_det, _Unicode(boundary_material)); bmat; ++bmat) {
    DD4hepDetectorHelper::xmlToProtoSurfaceMaterial(bmat, det_params, "boundary_material");
  }

  assembly.setAttributes(description, x_det.regionStr(), x_det.limitsStr(), x_det.visStr());
  sens.setType("tracker");

  const TrackingUnit unit =
      build_tracking_unit(description, sens, det_name, x_det.child(_U(module)));
  const SupportPlates plates = build_support_plates(description, x_det);

  // Plate mid-plane to module center, with the module's z_min side flush on the plate
  const double module_offset = plates.thickness / 2.0 - unit.z_min;

  xml_comp_t x_env          = x_det.child(_U(envelope));
  const double env_rmin_tol = x_env.attr<double>(_Unicode(rmin_tolerance));
  const double env_rmax_tol = x_env.attr<double>(_Unicode(rmax_tolerance));
  const double env_zmin_tol = x_env.attr<double>(_Unicode(zmin_tolerance));
  const double env_zmax_tol = x_env.attr<double>(_Unicode(zmax_tolerance));
  const std::string env_vis = getAttrOrDefault<std::string>(x_env, _Unicode(vis), "");

  for (xml_coll_t st(x_det, _Unicode(station)); st; ++st) {
    xml_comp_t x_station = st;
    const int station    = x_station.id();
    if (station < 1) {
      throw std::runtime_error(det_name + ": station ids start at 1");
    }
    xml_dim_t x_station_pos = x_station.position();
    const Position station_pos(x_station_pos.x(), x_station_pos.y(), x_station_pos.z());

    // The plates are passive: material mapping projects them onto the faces
    for (xml_coll_t sup(x_station, _Unicode(support)); sup; ++sup) {
      const std::string ref = xml_comp_t(sup).attr<std::string>(_Unicode(ref));
      const auto plate_it   = plates.slices.find(ref);
      if (plate_it == plates.slices.end()) {
        throw std::runtime_error(det_name + ": unknown <support ref=\"" + ref + "\">");
      }
      for (const auto& slice : plate_it->second) {
        assembly.placeVolume(slice.volume, station_pos + slice.position);
      }
    }

    const std::string layout_ref =
        xml_comp_t(x_station.child(_Unicode(layout))).attr<std::string>(_Unicode(ref));
    const xml_h x_layout = find_named(x_det, _Unicode(module_layout), layout_ref);
    if (!x_layout.ptr()) {
      throw std::runtime_error(det_name + ": unknown <layout ref=\"" + layout_ref + "\">");
    }

    for (xml_coll_t fc(x_layout, _Unicode(face)); fc; ++fc) {
      xml_comp_t x_face      = fc;
      const std::string side = x_face.attr<std::string>(_Unicode(side));
      if (side != "front" && side != "back") {
        throw std::runtime_error(det_name + ": <face side=\"" + side +
                                 "\">, expected \"front\" or \"back\"");
      }
      const bool is_front = side == "front";
      const int layer_id  = 2 * (station - 1) + (is_front ? 1 : 2);

      // ACTS centers the disc layer on the face origin, so put it at the unit
      // center, midway between the sensor planes
      const Position face_pos(station_pos.x(), station_pos.y(),
                              station_pos.z() + (is_front ? -module_offset : module_offset));

      const std::string face_name = det_name + "_station" + std::to_string(station) + "_" + side;
      Assembly face_vol(face_name);
      face_vol.setVisAttributes(description, env_vis);
      PlacedVolume face_pv = assembly.placeVolume(face_vol, face_pos);
      face_pv.addPhysVolID("layer", layer_id);

      DetElement face_de(sdet, face_name + "_P", layer_id);
      face_de.setPlacement(face_pv);

      // Front modules are flipped about y so they also sit on the plate; flipping about y
      // keeps the stave's orientation along y, so rotZ means the same on both faces
      int module_id = 1;
      for (xml_coll_t mp(x_face, _U(module)); mp; ++mp, ++module_id) {
        xml_comp_t x_mod = mp;
        RotationZYX mod_rot(getAttrOrDefault<double>(x_mod, _Unicode(rotZ), 0.0),
                            is_front ? M_PI : 0.0, 0.0);
        PlacedVolume mod_pv = face_vol.placeVolume(
            unit.assembly, Transform3D(mod_rot, Position(x_mod.x(), x_mod.y(), 0.0)));
        mod_pv.addPhysVolID("module", module_id);

        DetElement mod_de(face_de, _toString(module_id, "module%d"), module_id);
        mod_de.setPlacement(mod_pv);

        for (size_t i = 0; i < unit.sensors.size(); ++i) {
          const int sensor_id = static_cast<int>(i) + 1;
          DetElement sensor_de(mod_de, _toString(sensor_id, "sensor%d"), sensor_id);
          sensor_de.setPlacement(unit.sensors[i]);
          DD4hepDetectorHelper::ensureExtension<VariantParameters>(sensor_de).set<std::string>(
              "axis_definitions", "XYZ");
          volSurfaceList(sensor_de)->push_back(unit.surfaces[i]);
        }
      }
      face_vol->GetShape()->ComputeBBox();

      auto& face_params = DD4hepDetectorHelper::ensureExtension<VariantParameters>(face_de);
      face_params.set<double>("envelope_r_min", env_rmin_tol / dd4hep::mm);
      face_params.set<double>("envelope_r_max", env_rmax_tol / dd4hep::mm);
      face_params.set<double>("envelope_z_min", env_zmin_tol / dd4hep::mm);
      face_params.set<double>("envelope_z_max", env_zmax_tol / dd4hep::mm);
      for (xml_coll_t lmat(x_det, _Unicode(layer_material)); lmat; ++lmat) {
        DD4hepDetectorHelper::xmlToProtoSurfaceMaterial(lmat, face_params, "layer_material");
      }
    }
  }
  assembly->GetShape()->ComputeBBox();

  PlacedVolume pv = mother_vol.placeVolume(assembly, pos_and_rot);
  pv.addPhysVolID("system", det_id);
  sdet.setPlacement(pv);

  return sdet;
}

//@}
// clang-format off
DECLARE_DETELEMENT(ip6_B0Tracker, create_B0Tracker)
