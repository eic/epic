// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 - 2025 Whitney Armstrong, Jonathan Witte, Shujie Li

/** Curved silicon vertex barrel.
 *
 * One module is one physical RSU. Each RSU contains a full inactive-silicon
 * layer, a full 14 um silicon layer, and a full copper layer. Twelve silicon
 * daughters inside the 14 um layer mark the regions that record hits.
 *
 * @author Whitney Armstrong, Jonathan Witte, Shujie Li
 */

#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/Printout.h"
#include "DD4hep/Shapes.h"
#include "DDRec/DetectorData.h"
#include "DDRec/Surface.h"
#include "XML/Utilities.h"
#include "DD4hepDetectorHelper.h"

#include <algorithm>
#include <cmath>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

using namespace std;
using namespace dd4hep;
using namespace dd4hep::rec;
using namespace dd4hep::detail;

namespace {

using Placements = vector<PlacedVolume>;

// A prototype is built once per layer radius and reused at every phi-z position.
struct ModulePrototype {
  Volume volume;
  PlacedVolume active_layer;
  Placements tiles;
  vector<VolPlane> surfaces;
  double length{0.0};
};

} // namespace

static Ref_t create_CurvedBarrelTracker(Detector& description, xml_h e, SensitiveDetector sens) {
  xml_det_t x_det = e;
  Material air    = description.air();
  int det_id      = x_det.id();
  string det_name = x_det.nameStr();
  DetElement sdet(det_name, det_id);

  dd4hep::xml::setDetectorTypeFlag(x_det, sdet);
  auto& params = DD4hepDetectorHelper::ensureExtension<VariantParameters>(sdet);

  // These are measured gaps in the full 2-by-6 tile layout.
  xml_comp_t x_rsu_layout  = x_det.child(_Unicode(rsu_layout));
  int rsu_rows_rphi        = x_rsu_layout.attr<int>(_Unicode(rows_rphi));
  int rsu_halves_z         = x_rsu_layout.attr<int>(_Unicode(halves_z));
  int rsu_tiles_per_z_half = x_rsu_layout.attr<int>(_Unicode(tiles_per_z_half));
  double rsu_backbone_z    = x_rsu_layout.attr<double>(_Unicode(backbone_z));
  double rsu_biasing_rphi  = x_rsu_layout.attr<double>(_Unicode(biasing_rphi));
  double rsu_readout_rphi  = x_rsu_layout.attr<double>(_Unicode(readout_rphi));
  double rsu_power_z       = x_rsu_layout.attr<double>(_Unicode(power_z));
  double rsu_tile_rphi     = x_rsu_layout.attr<double>(_Unicode(tile_rphi));
  double rsu_tile_pitch_z  = x_rsu_layout.attr<double>(_Unicode(tile_pitch_z));
  double rsu_tile_z        = x_rsu_layout.attr<double>(_Unicode(tile_z));
  int rsu_tiles            = rsu_rows_rphi * rsu_halves_z * rsu_tiles_per_z_half;

  // The formulas below describe this specific EIC-LAS RSU. Four sensor-ID bits
  // allow values 1 through 12; zero is left unused.
  if (rsu_rows_rphi != 2 || rsu_halves_z != 2 || rsu_tiles_per_z_half != 3 || rsu_tiles != 12 ||
      rsu_backbone_z <= 0 || rsu_biasing_rphi <= 0 || rsu_readout_rphi <= 0 || rsu_power_z <= 0 ||
      rsu_tile_rphi <= 0 || rsu_tile_pitch_z <= 0 || rsu_tile_z <= 0) {
    printout(ERROR, "CurvedBarrelTracker", "Invalid RSU tile count or gap dimensions.");
    throw runtime_error("Invalid RSU layout in vertex_barrel.xml.");
  }

  for (xml_coll_t bmat(x_det, _Unicode(boundary_material)); bmat; ++bmat) {
    DD4hepDetectorHelper::xmlToProtoSurfaceMaterial(xml_comp_t(bmat), params, "boundary_material");
  }

  Assembly assembly(det_name);
  assembly.setVisAttributes(description.invisible());
  sens.setType("tracker");

  map<string, ModulePrototype> modules;

  // Build one complete RSU prototype for each layer radius.
  for (xml_coll_t mi(x_det, _U(module)); mi; ++mi) {
    xml_comp_t x_mod     = mi;
    string module_name   = x_mod.nameStr();
    double module_rmin   = x_mod.rmin();
    double module_width  = x_mod.width();
    double module_length = x_mod.length();

    if (modules.count(module_name) != 0) {
      printout(ERROR, "CurvedBarrelTracker", "Module '%s' is defined more than once.",
               module_name.c_str());
      throw runtime_error("Duplicate vertex-barrel module name.");
    }

    // Find the active layer before building any solids. Its middle radius is
    // the common angular reference for all three material layers.
    double total_thickness  = 0.0;
    double active_offset    = 0.0;
    double active_thickness = 0.0;
    int active_components   = 0;
    int component_count     = 0;
    for (xml_coll_t ci(x_mod, _U(module_component)); ci; ++ci) {
      xml_comp_t x_comp = ci;
      if (x_comp.thickness() <= 0) {
        throw runtime_error("Vertex-barrel component thickness must be positive.");
      }
      if (x_comp.isSensitive()) {
        active_offset    = total_thickness;
        active_thickness = x_comp.thickness();
        ++active_components;
      }
      total_thickness += x_comp.thickness();
      ++component_count;
    }

    if (component_count != 3 || active_components != 1 || module_width <= 0 || module_length <= 0 ||
        total_thickness <= 0) {
      printout(ERROR, "CurvedBarrelTracker",
               "Module '%s' needs three material layers and one active layer.",
               module_name.c_str());
      throw runtime_error("Invalid vertex-barrel material stack.");
    }

    double reference_radius = module_rmin + active_offset + active_thickness / 2.0;
    double module_dphi      = module_width / reference_radius;
    double half_rphi        = module_width / 2.0;
    double half_z           = module_length / 2.0;
    double length_tolerance = 1.0e-12 * max(module_width, module_length);

    // These equations check that the tile rows and columns fit exactly inside
    // the rectangular RSU. The small differences are the inactive gaps.
    if (abs(2.0 * (rsu_readout_rphi + rsu_tile_rphi + rsu_biasing_rphi) - module_width) >
            length_tolerance ||
        abs(2.0 * (rsu_backbone_z + rsu_tiles_per_z_half * rsu_tile_pitch_z) - module_length) >
            length_tolerance ||
        abs(rsu_tile_z + rsu_power_z - rsu_tile_pitch_z) > length_tolerance) {
      printout(ERROR, "CurvedBarrelTracker",
               "Module '%s' RSU dimensions do not add up to one rectangle.", module_name.c_str());
      throw runtime_error("Invalid full-RSU dimensions.");
    }

    Tube module_solid(module_rmin, module_rmin + total_thickness, module_length / 2.0, 0,
                      module_dphi);
    Volume module_volume(module_name, module_solid, air);
    module_volume.setVisAttributes(description.invisible());

    ModulePrototype prototype;
    prototype.volume = module_volume;
    prototype.length = module_length;

    double thickness_so_far = 0.0;
    for (xml_coll_t ci(x_mod, _U(module_component)); ci; ++ci) {
      xml_comp_t x_comp          = ci;
      string component_name      = module_name + "_" + x_comp.nameStr();
      double component_thickness = x_comp.thickness();
      double component_width     = getAttrOrDefault(x_comp, _U(width), module_width);
      double component_length    = getAttrOrDefault(x_comp, _U(length), module_length);
      double component_rmin      = module_rmin + thickness_so_far;
      double component_rmax      = component_rmin + component_thickness;

      if (abs(component_width - module_width) > length_tolerance ||
          abs(component_length - module_length) > length_tolerance) {
        printout(ERROR, "CurvedBarrelTracker",
                 "Module '%s' material layers must cover the complete RSU.", module_name.c_str());
        throw runtime_error("Partial full-RSU material layer.");
      }

      Tube component_solid(component_rmin, component_rmax, module_length / 2.0, 0, module_dphi);
      Volume component_volume(component_name, component_solid,
                              description.material(x_comp.materialStr()));
      component_volume.setRegion(description, x_comp.regionStr());
      component_volume.setLimitSet(description, x_comp.limitsStr());

      if (!x_comp.isSensitive()) {
        component_volume.setVisAttributes(description, x_comp.visStr());
        module_volume.placeVolume(component_volume);
        thickness_so_far += component_thickness;
        continue;
      }

      // Pattern for several sensitive areas in one module:
      // 1. Place one nonsensitive mother volume for the complete material layer.
      // 2. Place one daughter volume inside it for each sensitive area.
      // 3. Mark only the daughters sensitive and give each a sensor ID and surface.
      //
      // A daughter replaces the mother's material in its own region. Using the
      // same material therefore keeps one continuous silicon layer, while the
      // parts of the mother outside the daughters form the inactive gaps. The
      // daughters must be fully contained and must not overlap one another.
      // Here, sensitive="true" in the XML selects this component for tiling; it
      // does not make the full 14 um mother sensitive.
      component_volume.setVisAttributes(description.invisible());
      prototype.active_layer = module_volume.placeVolume(component_volume);

      double inner_thickness = thickness_so_far + component_thickness / 2.0;
      double outer_thickness = total_thickness - inner_thickness;
      int sensor_id          = 1;

      for (int row_id = 0; row_id < rsu_rows_rphi; ++row_id) {
        // From low to high r-phi: readout, tiles, two central bias bands,
        // tiles, readout. This is symmetric about the RSU center.
        double tile_rphi_start = row_id == 0 ? rsu_readout_rphi : half_rphi + rsu_biasing_rphi;

        for (int half_z_id = 0; half_z_id < rsu_halves_z; ++half_z_id) {
          double half_z_start = -half_z + half_z_id * half_z;

          for (int tile_id = 0; tile_id < rsu_tiles_per_z_half; ++tile_id, ++sensor_id) {
            double tile_z_start   = half_z_start + rsu_backbone_z + tile_id * rsu_tile_pitch_z;
            double tile_z_center  = tile_z_start + rsu_tile_z / 2.0;
            double tile_phi_start = tile_rphi_start / reference_radius;
            double tile_phi_end   = (tile_rphi_start + rsu_tile_rphi) / reference_radius;
            string tile_name      = module_name + _toString(sensor_id, "_tile%02d");

            Tube tile_solid(component_rmin, component_rmax, rsu_tile_z / 2.0, tile_phi_start,
                            tile_phi_end);
            Volume tile_volume(tile_name, tile_solid, description.material(x_comp.materialStr()));
            tile_volume.setRegion(description, x_comp.regionStr());
            tile_volume.setLimitSet(description, x_comp.limitsStr());
            tile_volume.setVisAttributes(description, x_comp.visStr());
            tile_volume.setSensitiveDetector(sens);

            PlacedVolume tile_pv =
                component_volume.placeVolume(tile_volume, Position(0, 0, tile_z_center));
            tile_pv.addPhysVolID("sensor", sensor_id);
            prototype.tiles.push_back(tile_pv);

            Vector3D u(-1., 0., 0.);
            Vector3D v(0., -1., 0.);
            Vector3D n(0., 0., 1.);
            SurfaceType type(SurfaceType::Sensitive);
            prototype.surfaces.emplace_back(tile_volume, type, inner_thickness, outer_thickness, u,
                                            v, n);
          }
        }
      }

      thickness_so_far += component_thickness;
    }

    if (!prototype.active_layer.isValid() || prototype.tiles.size() != 12 ||
        prototype.surfaces.size() != prototype.tiles.size()) {
      printout(ERROR, "CurvedBarrelTracker", "Module '%s' did not produce twelve tiles.",
               module_name.c_str());
      throw runtime_error("Incomplete vertex-barrel RSU prototype.");
    }

    modules.emplace(module_name, std::move(prototype));
  }

  // Place complete RSUs in each layer. The larger phi clearance appears only
  // at the two half-shell seams, never inside a physical RSU.
  for (xml_coll_t li(x_det, _U(layer)); li; ++li) {
    xml_comp_t x_layer  = li;
    xml_comp_t x_barrel = x_layer.child(_U(barrel_envelope));
    xml_comp_t x_layout = x_layer.child(_U(rphi_layout));
    xml_comp_t z_layout = x_layer.child(_U(z_layout));
    int layer_id        = x_layer.id();
    string module_name  = x_layer.moduleStr();
    string layer_name   = det_name + _toString(layer_id, "_layer%d");

    auto module_iter = modules.find(module_name);
    if (module_iter == modules.end()) {
      printout(ERROR, "CurvedBarrelTracker", "Layer '%s' uses unknown module '%s'.",
               layer_name.c_str(), module_name.c_str());
      throw runtime_error("Unknown vertex-barrel module.");
    }
    ModulePrototype& prototype = module_iter->second;

    Tube layer_solid(x_barrel.inner_r(), x_barrel.outer_r(), x_barrel.z_length() / 2.0);
    Volume layer_volume(layer_name, layer_solid, air);
    layer_volume.setVisAttributes(description.visAttributes(x_layer.visStr()));
    DetElement layer_element(sdet, layer_name, layer_id);

    auto& layer_params = DD4hepDetectorHelper::ensureExtension<VariantParameters>(layer_element);
    for (xml_coll_t lmat(x_layer, _Unicode(layer_material)); lmat; ++lmat) {
      DD4hepDetectorHelper::xmlToProtoSurfaceMaterial(xml_comp_t(lmat), layer_params,
                                                      "layer_material");
    }

    int nphi        = x_layout.nphi();
    int nz          = z_layout.nz();
    double phi0     = x_layout.phi0();
    double phi_step = (2.0 * M_PI - 2.0 * phi0) / nphi;
    double first_z  = z_layout.z0();
    if (nphi <= 0 || nphi % 2 != 0 || nz <= 0) {
      throw runtime_error("Vertex-barrel nphi must be positive and even; nz must be positive.");
    }

    int module_id = 1;
    for (int iphi = 0; iphi < nphi; ++iphi) {
      double seam_offset = iphi < nphi / 2 ? phi0 : 2.0 * phi0;
      double module_phi  = seam_offset + iphi * phi_step;
      double module_z    = first_z;

      for (int iz = 0; iz < nz; ++iz, module_z += prototype.length, ++module_id) {
        string placed_name = _toString(module_id, "module%d");
        DetElement module_element(layer_element, placed_name, module_id);
        Transform3D transform(RotationZYX(module_phi, 0, 0), Position(0, 0, module_z));
        PlacedVolume module_pv = layer_volume.placeVolume(prototype.volume, transform);
        module_pv.addPhysVolID("module", module_id);
        module_element.setPlacement(module_pv);

        // The DetElement tree must match the physical-volume tree:
        // module -> active mother -> sensitive tile. This lets DD4hep combine
        // the layer, module, and sensor IDs, and gives ACTS the tile's correct
        // world transform. Attach each VolPlane to its tile DetElement, not to
        // the common active mother.
        DetElement active_element(module_element, "active_si", 0);
        active_element.setPlacement(prototype.active_layer);
        for (size_t tile_index = 0; tile_index < prototype.tiles.size(); ++tile_index) {
          int sensor_id        = static_cast<int>(tile_index) + 1;
          PlacedVolume tile_pv = prototype.tiles[tile_index];
          DetElement tile_element(active_element, tile_pv.volume().name(), sensor_id);
          tile_element.setPlacement(tile_pv);
          auto& tile_params =
              DD4hepDetectorHelper::ensureExtension<VariantParameters>(tile_element);
          tile_params.set<string>("axis_definitions", "XYZ");
          volSurfaceList(tile_element)->push_back(prototype.surfaces[tile_index]);
        }
      }
    }

    Position layer_position(0, 0, getAttrOrDefault(x_barrel, _U(z0), 0.));
    PlacedVolume layer_pv = assembly.placeVolume(layer_volume, layer_position);
    layer_pv.addPhysVolID("layer", layer_id);
    layer_element.setAttributes(description, layer_volume, x_layer.regionStr(), x_layer.limitsStr(),
                                x_layer.visStr());
    layer_element.setPlacement(layer_pv);
  }

  sdet.setAttributes(description, assembly, x_det.regionStr(), x_det.limitsStr(), x_det.visStr());
  PlacedVolume detector_pv = description.pickMotherVolume(sdet).placeVolume(assembly);
  detector_pv.addPhysVolID("system", det_id);
  sdet.setPlacement(detector_pv);
  return sdet;
}

//@}
// clang-format off
DECLARE_DETELEMENT(epic_CylinderSVTBarrel,create_CurvedBarrelTracker)
