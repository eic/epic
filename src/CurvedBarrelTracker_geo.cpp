// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 - 2025 Whitney Armstrong, Jonathan Witte, Shujie Li

/** Curved Barrel tracker
 * - Derived from "BarrelTrackerWithFrame_geo.cpp".
 *
 * - Designed to process "vertex_barrel.xml":
 *       - When the upper and lower modules are provided, builds the four-section
 *         EIC-LAS RSU with three sensitive tiles and inactive gaps per section
 *
 *
 * \code
 * \endcode
 *
 *
 * @author Whitney Armstrong, Jonathan Witte, Shujie Li
 */

#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/Printout.h"
#include "DD4hep/Shapes.h"
#include "DDRec/DetectorData.h"
#include "DDRec/Surface.h"
#include "XML/Layering.h"
#include "XML/Utilities.h"
#include <algorithm>
#include <array>
#include <cmath>
#include "DD4hepDetectorHelper.h"

using namespace std;
using namespace dd4hep;
using namespace dd4hep::rec;
using namespace dd4hep::detail;

static Ref_t create_CurvedBarrelTracker(Detector& description, xml_h e, SensitiveDetector sens) {

  xml_det_t x_det = e;
  Material air    = description.air();
  int det_id      = x_det.id();
  string det_name = x_det.nameStr();
  DetElement sdet(det_name, det_id);

  typedef vector<PlacedVolume> Placements;
  map<string, Volume> volumes;
  map<string, Placements> sensitives;
  map<string, std::vector<VolPlane>> volplane_surfaces;
  map<string, std::vector<double>> module_length;

  PlacedVolume pv;

  // Set detector type flag
  dd4hep::xml::setDetectorTypeFlag(x_det, sdet);
  auto& params = DD4hepDetectorHelper::ensureExtension<dd4hep::rec::VariantParameters>(sdet);

  // The XML describes one quarter-RSU section. Upper/lower module prototypes mirror
  // the r-phi ordering, while consecutive z placements provide the other two sections.
  xml_comp_t x_rsu_layout   = x_det.child(_Unicode(rsu_layout));
  int rsu_tiles_per_section = x_rsu_layout.attr<int>(_Unicode(tiles_per_section));
  double rsu_backbone_z     = x_rsu_layout.attr<double>(_Unicode(backbone_z));
  double rsu_biasing_rphi   = x_rsu_layout.attr<double>(_Unicode(biasing_rphi));
  double rsu_readout_rphi   = x_rsu_layout.attr<double>(_Unicode(readout_rphi));
  double rsu_power_z        = x_rsu_layout.attr<double>(_Unicode(power_z));
  double rsu_tile_rphi      = x_rsu_layout.attr<double>(_Unicode(tile_rphi));
  double rsu_tile_pitch_z   = x_rsu_layout.attr<double>(_Unicode(tile_pitch_z));
  double rsu_tile_z         = x_rsu_layout.attr<double>(_Unicode(tile_z));

  // VertexBarrelHits allocates two bits to the sensor field, with zero reserved.
  if (rsu_tiles_per_section < 1 || rsu_tiles_per_section > 3 || rsu_backbone_z <= 0 ||
      rsu_biasing_rphi <= 0 || rsu_readout_rphi <= 0 || rsu_power_z <= 0 || rsu_tile_rphi <= 0 ||
      rsu_tile_pitch_z <= 0 || rsu_tile_z <= 0) {
    printout(ERROR, "CurvedBarrelTracker", "Invalid RSU layout dimensions or tile count.");
    throw runtime_error("Invalid RSU layout in vertex_barrel.xml.");
  }

  // Add the volume boundary material if configured
  for (xml_coll_t bmat(x_det, _Unicode(boundary_material)); bmat; ++bmat) {
    xml_comp_t x_boundary_material = bmat;
    DD4hepDetectorHelper::xmlToProtoSurfaceMaterial(x_boundary_material, params,
                                                    "boundary_material");
  }

  Assembly assembly(det_name);

  sens.setType("tracker");

  // loop over the modules
  for (xml_coll_t mi(x_det, _U(module)); mi; ++mi) {
    xml_comp_t x_mod = mi;
    string m_nam     = x_mod.nameStr();
    double m_rmin    = x_mod.rmin();
    double m_length  = x_mod.length();
    double m_width   = x_mod.width();
    module_length[m_nam].push_back(m_length);

    if (volumes.find(m_nam) != volumes.end()) {
      printout(ERROR, "CurvedBarrelTracker",
               string((string("Module with named ") + m_nam + string(" already exists."))).c_str());
      throw runtime_error("Logics error in building modules.");
    }
    int ncomponents        = 0;
    int sensor_number      = 1;
    double total_thickness = 0;

    // Compute module total thickness from components
    xml_coll_t ci(x_mod, _U(module_component));
    for (ci.reset(), total_thickness = 0.0; ci; ++ci) {
      total_thickness += xml_comp_t(ci).thickness();
    }
    // the module assembly volume
    Assembly m_vol(m_nam);
    volumes[m_nam] = m_vol;
    m_vol.setVisAttributes(description.visAttributes(x_mod.visStr()));

    double thickness_so_far = 0.0;
    for (xml_coll_t mci(x_mod, _U(module_component)); mci; ++mci, ++ncomponents) {
      xml_comp_t x_comp  = mci;
      const string c_nam = x_comp.nameStr(); // _toString(ncomponents, "component%d");
      const string c_mat = x_comp.materialStr();
      double c_thickness = x_comp.thickness();
      double c_width     = getAttrOrDefault(x_comp, _U(width), m_width);
      double c_length    = getAttrOrDefault(x_comp, _U(length), m_length);
      double c_rmin      = m_rmin + thickness_so_far;
      double c_dphi      = c_width / c_rmin;

      auto attach_sensitive = [&](Volume& component_volume, PlacedVolume& component_pv) {
        component_pv.addPhysVolID("sensor", sensor_number++);
        component_volume.setSensitiveDetector(sens);
        sensitives[m_nam].push_back(component_pv);
        // -------- create a measurement plane for the tracking surface attched to the sensitive volume -----
        Vector3D u(-1., 0., 0.);
        Vector3D v(0., -1., 0.);
        Vector3D n(0., 0., 1.);
        // compute the inner and outer thicknesses that need to be assigned to the tracking surface
        // depending on wether the support is above or below the sensor
        double inner_thickness = thickness_so_far + c_thickness / 2.0;
        double outer_thickness = total_thickness - inner_thickness;

        SurfaceType type(SurfaceType::Sensitive);
        VolPlane surf(component_volume, type, inner_thickness, outer_thickness, u, v, n);
        volplane_surfaces[m_nam].push_back(surf);
      };

      if (c_nam == "RSU") {
        const string c_type = x_comp.typeStr();
        if (c_type != "upper" && c_type != "lower") {
          printout(ERROR, "CurvedBarrelTracker",
                   string((string("Module ") + m_nam + string(": invalid RSU component type [") +
                           c_type + string("], should be upper or lower")))
                       .c_str());
          throw runtime_error("Logics error in building modules.");
        }

        const double partition_area =
            c_width * rsu_backbone_z +
            (rsu_biasing_rphi + rsu_readout_rphi) * (c_length - rsu_backbone_z) +
            rsu_tile_rphi * rsu_tiles_per_section * (rsu_tile_z + rsu_power_z);
        const double length_tolerance = 1.0e-12 * std::max(c_width, c_length);
        const double area_tolerance   = 1.0e-12 * c_width * c_length;
        if (std::abs(rsu_biasing_rphi + rsu_tile_rphi + rsu_readout_rphi - c_width) >
                length_tolerance ||
            std::abs(rsu_backbone_z + rsu_tiles_per_section * rsu_tile_pitch_z - c_length) >
                length_tolerance ||
            std::abs(rsu_tile_z + rsu_power_z - rsu_tile_pitch_z) > length_tolerance ||
            std::abs(partition_area - c_width * c_length) > area_tolerance) {
          printout(ERROR, "CurvedBarrelTracker",
                   string((string("Module ") + m_nam +
                           string(": RSU pieces do not form a positive rectangular section.")))
                       .c_str());
          throw runtime_error("Invalid RSU section geometry.");
        }

        const string component_base = m_nam + _toString(ncomponents, "_component%d_") + c_type;
        auto place_rsu_piece = [&](const string& piece_name, double rphi_start, double rphi_width,
                                   double z_start, double z_width, const string& vis,
                                   bool sensitive) {
          Tube piece_solid(c_rmin, c_rmin + c_thickness, z_width / 2.0, rphi_start / c_rmin,
                           (rphi_start + rphi_width) / c_rmin);
          Volume piece_volume(component_base + "_" + piece_name, piece_solid,
                              description.material(c_mat));
          piece_volume.setRegion(description, x_comp.regionStr());
          piece_volume.setLimitSet(description, x_comp.limitsStr());
          piece_volume.setVisAttributes(description, vis);
          PlacedVolume piece_pv = m_vol.placeVolume(
              piece_volume, Position(0, 0, -c_length / 2.0 + z_start + z_width / 2.0));
          if (sensitive) {
            attach_sensitive(piece_volume, piece_pv);
          }
        };

        if (x_comp.isSensitive()) {
          // The backbone owns the two passive-band corner areas, so all pieces are disjoint.
          place_rsu_piece("backbone", 0, c_width, 0, rsu_backbone_z, "SVTReadoutVis", false);

          const double first_passive_rphi = c_type == "upper" ? rsu_biasing_rphi : rsu_readout_rphi;
          const double second_passive_rphi =
              c_type == "upper" ? rsu_readout_rphi : rsu_biasing_rphi;
          const string first_passive_name  = c_type == "upper" ? "biasing" : "readout";
          const string second_passive_name = c_type == "upper" ? "readout" : "biasing";
          const double z_after_backbone    = c_length - rsu_backbone_z;
          place_rsu_piece(first_passive_name, 0, first_passive_rphi, rsu_backbone_z,
                          z_after_backbone, "SVTReadoutVis", false);
          place_rsu_piece(second_passive_name, first_passive_rphi + rsu_tile_rphi,
                          second_passive_rphi, rsu_backbone_z, z_after_backbone, "SVTReadoutVis",
                          false);

          for (int tile = 0; tile < rsu_tiles_per_section; ++tile) {
            const double tile_z_start = rsu_backbone_z + tile * rsu_tile_pitch_z;
            place_rsu_piece(_toString(tile + 1, "tile%d"), first_passive_rphi, rsu_tile_rphi,
                            tile_z_start, rsu_tile_z, x_comp.visStr(), true);
            place_rsu_piece(_toString(tile + 1, "power_gap%d"), first_passive_rphi, rsu_tile_rphi,
                            tile_z_start + rsu_tile_z, rsu_power_z, "SVTReadoutVis", false);
          }
        } else {
          // Inactive silicon and copper remain continuous across the complete section.
          place_rsu_piece("full", 0, c_width, 0, c_length, x_comp.visStr(), false);
        }

        thickness_so_far += c_thickness;
        continue;
      }

      // Regular components retain the original single-volume behavior.
      Tube c_tube(c_rmin, c_rmin + c_thickness, c_length / 2, 0, c_dphi);
      Volume c_vol(c_nam, c_tube, description.material(c_mat));
      pv = m_vol.placeVolume(c_vol, Position(0, 0, 0));
      c_vol.setRegion(description, x_comp.regionStr());
      c_vol.setLimitSet(description, x_comp.limitsStr());
      c_vol.setVisAttributes(description, x_comp.visStr());
      if (x_comp.isSensitive()) {
        attach_sensitive(c_vol, pv);
      }
      thickness_so_far += c_thickness;
    }
  }

  // now build the layers by alternating upper and lower modules.
  // each layer has a gap b/w the upper and lower half (see phi0 and phi1)
  for (xml_coll_t li(x_det, _U(layer)); li; ++li) {
    xml_comp_t x_layer  = li;
    xml_comp_t x_barrel = x_layer.child(_U(barrel_envelope));
    xml_comp_t x_layout = x_layer.child(_U(rphi_layout));
    xml_comp_t z_layout = x_layer.child(_U(z_layout)); // Get the <z_layout> element.
    int lay_id          = x_layer.id();
    string m_nam        = x_layer.moduleStr();
    string lay_nam      = det_name + _toString(x_layer.id(), "_layer%d");
    Tube lay_tub(x_barrel.inner_r(), x_barrel.outer_r(), x_barrel.z_length() / 2.0);

    Volume lay_vol(lay_nam, lay_tub, air); // Create the layer envelope volume.
    Position lay_pos(0, 0, getAttrOrDefault(x_barrel, _U(z0), 0.));
    lay_vol.setVisAttributes(description.visAttributes(x_layer.visStr()));

    int l_nphi      = x_layout.nphi(); // Number of modules in phi.
    double phi0     = x_layout.phi0(); // Starting phi of first module.
    int nphi[2]     = {int((l_nphi + 1) / 2),
                       int(l_nphi / 2)}; // number of modules in uppper and lower modules.
    double phi_incr = (M_PI * 2 - phi0 * 2) / l_nphi; // Phi increment for one module.
    double z0       = z_layout.z0();                  // Z position of first module in phi.
    double nz       = z_layout.nz();                  // Number of modules to place in z.

    Volume module_env[2];
    Placements sensVols[2];
    string m_nams[2] = {m_nam + "_upper", m_nam + "_lower"};
    // if both upper and lower modules are provided (for RSU tiles)
    if ((volumes.find(m_nams[0]) != volumes.end()) && (volumes.find(m_nams[1]) != volumes.end())) {
      if (nphi[0] != nphi[1]) {
        printout(ERROR, "CurvedBarrelTracker",
                 string((string("Layer ") + lay_nam +
                         string(": nphi must be even number to allow upper and lower modules"))
                            .c_str()));
        throw runtime_error("Logics error in building modules.");
      }
      module_env[0] = volumes[m_nams[0]];
      module_env[1] = volumes[m_nams[1]];
      sensVols[0]   = sensitives[m_nams[0]];
      sensVols[1]   = sensitives[m_nams[1]];
    }
    // for other regular modules
    else if (volumes.find(m_nam) != volumes.end()) {
      module_env[0] = volumes[m_nam];
      module_env[1] = volumes[m_nam];
      sensVols[0]   = sensitives[m_nam];
      sensVols[1]   = sensitives[m_nam];
      m_nams[0]     = m_nam;
      m_nams[1]     = m_nam;
    }
    DetElement lay_elt(sdet, lay_nam, lay_id);

    // the local coordinate systems of modules in dd4hep and acts differ
    // see http://acts.web.cern.ch/ACTS/latest/doc/group__DD4hepPlugins.html
    auto& layerParams =
        DD4hepDetectorHelper::ensureExtension<dd4hep::rec::VariantParameters>(lay_elt);

    for (xml_coll_t lmat(x_layer, _Unicode(layer_material)); lmat; ++lmat) {
      xml_comp_t x_layer_material = lmat;
      DD4hepDetectorHelper::xmlToProtoSurfaceMaterial(x_layer_material, layerParams,
                                                      "layer_material");
    }

    // Z increment for module placement along Z axis.
    // Adjust for z0 at center of module rather than
    // the end of cylindrical envelope.
    // double z_incr = nz > 1 ? (2.0 * abs(z0)) / (nz - 1) : 0.0;
    // Starting z for module placement along Z axis.
    int module = 1;
    // Loop over the number of modules in phi.
    for (int kk = 0; kk < 2; kk++) {
      int iphi      = nphi[kk];
      double z_incr = module_length[m_nams[kk]][0];
      for (int ii = 0; ii < iphi; ii++) {
        double dphi = phi0;
        if (ii > (iphi / 2) - 1)
          dphi = 2 * phi0;
        // Loop over the number of modules in z.
        double module_z = z0;
        for (int j = 0; j < nz; j++, module_z += z_incr) {
          string module_name = _toString(module, "module%d");
          DetElement mod_elt(lay_elt, module_name, module);
          Transform3D tr(
              RotationZYX(dphi + phi_incr * kk + phi_incr * ii * 2, 0, 0),
              Position(0, 0, module_z)); // altering upper and lower module to fill every other row
          pv = lay_vol.placeVolume(module_env[kk], tr);
          pv.addPhysVolID("module", module);
          mod_elt.setPlacement(pv);
          for (size_t ic = 0; ic < sensVols[kk].size(); ++ic) {
            PlacedVolume sens_pv = sensVols[kk][ic];
            DetElement comp_de(mod_elt, std::string("de_") + sens_pv.volume().name(), module);
            comp_de.setPlacement(sens_pv);
            auto& comp_de_params =
                DD4hepDetectorHelper::ensureExtension<dd4hep::rec::VariantParameters>(comp_de);
            comp_de_params.set<string>("axis_definitions", "XYZ");
            // note from Wouter: Region allows setting specific physics in regions of space. Limits similar but with production limits. Currently unused in simulations. But we'd sort of tried to make sure we don't prevent ourselves from using it entirely.
            // comp_de.setAttributes(description, sens_pv.volume(), x_layer.regionStr(), x_layer.limitsStr(),
            //                       xml_det_t(xmleles[m_nam]).visStr());
            //
            volSurfaceList(comp_de)->push_back(volplane_surfaces[m_nams[kk]][ic]);
          }

          /// Increase counters etc.
          module++;
        }
      }
    }
    // Create the PhysicalVolume for the layer.
    pv = assembly.placeVolume(lay_vol, lay_pos); // Place layer in mother
    pv.addPhysVolID("layer", lay_id);            // Set the layer ID.
    lay_elt.setAttributes(description, lay_vol, x_layer.regionStr(), x_layer.limitsStr(),
                          x_layer.visStr());
    lay_elt.setPlacement(pv);
  }
  sdet.setAttributes(description, assembly, x_det.regionStr(), x_det.limitsStr(), x_det.visStr());
  assembly.setVisAttributes(description.invisible());
  pv = description.pickMotherVolume(sdet).placeVolume(assembly);
  pv.addPhysVolID("system", det_id); // Set the subdetector system ID.
  sdet.setPlacement(pv);
  return sdet;
}

//@}
// clang-format off
DECLARE_DETELEMENT(epic_CylinderSVTBarrel,create_CurvedBarrelTracker)
