// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 - 2024, Nicolas Schmidt, Chun Yuen Tsang

/** \addtogroup Trackers Trackers
 * \brief Type: **Endcap Tracker with TOF**.
 *
 * \ingroup trackers
 *
 * @{
 */
#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/Printout.h"
#include "DD4hep/Shapes.h"
#include "DD4hepDetectorHelper.h"
#include "DDRec/DetectorData.h"
#include "DDRec/Surface.h"
#include "XML/Layering.h"
#include "XML/Utilities.h"
#include <array>
#include <map>
#include <unordered_set>

using namespace std;
using namespace dd4hep;
using namespace dd4hep::rec;
using namespace dd4hep::detail;

static Ref_t create_detector(Detector& description, xml_h e, SensitiveDetector sens) {
  typedef vector<PlacedVolume> Placements;
  xml_det_t x_det      = e;
  int det_id           = x_det.id();
  std::string det_name = x_det.nameStr();
  DetElement sdet(det_name, det_id);
  Material air = description.material("Air");
  PlacedVolume pv;

  map<string, std::array<double, 2>> module_thicknesses;

  // Set detector type flag
  dd4hep::xml::setDetectorTypeFlag(x_det, sdet);
  auto& params = DD4hepDetectorHelper::ensureExtension<dd4hep::rec::VariantParameters>(sdet);
  // Add the volume boundary material if configured
  for (xml_coll_t bmat(x_det, _Unicode(boundary_material)); bmat; ++bmat) {
    xml_comp_t x_boundary_material = bmat;
    DD4hepDetectorHelper::xmlToProtoSurfaceMaterial(x_boundary_material, params,
                                                    "boundary_material");
  }

  Assembly assembly(det_name);
  assembly.setVisAttributes(description.invisible());
  sens.setType("tracker");

  // dimensions of the modules (2x2 sensors)
  xml_comp_t x_modsz = x_det.child(_Unicode(modsize));

  double module_x       = x_modsz.length();
  double module_y       = x_modsz.width();
  double module_overlap = getAttrOrDefault(x_modsz, _Unicode(overlap), 0.); // x_modsz.overlap();
  double module_spacing = getAttrOrDefault(x_modsz, _Unicode(spacing), 0.); // x_modsz.overlap();
  double board_gap      = getAttrOrDefault(x_modsz, _Unicode(board_gap), 0.);

  map<string, Assembly> modules;
  map<string, Placements> sensitives;
  map<string, std::vector<VolPlane>> volplane_surfaces;
  map<string, double> mod_thickness;

  for (xml_coll_t mi(x_det, _U(module)); mi; ++mi) {
    xml_comp_t x_mod       = mi;
    double total_thickness = 0;
    // Compute module total thickness from components
    xml_coll_t ci(x_mod, _U(module_component));
    for (ci.reset(); ci; ++ci) {
      xml_comp_t x_comp    = ci;
      bool keep_same_layer = getAttrOrDefault<bool>(x_comp, _Unicode(keep_layer), false);
      if (!keep_same_layer)
        total_thickness += x_comp.thickness();
    }

    double thickness_so_far = 0.0;
    //double thickness_sum        = -total_thickness / 2.0;
    int sensitive_id        = 0;
    const std::string m_nam = x_mod.nameStr();
    mod_thickness[m_nam]    = total_thickness;

    // the module assembly volume
    Assembly m_vol(m_nam);
    m_vol.setVisAttributes(description.visAttributes(x_mod.visStr()));

    int ncomponents = 0;
    for (xml_coll_t mci(x_mod, _U(module_component)); mci; ++mci, ++ncomponents) {
      xml_comp_t x_comp  = mci;
      xml_comp_t x_pos   = x_comp.position(false);
      xml_comp_t x_rot   = x_comp.rotation(false);
      const string c_nam = x_comp.nameStr();

      Box c_box(x_comp.width() / 2, x_comp.length() / 2, x_comp.thickness() / 2);
      Volume c_vol(c_nam, c_box, description.material(x_comp.materialStr()));
      // Utility variable for the relative z-offset based off the previous components
      const double zoff = thickness_so_far + x_comp.thickness() / 2.0;
      if (x_pos && x_rot) {
        Position c_pos(x_pos.x(0), x_pos.y(0), x_pos.z(0) + zoff);
        RotationZYX c_rot(x_rot.z(0), x_rot.y(0), x_rot.x(0));
        pv = m_vol.placeVolume(c_vol, Transform3D(c_rot, c_pos));
      } else if (x_rot) {
        Position c_pos(0, 0, zoff);
        pv = m_vol.placeVolume(c_vol,
                               Transform3D(RotationZYX(x_rot.z(0), x_rot.y(0), x_rot.x(0)), c_pos));
      } else if (x_pos) {
        pv = m_vol.placeVolume(c_vol, Position(x_pos.x(0), x_pos.y(0), x_pos.z(0) + zoff));
      } else {
        pv = m_vol.placeVolume(c_vol, Position(0, 0, zoff));
      }
      c_vol.setRegion(description, x_comp.regionStr());
      c_vol.setLimitSet(description, x_comp.limitsStr());
      c_vol.setVisAttributes(description, x_comp.visStr());
      if (x_comp.isSensitive()) {
        pv.addPhysVolID("ids", sensitive_id);
        ++sensitive_id;

        c_vol.setSensitiveDetector(sens);
        module_thicknesses[m_nam] = {thickness_so_far + x_comp.thickness() / 2.0,
                                     total_thickness - thickness_so_far - x_comp.thickness() / 2.0};
        sensitives[m_nam].push_back(pv);

        // -------- create a measurement plane for the tracking surface attched to the sensitive volume -----
        Vector3D u(-1., 0., 0.);
        Vector3D v(0., -1., 0.);
        Vector3D n(0., 0., 1.);

        // compute the inner and outer thicknesses that need to be assigned to the tracking surface
        // depending on wether the support is above or below the sensor
        double inner_thickness = module_thicknesses[m_nam][0];
        double outer_thickness = module_thicknesses[m_nam][1];

        SurfaceType type(SurfaceType::Sensitive);

        VolPlane surf(c_vol, type, inner_thickness, outer_thickness, u, v, n);
        volplane_surfaces[m_nam].push_back(surf);

        //--------------------------------------------
      }
      bool keep_same_layer = getAttrOrDefault<bool>(x_comp, _Unicode(keep_layer), false);
      if (!keep_same_layer) {
        thickness_so_far += x_comp.thickness();
        // apply relative offsets in z-position used to stack components side-by-side
        if (x_pos) {
          thickness_so_far += x_pos.z(0);
        }
      }
    }
    modules[m_nam] = m_vol;
  }

  int module = 0;
  for (xml_coll_t li(x_det, _U(layer)); li; ++li) {
    xml_comp_t x_layer = li;
    bool front         = x_layer.attr<bool>(_Unicode(front));

    const std::string locStr = x_layer.nameStr();

    // now build the envelope for the detector
    xml_comp_t envelope = x_layer.child(_Unicode(envelope), false);
    int lay_id          = x_layer.id();
    string lay_nam      = det_name + _toString(lay_id, "_layer%d");
    double phimin       = dd4hep::getAttrOrDefault<double>(envelope, _Unicode(phimin), 0.);
    double phimax       = dd4hep::getAttrOrDefault<double>(envelope, _Unicode(phimax), 2 * M_PI);
    double xoffset      = getAttrOrDefault<double>(envelope, _Unicode(xoffset), 0);
    //int mod_num         = 0;

    // envelope thickness is the max layer thickness to accomodate all layers
    double envelope_length = 0;
    for (xml_coll_t llayout(x_layer, _Unicode(layout)); llayout; ++llayout) {
      xml_comp_t x_layout = llayout;
      string m_nam        = x_layout.moduleStr();
      envelope_length     = std::max(envelope_length, mod_thickness[m_nam]);
    }

    Tube lay_tub(envelope.rmin(), envelope.rmax(), envelope_length / 2.0, phimin, phimax);
    Volume lay_vol(lay_nam, lay_tub, air); // Create the layer envelope volume.
    Position lay_pos(xoffset, 0,
                     envelope.zstart() + (front ? -0.5 * envelope_length : 0.5 * envelope_length));
    lay_vol.setVisAttributes(description.visAttributes(x_layer.visStr()));

    DetElement lay_elt(sdet, lay_nam, lay_id);
    // Create the PhysicalVolume for the layer.
    pv = assembly.placeVolume(lay_vol, lay_pos); // Place layer in mother
    pv.addPhysVolID("layer", lay_id);            // Set the layer ID.
    lay_elt.setAttributes(description, lay_vol, x_layer.regionStr(), x_layer.limitsStr(),
                          x_layer.visStr());
    lay_elt.setPlacement(pv);

    // the local coordinate systems of modules in dd4hep and acts differ
    // see http://acts.web.cern.ch/ACTS/latest/doc/group__DD4hepPlugins.html
    auto& layerParams =
        DD4hepDetectorHelper::ensureExtension<dd4hep::rec::VariantParameters>(lay_elt);

    for (xml_coll_t lmat(x_layer, _Unicode(layer_material)); lmat; ++lmat) {
      xml_comp_t x_layer_material = lmat;
      DD4hepDetectorHelper::xmlToProtoSurfaceMaterial(x_layer_material, layerParams,
                                                      "layer_material");
    }

    for (xml_coll_t llayout(x_layer, _Unicode(layout)); llayout; ++llayout, ++module) {
      xml_comp_t x_layout    = llayout;
      bool left              = x_layout.attr<bool>(_Unicode(left));
      string m_nam           = x_layout.moduleStr();
      Volume m_vol           = modules[m_nam];
      double total_thickness = mod_thickness[m_nam];
      Placements& sensVols   = sensitives[m_nam];

      float ycoord = envelope.rmax() -
                     module_y / 2.; // y-center-coord of the top sensor. Start from the top row
      int iy = 0;

      for (xml_coll_t lrow(x_layout, _Unicode(row)); lrow; ++lrow) {
        xml_comp_t x_row = lrow;
        double deadspace = getAttrOrDefault<double>(x_row, _Unicode(deadspace), 0);
        if (deadspace > 0) {
          ycoord -= deadspace;
          continue;
        }
        double x_offset = getAttrOrDefault<double>(x_row, _Unicode(x_offset), 0);
        int nsensors    = getAttrOrDefault<int>(x_row, _Unicode(nsensors), 0);

        // find the sensor id that corrsponds to the rightmost sensor in a board
        // we need to know where to apply additional spaces between neighboring board
        std::unordered_set<int> sensors_id_board_edge;
        int curr_ix = nsensors; // the first sensor to the right of center has ix of nsensors
        for (xml_coll_t lboard(x_row, _Unicode(board)); lboard; ++lboard) {
          xml_comp_t x_board = lboard;
          int nboard_sensors = getAttrOrDefault<int>(x_board, _Unicode(nsensors), 1);
          curr_ix += nboard_sensors;
          sensors_id_board_edge.insert(curr_ix);
          sensors_id_board_edge.insert(2 * nsensors - curr_ix -
                                       1); // reflected to sensor id on the left
        }

        double accum_xoffset = x_offset;

        for (int ix = (left ? nsensors - 1 : nsensors); (ix >= 0) && (ix < 2 * nsensors);
             ix     = ix + (left ? -1 : 1)) {
          // add board spacing
          if (sensors_id_board_edge.find(ix) != sensors_id_board_edge.end())
            accum_xoffset = accum_xoffset + board_gap;

          // there is a hole in the middle, with radius = x_offset
          float xcoord = (ix - nsensors + 0.5) * (module_x + module_spacing) +
                         +(left ? -accum_xoffset : accum_xoffset);
          //! Note the module ordering is different for front and back side

          double module_z = -0.5 * envelope_length;
          if (front)
            module_z = 0.5 * envelope_length - total_thickness;

          // module built!
          Transform3D tr(RotationZYX(M_PI / 2, 0, 0), Position(xcoord, ycoord, module_z));

          pv = lay_vol.placeVolume(m_vol, tr);
          pv.addPhysVolID("idx", ix);
          pv.addPhysVolID("idy", iy);
          pv.addPhysVolID("module", module);

          string comp_nam = Form("%s_%s_ix%d_iy%d", lay_nam.c_str(), m_nam.c_str(), ix, iy);
          DetElement comp_elt(lay_elt, comp_nam, det_id);
          comp_elt.setPlacement(pv);

          for (size_t ic = 0; ic < sensVols.size(); ++ic) {
            PlacedVolume sens_pv = sensVols[ic];
            DetElement sensor_elt(comp_elt, sens_pv.volume().name(), module);
            sensor_elt.setPlacement(sens_pv);
            auto& sensor_elt_params =
                DD4hepDetectorHelper::ensureExtension<dd4hep::rec::VariantParameters>(sensor_elt);
            sensor_elt_params.set<string>("axis_definitions", "XYZ");
            volSurfaceList(sensor_elt)->push_back(volplane_surfaces[m_nam][ic]);
          }
          //++mod_num;
        }
        ycoord -= (module_y - module_overlap);
        ++iy;
      }
    }
  }
  pv = description.pickMotherVolume(sdet).placeVolume(assembly, Position(0, 0, 0));
  pv.addPhysVolID("system", det_id);
  sdet.setPlacement(pv);

  return sdet;
}

DECLARE_DETELEMENT(epic_TOFEndcap, create_detector)
