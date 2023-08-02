// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Chao Peng, Maria Zurek, Whitney Armstrong

// Detector plugin to support a hybrid central barrel calorimeter
// The detector consists of interlayers of Pb/ScFi (segmentation in global r, phi) and W/Si (segmentation in local x, y)
// Assembly is used as the envelope so two different detectors can be interlayered with each other
//
//
// 06/19/2021: Implementation of the Sci Fiber geometry. M. Żurek
// 07/09/2021: Support interlayers between multiple detectors. C. Peng
// 07/23/2021: Add assemblies as mother volumes of fibers to reduce the number of daughter volumes. C. Peng, M. Żurek
//     Reference: TGeo performance issue with large number of daughter volumes
//     https://indico.cern.ch/event/967418/contributions/4075358/attachments/2128099/3583278/201009_shKo_dd4hep.pdf
// 07/24/2021: Changed support implementation to avoid too many uses of boolean geometries. DAWN view seems to have
//     issue dealing with it. C. Peng

#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/Printout.h"
#include "DDRec/Surface.h"
#include "Math/Point2D.h"
#include "TGeoPolygon.h"
#include "XML/Layering.h"
#include <functional>

using namespace dd4hep;

typedef ROOT::Math::XYPoint Point;

// barrel ecal layers contained in an assembly
static Ref_t create_detector(Detector& desc, xml_h e, SensitiveDetector sens)
{
  xml_det_t   x_detector    = e;
  int         detector_id   = x_detector.id();
  std::string detector_name = x_detector.nameStr();
  double      offset        = x_detector.attr<double>(_Unicode(offset));
  xml_comp_t  x_dimensions  = x_detector.dimensions();
  int         nsides        = x_dimensions.numsides();
  double      inner_r       = x_dimensions.rmin();
  double      dphi          = (2 * M_PI / nsides);
  double      half_dphi     = dphi / 2;

  DetElement sdet(detector_name, detector_id);
  Volume     mother_volume = desc.pickMotherVolume(sdet);

  Assembly detector_volume(detector_name);
  detector_volume.setAttributes(desc, x_detector.regionStr(), x_detector.limitsStr(), x_detector.visStr());

  Transform3D  detector_tr = Translation3D(0, 0, offset) * RotationZ(half_dphi);
  PlacedVolume detector_pv = mother_volume.placeVolume(detector_volume, detector_tr);
  sens.setType("calorimeter");

  detector_pv.addPhysVolID("system", detector_id);
  sdet.setPlacement(detector_pv);

  // build a single sector
  DetElement sector_element("sector0", detector_id);
  Assembly   sector_volume("sector");
  if (x_detector.hasChild(_Unicode(sectors))) {
    xml_comp_t x_sectors = x_detector.child(_Unicode(sectors));
    sector_volume.setVisAttributes(desc.visAttributes(x_sectors.visStr()));
  }

  // keep tracking of the total thickness
  double layer_pos_z = inner_r;

  // Parameters for computing the layer X dimension:
  double tan_half_dphi = std::tan(half_dphi);
  double layer_dim_y   = x_dimensions.z() / 2.;

  // Loop over the modules
  PlacedVolume pv;
  using dd4hep::rec::VolPlane;
  std::map<std::string, Volume>                    volumes;
  std::map<std::string, std::vector<PlacedVolume>> sensitives;
  std::map<std::string, std::vector<VolPlane>>     volplane_surfaces;
  std::map<std::string, std::array<double, 2>>     module_thicknesses;
  for (xml_coll_t i_module(x_detector, _U(module)); i_module; ++i_module) {
    xml_comp_t  x_mod = i_module;
    std::string m_nam = x_mod.nameStr();

    if (volumes.find(m_nam) != volumes.end()) {
      printout(ERROR, "BarrelCalorimeterImaging", "Module with name %s already exists", m_nam.c_str());
      throw std::runtime_error("Logics error in building modules.");
    }

    // Compute module total thickness from components
    double     total_thickness = 0;
    xml_coll_t ci(x_mod, _U(module_component));
    for (ci.reset(), total_thickness = 0.0; ci; ++ci) {
      total_thickness += xml_comp_t(ci).thickness();
    }

    // Create the module assembly volume
    Assembly m_vol(m_nam);
    volumes[m_nam] = m_vol;
    m_vol.setVisAttributes(desc.visAttributes(x_mod.visStr()));

    // Add components
    int    sensor_number    = 1;
    int    ncomponents      = 0;
    double thickness_so_far = 0.0;
    double thickness_sum    = -total_thickness / 2.0;
    for (xml_coll_t mci(x_mod, _U(module_component)); mci; ++mci, ++ncomponents) {
      xml_comp_t        x_comp = mci;
      xml_comp_t        x_pos  = x_comp.position(false);
      xml_comp_t        x_rot  = x_comp.rotation(false);
      const std::string c_nam  = _toString(ncomponents, "component%d");
      Box               c_box(x_comp.width() / 2, x_comp.length() / 2, x_comp.thickness() / 2);
      Volume            c_vol(c_nam, c_box, desc.material(x_comp.materialStr()));

      // Utility variable for the relative z-offset based off the previous components
      const double zoff = thickness_sum + x_comp.thickness() / 2.0;
      if (x_pos && x_rot) {
        Position    c_pos(x_pos.x(0), x_pos.y(0), x_pos.z(0) + zoff);
        RotationZYX c_rot(x_rot.z(0), x_rot.y(0), x_rot.x(0));
        pv = m_vol.placeVolume(c_vol, Transform3D(c_rot, c_pos));
      } else if (x_rot) {
        Position c_pos(0, 0, zoff);
        pv = m_vol.placeVolume(c_vol, Transform3D(RotationZYX(x_rot.z(0), x_rot.y(0), x_rot.x(0)), c_pos));
      } else if (x_pos) {
        pv = m_vol.placeVolume(c_vol, Position(x_pos.x(0), x_pos.y(0), x_pos.z(0) + zoff));
      } else {
        pv = m_vol.placeVolume(c_vol, Position(0, 0, zoff));
      }
      c_vol.setRegion(desc, x_comp.regionStr());
      c_vol.setLimitSet(desc, x_comp.limitsStr());
      c_vol.setVisAttributes(desc, x_comp.visStr());
      if (x_comp.isSensitive()) {
        pv.addPhysVolID("sensor", sensor_number++);
        c_vol.setSensitiveDetector(sens);
        sensitives[m_nam].push_back(pv);
        module_thicknesses[m_nam] = {thickness_so_far + x_comp.thickness() / 2.0,
                                     total_thickness - thickness_so_far - x_comp.thickness() / 2.0};
      }
      thickness_sum += x_comp.thickness();
      thickness_so_far += x_comp.thickness();
      // apply relative offsets in z-position used to stack components side-by-side
      if (x_pos) {
        thickness_sum += x_pos.z(0);
        thickness_so_far += x_pos.z(0);
      }
    }
  }

  // Loop over the sets of layer elements in the detector.
  int layer_num = 1;
  for (xml_coll_t i_layer(x_detector, _U(layer)); i_layer; ++i_layer) {
    xml_comp_t x_layer             = i_layer;
    int        layer_repeat        = x_layer.repeat();
    double     layer_thickness     = x_layer.thickness();
    bool       layer_has_frame     = x_layer.hasChild(_Unicode(frame));
    double     layer_space_between = getAttrOrDefault(x_layer, _Unicode(space_between), 0.);
    double     layer_space_before  = getAttrOrDefault(x_layer, _Unicode(space_before), 0.);
    layer_pos_z += layer_space_before;

    // Loop over number of repeats for this layer.
    for (int layer_j = 0; layer_j < layer_repeat; layer_j++) {

      // Make an envelope for this layer
      std::string layer_name  = Form("layer%d", layer_num);
      double      layer_dim_x = tan_half_dphi * layer_pos_z;
      auto        layer_mat   = desc.air();

      Position   layer_pos(0, 0, layer_pos_z + layer_thickness / 2.);
      double     layer_trd_x1 = layer_dim_x;
      double     layer_trd_x2 = layer_dim_x + layer_thickness * tan_half_dphi;
      double     layer_trd_y1 = layer_dim_y;
      double     layer_trd_y2 = layer_trd_y1;
      double     layer_trd_z  = layer_thickness / 2; // account for frame
      Trapezoid  layer_shape(layer_trd_x1, layer_trd_x2, layer_trd_y1, layer_trd_y2, layer_trd_z);
      Volume     layer_volume(layer_name, layer_shape, layer_mat);
      DetElement layer_element(sector_element, layer_name, detector_id);

      // Set region, limitset, and vis of layer.
      layer_volume.setAttributes(desc, x_layer.regionStr(), x_layer.limitsStr(), x_layer.visStr());

      // Loop over the staves for this layer.
      int stave_num = 1;
      for (xml_coll_t i_stave(x_layer, _U(stave)); i_stave; ++i_stave) {
        xml_comp_t x_stave      = i_stave;
        int        stave_repeat = x_stave.repeat();
        double     stave_thick  = x_stave.thickness();
        double     stave_dim_x  = x_stave.width() / 2.0;
        double     stave_dim_y  = x_stave.length() / 2.0;
        double     stave_dim_z  = stave_thick / 2.0;
        double     stave_rot_y  = getAttrOrDefault(x_stave, _Unicode(angle), 0.);
        double     stave_off_z  = getAttrOrDefault(x_stave, _Unicode(offset), 0.);

        // Optional x_layout and z_layout specification
        [[maybe_unused]] bool has_x_layout = x_stave.hasChild(_U(z_layout));
        [[maybe_unused]] bool has_z_layout = x_stave.hasChild(_U(z_layout));

        // Arrange staves symmetrically around center of layer
        double stave_pos_x = -layer_dim_x + stave_dim_x;
        double stave_pitch = -2.0 * stave_pos_x / (stave_repeat - 1);

        // Make one stave
        std::string stave_name     = Form("stave%d", stave_num);
        auto        stave_material = desc.air();
        Box         stave_shape(stave_dim_x, stave_dim_y, stave_dim_z);
        Volume      stave_volume(stave_name, stave_shape, stave_material);
        stave_volume.setAttributes(desc, x_stave.regionStr(), x_stave.limitsStr(), x_stave.visStr());
        DetElement stave_element(layer_element, stave_name, detector_id);

        // Loop over the slices for this stave
        double slice_pos_z = -(stave_thick / 2.);

        // If x_layout and z_layout
        if (has_x_layout && has_z_layout) {
          auto x_layout      = x_stave.child(_Unicode(x_layout));
          auto z_layout      = x_stave.child(_Unicode(x_layout));
          auto module_str    = x_stave.moduleStr();
          auto module_volume = volumes[module_str];
          int  i_module      = 0;
          for (auto i_x = 0; i_x < x_layout.nx(); ++i_x) {
            for (auto i_z = 0; i_z < z_layout.nz(); ++i_z) {
              std::string module_name = _toString(i_module, "module%d");
              DetElement  module_element(stave_element, module_name, i_module);

              // Place module
              double x = x_layout.x0() + i_x * x_layout.dx();
              double z = z_layout.z0() + i_z * z_layout.dz();

              PlacedVolume module_pv = stave_volume.placeVolume(module_volume, Position(x, 0, z));
              module_pv.addPhysVolID("module", module);
              module.setPlacement(module_pv);

              i_module++;
            }
          }
          slice_pos_z += module_thicknesses[module_name].first;
        }

        // Loop over the slices for this stave
        int slice_num = 1;
        for (xml_coll_t i_slice(x_stave, _U(slice)); i_slice; ++i_slice) {
          xml_comp_t  x_slice     = i_slice;
          std::string slice_name  = Form("slice%d", slice_num);
          double      slice_thick = x_slice.thickness();
          double      slice_dim_x = stave_dim_x;
          double      slice_dim_y = stave_dim_y;
          double      slice_dim_z = slice_thick / 2.;
          Box         slice_shape(slice_dim_x, slice_dim_y, slice_dim_z);
          Volume      slice_volume(slice_name, slice_shape, desc.material(x_slice.materialStr()));
          slice_volume.setAttributes(desc, x_slice.regionStr(), x_slice.limitsStr(), x_slice.visStr());
          DetElement slice_element(stave_element, slice_name, detector_id);

          // Set sensitive
          if (x_slice.isSensitive()) {
            slice_volume.setSensitiveDetector(sens);
          }

          // Place slice
          PlacedVolume slice_pv = stave_volume.placeVolume(slice_volume, Position(0, 0, slice_pos_z + slice_thick / 2));
          slice_pv.addPhysVolID("slice", slice_num);
          slice.setPlacement(slice_pv);

          // Increment Z position of slice
          slice_pos_z += slice_thick;
          ++slice_num;
        }

        // Loop over number of repeats for this stave
        for (int stave_j = 0; stave_j < stave_repeat; stave_j++) {

          // Stave placement
          Position     stave_pos(stave_pos_x, 0, layer_j % 2 == 0 ? +stave_off_z : -stave_off_z);
          RotationY    stave_rot(layer_j % 2 == 0 ? +stave_rot_y : -stave_rot_y);
          Transform3D  stave_tr(stave_rot, stave_pos);
          PlacedVolume stave_pv = layer_volume.placeVolume(stave_volume, stave_tr);
          stave_pv.addPhysVolID("stave", stave_num);
          stave.setPlacement(stave_pv);

          // Increment X position of stave
          stave_pos_x += stave_pitch;
          ++stave_num;
        }
      }

      // Place frame
      if (layer_has_frame) {
        auto   x_frame         = x_layer.child(_Unicode(frame));
        double frame_height    = x_frame.attr<double>(_Unicode(height));
        double frame_thickness = x_frame.attr<double>(_Unicode(thickness));
        auto   frame_material  = desc.material(x_frame.attr<std::string>(_Unicode(material)));

        std::string frame1_name   = Form("frame_inner%d", layer_num);
        double      frame1_thick  = frame_thickness;
        double      frame1_trd_x1 = layer_dim_x + (layer_thickness / 2 - frame_height / 2) * tan_half_dphi;
        double    frame1_trd_x2 = layer_dim_x + (layer_thickness / 2 - frame_height / 2 + frame1_thick) * tan_half_dphi;
        double    frame1_trd_y1 = layer_trd_y1;
        double    frame1_trd_y2 = frame1_trd_y1;
        double    frame1_trd_z  = frame1_thick / 2.;
        Trapezoid frame1_shape(frame1_trd_x1, frame1_trd_x2, frame1_trd_y1, frame1_trd_y2, frame1_trd_z);
        Volume    frame1_volume(frame1_name, frame1_shape, frame_material);
        layer_volume.placeVolume(frame1_volume, Position(0, 0, -frame_height / 2.0 - frame_thickness / 2.0));

        std::string frame2_name   = Form("frame_outer%d", layer_num);
        double      frame2_thick  = frame_thickness;
        double      frame2_trd_x1 = layer_dim_x + (layer_thickness / 2 - frame_height / 2) * tan_half_dphi;
        double    frame2_trd_x2 = layer_dim_x + (layer_thickness / 2 - frame_height / 2 + frame2_thick) * tan_half_dphi;
        double    frame2_trd_y1 = layer_trd_y1;
        double    frame2_trd_y2 = frame2_trd_y1;
        double    frame2_trd_z  = frame2_thick / 2.;
        Trapezoid frame2_shape(frame2_trd_x1, frame2_trd_x2, frame2_trd_y1, frame2_trd_y2, frame2_trd_z);
        Volume    frame2_volume(frame2_name, frame2_shape, frame_material);
        layer_volume.placeVolume(frame2_volume, Position(0, 0, +frame_height / 2.0 + frame_thickness / 2.0));
      }

      // Place layer into sector
      PlacedVolume layer_pv = sector_volume.placeVolume(layer_volume, layer_pos);
      layer_pv.addPhysVolID("layer", layer_num);
      layer.setPlacement(layer_pv);

      // Increment to next layer Z position. Do not add space_between for the last layer
      layer_pos_z += layer_thickness;
      if (layer_j < layer_repeat - 1) {
        layer_pos_z += layer_space_between;
      }
      ++layer_num;
    }
  }

  // Phi start for a sector.
  double phi = M_PI / nsides;
  // Create nsides sectors.
  for (int i = 0; i < nsides; i++, phi -= dphi) { // i is sector number
    // Compute the sector position
    Transform3D  tr(RotationZYX(0, phi, M_PI * 0.5), Translation3D(0, 0, 0));
    PlacedVolume detector_pv = detector_volume.placeVolume(sector_volume, tr);
    detector_pv.addPhysVolID("sector", i + 1);
    DetElement sd = (i == 0) ? sector_element : sector_element.clone(Form("sector%d", i));
    sd.setPlacement(detector_pv);
    sdet.add(sd);
  }

  return sdet;
}

DECLARE_DETELEMENT(epic_EcalBarrelImaging, create_detector)
