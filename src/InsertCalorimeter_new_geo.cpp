// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Ryan Milton
// Copyright (C) 2026 Aiden Wu
// Modified in 2026 by Aiden Wu from the original implementation by Ryan Milton.

//==========================================================================
//  Implementation of forward insert calorimeter
//--------------------------------------------------------------------------
//  Author: Ryan Milton (UCR)
//  Insert Redesign: Aiden Wu (ORNL)
//==========================================================================
//  TODO: Add an LFHCAL-style module casing around the insert, add 5x5 cm tiles, add actual pcb boards.

#include "DD4hep/DetFactoryHelper.h"
#include <XML/Helper.h>
#include <XML/Layering.h>
#include <XML/Utilities.h>

using namespace dd4hep;

static Ref_t createDetector(Detector& desc, xml_h handle, SensitiveDetector sens) {
  xml_det_t detElem   = handle;
  std::string detName = detElem.nameStr();
  int detID           = detElem.id();

  xml_dim_t dim = detElem.dimensions();
  double width  = dim.x(); // Size along x-axis
  double height = dim.y(); // Size along y-axis
  double length = dim.z(); // Size along z-axis

  xml_dim_t pos = detElem.position(); // Position in global coordinates

  Material air = desc.material("Air");

  // Read the dimensions and positions that define the central cutout.
  const xml::Component& cutout_xml = detElem.child(_Unicode(beampipe_hole));

  // Shared circle center in the insert's local coordinates.
  const double center_x = dd4hep::getAttrOrDefault<double>(cutout_xml, _Unicode(center_x), -10.0 * cm) - pos.x();
  const double center_y = dd4hep::getAttrOrDefault<double>(cutout_xml, _Unicode(center_y), 0. * cm) - pos.y();
  // Circle radius for each side of the insert.
  const double left_radius = dd4hep::getAttrOrDefault<double>(cutout_xml, _Unicode(left_radius), 16.781 * cm);
  const double right_radius = dd4hep::getAttrOrDefault<double>(cutout_xml, _Unicode(right_radius), 16.786 * cm);

  // Horizontal separation between the independently constructed left and right halves.
  const double left_right_gap = dd4hep::getAttrOrDefault<double>(detElem, _Unicode(left_right_gap), 0.38 * cm);

  // Width, height, and center of the rectangular section joined to the right circle.
  const double rect_width = dd4hep::getAttrOrDefault<double>(cutout_xml, _Unicode(rect_width), 18.84 * cm);
  const double rect_height = dd4hep::getAttrOrDefault<double>(cutout_xml, _Unicode(rect_height), 33.41 * cm);
  const double rect_x = -left_right_gap / 2. - pos.x() - center_x - rect_width / 2.;

  // Width, height, and vertical center of the horizontal PCB slots.
  const double h_pcb_slot_width = dd4hep::getAttrOrDefault<double>(cutout_xml, _Unicode(h_pcb_slot_width), 15.0 * cm);
  const double h_pcb_slot_height = dd4hep::getAttrOrDefault<double>(cutout_xml, _Unicode(h_pcb_slot_height), 0.5 * cm);
  const double h_pcb_slot_y = (rect_height + h_pcb_slot_height) / 2.;

  // Width, side-specific heights, and edge clearance of the vertical PCB slots.
  const double v_pcb_slot_width = dd4hep::getAttrOrDefault<double>(cutout_xml, _Unicode(v_pcb_slot_width), 0.5 * cm);
  const double right_v_pcb_slot_height = dd4hep::getAttrOrDefault<double>(cutout_xml, _Unicode(right_v_pcb_slot_height), 11.0 * cm);
  const double left_v_pcb_slot_height = dd4hep::getAttrOrDefault<double>(cutout_xml, _Unicode(left_v_pcb_slot_height), 14.5 * cm);
  const double v_pcb_slot_offset_y = dd4hep::getAttrOrDefault<double>(cutout_xml, _Unicode(v_pcb_slot_offset_y), 1.0 * cm);
  const double right_v_pcb_slot_y = height / 2. - v_pcb_slot_offset_y - right_v_pcb_slot_height / 2.;
  const double left_v_pcb_slot_y = height / 2. - v_pcb_slot_offset_y - left_v_pcb_slot_height / 2.;

  // Horizontal centers of the vertical PCB slots relative to the shared circle center.
  const double right_v_pcb_slot_x = -left_right_gap / 2. - pos.x() - v_pcb_slot_width / 2. - center_x;
  const double left_v_pcb_slot_x = left_right_gap / 2. - pos.x() + v_pcb_slot_width / 2. - center_x;

  // Size the longitudinal PCBs with the same face and edge clearances as LFHCAL.
  const double pcb_length = dd4hep::getAttrOrDefault<double>(cutout_xml, _Unicode(pcb_length), 125.0 * cm);
  const double pcb_face_clearance = dd4hep::getAttrOrDefault<double>(cutout_xml, _Unicode(pcb_face_clearance), 0.05 * cm);
  const double pcb_edge_clearance = dd4hep::getAttrOrDefault<double>(cutout_xml, _Unicode(pcb_edge_clearance), 0.4 * cm);
  const double h_pcb_thickness = h_pcb_slot_height - 2. * pcb_face_clearance;
  const double v_pcb_thickness = v_pcb_slot_width - 2. * pcb_face_clearance;
  const double h_pcb_width = h_pcb_slot_width - 2. * pcb_edge_clearance;
  const double right_v_pcb_height = right_v_pcb_slot_height - 2. * pcb_edge_clearance;
  const double left_v_pcb_height = left_v_pcb_slot_height - 2. * pcb_edge_clearance;
  // Align both longitudinal PCB ends with the active layer stack.
  const double pcb_z = -(length - pcb_length) / 2.;

  // Assembly that will contain all the layers
  Assembly assembly(detName);
  // FIXME Workaround for https://github.com/eic/epic/issues/411
  assembly.setVisAttributes(desc.visAttributes("InvisibleWithDaughters"));
  PlacedVolume pv;
  for (int side_num = 0; side_num < 2; side_num++) { // 0 = right, 1 = left
    std::string side_name = side_num == 1 ? "L" : "R";
    // Keeps track of the z location as we move longiduinally through the insert
    double z_distance_traversed = 0.;

    int layer_num = 1;

    // Looping through all the different layer sections (W/Sc, Steel/Sc, backplate)
    for (xml_coll_t c(detElem, _U(layer)); c; c++) {
      xml_comp_t x_layer     = c;
      int repeat             = x_layer.repeat();
      double layer_thickness = x_layer.thickness();

      // Looping through the number of repeated layers in each section
      for (int i = 0; i < repeat; i++) {
        std::string layer_name = detName + _toString(layer_num, "_layer%d") + "_" + side_name;
        Box layer(width / 2., height / 2., layer_thickness / 2.);

        // Removing beampipe shape from each layer
        Solid layer_cutout = Tube(0., left_radius, layer_thickness / 2.);

        // Right side.
        if (side_num == 0) {
          // Beampipe cutout.
          Box rect_cutout(rect_width / 2., rect_height / 2., layer_thickness / 2.);
          Tube right_circle(0., right_radius, layer_thickness / 2.);
          layer_cutout = UnionSolid(right_circle, rect_cutout, Position(rect_x, 0., 0.));

          // Horizontal PCB slots.
          Box h_pcb_slot(h_pcb_slot_width / 2., h_pcb_slot_height / 2., layer_thickness / 2.);
          layer_cutout = UnionSolid(layer_cutout, h_pcb_slot, Position(rect_x, h_pcb_slot_y, 0.));
          layer_cutout = UnionSolid(layer_cutout, h_pcb_slot, Position(rect_x, -h_pcb_slot_y, 0.));

          // Vertical PCB slots.
          Box v_pcb_slot(v_pcb_slot_width / 2., right_v_pcb_slot_height / 2., layer_thickness / 2.);
          layer_cutout = UnionSolid(layer_cutout, v_pcb_slot, Position(right_v_pcb_slot_x, right_v_pcb_slot_y, 0.));
          layer_cutout = UnionSolid(layer_cutout, v_pcb_slot, Position(right_v_pcb_slot_x, -right_v_pcb_slot_y, 0.));
        }

        // Left side.
        else {
          // Vertical PCB slots.
          Box v_pcb_slot(v_pcb_slot_width / 2., left_v_pcb_slot_height / 2., layer_thickness / 2.);
          layer_cutout = UnionSolid(layer_cutout, v_pcb_slot, Position(left_v_pcb_slot_x, left_v_pcb_slot_y, 0.));
          layer_cutout = UnionSolid(layer_cutout, v_pcb_slot, Position(left_v_pcb_slot_x, -left_v_pcb_slot_y, 0.));
        }
        SubtractionSolid layer_with_cutout(layer, layer_cutout, Position(center_x, center_y, 0.));

        // Only select the left or right side of the layer
        Box side_cut(width / 2., height, layer_thickness);
        Position side_cut_position((width / 2 - left_right_gap / 2) * (1 - 2 * side_num) - pos.x(), 0, 0);
        SubtractionSolid layer_side_with_cutout(layer_with_cutout, side_cut, side_cut_position);
        Volume layer_vol(layer_name, layer_side_with_cutout, air);

        int slice_num  = 1;
        double slice_z = -layer_thickness / 2.; // Keeps track of slices' z locations in each layer

        // Looping over each layer's slices
        for (xml_coll_t l(x_layer, _U(slice)); l; l++) {
          xml_comp_t x_slice     = l;
          double slice_thickness = x_slice.thickness();
          std::string slice_name = layer_name + _toString(slice_num, "slice%d");
          Material slice_mat     = desc.material(x_slice.materialStr());
          slice_z += slice_thickness / 2.; // Going to slice halfway point

          // Each slice within a layer has the same cutout dimensions and position
          Box slice(width / 2., height / 2., slice_thickness / 2.);
          Solid slice_cutout = Tube(0., left_radius, slice_thickness / 2.);

          // Right side.
          if (side_num == 0) {
            // Beampipe cutout.
            Box rect_cutout(rect_width / 2., rect_height / 2., slice_thickness / 2.);
            Tube right_circle(0., right_radius, slice_thickness / 2.);
            slice_cutout = UnionSolid(right_circle, rect_cutout, Position(rect_x, 0., 0.));

            // Horizontal PCB slots.
            Box h_pcb_slot(h_pcb_slot_width / 2., h_pcb_slot_height / 2., slice_thickness / 2.);
            slice_cutout = UnionSolid(slice_cutout, h_pcb_slot, Position(rect_x, h_pcb_slot_y, 0.));
            slice_cutout = UnionSolid(slice_cutout, h_pcb_slot, Position(rect_x, -h_pcb_slot_y, 0.));

            // Vertical PCB slots.
            Box v_pcb_slot(v_pcb_slot_width / 2., right_v_pcb_slot_height / 2., slice_thickness / 2.);
            slice_cutout = UnionSolid(slice_cutout, v_pcb_slot, Position(right_v_pcb_slot_x, right_v_pcb_slot_y, 0.));
            slice_cutout = UnionSolid(slice_cutout, v_pcb_slot, Position(right_v_pcb_slot_x, -right_v_pcb_slot_y, 0.));
          }

          // Left side.
          else {
            // Vertical PCB slots.
            Box v_pcb_slot(v_pcb_slot_width / 2., left_v_pcb_slot_height / 2., slice_thickness / 2.);
            slice_cutout = UnionSolid(slice_cutout, v_pcb_slot, Position(left_v_pcb_slot_x, left_v_pcb_slot_y, 0.));
            slice_cutout = UnionSolid(slice_cutout, v_pcb_slot, Position(left_v_pcb_slot_x, -left_v_pcb_slot_y, 0.));
          }

          SubtractionSolid slice_with_cutout(slice, slice_cutout, Position(center_x, center_y, 0.));
          Box side_cut_slice(width / 2., height, layer_thickness);
          Position side_cut_position_slice((width / 2 - left_right_gap / 2) * (1 - 2 * side_num) - pos.x(), 0, 0);
          SubtractionSolid slice_side_with_cutout(slice_with_cutout, side_cut_slice, side_cut_position_slice);
          Volume slice_vol(slice_name, slice_side_with_cutout, slice_mat);

          // Setting appropriate slices as sensitive
          if (x_slice.isSensitive()) {
            sens.setType("calorimeter");
            slice_vol.setSensitiveDetector(sens);
          }

          // Setting slice attributes
          slice_vol.setAttributes(desc, x_slice.regionStr(), x_slice.limitsStr(), x_slice.visStr());

          // Placing slice within layer
          pv = layer_vol.placeVolume(slice_vol,
                                     Transform3D(RotationZYX(0, 0, 0), Position(0., 0., slice_z)));
          pv.addPhysVolID("slice", slice_num);
          pv.addPhysVolID("side", side_num);
          slice_z += slice_thickness / 2.;
          z_distance_traversed += slice_thickness;
          slice_num++;
        }

        // Setting layer attributes
        layer_vol.setAttributes(desc, x_layer.regionStr(), x_layer.limitsStr(), x_layer.visStr());
        /*
          Placing each layer inside assembly
          -length/2. is front of detector in global coordinate system
          + (z_distance_traversed - layer_thickness) goes to the front of each layer
          + layer_thickness/2. places layer in correct spot
          Example: After placement of slices in first layer, z_distance_traversed = layer_thickness
                   Subtracting layer_thickness goes back to the front of the first slice (Now, z = -length/2)
                   Adding layer_thickness / 2. goes to half the first layer thickness (proper place to put layer)
                   Each loop over repeat will increases z_distance_traversed by layer_thickness
        */
        pv = assembly.placeVolume(
            layer_vol,
            Transform3D(RotationZYX(0, 0, 0),
                        Position(0., 0.,
                                 -length / 2. + (z_distance_traversed - layer_thickness) +
                                     layer_thickness / 2.)));

        pv.addPhysVolID("layer", layer_num);
        pv.addPhysVolID("side", side_num);
        layer_num++;
      }
    }
  }

  // Place two horizontal and four vertical readout PCBs inside their longitudinal slots.
  Material pcb_material = desc.material("Fr4");
  Box h_pcb(h_pcb_width / 2., h_pcb_thickness / 2., pcb_length / 2.);
  Volume h_pcb_vol(detName + "_HorizontalPCB", h_pcb, pcb_material);
  h_pcb_vol.setVisAttributes(desc.visAttributes("AnlDarkGreen"));
  assembly.placeVolume(h_pcb_vol, Position(center_x + rect_x, center_y + h_pcb_slot_y, pcb_z));
  assembly.placeVolume(h_pcb_vol, Position(center_x + rect_x, center_y - h_pcb_slot_y, pcb_z));

  Box right_v_pcb(v_pcb_thickness / 2., right_v_pcb_height / 2., pcb_length / 2.);
  Volume right_v_pcb_vol(detName + "_RightVerticalPCB", right_v_pcb, pcb_material);
  right_v_pcb_vol.setVisAttributes(desc.visAttributes("AnlDarkGreen"));
  assembly.placeVolume(
      right_v_pcb_vol,
      Position(center_x + right_v_pcb_slot_x, center_y + right_v_pcb_slot_y, pcb_z));
  assembly.placeVolume(
      right_v_pcb_vol,
      Position(center_x + right_v_pcb_slot_x, center_y - right_v_pcb_slot_y, pcb_z));

  Box left_v_pcb(v_pcb_thickness / 2., left_v_pcb_height / 2., pcb_length / 2.);
  Volume left_v_pcb_vol(detName + "_LeftVerticalPCB", left_v_pcb, pcb_material);
  left_v_pcb_vol.setVisAttributes(desc.visAttributes("AnlDarkGreen"));
  assembly.placeVolume(
      left_v_pcb_vol,
      Position(center_x + left_v_pcb_slot_x, center_y + left_v_pcb_slot_y, pcb_z));
  assembly.placeVolume(
      left_v_pcb_vol,
      Position(center_x + left_v_pcb_slot_x, center_y - left_v_pcb_slot_y, pcb_z));

  DetElement det(detName, detID);
  Volume motherVol = desc.pickMotherVolume(det);

  // apply any detector type flags set in XML
  dd4hep::xml::setDetectorTypeFlag(detElem, det);

  // Placing insert in world volume
  auto tr          = Transform3D(Position(pos.x(), pos.y(), pos.z() + length / 2.));
  PlacedVolume phv = motherVol.placeVolume(assembly, tr);
  phv.addPhysVolID("system", detID);
  det.setPlacement(phv);

  return det;
}
DECLARE_DETELEMENT(epic_InsertCalorimeterNew, createDetector)
