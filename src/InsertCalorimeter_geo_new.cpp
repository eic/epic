// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Ryan Milton
// Copyright (C) 2026 Aiden Wu
// Updated geometry dimensions provided by Eliott Fountain.
/*
==========================================================================
 Implementation of forward insert calorimeter
--------------------------------------------------------------------------
 Author: Ryan Milton (UCR)
 Insert Redesign: Aiden Wu (ORNL (NGSI))
==========================================================================
Cell ID bit layout:
  [7:0]   system   [8]     side     [16:9]  layer
  [23:17] slice    [31:24] unused   [47:32] x
  [63:48] y
column_id is decoded from x; row_id is decoded from y.
Columns increase from each outer edge toward the center; rows increase top to bottom.
Insert-local tile-center coordinate mapping:
  Due to the Insert's orientation along z, the front face's left side lies on +x
  and its right side lies on -x.
  right_x0 and left_x0 are the first outer tile-center x positions.
  top_y0 is the first row's tile-center y position.
  Both sides start at column x=0 and move inward.

  Using that logic:
    right_x0 = -width/2 + right_tile_outer_margin + tile_size/2
    left_x0 = width/2 - tile_gap - tile_size/2
    top_y0 = 5.5*tile_pitch

    x_center (right side) = right_x0 + column_id*tile_pitch
    x_center (left side) = left_x0 - column_id*tile_pitch
    y_center = top_y0 - row_id*tile_pitch
==========================================================================
*/

#include "DD4hep/DetFactoryHelper.h"
#include <XML/Helper.h>
#include <XML/Utilities.h>
#include <array>
#include <cmath>

using namespace dd4hep;

static Ref_t createDetector(Detector& desc, xml_h handle, SensitiveDetector sens) {
  // Detector identity.
  xml_det_t detElem   = handle;
  std::string detName = detElem.nameStr();
  int detID           = detElem.id();
  sens.setType("calorimeter");

  // Detector envelope.
  xml_dim_t dim = detElem.dimensions();
  double width  = dim.x();
  double height = dim.y();
  double length = dim.z();

  // World position and container material.
  xml_dim_t pos = detElem.position();
  Material air  = desc.material("Air");

  // Insert component configuration.
  const xml::Component& beam_xml         = detElem.child(_Unicode(beampipe_cutout));
  const xml::Component& slots_xml        = detElem.child(_Unicode(pcb_slots));
  const xml::Component& pcb_xml          = detElem.child(_Unicode(pcb_boards));
  const xml::Component& tiles_xml        = detElem.child(_Unicode(tiles));
  const xml::Component& casing_xml       = detElem.child(_Unicode(casing));
  const xml::Component& cover_xml        = detElem.child(_Unicode(pcb_covers));
  const xml::Component& back_cutouts_xml = detElem.child(_Unicode(backplate_cutouts));

  // Circle center in insert coordinates.
  const double center_x =
      dd4hep::getAttrOrDefault<double>(beam_xml, _Unicode(center_x), -10.0 * cm) - pos.x();
  const double center_y =
      dd4hep::getAttrOrDefault<double>(beam_xml, _Unicode(center_y), 0. * cm) - pos.y();

  // Circle radii.
  const double left_radius =
      dd4hep::getAttrOrDefault<double>(beam_xml, _Unicode(left_radius), 16.781 * cm);
  const double right_radius =
      dd4hep::getAttrOrDefault<double>(beam_xml, _Unicode(right_radius), 16.786 * cm);

  // Gap between insert sides.
  const double left_right_gap =
      dd4hep::getAttrOrDefault<double>(detElem, _Unicode(left_right_gap), 0.38 * cm);
  const double right_split_x = -left_right_gap / 2. - pos.x();
  const double left_split_x  = left_right_gap / 2. - pos.x();

  // Right-side rectangular cutout.
  const double rect_width =
      dd4hep::getAttrOrDefault<double>(beam_xml, _Unicode(rect_width), 18.84 * cm);
  const double rect_height =
      dd4hep::getAttrOrDefault<double>(beam_xml, _Unicode(rect_height), 33.41 * cm);
  const double rect_x = right_split_x - center_x - rect_width / 2.;

  // Horizontal PCB slots.
  const double h_pcb_slot_width =
      dd4hep::getAttrOrDefault<double>(slots_xml, _Unicode(horizontal_width), 15.0 * cm);
  const double h_pcb_slot_height =
      dd4hep::getAttrOrDefault<double>(slots_xml, _Unicode(horizontal_height), 0.5 * cm);
  const double h_pcb_slot_y = (rect_height + h_pcb_slot_height) / 2.;

  // Vertical PCB slots.
  const double v_pcb_slot_width =
      dd4hep::getAttrOrDefault<double>(slots_xml, _Unicode(vertical_width), 0.5 * cm);
  const double right_v_pcb_slot_height =
      dd4hep::getAttrOrDefault<double>(slots_xml, _Unicode(right_vertical_height), 11.0 * cm);
  const double left_v_pcb_slot_height =
      dd4hep::getAttrOrDefault<double>(slots_xml, _Unicode(left_vertical_height), 14.5 * cm);
  const double v_pcb_slot_offset_y =
      dd4hep::getAttrOrDefault<double>(slots_xml, _Unicode(vertical_offset_y), 1.0 * cm);
  const double right_v_pcb_slot_y =
      height / 2. - v_pcb_slot_offset_y - right_v_pcb_slot_height / 2.;
  const double left_v_pcb_slot_y = height / 2. - v_pcb_slot_offset_y - left_v_pcb_slot_height / 2.;

  // Vertical PCB slot x positions.
  const double right_v_pcb_slot_x = right_split_x - v_pcb_slot_width / 2. - center_x;
  const double left_v_pcb_slot_x  = left_split_x + v_pcb_slot_width / 2. - center_x;

  // PCB length, thickness, and slot clearance.
  const double pcb_length = dd4hep::getAttrOrDefault<double>(pcb_xml, _Unicode(length), 125.0 * cm);
  const double pcb_thickness =
      dd4hep::getAttrOrDefault<double>(pcb_xml, _Unicode(thickness), 0.16 * cm);
  const double pcb_edge_clearance =
      dd4hep::getAttrOrDefault<double>(pcb_xml, _Unicode(edge_clearance), 0.4 * cm);

  // Steel casing dimensions.
  const double casing_length =
      dd4hep::getAttrOrDefault<double>(casing_xml, _Unicode(length), 132.0 * cm);
  const double casing_thickness =
      dd4hep::getAttrOrDefault<double>(casing_xml, _Unicode(thickness), 0.34163 * cm);
  const double casing_left_radius =
      dd4hep::getAttrOrDefault<double>(casing_xml, _Unicode(left_radius), 16.781 * cm);
  const double cover_thickness =
      dd4hep::getAttrOrDefault<double>(cover_xml, _Unicode(thickness), 0.19844 * cm);
  const double right_back_cutout_margin =
      dd4hep::getAttrOrDefault<double>(back_cutouts_xml, _Unicode(right_margin), 4.0 * cm);
  const double left_back_cutout_margin =
      dd4hep::getAttrOrDefault<double>(back_cutouts_xml, _Unicode(left_margin), 3.5 * cm);
  const double back_cutout_vertical_margin =
      dd4hep::getAttrOrDefault<double>(back_cutouts_xml, _Unicode(vertical_margin), 3.5 * cm);
  const double right_back_cutout_width =
      width / 2. - left_right_gap / 2. - pos.x() - 2. * right_back_cutout_margin;
  const double left_back_cutout_width =
      width / 2. - left_right_gap / 2. + pos.x() - 2. * left_back_cutout_margin;
  const double back_cutout_height =
      height / 2. - center_y - rect_height / 2. - 2. * back_cutout_vertical_margin;

  // Shorten each PCB along its slot by the edge clearance at both ends.
  const double h_pcb_width        = h_pcb_slot_width - 2. * pcb_edge_clearance;
  const double right_v_pcb_height = right_v_pcb_slot_height - 2. * pcb_edge_clearance;
  const double left_v_pcb_height  = left_v_pcb_slot_height - 2. * pcb_edge_clearance;

  // Align the PCB front with the casing's inner face.
  const double pcb_z = -(length - pcb_length) / 2. + casing_thickness;

  // Tile size, gap, and pitch.
  const double tile_size  = dd4hep::getAttrOrDefault<double>(tiles_xml, _Unicode(size), 4.7 * cm);
  const double tile_gap   = dd4hep::getAttrOrDefault<double>(tiles_xml, _Unicode(gap), 0.02 * cm);
  const double tile_pitch = tile_size + tile_gap;

  // Make the gap from the right tile grid to the central PCB slot four times its outer-edge gap.
  const double right_tile_outer_margin =
      (width / 2. + center_x + right_v_pcb_slot_x - v_pcb_slot_width / 2. -
       (8. * tile_size + 7. * tile_gap)) /
      5.;
  const double left_x0  = width / 2. - tile_gap - tile_size / 2.;
  const double right_x0 = -width / 2. + right_tile_outer_margin + tile_size / 2.;
  const double top_y0   = 5.5 * tile_pitch;

  // Tiles per row, top to bottom.
  const std::array<int, 12> left_tiles_per_row  = {4, 4, 4, 3, 3, 2, 2, 3, 3, 4, 4, 4};
  const std::array<int, 12> right_tiles_per_row = {8, 8, 4, 3, 2, 2, 2, 2, 3, 4, 8, 8};

  // Place left-side tiles from the +x outer edge inward and from top to bottom.
  auto left_tile_position = [left_x0, top_y0, tile_pitch](int row, int column) {
    double x = left_x0 - column * tile_pitch;
    double y = top_y0 - row * tile_pitch;
    return Position(x, y, 0.);
  };

  // Place right-side tiles from the -x outer edge inward and from top to bottom.
  auto right_tile_position = [right_x0, top_y0, tile_pitch](int row, int column) {
    double x = right_x0 + column * tile_pitch;
    double y = top_y0 - row * tile_pitch;
    return Position(x, y, 0.);
  };

  // Build the beampipe hole and PCB slots for the layer container and its passive slices.
  auto make_cutout = [&](double thickness, int side) {
    Solid cutout = Tube(0., left_radius, thickness / 2.);

    // Add the right-side rectangle and PCB slots to its circular cutout.
    if (side == 0) {
      Box rect_cutout(rect_width / 2., rect_height / 2., thickness / 2.);
      Tube right_circle(0., right_radius, thickness / 2.);
      cutout = UnionSolid(right_circle, rect_cutout, Position(rect_x, 0., 0.));

      Box h_pcb_slot(h_pcb_slot_width / 2., h_pcb_slot_height / 2., thickness / 2.);
      cutout = UnionSolid(cutout, h_pcb_slot, Position(rect_x, h_pcb_slot_y, 0.));
      cutout = UnionSolid(cutout, h_pcb_slot, Position(rect_x, -h_pcb_slot_y, 0.));

      Box v_pcb_slot(v_pcb_slot_width / 2., right_v_pcb_slot_height / 2., thickness / 2.);
      cutout = UnionSolid(cutout, v_pcb_slot, Position(right_v_pcb_slot_x, right_v_pcb_slot_y, 0.));
      cutout =
          UnionSolid(cutout, v_pcb_slot, Position(right_v_pcb_slot_x, -right_v_pcb_slot_y, 0.));
    }

    // Add the left-side PCB slots to its circular cutout.
    else {
      Box v_pcb_slot(v_pcb_slot_width / 2., left_v_pcb_slot_height / 2., thickness / 2.);
      cutout = UnionSolid(cutout, v_pcb_slot, Position(left_v_pcb_slot_x, left_v_pcb_slot_y, 0.));
      cutout = UnionSolid(cutout, v_pcb_slot, Position(left_v_pcb_slot_x, -left_v_pcb_slot_y, 0.));
    }

    return cutout;
  };

  // Insert assembly.
  Assembly assembly(detName);
  // FIXME Workaround for https://github.com/eic/epic/issues/411
  assembly.setVisAttributes(desc.visAttributes("InvisibleWithDaughters"));
  PlacedVolume pv;

  // Loop through the front-face right (-x) and left (+x) sides.
  for (int side = 0; side < 2; side++) { // 0 = right (-x), 1 = left (+x)
    std::string side_name = side == 1 ? "L" : "R";

    double z_distance_traversed = 0.;
    int layer_num               = 1;

    // Loop through layer sections.
    for (xml_coll_t c(detElem, _U(layer)); c; c++) {
      xml_comp_t x_layer     = c;
      int repeat             = x_layer.repeat();
      double layer_thickness = x_layer.thickness();

      // Loop through each layer in this section.
      for (int i = 0; i < repeat; i++) {
        std::string layer_name = detName + _toString(layer_num, "_layer%d") + "_" + side_name;
        Box layer(width / 2., height / 2., layer_thickness / 2.);

        // Build the layer container with the side-specific cutout.
        Solid layer_cutout = make_cutout(layer_thickness, side);
        Solid layer_with_cutout =
            SubtractionSolid(layer, layer_cutout, Position(center_x, center_y, 0.));

        // Extend the left-side Air container around two tiles that cross the circular cutout.
        if (side == 1) {
          Box tile_pocket(tile_size / 2., tile_size / 2., layer_thickness / 2.);
          layer_with_cutout = UnionSolid(layer_with_cutout, tile_pocket, left_tile_position(4, 2));
          layer_with_cutout = UnionSolid(layer_with_cutout, tile_pocket, left_tile_position(7, 2));
        }

        // Side boundary.
        Box side_cut(width / 2., height, layer_thickness);
        Position side_cut_position(
            side == 0 ? right_split_x + width / 2. : left_split_x - width / 2., 0., 0.);
        SubtractionSolid layer_side_with_cutout(layer_with_cutout, side_cut, side_cut_position);
        Volume layer_vol(layer_name, layer_side_with_cutout, air);

        int slice_num  = 1;
        double slice_z = -layer_thickness / 2.;

        // Loop through layer slices.
        for (xml_coll_t l(x_layer, _U(slice)); l; l++) {
          xml_comp_t x_slice     = l;
          double slice_thickness = x_slice.thickness();
          std::string slice_name = layer_name + _toString(slice_num, "slice%d");
          Material slice_mat     = desc.material(x_slice.materialStr());
          slice_z += slice_thickness / 2.;

          // Build a sensitive slice from physical tiles.
          if (x_slice.isSensitive()) {
            Assembly tile_assembly(slice_name + "_tiles");
            Box tile(tile_size / 2., tile_size / 2., slice_thickness / 2.);
            Volume tile_vol(slice_name + "_tile", tile, slice_mat);

            // Tile attributes.
            tile_vol.setSensitiveDetector(sens);
            tile_vol.setAttributes(desc, x_slice.regionStr(), x_slice.limitsStr(),
                                   x_slice.visStr());

            // Loop through all rows of tiles.
            for (std::size_t row = 0; row < left_tiles_per_row.size(); row++) {

              int tiles_in_row = side == 0 ? right_tiles_per_row[row] : left_tiles_per_row[row];

              // Per row, loop through each column.
              for (int column = 0; column < tiles_in_row; column++) {

                // Place tile volume depending on side.
                Position tile_position =
                    side == 0 ? right_tile_position(row, column) : left_tile_position(row, column);
                PlacedVolume tile_pv = tile_assembly.placeVolume(tile_vol, tile_position);

                // Tile IDs.
                tile_pv.addPhysVolID("x", column).addPhysVolID("y", row);
              }
            }

            // Tile assembly.
            pv = layer_vol.placeVolume(
                tile_assembly, Transform3D(RotationZYX(0, 0, 0), Position(0., 0., slice_z)));
            pv.addPhysVolID("slice", slice_num);
            pv.addPhysVolID("side", side);
            slice_z += slice_thickness / 2.;
            z_distance_traversed += slice_thickness;
            slice_num++;
            continue;
          }

          // Passive slice.
          Box slice(width / 2., height / 2., slice_thickness / 2.);
          Solid slice_cutout = make_cutout(slice_thickness, side);
          SubtractionSolid slice_with_cutout(slice, slice_cutout, Position(center_x, center_y, 0.));
          Box side_cut_slice(width / 2., height, layer_thickness);
          Position side_cut_position_slice(
              side == 0 ? right_split_x + width / 2. : left_split_x - width / 2., 0, 0.);
          SubtractionSolid slice_side_with_cutout(slice_with_cutout, side_cut_slice,
                                                  side_cut_position_slice);
          Volume slice_vol(slice_name, slice_side_with_cutout, slice_mat);

          // Slice attributes.
          slice_vol.setAttributes(desc, x_slice.regionStr(), x_slice.limitsStr(), x_slice.visStr());

          // Slice placement.
          pv = layer_vol.placeVolume(slice_vol,
                                     Transform3D(RotationZYX(0, 0, 0), Position(0., 0., slice_z)));
          pv.addPhysVolID("slice", slice_num);
          pv.addPhysVolID("side", side);
          slice_z += slice_thickness / 2.;
          z_distance_traversed += slice_thickness;
          slice_num++;
        }

        // Layer attributes.
        layer_vol.setAttributes(desc, x_layer.regionStr(), x_layer.limitsStr(), x_layer.visStr());

        // Layer placement.
        pv = assembly.placeVolume(
            layer_vol, Transform3D(RotationZYX(0, 0, 0),
                                   Position(0., 0.,
                                            -length / 2. + casing_thickness +
                                                (z_distance_traversed - layer_thickness) +
                                                layer_thickness / 2.)));

        pv.addPhysVolID("layer", layer_num);
        pv.addPhysVolID("side", side);
        layer_num++;
      }
    }
  }

  // Build one casing box and divide it into right and left halves.
  Material casing_material   = desc.material(casing_xml.materialStr());
  const double casing_width  = width + 2. * casing_thickness;
  const double casing_height = height + 2. * casing_thickness;
  const double casing_z      = -(length - casing_length) / 2.;

  Box casing_outer(casing_width / 2., casing_height / 2., casing_length / 2.);
  Box casing_inner(width / 2., height / 2., casing_length / 2. - casing_thickness);
  SubtractionSolid casing_box(casing_outer, casing_inner);
  const double back_cutout_y = (height / 2. + center_y + rect_height / 2.) / 2.;
  const double back_cutout_z = casing_length / 2. - casing_thickness / 2.;

  // Loop through both sides.
  for (int side = 0; side < 2; side++) {
    Box side_cut(casing_width / 2., casing_height, casing_length);
    Position side_cut_position(
        side == 0 ? right_split_x + casing_width / 2. : left_split_x - casing_width / 2., 0., 0.);
    Solid casing_side = SubtractionSolid(casing_box, side_cut, side_cut_position);

    // Cut the side's beampipe opening through the casing.
    Solid casing_cutout = Tube(0., casing_left_radius, casing_length / 2. + casing_thickness);
    if (side == 0) {
      const double right_rect_width = right_split_x - center_x;
      Tube right_cap(0., rect_height / 2., casing_length / 2. + casing_thickness, M_PI / 2.,
                     3. * M_PI / 2.);
      Box right_rect(right_rect_width / 2., rect_height / 2.,
                     casing_length / 2. + casing_thickness);
      casing_cutout = UnionSolid(right_cap, right_rect, Position(right_rect_width / 2., 0., 0.));
    }

    Solid casing_with_cutout =
        SubtractionSolid(casing_side, casing_cutout, Position(center_x, center_y, 0.));

    // Cut mirrored openings through each rear plate.
    const double cutout_width = side == 0 ? right_back_cutout_width : left_back_cutout_width;
    const double cutout_x     = side == 0 ? (-casing_width / 2. + right_split_x) / 2.
                                          : (left_split_x + casing_width / 2.) / 2.;
    Box back_cutout(cutout_width / 2., back_cutout_height / 2., casing_thickness);
    casing_with_cutout      = SubtractionSolid(casing_with_cutout, back_cutout,
                                               Position(cutout_x, back_cutout_y, back_cutout_z));
    casing_with_cutout      = SubtractionSolid(casing_with_cutout, back_cutout,
                                               Position(cutout_x, -back_cutout_y, back_cutout_z));
    std::string casing_name = detName + (side == 0 ? "_RightCasing" : "_LeftCasing");
    Volume casing_vol(casing_name, casing_with_cutout, casing_material);
    casing_vol.setVisAttributes(desc.visAttributes(casing_xml.visStr()));
    assembly.placeVolume(casing_vol, Position(0., 0., casing_z));
  }

  // Keep shared cover geometry inputs available for the remaining PCB covers.
  const double left_split_from_center = left_split_x - center_x;
  Material cover_material             = desc.material(cover_xml.materialStr());

  // (maybe) FIXME: Restore the left-side conical cover after figuring out what to do with the tiles that stick out.
  /*
  const double cover_phi = std::acos(left_split_from_center / casing_left_radius);
  ConeSegment conical_cover(pcb_length / 2., casing_left_radius - cover_thickness,
                            casing_left_radius,
                            casing_left_radius - cover_thickness, casing_left_radius,
                            -cover_phi, cover_phi);
  Volume cover_vol(detName + "_LeftPCBConicalCover", conical_cover, cover_material);
  cover_vol.setVisAttributes(desc.visAttributes(cover_xml.visStr()));
  assembly.placeVolume(cover_vol, Position(center_x, center_y, pcb_z));
  */

  // Place Horizontal PCBs.
  Material pcb_material = desc.material("Fr4");
  Box h_pcb(h_pcb_width / 2., pcb_thickness / 2., pcb_length / 2.);
  Volume h_pcb_vol(detName + "_HorizontalPCB", h_pcb, pcb_material);
  h_pcb_vol.setVisAttributes(desc.visAttributes("AnlDarkGreen"));
  assembly.placeVolume(h_pcb_vol, Position(center_x + rect_x, center_y + h_pcb_slot_y, pcb_z));
  assembly.placeVolume(h_pcb_vol, Position(center_x + rect_x, center_y - h_pcb_slot_y, pcb_z));

  // Line the top and bottom edges of the rectangular beampipe opening.
  Box h_pcb_cover(rect_width / 2., cover_thickness / 2., pcb_length / 2.);
  Volume h_pcb_cover_vol(detName + "_HorizontalPCBCover", h_pcb_cover, cover_material);
  h_pcb_cover_vol.setVisAttributes(desc.visAttributes(cover_xml.visStr()));
  assembly.placeVolume(
      h_pcb_cover_vol,
      Position(center_x + rect_x, center_y + rect_height / 2. - cover_thickness / 2., pcb_z));
  assembly.placeVolume(
      h_pcb_cover_vol,
      Position(center_x + rect_x, center_y - rect_height / 2. + cover_thickness / 2., pcb_z));

  // Place vertical PCBs for the right side.
  Box right_v_pcb(pcb_thickness / 2., right_v_pcb_height / 2., pcb_length / 2.);
  Volume right_v_pcb_vol(detName + "_RightVerticalPCB", right_v_pcb, pcb_material);
  right_v_pcb_vol.setVisAttributes(desc.visAttributes("AnlDarkGreen"));
  assembly.placeVolume(right_v_pcb_vol, Position(center_x + right_v_pcb_slot_x,
                                                 center_y + right_v_pcb_slot_y, pcb_z));
  assembly.placeVolume(right_v_pcb_vol, Position(center_x + right_v_pcb_slot_x,
                                                 center_y - right_v_pcb_slot_y, pcb_z));

  // Then place the covers for the right vertical PCBs.
  const double right_v_pcb_cover_inner_y = rect_height / 2. - cover_thickness;
  const double right_v_pcb_cover_outer_y = height / 2. + casing_thickness;
  const double right_v_pcb_cover_height  = right_v_pcb_cover_outer_y - right_v_pcb_cover_inner_y;
  const double right_v_pcb_cover_y = (right_v_pcb_cover_inner_y + right_v_pcb_cover_outer_y) / 2.;
  Box right_v_pcb_cover(cover_thickness / 2., right_v_pcb_cover_height / 2., pcb_length / 2.);
  Volume right_v_pcb_cover_vol(detName + "_RightVerticalPCBCover", right_v_pcb_cover,
                               cover_material);
  right_v_pcb_cover_vol.setVisAttributes(desc.visAttributes(cover_xml.visStr()));
  assembly.placeVolume(
      right_v_pcb_cover_vol,
      Position(center_x + right_v_pcb_slot_x + v_pcb_slot_width / 2. + cover_thickness / 2.,
               center_y + right_v_pcb_cover_y, pcb_z));
  assembly.placeVolume(
      right_v_pcb_cover_vol,
      Position(center_x + right_v_pcb_slot_x + v_pcb_slot_width / 2. + cover_thickness / 2.,
               center_y - right_v_pcb_cover_y, pcb_z));

  // Place vertical PCBs for the left side.
  Box left_v_pcb(pcb_thickness / 2., left_v_pcb_height / 2., pcb_length / 2.);
  Volume left_v_pcb_vol(detName + "_LeftVerticalPCB", left_v_pcb, pcb_material);
  left_v_pcb_vol.setVisAttributes(desc.visAttributes("AnlDarkGreen"));
  assembly.placeVolume(left_v_pcb_vol,
                       Position(center_x + left_v_pcb_slot_x, center_y + left_v_pcb_slot_y, pcb_z));
  assembly.placeVolume(left_v_pcb_vol,
                       Position(center_x + left_v_pcb_slot_x, center_y - left_v_pcb_slot_y, pcb_z));

  // Then place the covers for the left vertical PCBs.
  const double left_v_pcb_cover_inner_y = std::sqrt(
      casing_left_radius * casing_left_radius -
      left_split_from_center * left_split_from_center);
  const double left_v_pcb_cover_outer_y = height / 2.;
  const double left_v_pcb_cover_height  = left_v_pcb_cover_outer_y - left_v_pcb_cover_inner_y;
  const double left_v_pcb_cover_y = (left_v_pcb_cover_inner_y + left_v_pcb_cover_outer_y) / 2.;
  Box left_v_pcb_cover(cover_thickness / 2., left_v_pcb_cover_height / 2., pcb_length / 2.);
  Volume left_v_pcb_cover_vol(detName + "_LeftVerticalPCBCover", left_v_pcb_cover, cover_material);
  left_v_pcb_cover_vol.setVisAttributes(desc.visAttributes(cover_xml.visStr()));
  assembly.placeVolume(
      left_v_pcb_cover_vol,
      Position(center_x + left_v_pcb_slot_x - v_pcb_slot_width / 2. - cover_thickness / 2.,
               center_y + left_v_pcb_cover_y, pcb_z));
  assembly.placeVolume(
      left_v_pcb_cover_vol,
      Position(center_x + left_v_pcb_slot_x - v_pcb_slot_width / 2. - cover_thickness / 2.,
               center_y - left_v_pcb_cover_y, pcb_z));

  DetElement det(detName, detID);
  Volume motherVol = desc.pickMotherVolume(det);

  // Detector flags.
  dd4hep::xml::setDetectorTypeFlag(detElem, det);

  // World placement.
  auto tr          = Transform3D(Position(pos.x(), pos.y(), pos.z() + length / 2.));
  PlacedVolume phv = motherVol.placeVolume(assembly, tr);
  phv.addPhysVolID("system", detID);
  det.setPlacement(phv);

  return det;
}
DECLARE_DETELEMENT(epic_InsertCalorimeter_new, createDetector)
