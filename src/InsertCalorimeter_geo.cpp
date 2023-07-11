// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Ryan Milton

//==========================================================================
//  Implementation of forward insert calorimeter
//--------------------------------------------------------------------------
//  Author: Ryan Milton (UCR)
//==========================================================================

#include "DD4hep/DetFactoryHelper.h"
#include <XML/Helper.h>
#include <XML/Layering.h>
#include <tuple>
#include <vector>

using namespace dd4hep;

static Ref_t createDetector(Detector& desc, xml_h handle, SensitiveDetector sens)
{
  xml_det_t   detElem = handle;
  std::string detName = detElem.nameStr();
  int         detID   = detElem.id();

  xml_dim_t dim    = detElem.dimensions();
  double    width  = dim.x();         // Size along x-axis
  double    height = dim.y();         // Size along y-axis
  double    length = dim.z();         // Size along z-axis

  xml_dim_t pos = detElem.position(); // Position in global coordinates

  Material air = desc.material("Air");

  // Getting beampipe hole dimensions
  const xml::Component& beampipe_hole_xml = detElem.child(_Unicode(beampipe_hole));
  const double          hole_radius_initial =
      dd4hep::getAttrOrDefault<double>(beampipe_hole_xml, _Unicode(initial_hole_radius), 14.61 * cm);
  const double hole_radius_final =
      dd4hep::getAttrOrDefault<double>(beampipe_hole_xml, _Unicode(final_hole_radius), 17.04 * cm);
  const std::pair<double, double> hole_radii_parameters(hole_radius_initial, hole_radius_final);

  // Subtract by pos.x() and pos.y() to convert from global to local coordinates
  const double hole_x_initial =
      dd4hep::getAttrOrDefault<double>(beampipe_hole_xml, _Unicode(initial_hole_x), -7.20 * cm) - pos.x();
  const double hole_x_final =
      dd4hep::getAttrOrDefault<double>(beampipe_hole_xml, _Unicode(final_hole_x), -10.27 * cm) - pos.x();
  const std::pair<double, double> hole_x_parameters(hole_x_initial, hole_x_final);

  const double hole_y_initial =
      dd4hep::getAttrOrDefault<double>(beampipe_hole_xml, _Unicode(initial_hole_y), 0. * cm) - pos.y();
  const double hole_y_final =
      dd4hep::getAttrOrDefault<double>(beampipe_hole_xml, _Unicode(final_hole_y), 0. * cm) - pos.y();
  const std::pair<double, double> hole_y_parameters(hole_y_initial, hole_y_final);

  // Getting thickness of backplate
  /*
    The hole radius & position is determined by liner interpolation
    The HCal insert has multiple layers; interpolation goes from front of insert to front of final layer
    So need final layer thickness (backplate thickness) for interpolation

    For the ECal insert, the hole radius & position is constant
    Also has only one layer so don't have a backplate_thickness there (so set to 0)
  */
  auto backplate_thickness =
      detElem.hasChild(_Unicode(backplate)) ? detElem.child(_Unicode(backplate)).attr<double>(_Unicode(thickness)) : 0.;

  // Function that returns a linearly interpolated hole radius, x-position, and y-position at a given z
  auto get_hole_rxy = [hole_radii_parameters, hole_x_parameters, hole_y_parameters, length,
                       backplate_thickness](double z_pos) {
    /*
      radius = (hole_radius_final - hole_radius_initial)/(length) * z_pos +  hole_radius_initial
      Treats the beginning of the insert as z = 0
      At z = 0 (beginning of first layer), hole radius is hole_radius_initial
      The radius is hole_radius_final at the beginning of the backplate,
        i.e. z = length - backplate_thickness
    */
    double hole_radius_slope =
        (hole_radii_parameters.second - hole_radii_parameters.first) / (length - backplate_thickness);
    double hole_radius_at_z = hole_radius_slope * z_pos + hole_radii_parameters.first;

    double hole_xpos_slope = (hole_x_parameters.second - hole_x_parameters.first) / (length - backplate_thickness);
    double hole_xpos_at_z  = hole_xpos_slope * z_pos + hole_x_parameters.first;

    double hole_ypos_slope = (hole_y_parameters.second - hole_y_parameters.first) / (length - backplate_thickness);
    double hole_ypos_at_z  = hole_ypos_slope * z_pos + hole_y_parameters.first;

    std::tuple<double, double, double> hole_rxy = std::make_tuple(hole_radius_at_z, hole_xpos_at_z, hole_ypos_at_z);
    return hole_rxy;
  };

  // Assembly that will contain all the layers
  Assembly assembly(detName);
  // FIXME Workaround for https://github.com/eic/epic/issues/411
  assembly.setVisAttributes(desc.visAttributes("InvisibleWithDaughters"));
  PlacedVolume pv;

  // Keeps track of the z location as we move longiduinally through the insert
  // Will use this tracking variable as input to get_hole_rxy
  double z_distance_traversed = 0.;

  int layer_num = 1;

  // Looping through all the different layer sections (W/Sc, Steel/Sc, backplate)
  for (xml_coll_t c(detElem, _U(layer)); c; c++) {
    xml_comp_t x_layer         = c;
    int        repeat          = x_layer.repeat();
    double     layer_thickness = x_layer.thickness();

    // Looping through the number of repeated layers in each section
    for (int i = 0; i < repeat; i++) {
      std::string layer_name = detName + _toString(layer_num, "_layer%d");
      Box         layer(width / 2., height / 2., layer_thickness / 2.);

      // Hole radius and position for each layer is determined from z position at the front of the layer
      const auto hole_rxy = get_hole_rxy(z_distance_traversed);
      double     hole_r   = std::get<0>(hole_rxy);
      double     hole_x   = std::get<1>(hole_rxy);
      double     hole_y   = std::get<2>(hole_rxy);

      // Removing beampipe shape from each layer
      Tube             layer_hole(0., hole_r, layer_thickness / 2.);
      SubtractionSolid layer_with_hole(layer, layer_hole, Position(hole_x, hole_y, 0.));
      Volume           layer_vol(layer_name, layer_with_hole, air);

      int    slice_num = 1;
      double slice_z   = -layer_thickness / 2.; // Keeps track of slices' z locations in each layer

      // Looping over each layer's slices
      for (xml_coll_t l(x_layer, _U(slice)); l; l++) {
        xml_comp_t  x_slice         = l;
        double      slice_thickness = x_slice.thickness();
        std::string slice_name      = layer_name + _toString(slice_num, "slice%d");
        Material    slice_mat       = desc.material(x_slice.materialStr());
        slice_z += slice_thickness / 2.; // Going to slice halfway point

        // Each slice within a layer has the same hole radius and x-y position
        Box              slice(width / 2., height / 2., slice_thickness / 2.);
        Tube             slice_hole(0., hole_r, slice_thickness / 2.);
        SubtractionSolid slice_with_hole(slice, slice_hole, Position(hole_x, hole_y, 0.));
        Volume           slice_vol(slice_name, slice_with_hole, slice_mat);

        // Setting appropriate slices as sensitive
        if (x_slice.isSensitive()) {
          sens.setType("calorimeter");
          slice_vol.setSensitiveDetector(sens);
        }

        // Setting slice attributes
        slice_vol.setAttributes(desc, x_slice.regionStr(), x_slice.limitsStr(), x_slice.visStr());

        // Placing slice within layer
        pv = layer_vol.placeVolume(slice_vol, Transform3D(RotationZYX(0, 0, 0), Position(0., 0., slice_z)));
        pv.addPhysVolID("slice", slice_num);
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
          Transform3D(
              RotationZYX(0, 0, 0),
              Position(0., 0., -length / 2. + (z_distance_traversed - layer_thickness) + layer_thickness / 2.)));

      pv.addPhysVolID("layer", layer_num);
      layer_num++;
    }
  }

  DetElement det(detName, detID);
  Volume     motherVol = desc.pickMotherVolume(det);

  // Placing insert in world volume
  auto         tr  = Transform3D(Position(pos.x(), pos.y(), pos.z() + length / 2.));
  PlacedVolume phv = motherVol.placeVolume(assembly, tr);
  phv.addPhysVolID("system", detID);
  det.setPlacement(phv);

  return det;
}
DECLARE_DETELEMENT(epic_InsertCalorimeter, createDetector)
