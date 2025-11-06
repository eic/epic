// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2023-2025 Justin Chan, Simon Gardner

#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/Printout.h"
#include "TMath.h"
#include <XML/Helper.h>

using namespace std;
using namespace dd4hep;

static Ref_t create_detector(Detector& det, xml_h e, SensitiveDetector /* sens */) {

  using namespace ROOT::Math;
  xml_det_t x_det = e;
  string det_name = x_det.nameStr();
  DetElement sdet(det_name, x_det.id());
  Assembly assembly(det_name + "_assembly");
  Material m_Iron   = det.material("Iron");
  // Material m_Copper = det.material("Copper");
  const string vis1 = getAttrOrDefault<string>(x_det, _Unicode(vis), "AnlGreen");
  // auto vis_air  = det.visAttributes("BackwardsVac");
  // auto vis_iron = det.visAttributes(vis1);

  // Creates the outer box for the main body
  xml::Component outer_box_dim = x_det.child(_Unicode(dimensions_mainbody_outer));
  double outer_height          = outer_box_dim.attr<double>(_Unicode(y));
  double outer_width           = outer_box_dim.attr<double>(_Unicode(x));
  double outer_depth           = outer_box_dim.attr<double>(_Unicode(z));

  Box box_outer(outer_width / 2., outer_height / 2., outer_depth / 2.);

  // Creates the innner box for the main body
  xml::Component inner_box_dim = x_det.child(_Unicode(dimensions_mainbody_inner));
  double inner_height          = inner_box_dim.attr<double>(_Unicode(y));
  double inner_width           = inner_box_dim.attr<double>(_Unicode(x));

  // Yoke/coil parameters
  double pole_gap = getAttrOrDefault<double>(x_det, _Unicode(pole_gap), inner_height/2.0);
  double yoke_height = inner_height / 2.0 - pole_gap/2;
  double yoke_outer_dim = getAttrOrDefault<double>(x_det, _Unicode(yoke_outer_width), inner_width);
  double yoke_inner_dim = getAttrOrDefault<double>(x_det, _Unicode(yoke_inner_width), inner_width);


  // // Creates the outer box for the shape of the coils
  // xml::Component box_dim_3 = x_det.child(_Unicode(dimensions_coils_outer));
  // double height_3          = box_dim_3.attr<double>(_Unicode(y));
  // double width_3           = box_dim_3.attr<double>(_Unicode(x));
  // double depth_3           = box_dim_3.attr<double>(_Unicode(z));

  // Box coils_outer(width_3 / 2., height_3 / 2., depth_3 / 2.);

  // // Creates the first inner box for the shape of the coils
  // xml::Component box_dim_4 = x_det.child(_Unicode(dimensions_coils_inner_1));
  // double height_4          = box_dim_4.attr<double>(_Unicode(y));
  // double width_4           = box_dim_4.attr<double>(_Unicode(x));
  // double depth_4           = box_dim_4.attr<double>(_Unicode(z));

  // Box coils_inner_1(width_4 / 2., height_4 / 2., depth_4 / 2.);

  // // Creates the second inner box for the shape of the coils
  // xml::Component box_dim_5 = x_det.child(_Unicode(dimensions_coils_inner_2));
  // double height_5          = box_dim_5.attr<double>(_Unicode(y));
  // double width_5           = box_dim_5.attr<double>(_Unicode(x));
  // double depth_5           = box_dim_5.attr<double>(_Unicode(z));

  // Box coils_inner_2(width_5 / 2., height_5 / 2., depth_5 / 2.);

  // // Creates the outer box for the shape of the yoke
  // xml::Component box_dim_6 = x_det.child(_Unicode(dimensions_yoke_outer));
  // double height_6          = box_dim_6.attr<double>(_Unicode(y));
  // double width_6           = box_dim_6.attr<double>(_Unicode(x));
  // double depth_6           = box_dim_6.attr<double>(_Unicode(z));

  // Box yoke_outer(width_6 / 2., height_6 / 2., depth_6 / 2.);

  // // Creates the first inner box for the shape of the coils
  // xml::Component box_dim_7 = x_det.child(_Unicode(dimensions_yoke_inner));
  // double height_7          = box_dim_7.attr<double>(_Unicode(y));
  // double width_7           = box_dim_7.attr<double>(_Unicode(x));
  // double depth_7           = box_dim_7.attr<double>(_Unicode(z));

  // Box yoke_inner(width_7 / 2., height_7 / 2., depth_7 / 2.);

  // // Creates the outer box for the shape of the legs
  // xml::Component box_dim_8 = x_det.child(_Unicode(dimensions_leg_outer));
  // double height_8          = box_dim_8.attr<double>(_Unicode(y));
  // double width_8           = box_dim_8.attr<double>(_Unicode(x));
  // double depth_8           = box_dim_8.attr<double>(_Unicode(z));

  // Box leg_outer(width_8 / 2., height_8 / 2., depth_8 / 2.);

  // // Creates the inner box for the shape of the legs
  // xml::Component box_dim_9 = x_det.child(_Unicode(dimensions_leg_inner));
  // double height_9          = box_dim_9.attr<double>(_Unicode(y));
  // double width_9           = box_dim_9.attr<double>(_Unicode(x));
  // double depth_9           = box_dim_9.attr<double>(_Unicode(z));

  // Box leg_inner(width_9 / 2., height_9 / 2., depth_9 / 2.);

  // Sets box position
  xml::Component box_pos = x_det.child(_Unicode(position));
  double x               = box_pos.attr<double>(_Unicode(x));
  double y               = box_pos.attr<double>(_Unicode(y));
  double z               = box_pos.attr<double>(_Unicode(z));

  // Calculate wall thickness
  double thickness_x = (outer_width - inner_width) / 2.0;
  double thickness_y = (outer_height - inner_height) / 2.0;

  // Create 4 separate box volumes for the sides of the rectangular tube
  // Top wall (full width, thickness in y, full depth)
  Box box_top(outer_width / 2.0, thickness_y / 2.0, outer_depth / 2.0);
  Volume vol_top(det_name + "_vol_top", box_top, m_Iron);
  vol_top.setAttributes(det, x_det.regionStr(), x_det.limitsStr(), vis1);

  // Bottom wall (full width, thickness in y, full depth)
  Volume vol_bottom(det_name + "_vol_bottom", box_top, m_Iron);
  vol_bottom.setAttributes(det, x_det.regionStr(), x_det.limitsStr(), vis1);

  // Left wall (thickness in x, inner height, full depth)
  Box box_side(thickness_x / 2.0, inner_height / 2.0, outer_depth / 2.0);
  Volume vol_left(det_name + "_vol_left", box_side, m_Iron);
  vol_left.setAttributes(det, x_det.regionStr(), x_det.limitsStr(), vis1);

  // Right wall (thickness in x, inner height, full depth)
  Volume vol_right(det_name + "_vol_right", box_side, m_Iron);
  vol_right.setAttributes(det, x_det.regionStr(), x_det.limitsStr(), vis1);

  // Upper coils/yoke - This is not entirely accurate as CAD has curved edges
  Trd1 yoke_trap(yoke_outer_dim / 2.0, yoke_inner_dim / 2.0,
                  outer_depth / 2.0, yoke_height / 2.0);
  Volume vol_yoke(det_name + "_vol_yoke", yoke_trap, m_Iron);
  vol_yoke.setAttributes(det, x_det.regionStr(), x_det.limitsStr(), vis1);

  // Create assembly for the magnet
  Assembly magnet_assembly(det_name + "_magnet_assembly");

  // Place the 4 walls
  // Top wall at +y
  magnet_assembly.placeVolume(vol_top, Position(0, outer_height / 2.0 - thickness_y / 2.0, 0));
  
  // Bottom wall at -y
  magnet_assembly.placeVolume(vol_bottom, Position(0, -outer_height / 2.0 + thickness_y / 2.0, 0));
  
  // Left wall at -x
  magnet_assembly.placeVolume(vol_left, Position(-outer_width / 2.0 + thickness_x / 2.0, 0, 0));
  
  // Right wall at +x
  magnet_assembly.placeVolume(vol_right, Position(outer_width / 2.0 - thickness_x / 2.0, 0, 0));

  //Translation for top and bottom yoke
  

  // Place the top yoke with position and rotation around X axis
  magnet_assembly.placeVolume(
      vol_yoke, Transform3D(RotationX(TMath::Pi()/2), Position(0, inner_height / 2.0 - yoke_height / 2.0, 0)));

  // Place the bottom yoke
  magnet_assembly.placeVolume(
      vol_yoke, Transform3D(RotationX(-TMath::Pi()/2), Position(0, -inner_height / 2.0 + yoke_height / 2.0, 0)));

  // Bar parameters
 for (xml_coll_t bar_coll(x_det, _Unicode(bar)); bar_coll; ++bar_coll) {
    double bar_width  = bar_coll.attr<double>(_Unicode(width));
    double bar_height = bar_coll.attr<double>(_Unicode(height));
    double bar_length = bar_coll.attr<double>(_Unicode(length));
    // Create box for bar
    Box bar_box(bar_width / 2., bar_height / 2., bar_length / 2.);
    Volume vol_bar(det_name + "_vol_bar", bar_box, m_Iron);
    vol_bar.setAttributes(det, x_det.regionStr(), x_det.limitsStr(), vis1);
    for(xml_coll_t pos(bar_coll, _Unicode(position)); pos; ++pos) {
      double posX = pos.attr<double>(_Unicode(x));
      double posY = pos.attr<double>(_Unicode(y));
      double posZ = pos.attr<double>(_Unicode(z));
      // Place bar at specified position
      magnet_assembly.placeVolume(vol_bar, Position(posX, posY, posZ));
    }
  }


  // Subtractes the volume of the inner box from the outer box for the main body
  // BooleanSolid main_body = SubtractionSolid(box_outer, box_inner);
  // Volume v_main_body(det_name + "_vol_main_body", main_body, m_Iron);
  // sdet.setAttributes(det, v_main_body, x_det.regionStr(), x_det.limitsStr(), vis1);
  // assembly.placeVolume(v_main_body, Position(x, y, z));

  // // Creates panels by subtracting the inner boxes of the coils from the outer box of the coilss
  // BooleanSolid coils_1 = SubtractionSolid(coils_outer, coils_inner_1);
  // BooleanSolid coils   = SubtractionSolid(coils_1, coils_inner_2);
  // Volume v_coils(det_name + "_vol_coils", coils, m_Copper);
  // sdet.setAttributes(det, v_coils, x_det.regionStr(), x_det.limitsStr(), vis2);
  // assembly.placeVolume(v_coils, Position(x, y, z));

  // // Creates coils by subtracting the inner box of the yoke from the outer box of the yoke
  // BooleanSolid yoke = SubtractionSolid(yoke_outer, yoke_inner);
  // Volume v_yoke(det_name + "_vol_yoke", yoke, m_Iron);
  // sdet.setAttributes(det, v_yoke, x_det.regionStr(), x_det.limitsStr(), vis1);
  // assembly.placeVolume(v_yoke, Position(x, y, z));

  // Creates the legs by subtracting the inner box of the legs from the outer box of the legs
  // BooleanSolid legs = SubtractionSolid(leg_outer, leg_inner);
  // Volume v_legs(det_name + "_vol_legs", legs, m_Iron);
  // sdet.setAttributes(det, v_legs, x_det.regionStr(), x_det.limitsStr(), vis3);
  // assembly.placeVolume(v_legs, Position(x, y - height_9 / 2. - height / 2., z));

  // Final placement
  auto pv_assembly = det.pickMotherVolume(sdet).placeVolume(
      magnet_assembly, Transform3D(RotationZYX(0.0, 0.0, 0.0), Position(x, y, z)));

  sdet.setPlacement(pv_assembly);

  // assembly->GetShape()->ComputeBBox();

  return sdet;
}

DECLARE_DETELEMENT(BeamlineDipoleMagnet, create_detector)
