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
  const string vis1 = getAttrOrDefault<string>(x_det, _Unicode(vis), "AnlGreen");

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
  double pole_gap       = getAttrOrDefault<double>(x_det, _Unicode(pole_gap), inner_height / 2.0);
  double yoke_height    = inner_height / 2.0 - pole_gap / 2;
  double yoke_outer_dim = getAttrOrDefault<double>(x_det, _Unicode(yoke_outer_width), inner_width);
  double yoke_inner_dim = getAttrOrDefault<double>(x_det, _Unicode(yoke_inner_width), inner_width);

  // Sets box position
  xml::Component box_pos = x_det.child(_Unicode(position));
  double x               = box_pos.attr<double>(_Unicode(x));
  double y               = box_pos.attr<double>(_Unicode(y));
  double z               = box_pos.attr<double>(_Unicode(z));

  // Sets box rotation
  xml::Component box_rot = x_det.child(_Unicode(rotation));
  double rot_x           = box_rot.attr<double>(_Unicode(x));
  double rot_y           = box_rot.attr<double>(_Unicode(y));
  double rot_z           = box_rot.attr<double>(_Unicode(z));

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
  Trd1 yoke_trap(yoke_outer_dim / 2.0, yoke_inner_dim / 2.0, outer_depth / 2.0, yoke_height / 2.0);
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
  magnet_assembly.placeVolume(vol_yoke,
                              Transform3D(RotationX(TMath::Pi() / 2),
                                          Position(0, inner_height / 2.0 - yoke_height / 2.0, 0)));

  // Place the bottom yoke
  magnet_assembly.placeVolume(vol_yoke,
                              Transform3D(RotationX(-TMath::Pi() / 2),
                                          Position(0, -inner_height / 2.0 + yoke_height / 2.0, 0)));

  // Coil parameters
  for (xml_coll_t coil_coll(x_det, _Unicode(coil)); coil_coll; ++coil_coll) {
    double coil_width  = coil_coll.attr<double>(_Unicode(width));
    double coil_height = coil_coll.attr<double>(_Unicode(height));
    double coil_length = coil_coll.attr<double>(_Unicode(length));
    // Create box for coil
    Box coil_box(coil_width / 2., coil_height / 2., coil_length / 2.);
    Volume vol_coil(det_name + "_vol_coil", coil_box, m_Iron);
    vol_coil.setAttributes(det, x_det.regionStr(), x_det.limitsStr(), vis1);
    for (xml_coll_t pos(coil_coll, _Unicode(position)); pos; ++pos) {
      double posX = pos.attr<double>(_Unicode(x));
      double posY = pos.attr<double>(_Unicode(y));
      double posZ = pos.attr<double>(_Unicode(z));
      // Place coil at specified position
      magnet_assembly.placeVolume(vol_coil, Position(posX, posY, posZ));
    }
  }

  // Final placement
  auto pv_assembly = det.pickMotherVolume(sdet).placeVolume(
      magnet_assembly, Transform3D(RotationZYX(rot_z, rot_y, rot_x), Position(x, y, z)));

  sdet.setPlacement(pv_assembly);

  // assembly->GetShape()->ComputeBBox();

  return sdet;
}

DECLARE_DETELEMENT(BeamlineDipoleMagnet, create_detector)
