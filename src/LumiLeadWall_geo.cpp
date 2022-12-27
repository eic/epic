// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2023 Dhevan Gangadharan

//==========================================================================
//
// Lead wall to block secondaries produced upstream of lumi detectors.
// Cutout in center for photon beam clearance.
//
//==========================================================================

#include "DD4hep/DetFactoryHelper.h"
#include <XML/Helper.h>

using namespace std;
using namespace dd4hep;

static Ref_t create_detector(Detector& description, xml_h e, SensitiveDetector /* sens */)
{

  using namespace ROOT::Math;
  xml_det_t  x_det    = e;
  string     det_name = x_det.nameStr();
  DetElement sdet(det_name, x_det.id());
  Assembly   assembly(det_name + "_assembly");
  Material   m_Steel  = description.material("Lead");
  string     vis_name = x_det.visStr();

  // Create outer box
  xml::Component box_dim = x_det.child(_Unicode(dimensions_outer));
  double         height  = box_dim.attr<double>(_Unicode(y));;
  double         width   = box_dim.attr<double>(_Unicode(x));;
  double         depth   = box_dim.attr<double>(_Unicode(z));;

  Box box_outer(width, height, depth);

  // Create inner box
  xml::Component box_dim_2 = x_det.child(_Unicode(dimensions_inner));
  double         height_2  = box_dim_2.attr<double>(_Unicode(y));;
  double         width_2   = box_dim_2.attr<double>(_Unicode(x));;
  double         depth_2   = box_dim_2.attr<double>(_Unicode(z));;

  Box box_inner(width_2, height_2, depth_2);

  // Sets box positions
  xml::Component box_pos = x_det.child(_Unicode(position));
  double         x  = box_pos.attr<double>(_Unicode(x));
  double         y  = box_pos.attr<double>(_Unicode(y));
  double         z  = box_pos.attr<double>(_Unicode(z));

  // Subtractes the volume of the inner box from the outer box
  BooleanSolid collimator =  SubtractionSolid(box_outer, box_inner);

  // Assembles the collimator and sets its material
  Volume v_collimator("v_leadwall", collimator , m_Steel);

  sdet.setAttributes(description, v_collimator, x_det.regionStr(), x_det.limitsStr(), vis_name);

  assembly.placeVolume(v_collimator, Position(x, y, z));

  // Final placement
  auto pv_assembly = description.pickMotherVolume(sdet).placeVolume(
      assembly, Transform3D(RotationZYX(0.0, 0.0, 0.0), Position(0, 0, 0)));
  
  sdet.setPlacement(pv_assembly);
  
  assembly->GetShape()->ComputeBBox();

  return sdet;
}

DECLARE_DETELEMENT(LumiLeadWall, create_detector)
