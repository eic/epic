// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Justin Chan

#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/Printout.h"
#include "TMath.h"
#include <XML/Helper.h>

using namespace std;
using namespace dd4hep;

static Ref_t create_detector(Detector& det, xml_h e, SensitiveDetector /* sens */)
{

  using namespace ROOT::Math;
  xml_det_t  x_det    = e;
  string     det_name = x_det.nameStr();
  DetElement sdet(det_name, x_det.id());
  Assembly   assembly(det_name + "_assembly");
  Material   m_Iron  = det.material("Iron");
  Material   m_Copper  = det.material("Copper");
  const string vis1 = getAttrOrDefault<string>(x_det, _Unicode(vis_name1), "AnlGreen");
  const string vis2 = getAttrOrDefault<string>(x_det, _Unicode(vis_name2), "AnlRed");
  const string vis3 = getAttrOrDefault<string>(x_det, _Unicode(vis_name3), "AnlGray");

  // Creates the outer box for the main body
  xml::Component box_dim = x_det.child(_Unicode(dimensions_mainbody_outer));
  double         height  = box_dim.attr<double>(_Unicode(y));
  double         width   = box_dim.attr<double>(_Unicode(x));
  double         depth   = box_dim.attr<double>(_Unicode(z));

  Box box_outer(width/2., height/2., depth/2.);

  // Creates the innner box for the main body
  xml::Component box_dim_2 = x_det.child(_Unicode(dimensions_mainbody_inner));
  double         height_2  = box_dim_2.attr<double>(_Unicode(y));
  double         width_2   = box_dim_2.attr<double>(_Unicode(x));
  double         depth_2   = box_dim_2.attr<double>(_Unicode(z));

  Box box_inner(width_2/2., height_2/2., depth_2/2.);

  // Creates the outer box for the shape of the panels
  xml::Component box_dim_3 = x_det.child(_Unicode(dimensions_panel_outer));
  double         height_3  = box_dim_3.attr<double>(_Unicode(y));
  double         width_3   = box_dim_3.attr<double>(_Unicode(x));
  double         depth_3   = box_dim_3.attr<double>(_Unicode(z));

  Box panel_outer(width_3/2., height_3/2., depth_3/2.);

  // Creates the inner box for the shape of the panels
  xml::Component box_dim_4 = x_det.child(_Unicode(dimensions_panel_inner));
  double         height_4  = box_dim_4.attr<double>(_Unicode(y));
  double         width_4   = box_dim_4.attr<double>(_Unicode(x));
  double         depth_4   = box_dim_4.attr<double>(_Unicode(z));

  Box panel_inner(width_4/2., height_4/2., depth_4/2.);

  // Creates the outer box for the shape of the legs
  xml::Component box_dim_5 = x_det.child(_Unicode(dimensions_leg_outer));
  double         height_5  = box_dim_5.attr<double>(_Unicode(y));
  double         width_5   = box_dim_5.attr<double>(_Unicode(x));
  double         depth_5   = box_dim_5.attr<double>(_Unicode(z));

  Box leg_outer(width_5/2., height_5/2., depth_5/2.);

  // Creates the inner box for the shape of the legs
  xml::Component box_dim_6 = x_det.child(_Unicode(dimensions_leg_inner));
  double         height_6  = box_dim_6.attr<double>(_Unicode(y));
  double         width_6   = box_dim_6.attr<double>(_Unicode(x));
  double         depth_6   = box_dim_6.attr<double>(_Unicode(z));

  Box leg_inner(width_6/2., height_6/2., depth_6/2.);

  // Sets box position
  xml::Component box_pos = x_det.child(_Unicode(position));
  double         x  = box_pos.attr<double>(_Unicode(x));
  double         y  = box_pos.attr<double>(_Unicode(y));
  double         z  = box_pos.attr<double>(_Unicode(z));

  // Subtractes the volume of the inner box from the outer box for the main body
  BooleanSolid main_body =  SubtractionSolid(box_outer, box_inner, Transform3D(Translation3D(0, (height - height_2)/2. , 0)));

  Volume v_main_body( det_name + "_vol_main_body", main_body, m_Iron);

  sdet.setAttributes(det, v_main_body, x_det.regionStr(), x_det.limitsStr(), vis1);

  assembly.placeVolume(v_main_body, Position(x, y, z));

  // Creates panels by subtracting the inner box of the panels from the outer box of the panels
  BooleanSolid panels =  SubtractionSolid(panel_outer, panel_inner);

  Volume v_panels( det_name + "_vol_panels", panels, m_Copper);

  sdet.setAttributes(det, v_panels, x_det.regionStr(), x_det.limitsStr(), vis2);

  assembly.placeVolume(v_panels, Position(x, y + (height - height_2) - (height - height_3)/2., z));

  // Creates the legs by subtracting the inner box of the legs from the outer box of the legs
  BooleanSolid legs =  SubtractionSolid(leg_outer, leg_inner);

  Volume v_legs( det_name + "_vol_legs", legs, m_Iron);

  sdet.setAttributes(det, v_legs, x_det.regionStr(), x_det.limitsStr(), vis3);

  assembly.placeVolume(v_legs, Position(x, y - (height_5 + height)/2., z));

  // Final placement
  auto pv_assembly = det.pickMotherVolume(sdet).placeVolume(
      assembly, Transform3D(RotationZYX(0.0, 0.0, 0.0), Position(0, 0, 0)));

  sdet.setPlacement(pv_assembly);

  assembly->GetShape()->ComputeBBox();

  return sdet;
}

DECLARE_DETELEMENT(LumiMagnets, create_detector)
