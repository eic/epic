// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Alex Jentsch, Whitney Armstrong, Wouter Deconinck

//==========================================================================
//
//      <detector name ="DetName" type="Beampipe" >
//      <layer id="#(int)" inner_r="#(double)" outer_z="#(double)" >
//      <slice material="string" thickness="#(double)" >
//      </layer>
//      </detector>
//==========================================================================
#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/Printout.h"
#include "TMath.h"
#include <XML/Helper.h>

using namespace std;
using namespace dd4hep;

/** \addtogroup beamline Beamline Instrumentation
 */

/** \addtogroup IRChamber Interaction Region Vacuum Chamber.
 * \brief Type: **IRChamber**.
 * \ingroup beamline
 *
 *
 * \code
 *   <detector>
 *   </detector>
 * \endcode
 *
 */
static Ref_t create_detector(Detector& det, xml_h e, SensitiveDetector /* sens */)
{

  using namespace ROOT::Math;
  xml_det_t x_det    = e;
  string    det_name = x_det.nameStr();
  // Material   air       = det.air();
  DetElement sdet(det_name, x_det.id());
  Assembly   assembly(det_name + "_assembly");
  Material   m_Vac    = det.material("Vacuum");
  string     vis_name = x_det.visStr();

  PlacedVolume pv_assembly;


  /// hard-code defintion here, then refine and make more general


  double radius_b0pf = 2.9*dd4hep::cm;
  double length_b0pf = 250.028*dd4hep::cm;   //cm -- shorten by 6mm
  double rotation_b0pf = -0.025;  // radians
  double x_b0pf = -14.8832*dd4hep::cm;
  double y_b0pf = 0.0*dd4hep::cm;
  double z_b0pf = 575.314*dd4hep::cm;

  double radius_b0apf = 4.3*dd4hep::cm;   
  double length_b0apf = 149.352*dd4hep::cm;   //cm -- shorten by 6mm
  double rotation_b0apf = -0.025;  // radians
  double x_b0apf = -19.9248*dd4hep::cm;
  double y_b0apf = 0.0*dd4hep::cm;
  double z_b0apf = 774.957*dd4hep::cm;

  double radius_q1apf = 5.6*dd4hep::cm;   
  double length_q1apf = 185.699*dd4hep::cm; //cm -- shorten by 3mm  
  double rotation_q1apf = -0.0195;  // radians
  double x_q1apf = -25.0455*dd4hep::cm;
  double y_q1apf = 0.0*dd4hep::cm;
  double z_q1apf = 942.885*dd4hep::cm;

  double radius_q1bpf = 7.8*dd4hep::cm;   
  double length_q1bpf = 200.698*dd4hep::cm; //cm -- shorten by 3mm
  double rotation_q1bpf = -0.015;  // radians
  double x_q1bpf = -30.9852*dd4hep::cm;
  double y_q1bpf = 0.0*dd4hep::cm;
  double z_q1bpf = 1136.31*dd4hep::cm;

  double radius_q2pf = 13.15*dd4hep::cm;
  double length_q2pf = 419.495*dd4hep::cm;   //cm -- shorten by 5mm
  double rotation_q2pf = -0.0148;  // radians
  double x_q2pf = -40.4423*dd4hep::cm;
  double y_q2pf = 0.0*dd4hep::cm;
  double z_q2pf = 1446.73*dd4hep::cm;

  double radius_b1pf = 13.5*dd4hep::cm;
  double length_b1pf = 349.718*dd4hep::cm;   //cm -- shorten by 3mm
  double rotation_b1pf = -0.034;  // radians
  double x_b1pf = -49.4721*dd4hep::cm;
  double y_b1pf = 0.0*dd4hep::cm;
  double z_b1pf = 1831.59*dd4hep::cm;

  double radius_b1apf = 16.8*dd4hep::cm;   
  double length_b1apf = 199.725*dd4hep::cm; //cm -- shorten by 3mm  
  double rotation_b1apf = -0.025;  // radians
  double x_b1apf = -60.6686*dd4hep::cm;
  double y_b1apf = 0.0*dd4hep::cm;
  double z_b1apf = 2106.41*dd4hep::cm;



  // define shapes here

  Cone   b0pf_vacuum(length_b0pf / 2.0, 0.0, radius_b0pf, 0.0, radius_b0pf);
  Volume v_b0pf_vacuum("v_b0pf_vacuum", b0pf_vacuum, m_Vac);
  sdet.setAttributes(det, v_b0pf_vacuum, x_det.regionStr(), x_det.limitsStr(), vis_name);

  Cone   b0apf_vacuum(length_b0apf / 2.0, 0.0, radius_b0apf, 0.0, radius_b0apf);
  Volume v_b0apf_vacuum("v_b0apf_vacuum", b0apf_vacuum, m_Vac);
  sdet.setAttributes(det, v_b0apf_vacuum, x_det.regionStr(), x_det.limitsStr(), vis_name);

  Cone   q1apf_vacuum(length_q1apf / 2.0, 0.0, radius_q1apf, 0.0, radius_q1apf);
  Volume v_q1apf_vacuum("v_q1apf_vacuum", q1apf_vacuum, m_Vac);
  sdet.setAttributes(det, v_q1apf_vacuum, x_det.regionStr(), x_det.limitsStr(), vis_name);

  Cone   q1bpf_vacuum(length_q1bpf / 2.0, 0.0, radius_q1bpf, 0.0, radius_q1bpf);
  Volume v_q1bpf_vacuum("v_q1bpf_vacuum", q1bpf_vacuum, m_Vac);
  sdet.setAttributes(det, v_q1bpf_vacuum, x_det.regionStr(), x_det.limitsStr(), vis_name);

  Cone   q2pf_vacuum(length_q2pf / 2.0, 0.0, radius_q2pf, 0.0, radius_q2pf);
  Volume v_q2pf_vacuum("v_q2pf_vacuum", q2pf_vacuum, m_Vac);
  sdet.setAttributes(det, v_q2pf_vacuum, x_det.regionStr(), x_det.limitsStr(), vis_name);

  Cone   b1pf_vacuum(length_b1pf / 2.0, 0.0, radius_b1pf, 0.0, radius_b1pf);
  Volume v_b1pf_vacuum("v_b1pf_vacuum", b1pf_vacuum, m_Vac);
  sdet.setAttributes(det, v_b1pf_vacuum, x_det.regionStr(), x_det.limitsStr(), vis_name);

  Cone   b1apf_vacuum(length_b1apf / 2.0, 0.0, radius_b1apf, 0.0, radius_b1apf);
  Volume v_b1apf_vacuum("v_b1apf_vacuum", b1apf_vacuum, m_Vac);
  sdet.setAttributes(det, v_b1apf_vacuum, x_det.regionStr(), x_det.limitsStr(), vis_name);

  //----------------------------//

  auto pv_b0pf_vacuum =
      assembly.placeVolume(v_b0pf_vacuum, Transform3D(RotationY(rotation_b0pf), Position(x_b0pf, y_b0pf, z_b0pf)));
  pv_b0pf_vacuum.addPhysVolID("sector", 1);
  DetElement tube_de_1(sdet, "sector1_de", 1);
  tube_de_1.setPlacement(pv_b0pf_vacuum);

  auto pv_b0apf_vacuum =
      assembly.placeVolume(v_b0apf_vacuum, Transform3D(RotationY(rotation_b0apf), Position(x_b0apf, y_b0apf, z_b0apf)));
  pv_b0apf_vacuum.addPhysVolID("sector", 1);
  DetElement tube_de_2(sdet, "sector2_de", 1);
  tube_de_2.setPlacement(pv_b0apf_vacuum);

  auto pv_q1apf_vacuum =
      assembly.placeVolume(v_q1apf_vacuum, Transform3D(RotationY(rotation_q1apf), Position(x_q1apf, y_q1apf, z_q1apf)));
  pv_q1apf_vacuum.addPhysVolID("sector", 1);
  DetElement tube_de_3(sdet, "sector3_de", 1);
  tube_de_3.setPlacement(pv_q1apf_vacuum);

  auto pv_q1bpf_vacuum =
      assembly.placeVolume(v_q1bpf_vacuum, Transform3D(RotationY(rotation_q1bpf), Position(x_q1bpf, y_q1bpf, z_q1bpf)));
  pv_q1bpf_vacuum.addPhysVolID("sector", 1);
  DetElement tube_de_4(sdet, "sector4_de", 1);
  tube_de_4.setPlacement(pv_q1bpf_vacuum);

  auto pv_q2pf_vacuum =
      assembly.placeVolume(v_q2pf_vacuum, Transform3D(RotationY(rotation_q2pf), Position(x_q2pf, y_q2pf, z_q2pf)));
  pv_q2pf_vacuum.addPhysVolID("sector", 1);
  DetElement tube_de_5(sdet, "sector5_de", 1);
  tube_de_5.setPlacement(pv_q2pf_vacuum);

  auto pv_b1pf_vacuum =
      assembly.placeVolume(v_b1pf_vacuum, Transform3D(RotationY(rotation_b1pf), Position(x_b1pf, y_b1pf, z_b1pf)));
  pv_b1pf_vacuum.addPhysVolID("sector", 1);
  DetElement tube_de_6(sdet, "sector6_de", 1);
  tube_de_6.setPlacement(pv_q1apf_vacuum);

  auto pv_b1apf_vacuum =
      assembly.placeVolume(v_b1apf_vacuum, Transform3D(RotationY(rotation_b1apf), Position(x_b1apf, y_b1apf, z_b1apf)));
  pv_b1apf_vacuum.addPhysVolID("sector", 1);
  DetElement tube_de_7(sdet, "sector7_de", 1);
  tube_de_7.setPlacement(pv_b1apf_vacuum);


  pv_assembly = det.pickMotherVolume(sdet).placeVolume(assembly); //, posAndRot);
  pv_assembly.addPhysVolID("system", x_det.id()).addPhysVolID("barrel", 1);
  sdet.setPlacement(pv_assembly);
  assembly->GetShape()->ComputeBBox();
  return sdet;
}

DECLARE_DETELEMENT(magnetElementInnerVacuum, create_detector)
