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
static Ref_t create_detector(Detector& det, xml_h e, SensitiveDetector sens)  {

  using namespace ROOT::Math;
  xml_det_t  x_det     = e;
  string     det_name  = x_det.nameStr();
  Material   air       = det.air();
  DetElement sdet        (det_name,x_det.id());
  Assembly   assembly    (det_name+"_assembly");
  Material   m_Cu    = det.material("Copper");
  Material   m_Al    = det.material("Aluminum");
  Material   m_Be    = det.material("Beryllium");
  Material   m_SS    = det.material("StainlessSteel");
  string     vis_name  = x_det.visStr();

  PlacedVolume pv_assembly;

  xml::Component pos   = x_det.position();
  xml::Component rot   = x_det.rotation();

  /// hard-code defintion here, then refine and make more general

  double drift_beam_pipe_angle = -0.0475492667;

  double b0_hadron_tube_inner_r = 2.9; // cm
  double b0_hadron_tube_outer_r = 3.1; //cm
  double b0_hadron_tube_length  = 120.0; //cm

  double drift_hadron_section_1_inner_r = 19.5;
  double drift_hadron_section_1_outer_r = 20.5;
  double drift_hadron_section_1_length  = 392.0; //393.4334363;

  double drift_hadron_section_2_inner_r = 19.5;
  double drift_hadron_section_2_outer_r = 20.5;
  double drift_hadron_section_2_length  = 300.0;

  double drift_hadron_section_3_inner_r_ent = 19.5;
  double drift_hadron_section_3_outer_r_ent = 20.5;
  double drift_hadron_section_3_inner_r_ex  = 5.0;
  double drift_hadron_section_3_outer_r_ex  = 5.2;
  double drift_hadron_section_3_length  = 250.0;

  double drift_hadron_section_4_inner_r = 5.0;
  double drift_hadron_section_4_outer_r = 5.2;
  double drift_hadron_section_4_length  = 900.0;

  //This is the beam tube in the B0 magnet for the hadron beam

  Tube b0_hadron_tube(b0_hadron_tube_inner_r, b0_hadron_tube_outer_r, b0_hadron_tube_length/2.0);
  Volume v_b0_hadron_tube("v_b0_hadron_tube", b0_hadron_tube, m_Be);
  sdet.setAttributes(det, v_b0_hadron_tube  , x_det.regionStr(), x_det.limitsStr(), vis_name);

  //The tube that goes from B0pf to the start of the RP

  Tube drift_tube_section_1(drift_hadron_section_1_inner_r, drift_hadron_section_1_outer_r, drift_hadron_section_1_length/2.0);
  Volume v_drift_tube_section_1("v_drift_tube_section_1", drift_tube_section_1, m_SS);
  sdet.setAttributes(det, v_drift_tube_section_1  , x_det.regionStr(), x_det.limitsStr(), vis_name);

  //The tube that serves as a scattering chamber

  Tube drift_tube_section_2(drift_hadron_section_2_inner_r, drift_hadron_section_2_outer_r, drift_hadron_section_2_length/2.0);
  Volume v_drift_tube_section_2("v_drift_tube_section_2", drift_tube_section_2, m_SS);
  sdet.setAttributes(det, v_drift_tube_section_2  , x_det.regionStr(), x_det.limitsStr(), vis_name);

  //The taper from the RP to last straight section

  Cone drift_tube_section_3(drift_hadron_section_3_length/2.0, drift_hadron_section_3_inner_r_ent, drift_hadron_section_3_outer_r_ent,drift_hadron_section_3_inner_r_ex, drift_hadron_section_3_outer_r_ex);
  Volume v_drift_tube_section_3("v_drift_tube_section_3", drift_tube_section_3, m_SS);
  sdet.setAttributes(det, v_drift_tube_section_3  , x_det.regionStr(), x_det.limitsStr(), vis_name);

  //Final tube from taper to B2pf magnet

  Tube drift_tube_section_4(drift_hadron_section_4_inner_r, drift_hadron_section_4_outer_r, drift_hadron_section_4_length/2.0);
  Volume v_drift_tube_section_4("v_drift_tube_section_4", drift_tube_section_4, m_SS);
  sdet.setAttributes(det, v_drift_tube_section_4  , x_det.regionStr(), x_det.limitsStr(), vis_name);

  //----------------------------//

  auto pv_b0_hadron_tube = assembly.placeVolume( v_b0_hadron_tube, Transform3D(RotationY(-0.025), Position(pos.x(), pos.y(), pos.z())));
  pv_b0_hadron_tube.addPhysVolID("sector", 1);
  DetElement tube_de_1(sdet, "sector1_de", 1);
  tube_de_1.setPlacement(pv_b0_hadron_tube);

  auto pv_drift_tube_section_1 = assembly.placeVolume( v_drift_tube_section_1, Transform3D(RotationY(drift_beam_pipe_angle), Position(-71.269416, 0.0, 2353.7))); //2353.06094)));
  pv_drift_tube_section_1.addPhysVolID("sector", 1);
  DetElement tube_de_2(sdet, "sector2_de", 1);
  tube_de_2.setPlacement(pv_drift_tube_section_1);
  
  auto pv_drift_tube_section_2 = assembly.placeVolume( v_drift_tube_section_2, Transform3D(RotationY(drift_beam_pipe_angle), Position(-87.74933, 0.0, 2699.38578)));
  pv_drift_tube_section_2.addPhysVolID("sector", 1);
  DetElement tube_de_3(sdet, "sector3_de", 1);
  tube_de_3.setPlacement(pv_drift_tube_section_2);

  auto pv_drift_tube_section_3 = assembly.placeVolume( v_drift_tube_section_3, Transform3D(RotationY(drift_beam_pipe_angle), Position(-100.820452, 0.0, 2974.07496)));
  pv_drift_tube_section_3.addPhysVolID("sector", 1);
  DetElement tube_de_4(sdet, "sector4_de", 1);
  tube_de_4.setPlacement(pv_drift_tube_section_3);

  auto pv_drift_tube_section_4 = assembly.placeVolume( v_drift_tube_section_4, Transform3D(RotationY(drift_beam_pipe_angle), Position(-128.150979, 0.0, 3548.42507)));
  pv_drift_tube_section_4.addPhysVolID("sector", 1);
  DetElement tube_de_5(sdet, "sector5_de", 1);
  tube_de_5.setPlacement(pv_drift_tube_section_4);

  //Transform3D posAndRot(RotationZYX(rot.z(), rot.y(), rot.x()), Position(pos.x(), pos.y(), pos.z()));
  //Transform3D posAndRot(RotationZYX(rot.z(), rot.y(), rot.x()), Position(x_position, y_position, z_position));  

  pv_assembly = det.pickMotherVolume(sdet).placeVolume(assembly); //, posAndRot);
  pv_assembly.addPhysVolID("system",x_det.id()).addPhysVolID("barrel",1);
  sdet.setPlacement(pv_assembly);
  assembly->GetShape()->ComputeBBox() ;
  return sdet;
}

DECLARE_DETELEMENT(hadronDownstreamBeamPipe,create_detector)
