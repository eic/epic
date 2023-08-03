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
  // Material   m_Cu    = det.material("Copper");
  // Material   m_Al    = det.material("Aluminum");
  Material m_Be     = det.material("Beryllium");
  Material m_SS     = det.material("StainlessSteel");
  Material m_vac    = det.material("Vacuum");
  string   vis_name = x_det.visStr();

  PlacedVolume pv_assembly;

  //xml::Component pos = x_det.position();
  // xml::Component rot   = x_det.rotation();
  
  //int numPipePieces = 4; //number of individual pipe sections 


  //double drift_beam_pipe_angle = -0.047666638;
  // double zPosShift             = 50.0; //cm

  double b0_hadron_tube_inner_r = 2.9;   // cm
  double b0_hadron_tube_outer_r = 3.1;   // cm
  double b0_hadron_tube_length  = 120.0; // cm
  
  
  double drift_hadron_section_1_inner_r = 19.5;
  double drift_hadron_section_1_outer_r = 20.5;
  //double drift_hadron_section_1_length  = 342.225466; // 393.4334363;
  /*
  double drift_hadron_section_2_inner_r = 19.5;
  double drift_hadron_section_2_outer_r = 20.5;
  double drift_hadron_section_2_length  = 300.0;
  */
  double drift_hadron_section_3_inner_r_ent = 19.5;
  double drift_hadron_section_3_outer_r_ent = 20.5;
  double drift_hadron_section_3_inner_r_ex  = 5.0;
  double drift_hadron_section_3_outer_r_ex  = 5.2;
  double drift_hadron_section_3_length      = 150.0;

  double drift_hadron_section_3_x = 0.0;
  double drift_hadron_section_3_z = 0.0;

  double drift_hadron_section_4_inner_r = 5.0;
  double drift_hadron_section_4_outer_r = 5.2;
  double drift_hadron_section_4_length  = 850.0;
  
  double drift_hadron_section_4_x = 0.0;
  double drift_hadron_section_4_z = 0.0;

  //Calculatate full drift region from line formula for proton orbit
  					
					    //(z, x)
  //double orbit_start[2] = {22.0623828, -0.6543372}; //meters!
  //double orbit_end[2] = {38.5445361, -1.4039456}; //meters!

  double orbit_start[2] = {22.07774534, -0.650777226};
  double orbit_end[2] = {38.54362489, -1.436245325};

  //calculate straight line formula x = slope*z + intercept
  
  double slope = (orbit_end[1]-orbit_start[1])/(orbit_end[0]-orbit_start[0]);
  double intercept = orbit_start[1]-(slope*orbit_start[0]);
  
  
  // This is the beam tube in the B0 magnet for the hadron beam  
  // doesn't use the slope information calculated before - it stands alone

  Tube   b0_hadron_tube(b0_hadron_tube_inner_r, b0_hadron_tube_outer_r, b0_hadron_tube_length / 2.0);
  Volume v_b0_hadron_tube("v_b0_hadron_tube", b0_hadron_tube, m_Be);
  sdet.setAttributes(det, v_b0_hadron_tube, x_det.regionStr(), x_det.limitsStr(), vis_name);

  //----------------------------------
  //    build drift beam pipe here
  //----------------------------------

  double z_start_pipe[3]    = {orbit_start[0]+0.03, 30.000,  31.500 };
  double z_end_pipe[3]      = {30.000,              31.500,  40.000 };
  
  
  for(int iSection = 0; iSection < 3; iSection++){
  	
	  double z_endpoint   = z_end_pipe[iSection]; //meters
	  double x_endpoint   = (slope*z_endpoint) + intercept;
  	  double x_startpoint = (slope*z_start_pipe[iSection]) + intercept;
  
	  double length = sqrt(pow(z_endpoint - z_start_pipe[iSection],2) + pow(x_endpoint - x_startpoint,2));
	  double z_center = (0.5*length + z_start_pipe[iSection])*cos(slope);
	  double x_center = (slope*z_center) + intercept;
  
  	  double entrance_r_inner  = 0.0; //drift_hadron_section_1_inner_r;
  	  double exit_radius_inner = 0.0;
  	  double entrance_r_outer  = 0.0; //drift_hadron_section_1_inner_r;
  	  double exit_radius_outer = 0.0;
	  
	  if(iSection < 1){ 
		  entrance_r_inner   = drift_hadron_section_1_inner_r;
	  	  exit_radius_inner  = drift_hadron_section_1_inner_r;
		  entrance_r_outer   = drift_hadron_section_1_outer_r; 
		  exit_radius_outer  = drift_hadron_section_1_outer_r;
		  //length = drift_hadron_section_1_length + drift_hadron_section_2_length;
	  }
	  if(iSection == 1){ 
		  entrance_r_inner   = drift_hadron_section_3_inner_r_ent;
	  	  exit_radius_inner  = drift_hadron_section_3_inner_r_ex;
		  entrance_r_outer   = drift_hadron_section_3_outer_r_ent; 
		  exit_radius_outer  = drift_hadron_section_3_outer_r_ex;
		  drift_hadron_section_3_x = x_center;
		  drift_hadron_section_3_z = z_center;
		  drift_hadron_section_3_length = length;
		  //x_center = -99.25250431/100.0; 
		  //z_center = 2924.185347/100.0;
		  //length = drift_hadron_section_3_length/100.0;
	  }
	  if(iSection == 2){ 
		  entrance_r_inner   = drift_hadron_section_4_inner_r;
		  exit_radius_inner  = drift_hadron_section_4_inner_r;
		  entrance_r_outer   = drift_hadron_section_4_outer_r; 
		  exit_radius_outer  = drift_hadron_section_4_outer_r;
		  drift_hadron_section_4_x = x_center;
		  drift_hadron_section_4_z = z_center;
		  drift_hadron_section_4_length = length;
		 // x_center =  -123.076799/100.0; 
		 // z_center = 3423.617428/100.0;
		  //length = drift_hadron_section_4_length/100.0;
	  }
	  
	  Cone drift_pipe((length*100.0) / 2.0, entrance_r_inner, entrance_r_outer, exit_radius_inner, exit_radius_outer);
	  
	  Volume v_pipe(Form("v_drift_tube_pipe_%d", iSection), drift_pipe, m_SS);
	  sdet.setAttributes(det, v_pipe, x_det.regionStr(), x_det.limitsStr(), vis_name);
    
	  auto pv_pipe = assembly.placeVolume(v_pipe, Transform3D(RotationY(slope), Position(100.0*x_center, 0.0, 100.0*z_center))); // 2353.06094)));
	  pv_pipe.addPhysVolID("sector", 1);
	  DetElement pipe_de(sdet, Form("sector_pipe_%d_de", iSection), 1);
	  pipe_de.setPlacement(pv_pipe);
  }
  

  
  // The tube that goes from B1apf to the start of the RP

  /*

  Tube   drift_tube_section_1(drift_hadron_section_1_inner_r, drift_hadron_section_1_outer_r,
                              drift_hadron_section_1_length / 2.0);
  Volume v_drift_tube_section_1("v_drift_tube_section_1", drift_tube_section_1, m_SS);
  sdet.setAttributes(det, v_drift_tube_section_1, x_det.regionStr(), x_det.limitsStr(), vis_name);

  // The tube that serves as a scattering chamber

  Tube   drift_tube_section_2(drift_hadron_section_2_inner_r, drift_hadron_section_2_outer_r,
                              drift_hadron_section_2_length / 2.0);
  Volume v_drift_tube_section_2("v_drift_tube_section_2", drift_tube_section_2, m_SS);
  sdet.setAttributes(det, v_drift_tube_section_2, x_det.regionStr(), x_det.limitsStr(), vis_name);

  Tube   drift_tube_vacuum_2(0.0, drift_hadron_section_2_inner_r, drift_hadron_section_2_length / 2.0);
  Volume v_drift_tube_vacuum_2("v_drift_tube_vacuum_2", drift_tube_vacuum_2, m_vac);
  sdet.setAttributes(det, v_drift_tube_vacuum_2, x_det.regionStr(), x_det.limitsStr(), vis_name);

  // The taper from the RP to last straight section

  Cone   drift_tube_section_3(drift_hadron_section_3_length / 2.0, drift_hadron_section_3_inner_r_ent,
                              drift_hadron_section_3_outer_r_ent, drift_hadron_section_3_inner_r_ex,
                              drift_hadron_section_3_outer_r_ex);
  Volume v_drift_tube_section_3("v_drift_tube_section_3", drift_tube_section_3, m_SS);
  sdet.setAttributes(det, v_drift_tube_section_3, x_det.regionStr(), x_det.limitsStr(), vis_name);

  // Final tube from taper to B2pf magnet

  Tube   drift_tube_section_4(drift_hadron_section_4_inner_r, drift_hadron_section_4_outer_r,
                              drift_hadron_section_4_length / 2.0);
  Volume v_drift_tube_section_4("v_drift_tube_section_4", drift_tube_section_4, m_SS);
  sdet.setAttributes(det, v_drift_tube_section_4, x_det.regionStr(), x_det.limitsStr(), vis_name);

  //----------------------------//

  auto pv_b0_hadron_tube =
      assembly.placeVolume(v_b0_hadron_tube, Transform3D(RotationY(-0.025), Position(pos.x(), pos.y(), pos.z())));
  pv_b0_hadron_tube.addPhysVolID("sector", 1);
  DetElement tube_de_1(sdet, "sector1_de", 1);
  tube_de_1.setPlacement(pv_b0_hadron_tube);

  // first tube section - right after b1apf - has same size as RP chamber, but keeping separate.
  auto pv_drift_tube_section_1 = assembly.placeVolume(
      v_drift_tube_section_1,
      Transform3D(RotationY(drift_beam_pipe_angle), Position(-73.23100294, 0.0, 2378.69291))); // 2353.06094)));
  pv_drift_tube_section_1.addPhysVolID("sector", 1);
  DetElement tube_de_2(sdet, "sector2_de", 1);
  tube_de_2.setPlacement(pv_drift_tube_section_1);
  */
   

  /*
  // Second section - RP scattering chamber -- keeping separate for now.
  auto pv_drift_tube_section_2 = assembly.placeVolume(
      v_drift_tube_section_2, Transform3D(RotationY(drift_beam_pipe_angle), Position(-88.5315717, 0.0, 2699.440911)));
  pv_drift_tube_section_2.addPhysVolID("sector", 1);
  DetElement tube_de_3(sdet, "sector3_de", 1);
  tube_de_3.setPlacement(pv_drift_tube_section_2);

  auto pv_drift_tube_vacuum_2 = assembly.placeVolume(
      v_drift_tube_vacuum_2, Transform3D(RotationY(drift_beam_pipe_angle), Position(-88.5315717, 0.0, 2699.440911)));
  pv_drift_tube_vacuum_2.addPhysVolID("sector", 1);
  DetElement tube_de_4(sdet, "sector4_de", 1);
  tube_de_4.setPlacement(pv_drift_tube_vacuum_2);

  // Third section -- tapered section acting as poor man's universal exit window.
  auto pv_drift_tube_section_3 = assembly.placeVolume(
      v_drift_tube_section_3, Transform3D(RotationY(drift_beam_pipe_angle), Position(-99.25250431, 0.0, 2924.185347)));
  pv_drift_tube_section_3.addPhysVolID("sector", 1);
  DetElement tube_de_5(sdet, "sector5_de", 1);
  tube_de_5.setPlacement(pv_drift_tube_section_3);

  auto pv_drift_tube_section_4 = assembly.placeVolume(
      v_drift_tube_section_4, Transform3D(RotationY(drift_beam_pipe_angle), Position(-123.076799, 0.0, 3423.617428)));
  pv_drift_tube_section_4.addPhysVolID("sector", 1);
  DetElement tube_de_6(sdet, "sector6_de", 1);
  tube_de_6.setPlacement(pv_drift_tube_section_4);
  */
  //make vacuum volumes here
  
  //last two entries are dummy numbers for now
  double z_start_points[7]    = {orbit_start[0]+0.03, 22.590,  24.590,   26.040, 28.040 , 20.0, 20.0 };
  double z_endpoints_array[7] = {22.499,              24.499,  25.990,   27.990, z_start_pipe[1] , 25.0, 25.0 };
  
  
  for(int iVac = 0; iVac < 7; iVac++){
  	
	  double z_endpoint   = z_endpoints_array[iVac]; //meters
	  double x_endpoint   = (slope*z_endpoint) + intercept;
  	  double x_startpoint = (slope*z_start_points[iVac]) + intercept;
  
	  double length = sqrt(pow(z_endpoint - z_start_points[iVac],2) + pow(x_endpoint - x_startpoint,2));
	  double z_center = (0.5*length + z_start_points[iVac])*cos(slope);
	  double x_center = (slope*z_center) + intercept;
  
  	  double entrance_r_inner  = 0.0; //drift_hadron_section_1_inner_r;
  	  double exit_radius_inner = 0.0;
	  
	  if(iVac < 5){ 
		  entrance_r_inner   = drift_hadron_section_1_inner_r;
	  	  exit_radius_inner  = drift_hadron_section_1_inner_r;
	  }
	  if(iVac == 5){ 
		  entrance_r_inner   = drift_hadron_section_1_inner_r;
		  exit_radius_inner  = drift_hadron_section_3_inner_r_ex;
		  x_center = drift_hadron_section_3_x;//-99.25250431/100.0;
		  z_center = drift_hadron_section_3_z;//2924.185347/100.0;
		  length = drift_hadron_section_3_length; ///100.0;
	  }
	  if(iVac == 6){ 
		  entrance_r_inner   = drift_hadron_section_4_inner_r;
		  exit_radius_inner  = drift_hadron_section_4_inner_r;
		  //x_center =  -123.076799/100.0; 
		  //z_center = 3423.617428/100.0;
		  //length = drift_hadron_section_4_length/100.0;
		  x_center = drift_hadron_section_4_x;
		  z_center = drift_hadron_section_4_z;
		  length = drift_hadron_section_4_length;
	  }
	  
	  Cone drift_vacuum((length*100.0) / 2.0, 0.0, entrance_r_inner-0.5, 0.0, exit_radius_inner-0.5);
	  
	  Volume v_vacuum(Form("v_drift_tube_vacuum_%d", iVac), drift_vacuum, m_vac);
	  sdet.setAttributes(det, v_vacuum, x_det.regionStr(), x_det.limitsStr(), "AnlBlue");
    
	  auto pv_vacuum = assembly.placeVolume(v_vacuum, Transform3D(RotationY(slope), Position(100.0*x_center, 0.0, 100.0*z_center))); // 2353.06094)));
	  pv_vacuum.addPhysVolID("sector", 1);
	  DetElement vacuum_de(sdet, Form("sector_vac_%d_de", iVac), 1);
	  vacuum_de.setPlacement(pv_vacuum);
  }

  // Transform3D posAndRot(RotationZYX(rot.z(), rot.y(), rot.x()), Position(pos.x(), pos.y(), pos.z()));
  // Transform3D posAndRot(RotationZYX(rot.z(), rot.y(), rot.x()), Position(x_position, y_position, z_position));

  pv_assembly = det.pickMotherVolume(sdet).placeVolume(assembly); //, posAndRot);
  pv_assembly.addPhysVolID("system", x_det.id()).addPhysVolID("barrel", 1);
  sdet.setPlacement(pv_assembly);
  assembly->GetShape()->ComputeBBox();
  return sdet;
}

DECLARE_DETELEMENT(hadronDownstreamBeamPipe, create_detector)
