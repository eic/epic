// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2024 Alex Jentsch

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
	Material m_Al     = det.material("Aluminum");
	Material m_Be     = det.material("Beryllium");
	Material m_SS     = det.material("StainlessSteel");
	//Material m_vac    = det.material("Vacuum");
	string   vis_name = x_det.visStr();

	PlacedVolume pv_assembly;


	double b0_hadron_tube_inner_r = 2.9;   // cm
	double b0_hadron_tube_outer_r = 3.1;   // cm
	double b0_hadron_tube_length  = 120.0; // cm


	//double scatChamberRadius = 72.0;
	double pipeThickness     = 5.0;
  
  
	Double_t length[10];
	Double_t innerXRadius[10];  
	Double_t innerYRadius[10]; 
	Double_t outerXRadius[10];
	Double_t outerYRadius[10];
	Double_t rotationAngle[10];
	Double_t xCenter[10];
	Double_t yCenter[10];
	Double_t zCenter[10];

  
	//TGeoConeSeg *piece[10];
	//TGeoVolume  *vpiece[10];

	double globRotationAngle = -0.0454486856; //This is the angle of the proton orbit from the end of B1APF to the beginning of B2PF

	double b1APFEndPoint_z = 22.0623828*1000; //in mm
	double b1APFEndPoint_x = 0.6543372*1000; //in mm

	double tmp_endpoint_z = 0.0;
	double tmp_endpoint_x = 0.0;

	//forumula -> Z = b1APFEndPoint_z+((0.5*elementLengt)*Cos(globRotationAngle))
	//forumula -> X = b1APFEndPoint_z+((0.5*elementLengt)*Sin(globRotationAngle))

	//------------------------------------------------------------------------------------
	//primary pipe after B1APF, before neutral exit window + transition to smaller pipe
	//rectangular cross-section!!!!!
	//------------------------------------------------------------------------------------
  
  	length[0]            = 7615.486; // from VPC drawings
	innerXRadius[0]      = 275.0; //
	innerYRadius[0]      = 175.0; //
	outerXRadius[0]      = innerXRadius[0] + pipeThickness; //
	outerYRadius[0]      = innerYRadius[0] + pipeThickness; //

  	xCenter[0]           = -1*(b1APFEndPoint_x+((0.5*length[0])*TMath::Sin(-globRotationAngle))); //-745.20328;
  	yCenter[0]           = 0.0;
  	zCenter[0]           = (b1APFEndPoint_z+((0.5*length[0])*TMath::Cos(-globRotationAngle)));//24060.3176;
  	rotationAngle[0]     = globRotationAngle;
 
	cout << "first component: " << " z = " << zCenter[0] << " x = " << xCenter[0] << endl;

	//NEEDS TO BE TGeoBBox!!!!!!!!

	tmp_endpoint_z = zCenter[0]+((0.5*length[0])*TMath::Cos(-globRotationAngle));
	tmp_endpoint_x = -1*xCenter[0]+((0.5*length[0])*TMath::Sin(-globRotationAngle));
	
	cout << "large section ends at (z,x) = " << tmp_endpoint_z << ", " << tmp_endpoint_x << endl;

	double neutralExitWindowZ = tmp_endpoint_z - 5.0;
	double neutralExitWindowX = -0.025*neutralExitWindowZ;
	double windowRadius = neutralExitWindowZ*TMath::Tan(0.004);

  	//------------------------------------------------------------------------------------
  	//first small pipe section, between primary vessel and RP station 1
	//rectangular cross-section!!!!!
	//------------------------------------------------------------------------------------
  
  	length[1]            = 2780.273; // from VPC drawings
	innerXRadius[1]      = 150.0; //
	innerYRadius[1]      = 30.0; //
	outerXRadius[1]      = innerXRadius[1] + pipeThickness; //
	outerYRadius[1]      = innerYRadius[1] + pipeThickness; //

	xCenter[1]           = -1*(tmp_endpoint_x+((0.5*length[1])*TMath::Sin(-globRotationAngle)));//-972.36849;
	yCenter[1]           = 0.0;
	zCenter[1]           = tmp_endpoint_z+((0.5*length[1])*TMath::Cos(-globRotationAngle));//29055.1545;
	rotationAngle[1]     = globRotationAngle;
 
	cout << "second component: " << " z = " << zCenter[1] << " x = " << xCenter[1] << endl;

	tmp_endpoint_z = zCenter[1]+((0.5*length[1])*TMath::Cos(-globRotationAngle));
  	tmp_endpoint_x = -1*xCenter[1]+((0.5*length[1])*TMath::Sin(-globRotationAngle));   

  	//------------------------------------------------------------------------------------
  	//First roman pots scattering chamber
	//------------------------------------------------------------------------------------
  
  	length[2]            = 200; // from VPC drawings
	innerXRadius[2]      = 200.0; //
	innerYRadius[2]      = 125.0; //
	outerXRadius[2]      = innerXRadius[2] + pipeThickness; //
	outerYRadius[2]      = innerYRadius[2] + pipeThickness; //

  	xCenter[2]           = -1*(tmp_endpoint_x+((0.5*length[2])*TMath::Sin(-globRotationAngle)));//-972.36849;
  	yCenter[2]           = 0.0;
  	zCenter[2]           = tmp_endpoint_z+((0.5*length[2])*TMath::Cos(-globRotationAngle));//29055.1545;
  	rotationAngle[2]     = globRotationAngle;
 
	cout << "third component: " << " z = " << zCenter[2] << " x = " << xCenter[2] << endl;

	tmp_endpoint_z = zCenter[2]+((0.5*length[2])*TMath::Cos(-globRotationAngle));
  	tmp_endpoint_x = -1*xCenter[2]+((0.5*length[2])*TMath::Sin(-globRotationAngle));  


  	//------------------------------------------------------------------------------------
  	//pipe between RP 1 and RP 2 stations
	//rectangular cross-section!!!!!
	//------------------------------------------------------------------------------------
  
  	length[3]            = 1500.0; // from VPC drawings
	innerXRadius[3]      = 150.0; //
	innerYRadius[3]      = 30.0; //
	outerXRadius[3]      = innerXRadius[3] + pipeThickness; //
	outerYRadius[3]      = innerYRadius[3] + pipeThickness; //

  	xCenter[3]           = -1*(tmp_endpoint_x+((0.5*length[3])*TMath::Sin(-globRotationAngle)));//-972.36849;
  	yCenter[3]           = 0.0;
  	zCenter[3]           = tmp_endpoint_z+((0.5*length[3])*TMath::Cos(-globRotationAngle));//29055.1545;
  	rotationAngle[3]     = globRotationAngle;
 
	cout << "fourth component: " << " z = " << zCenter[3] << " x = " << xCenter[3] << endl;

	tmp_endpoint_z = zCenter[3]+((0.5*length[3])*TMath::Cos(-globRotationAngle));
  	tmp_endpoint_x = -1*xCenter[3]+((0.5*length[3])*TMath::Sin(-globRotationAngle));  


  	//------------------------------------------------------------------------------------
  	//second roman pots scattering chamber
	//------------------------------------------------------------------------------------
  
  	length[4]            = 200; // from VPC drawings
	innerXRadius[4]      = 200.0; //
	innerYRadius[4]      = 125.0; //
	outerXRadius[4]      = innerXRadius[4] + pipeThickness; //
	outerYRadius[4]      = innerYRadius[4] + pipeThickness; //

  	xCenter[4]           = -1*(tmp_endpoint_x+((0.5*length[4])*TMath::Sin(-globRotationAngle)));//-972.36849;
  	yCenter[4]           = 0.0;
  	zCenter[4]           = tmp_endpoint_z+((0.5*length[4])*TMath::Cos(-globRotationAngle));//29055.1545;
  	rotationAngle[4]     = globRotationAngle;
 
	cout << "fifth component: " << " z = " << zCenter[4] << " x = " << xCenter[4] << endl;


  	//------------------------------------------------------------------------------------
  	// Inner pipe to subtract from everything for proper entrance/exits
	//------------------------------------------------------------------------------------
  
  	length[5]            = 20000.0; // from VPC drawings
	innerXRadius[5]      = 150.0; //
	innerYRadius[5]      = 30.0; //
	outerXRadius[5]      = innerXRadius[5] + pipeThickness; //
	outerYRadius[5]      = innerYRadius[5] + pipeThickness; //

  	xCenter[5]           = -1*(b1APFEndPoint_x+((0.5*length[0])*TMath::Sin(-globRotationAngle))); //-745.20328;
  	yCenter[5]           = 0.0;
  	zCenter[5]           = (b1APFEndPoint_z+((0.5*length[0])*TMath::Cos(-globRotationAngle)));//24060.3176;
  	rotationAngle[5]     = globRotationAngle;
 

	//tmp_endpoint_z = zCenter[2]+((0.5*length[2])*TMath::Cos(-globRotationAngle));
  	//tmp_endpoint_x = -1*xCenter[2]+((0.5*length[2])*TMath::Sin(-globRotationAngle));  


  ///-------------------------------------------


  	int pieceIdx = 0;
	
	/*
  
    Cone drift_pipe((length*100.0) / 2.0, entrance_r_inner, entrance_r_outer, exit_radius_inner, exit_radius_outer);

    Volume v_pipe(Form("v_drift_tube_pipe_%d", iSection), drift_pipe, m_SS);
    sdet.setAttributes(det, v_pipe, x_det.regionStr(), x_det.limitsStr(), vis_name);

    auto pv_pipe = assembly.placeVolume(v_pipe, Transform3D(RotationY(slope), Position(100.0*x_center, 0.0, 100.0*z_center))); // 2353.06094)));
    pv_pipe.addPhysVolID("sector", 1);
    DetElement pipe_de(sdet, Form("sector_pipe_%d_de", iSection), 1);
    pipe_de.setPlacement(pv_pipe);
  
 
  	double origin[3] = {0.0, 0.0, 0.0};
	*/
	//TGeoBBox * inner_subtract_outer = new TGeoBBox("inner_subtract_outer_0", 0.1*outerXRadius[5], 0.1*outerYRadius[5], 0.1*length[5]/2, origin);
  	//TGeoBBox * inner_subtract_inner = new TGeoBBox("inner_subtract_inner_0",  0.1*innerXRadius[5], 0.1*innerYRadius[5], 0.1*(length[5]+5.0)/2, origin);

	//TGeoBBox * RP_subtract_outer = new TGeoBBox("RP_subtract_outer_0", 0.1*outerXRadius[1], 0.1*outerYRadius[1], 0.1*(length[2]+5.0)/2, origin);
  	//TGeoBBox * RP_subtract_inner = new TGeoBBox("RP_subtract_inner_0",  0.1*innerXRadius[1], 0.1*innerYRadius[1], 0.1*(length[2]+10.0)/2, origin);

	Box inner_subtract_outer(0.1*outerXRadius[5], 0.1*outerYRadius[5], 0.1*length[5]/2);
	Box inner_subtract_inner(0.1*innerXRadius[5], 0.1*innerYRadius[5], 0.1*(length[5]+5.0)/2);
	
	Box RP_subtract_outer(0.1*outerXRadius[1], 0.1*outerYRadius[1], 0.1*(length[2]+5.0)/2);
	Box RP_subtract_inner(0.1*innerXRadius[1], 0.1*innerYRadius[1], 0.1*(length[2]+10.0)/2);

	//Begin building volumes here

	pieceIdx = 0; //Larger, rectangular pipe transporting proton and neutral envelopes (neutral exit window and transfer to smaller proton line at the end)

	Box pipeAfterB1APF_outer(0.1*outerXRadius[pieceIdx], 0.1*outerYRadius[pieceIdx], 0.1*length[pieceIdx]/2);
	Box pipeAfterB1APF_inner(0.1*innerXRadius[pieceIdx], 0.1*innerYRadius[pieceIdx], 0.1*(length[pieceIdx]+5.0)/2);
	SubtractionSolid pipeAfterB1APF(pipeAfterB1APF_outer, pipeAfterB1APF_inner);
	
    Volume v_pipeAfterB1APF(Form("v_pipeAfterB1APF_%d", pieceIdx), pipeAfterB1APF, m_SS);
    sdet.setAttributes(det, v_pipeAfterB1APF, x_det.regionStr(), x_det.limitsStr(), vis_name);

    auto pv_pipe_0 = assembly.placeVolume(v_pipeAfterB1APF, Transform3D(RotationY(-0.025), Position(0.1 * xCenter[pieceIdx] + 5.0, 0.1 *  yCenter[pieceIdx], 0.1 * zCenter[pieceIdx]))); // 2353.06094)));
    pv_pipe_0.addPhysVolID("sector", 1);
    DetElement pipe_de_0(sdet, Form("sector_pipe_%d_de", pieceIdx), 1);
    pipe_de_0.setPlacement(pv_pipe_0);
	
	pieceIdx = 1; //smaller rectangular pipe for the protons

	Box pipe_1_outer(0.1*outerXRadius[pieceIdx], 0.1*outerYRadius[pieceIdx], 0.1*length[pieceIdx]/2);
	Box pipe_1_inner(0.1*innerXRadius[pieceIdx], 0.1*innerYRadius[pieceIdx], 0.1*(length[pieceIdx]+5.0)/2);
	SubtractionSolid pipe_1(pipe_1_outer, pipe_1_inner);
	
    Volume v_pipe_1(Form("v_pipe_1_%d", pieceIdx), pipe_1, m_SS);
    sdet.setAttributes(det, v_pipe_1, x_det.regionStr(), x_det.limitsStr(), vis_name);

    auto pv_pipe_1 = assembly.placeVolume(v_pipe_1, Transform3D(RotationY(rotationAngle[pieceIdx]), Position(0.1 * xCenter[pieceIdx], 0.1 *  yCenter[pieceIdx], 0.1 * zCenter[pieceIdx]))); // 2353.06094)));
    pv_pipe_1.addPhysVolID("sector", 1);
    DetElement pipe_de_1(sdet, Form("sector_pipe_%d_de", pieceIdx), 1);
    pipe_de_1.setPlacement(pv_pipe_1);

	//------------------------------------------------------------------------

	pieceIdx = 2; //first roman pots scattering chamber
	
	Box box_rp_station_1_outer(0.1*outerXRadius[pieceIdx], 0.1*outerYRadius[pieceIdx], 0.1*length[pieceIdx]/2);
	Box box_rp_station_1_inner(0.1*innerXRadius[pieceIdx], 0.1*innerYRadius[pieceIdx], 0.1*length[pieceIdx]/2);
	SubtractionSolid tmp(box_rp_station_1_outer, box_rp_station_1_inner);
	SubtractionSolid rpStation1(tmp, RP_subtract_outer);
	
    Volume v_rpStation1(Form("v_rpStation1_%d", pieceIdx), rpStation1, m_SS);
    sdet.setAttributes(det, v_rpStation1, x_det.regionStr(), x_det.limitsStr(), vis_name);

    auto pv_pipe_2 = assembly.placeVolume(v_rpStation1, Transform3D(RotationY(rotationAngle[pieceIdx]), Position(0.1 * xCenter[pieceIdx], 0.1 *  yCenter[pieceIdx], 0.1 * zCenter[pieceIdx]))); // 2353.06094)));
    pv_pipe_2.addPhysVolID("sector", 1);
    DetElement pipe_de_2(sdet, Form("sector_pipe_%d_de", pieceIdx), 1);
    pipe_de_2.setPlacement(pv_pipe_2);

	//---------------------------------------------------------------------
	
	pieceIdx = 3;

	Box pipe_3_outer(0.1*outerXRadius[pieceIdx], 0.1*outerYRadius[pieceIdx], 0.1*length[pieceIdx]/2);
	Box pipe_3_inner(0.1*innerXRadius[pieceIdx], 0.1*innerYRadius[pieceIdx], 0.1*(length[pieceIdx]+5.0)/2);
	SubtractionSolid pipe_3(pipe_1_outer, pipe_1_inner);
	
    Volume v_pipe_3(Form("v_pipe_3_%d", pieceIdx), pipe_3, m_SS);
    sdet.setAttributes(det, v_pipe_3, x_det.regionStr(), x_det.limitsStr(), vis_name);

    auto pv_pipe_3 = assembly.placeVolume(v_pipe_3, Transform3D(RotationY(rotationAngle[pieceIdx]), Position(0.1 * xCenter[pieceIdx], 0.1 *  yCenter[pieceIdx], 0.1 * zCenter[pieceIdx]))); // 2353.06094)));
    pv_pipe_3.addPhysVolID("sector", 1);
    DetElement pipe_de_3(sdet, Form("sector_pipe_%d_de", pieceIdx), 1);
    pipe_de_3.setPlacement(pv_pipe_3);

	//---------------------------------------------------------------------

	pieceIdx = 4; //Second roman pots scattering chamber

	Box box_rp_station_2_outer(0.1*outerXRadius[pieceIdx], 0.1*outerYRadius[pieceIdx], 0.1*length[pieceIdx]/2);
	Box box_rp_station_2_inner(0.1*innerXRadius[pieceIdx], 0.1*innerYRadius[pieceIdx], 0.1*length[pieceIdx]/2);
	SubtractionSolid tmp_2(box_rp_station_2_outer, box_rp_station_2_inner);
	SubtractionSolid rpStation2(tmp_2, RP_subtract_outer);
	
    Volume v_rpStation2(Form("v_rpStation2_%d", pieceIdx), rpStation2, m_SS);
    sdet.setAttributes(det, v_rpStation2, x_det.regionStr(), x_det.limitsStr(), vis_name);

    auto pv_pipe_4 = assembly.placeVolume(v_rpStation2, Transform3D(RotationY(rotationAngle[pieceIdx]), Position(0.1 * xCenter[pieceIdx], 0.1 *  yCenter[pieceIdx], 0.1 * zCenter[pieceIdx]))); // 2353.06094)));
    pv_pipe_4.addPhysVolID("sector", 1);
    DetElement pipe_de_4(sdet, Form("sector_pipe_%d_de", pieceIdx), 1);
    pipe_de_4.setPlacement(pv_pipe_4);
	

  	// This is the beam tube in the B0 magnet for the hadron beam
  	// doesn't use the slope information calculated before - it stands alone

	pieceIdx = 5;

  	Tube   b0_hadron_tube(b0_hadron_tube_inner_r, b0_hadron_tube_outer_r, b0_hadron_tube_length / 2.0);
  	Volume v_b0_hadron_tube("v_b0_hadron_tube", b0_hadron_tube, m_Be);
  	sdet.setAttributes(det, v_b0_hadron_tube, x_det.regionStr(), x_det.limitsStr(), vis_name);

   
    auto pv_pipe_5 = assembly.placeVolume(v_b0_hadron_tube, Transform3D(RotationY(-0.025), Position(-16.5, 0.0, 640.0))); // 2353.06094)));
    pv_pipe_5.addPhysVolID("sector", 1);
    DetElement pipe_de_5(sdet, Form("sector_pipe_%d_de", pieceIdx), 1);
    pipe_de_5.setPlacement(pv_pipe_5);
	
	pieceIdx = 6; //neutral exit window

  	Tube   neutral_exit_window(0.0, 0.1*windowRadius, 1.0); // 1.0cm thick
  	Volume v_neutral_exit_window("v_neutral_exit_window", neutral_exit_window, m_Al);
  	sdet.setAttributes(det, v_neutral_exit_window, x_det.regionStr(), x_det.limitsStr(), vis_name);

	
    auto pv_pipe_6 = assembly.placeVolume(v_neutral_exit_window, Transform3D(RotationY(-0.025), Position(0.1*neutralExitWindowX, 0.0, 0.1*neutralExitWindowZ))); // 2353.06094)));
    pv_pipe_6.addPhysVolID("sector", 1);
    DetElement pipe_de_6(sdet, Form("sector_pipe_%d_de", pieceIdx), 1);
    pipe_de_6.setPlacement(pv_pipe_6);

  	//----------------------------------
  	//    build drift beam pipe here
  	//----------------------------------

	/*

  	double z_start_pipe[3]    = {orbit_start[0], 30.000,  31.500 };
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
                  length = length - 0.04;
          }
          if(iSection == 1){
                  entrance_r_inner   = drift_hadron_section_3_inner_r_ent;
                  exit_radius_inner  = drift_hadron_section_3_inner_r_ex;
                  entrance_r_outer   = drift_hadron_section_3_outer_r_ent;
                  exit_radius_outer  = drift_hadron_section_3_outer_r_ex;
                  drift_hadron_section_3_x = x_center;
                  drift_hadron_section_3_z = z_center;
                  drift_hadron_section_3_length = length - 0.02;
                  //old numbers commented out for reference - A. Jentsch
                  //length = length - 0.02;
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
                  drift_hadron_section_4_length = length - 0.02;
                  length = length - 0.02;
                  //old numbers commented out for reference - A. Jentsch
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
	*/
  	//------------------------------
  	//  make vacuum volumes here
  	//------------------------------

  	/*
  	//last two entries are dummy numbers for now
  	double z_start_points[7]    = {orbit_start[0]+0.03, 22.590,  24.590,   26.055, 28.055 , 20.0, 20.0 };
  	double z_endpoints_array[7] = {22.499,              24.499,  25.980,   27.980, z_start_pipe[1] , 25.0, 25.0 };


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
                  length = drift_hadron_section_3_length - 0.02; ///100.0;
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
	*/
 	 // Transform3D posAndRot(RotationZYX(rot.z(), rot.y(), rot.x()), Position(pos.x(), pos.y(), pos.z()));
  	// Transform3D posAndRot(RotationZYX(rot.z(), rot.y(), rot.x()), Position(x_position, y_position, z_position));

  	pv_assembly = det.pickMotherVolume(sdet).placeVolume(assembly); //, posAndRot);
  	pv_assembly.addPhysVolID("system", x_det.id()).addPhysVolID("barrel", 1);
  	sdet.setPlacement(pv_assembly);
  	assembly->GetShape()->ComputeBBox();
  	return sdet;
}

DECLARE_DETELEMENT(forwardBeamPipeBrazil, create_detector)
