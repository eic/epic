// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2024 Alex Jentsch

#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/Printout.h"
#include "TMath.h"
#include <XML/Helper.h>

using namespace std;
using namespace dd4hep;

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
	Material m_vac    = det.material("Vacuum");
	string   vis_name = x_det.visStr();

	PlacedVolume pv_assembly;

	double b0_hadron_tube_inner_r = 2.9 * dd4hep::cm;   
	double b0_hadron_tube_outer_r = 3.1 * dd4hep::cm;   
	double b0_hadron_tube_length  = 120.0 * dd4hep::cm; 

	double pipeThickness     = 5.0;
  
	//Double_t length[10];
	//Double_t innerXRadius[10];  
	//Double_t innerYRadius[10]; 
	//Double_t outerXRadius[10];
	//Double_t outerYRadius[10];
	//Double_t rotationAngle[10];
	//Double_t xCenter[10];
	//Double_t yCenter[10];
	//Double_t zCenter[10];

	struct beampipe_dimensions_t {
    	Double_t length;
    	Double_t innerXRadius;
    	Double_t innerYRadius;
    	Double_t outerXRadius;
    	Double_t outerYRadius;
    	Double_t xCenter;
    	Double_t yCenter;
    	Double_t zCenter;
		Double_t rotationAngle;
	};
        
	std::vector<beampipe_dimensions_t> beampipe_dimensions;

	std::cout << "Empty vector of structs created..." << std::endl;
  
	double globRotationAngle = -0.0454486856; //This is the angle of the proton orbit from the end of B1APF to the beginning of B2PF

	double b1APFEndPoint_z = 22062.3828; //location of proton orbit at b1APF exit -- in mm
	double b1APFEndPoint_x = 654.3372; //location of proton orbit at b1APF exit -- in mm

	double tmp_endpoint_z = 0.0;
	double tmp_endpoint_x = 0.0;

	//forumula -> Z = b1APFEndPoint_z+((0.5*elementLengt)*Cos(globRotationAngle))
	//forumula -> X = b1APFEndPoint_z+((0.5*elementLengt)*Sin(globRotationAngle))

	//------------------------------------------------------------------------------------
	//primary pipe after B1APF, before neutral exit window + transition to smaller pipe
	//rectangular cross-section!!!!!
	//------------------------------------------------------------------------------------
  
  	//length[0]            = 7615.486; //from VPC drawings, in mm
	//innerXRadius[0]      = 275.0;    //
	//innerYRadius[0]      = 175.0;    //
	//outerXRadius[0]      = innerXRadius[0] + pipeThickness; //
	//outerYRadius[0]      = innerYRadius[0] + pipeThickness; //

  	//xCenter[0]           = -1*(b1APFEndPoint_x+((0.5*length[0])*TMath::Sin(-globRotationAngle))); //-745.20328;
	//xCenter[0]           = xCenter[0]; //shift of center location to account for offset
  	//yCenter[0]           = 0.0;
  	//zCenter[0]           = (b1APFEndPoint_z+((0.5*length[0])*TMath::Cos(-globRotationAngle)));//24060.3176;
  	//rotationAngle[0]     = globRotationAngle;

	beampipe_dimensions.push_back({

		.length            = 7615.486, //from VPC drawings, in mm
    	.innerXRadius      = 275.0,    //
    	.innerYRadius      = 175.0,    //
    	.outerXRadius      = 0.0, //
    	.outerYRadius      = 0.0, //
		.xCenter           = 0.0, //-745.20328;
    	.yCenter           = 0.0,
    	.zCenter           = 0.0,//24060.3176;
    	.rotationAngle     = globRotationAngle
	});
	
	beampipe_dimensions[0].outerXRadius      = beampipe_dimensions[0].innerXRadius + pipeThickness, //
    beampipe_dimensions[0].outerYRadius      = beampipe_dimensions[0].innerYRadius + pipeThickness, //
	beampipe_dimensions[0].xCenter           = -1*(b1APFEndPoint_x+((0.5*beampipe_dimensions[0].length)*TMath::Sin(-globRotationAngle))), //-745.20328;
    beampipe_dimensions[0].yCenter           = 0.0,
    beampipe_dimensions[0].zCenter           = (b1APFEndPoint_z+((0.5*beampipe_dimensions[0].length)*TMath::Cos(-globRotationAngle))),//24060.3176;
	 
	tmp_endpoint_z = beampipe_dimensions[0].zCenter+((0.5*beampipe_dimensions[0].length)*TMath::Cos(-globRotationAngle));
	tmp_endpoint_x = -1*beampipe_dimensions[0].xCenter+((0.5*beampipe_dimensions[0].length)*TMath::Sin(-globRotationAngle));

	double windowRadius = 110.0;


  	//------------------------------------------------------------------------------------
  	//first small pipe section, between primary vessel and RP station 1
	//rectangular cross-section!!!!!
	//------------------------------------------------------------------------------------
  
  	//length[1]            = 2780.273; // from VPC drawings
	//innerXRadius[1]      = 150.0; //
	//innerYRadius[1]      = 30.0; //
	//outerXRadius[1]      = innerXRadius[1] + pipeThickness; //
	//outerYRadius[1]      = innerYRadius[1] + pipeThickness; //

	//xCenter[1]           = -1*(tmp_endpoint_x+((0.5*length[1])*TMath::Sin(-globRotationAngle)));//-972.36849;
	//yCenter[1]           = 0.0;
	//zCenter[1]           = tmp_endpoint_z+((0.5*length[1])*TMath::Cos(-globRotationAngle));//29055.1545;
	//rotationAngle[1]     = globRotationAngle;

	beampipe_dimensions.push_back({ 

		.length            = 2780.273, // from VPC drawings
    	.innerXRadius      = 150.0, //
    	.innerYRadius      = 30.0, //
    	.outerXRadius      = 0.0, //
    	.outerYRadius      = 0.0, //
		.xCenter           = 0.0,//-972.36849;
    	.yCenter           = 0.0,
    	.zCenter           = 0.0,//29055.1545;
    	.rotationAngle     = globRotationAngle
	});

	beampipe_dimensions[1].outerXRadius      = beampipe_dimensions[1].innerXRadius + pipeThickness, //
	beampipe_dimensions[1].outerYRadius      = beampipe_dimensions[1].innerYRadius + pipeThickness, //
	beampipe_dimensions[1].xCenter           = -1*(tmp_endpoint_x+((0.5*beampipe_dimensions[1].length)*TMath::Sin(-globRotationAngle))),//-972.36849;
	beampipe_dimensions[1].yCenter           = 0.0,
	beampipe_dimensions[1].zCenter           = tmp_endpoint_z+((0.5*beampipe_dimensions[1].length)*TMath::Cos(-globRotationAngle)),//29055.1545;

	tmp_endpoint_z = beampipe_dimensions[1].zCenter+((0.5*beampipe_dimensions[1].length)*TMath::Cos(-globRotationAngle));
  	tmp_endpoint_x = -1*beampipe_dimensions[1].xCenter+((0.5*beampipe_dimensions[1].length)*TMath::Sin(-globRotationAngle));   

  	//------------------------------------------------------------------------------------
  	//First roman pots scattering chamber
	//------------------------------------------------------------------------------------
  
  	//length[2]            = 200; // from VPC drawings
	//innerXRadius[2]      = 200.0; //
	//innerYRadius[2]      = 125.0; //
	//outerXRadius[2]      = innerXRadius[2] + pipeThickness; //
	//outerYRadius[2]      = innerYRadius[2] + pipeThickness; //

  	//xCenter[2]           = -1*(tmp_endpoint_x+((0.5*length[2])*TMath::Sin(-globRotationAngle)));//-972.36849;
  	//yCenter[2]           = 0.0;
  	//zCenter[2]           = tmp_endpoint_z+((0.5*length[2])*TMath::Cos(-globRotationAngle));//29055.1545;
  	//rotationAngle[2]     = globRotationAngle;

	beampipe_dimensions.push_back({

		.length            = 200, // from VPC drawings
    	.innerXRadius      = 200.0, //
    	.innerYRadius      = 125.0, //
    	.outerXRadius      = 0.0, //
    	.outerYRadius      = 0.0, //
		.xCenter           = 0.0,//-972.36849;
    	.yCenter           = 0.0,
    	.zCenter           = 0.0,//29055.1545;
    	.rotationAngle     = globRotationAngle
	});
	
	beampipe_dimensions[2].outerXRadius      = beampipe_dimensions[2].innerXRadius + pipeThickness, //
	beampipe_dimensions[2].outerYRadius      = beampipe_dimensions[2].innerYRadius + pipeThickness, //
	beampipe_dimensions[2].xCenter           = -1*(tmp_endpoint_x+((0.5*beampipe_dimensions[2].length)*TMath::Sin(-globRotationAngle))),//-972.36849;
	beampipe_dimensions[2].yCenter           = 0.0,
	beampipe_dimensions[2].zCenter           = tmp_endpoint_z+((0.5*beampipe_dimensions[2].length)*TMath::Cos(-globRotationAngle)),//29055.1545;

	tmp_endpoint_z = beampipe_dimensions[2].zCenter+((0.5*beampipe_dimensions[2].length)*TMath::Cos(-globRotationAngle));
  	tmp_endpoint_x = -1*beampipe_dimensions[2].xCenter+((0.5*beampipe_dimensions[2].length)*TMath::Sin(-globRotationAngle));  

  	//------------------------------------------------------------------------------------
  	//pipe between RP 1 and RP 2 stations
	//rectangular cross-section!!!!!
	//------------------------------------------------------------------------------------
  
  	//length[3]            = 1500.0; // from VPC drawings
	//innerXRadius[3]      = 150.0; //
	//innerYRadius[3]      = 30.0; //
	//outerXRadius[3]      = innerXRadius[3] + pipeThickness; //
	//outerYRadius[3]      = innerYRadius[3] + pipeThickness; //

  	//xCenter[3]           = -1*(tmp_endpoint_x+((0.5*length[3])*TMath::Sin(-globRotationAngle)));//-972.36849;
  	//yCenter[3]           = 0.0;
  	//zCenter[3]           = tmp_endpoint_z+((0.5*length[3])*TMath::Cos(-globRotationAngle));//29055.1545;
  	//rotationAngle[3]     = globRotationAngle;

	beampipe_dimensions.push_back({

		.length            = 1500.0, // from VPC drawings
    	.innerXRadius      = 150.0, //
    	.innerYRadius      = 30.0, //
    	.outerXRadius      = 0.0, //
    	.outerYRadius      = 0.0, //
		.xCenter           = 0.0,//-972.36849;
    	.yCenter           = 0.0,
    	.zCenter           = 0.0,//29055.1545;
    	.rotationAngle     = globRotationAngle
	});
	
	beampipe_dimensions[3].outerXRadius      = beampipe_dimensions[3].innerXRadius + pipeThickness, //
	beampipe_dimensions[3].outerYRadius      = beampipe_dimensions[3].innerYRadius + pipeThickness, //
	beampipe_dimensions[3].xCenter           = -1*(tmp_endpoint_x+((0.5*beampipe_dimensions[3].length)*TMath::Sin(-globRotationAngle))),//-972.36849;
	beampipe_dimensions[3].yCenter           = 0.0,
	beampipe_dimensions[3].zCenter           = tmp_endpoint_z+((0.5*beampipe_dimensions[3].length)*TMath::Cos(-globRotationAngle)),//29055.1545;

	tmp_endpoint_z = beampipe_dimensions[3].zCenter+((0.5*beampipe_dimensions[3].length)*TMath::Cos(-globRotationAngle));
  	tmp_endpoint_x = -1*beampipe_dimensions[3].xCenter+((0.5*beampipe_dimensions[3].length)*TMath::Sin(-globRotationAngle));  

  	//------------------------------------------------------------------------------------
  	//second roman pots scattering chamber
	//------------------------------------------------------------------------------------
  
  	//length[4]            = 200; // from VPC drawings
	//innerXRadius[4]      = 200.0; //
	//innerYRadius[4]      = 125.0; //
	//outerXRadius[4]      = innerXRadius[4] + pipeThickness; //
	//outerYRadius[4]      = innerYRadius[4] + pipeThickness; //

  	//xCenter[4]           = -1*(tmp_endpoint_x+((0.5*length[4])*TMath::Sin(-globRotationAngle)));//-972.36849;
  	//yCenter[4]           = 0.0;
  	//zCenter[4]           = tmp_endpoint_z+((0.5*length[4])*TMath::Cos(-globRotationAngle));//29055.1545;
  	//rotationAngle[4]     = globRotationAngle;
 
	beampipe_dimensions.push_back({

		.length            = 200, // from VPC drawings
    	.innerXRadius      = 200.0, //
    	.innerYRadius      = 125.0, //
    	.outerXRadius      = 0.0, //
    	.outerYRadius      = 0.0, //
		.xCenter           = 0.0,//-972.36849;
    	.yCenter           = 0.0,
    	.zCenter           = 0.0,//29055.1545;
    	.rotationAngle     = globRotationAngle
	});
	
	beampipe_dimensions[4].outerXRadius      = beampipe_dimensions[4].innerXRadius + pipeThickness, //
	beampipe_dimensions[4].outerYRadius      = beampipe_dimensions[4].innerYRadius + pipeThickness, //
	beampipe_dimensions[4].xCenter           = -1*(tmp_endpoint_x+((0.5*beampipe_dimensions[4].length)*TMath::Sin(-globRotationAngle))),//-972.36849;
	beampipe_dimensions[4].yCenter           = 0.0,
	beampipe_dimensions[4].zCenter           = tmp_endpoint_z+((0.5*beampipe_dimensions[4].length)*TMath::Cos(-globRotationAngle)),//29055.1545;

	tmp_endpoint_z = beampipe_dimensions[4].zCenter+((0.5*beampipe_dimensions[4].length)*TMath::Cos(-globRotationAngle));
    tmp_endpoint_x = -1*beampipe_dimensions[4].xCenter+((0.5*beampipe_dimensions[4].length)*TMath::Sin(-globRotationAngle));

	//------------------------------------------------------------------------------------
    // Pipe from second RP chamber to taper 
    //------------------------------------------------------------------------------------

	//length[5]            = 100.0; // from VPC drawings
    //innerXRadius[5]      = 150.0; //
    //innerYRadius[5]      = 30.0; //
    //outerXRadius[5]      = innerXRadius[5] + pipeThickness; //
    //outerYRadius[5]      = innerYRadius[5] + pipeThickness; //

    //xCenter[5]           = -1*(tmp_endpoint_x+((0.5*length[5])*TMath::Sin(-globRotationAngle)));//-972.36849;
    //yCenter[5]           = 0.0;
    //zCenter[5]           = tmp_endpoint_z+((0.5*length[5])*TMath::Cos(-globRotationAngle));//29055.1545;
    //rotationAngle[5]     = globRotationAngle;

	beampipe_dimensions.push_back({

		.length            = 100.0, // from VPC drawings
    	.innerXRadius      = 150.0, //
    	.innerYRadius      = 30.0, //
    	.outerXRadius      = 0.0, //
    	.outerYRadius      = 0.0, //
		.xCenter           = 0.0,//-972.36849;
    	.yCenter           = 0.0,
    	.zCenter           = 0.0,//29055.1545;
    	.rotationAngle     = globRotationAngle
	});
	
	beampipe_dimensions[5].outerXRadius      = beampipe_dimensions[5].innerXRadius + pipeThickness, //
	beampipe_dimensions[5].outerYRadius      = beampipe_dimensions[5].innerYRadius + pipeThickness, //
	beampipe_dimensions[5].xCenter           = -1*(tmp_endpoint_x+((0.5*beampipe_dimensions[5].length)*TMath::Sin(-globRotationAngle))),//-972.36849;
	beampipe_dimensions[5].yCenter           = 0.0,
	beampipe_dimensions[5].zCenter           = tmp_endpoint_z+((0.5*beampipe_dimensions[5].length)*TMath::Cos(-globRotationAngle)),//29055.1545;

    tmp_endpoint_z = beampipe_dimensions[5].zCenter+((0.5*beampipe_dimensions[5].length)*TMath::Cos(-globRotationAngle));
    tmp_endpoint_x = -1*beampipe_dimensions[5].xCenter+((0.5*beampipe_dimensions[5].length)*TMath::Sin(-globRotationAngle));

	//------------------------------------------------------------------------------------
    // taper near ZDC
    //------------------------------------------------------------------------------------

	//numbers here are not really correct for the full taper, just for the opening

    //length[6]            = 599.692; // from VPC drawings
    //innerXRadius[6]      = 150.0; //
    //innerYRadius[6]      = 30.0; //
    //outerXRadius[6]      = innerXRadius[6] + pipeThickness; //
    //outerYRadius[6]      = innerYRadius[6] + pipeThickness; //

    //xCenter[6]           = -1*(tmp_endpoint_x+((0.5*length[6])*TMath::Sin(-globRotationAngle)));//-972.36849;
    //yCenter[6]           = 0.0;
    //zCenter[6]           = tmp_endpoint_z+((0.5*length[6])*TMath::Cos(-globRotationAngle));//29055.1545;
    //rotationAngle[6]     = globRotationAngle;

	beampipe_dimensions.push_back({

		.length            = 599.692, // from VPC drawings
    	.innerXRadius      = 150.0, //
    	.innerYRadius      = 30.0, //
    	.outerXRadius      = 0.0, //
    	.outerYRadius      = 0.0, //
		.xCenter           = 0.0,//-972.36849;
    	.yCenter           = 0.0,
    	.zCenter           = 0.0,//29055.1545;
    	.rotationAngle     = globRotationAngle

	});
	
	beampipe_dimensions[6].outerXRadius      = beampipe_dimensions[6].innerXRadius + pipeThickness, //
	beampipe_dimensions[6].outerYRadius      = beampipe_dimensions[6].innerYRadius + pipeThickness, //
	beampipe_dimensions[6].xCenter           = -1*(tmp_endpoint_x+((0.5*beampipe_dimensions[6].length)*TMath::Sin(-globRotationAngle))),//-972.36849;
	beampipe_dimensions[6].yCenter           = 0.0,
	beampipe_dimensions[6].zCenter           = tmp_endpoint_z+((0.5*beampipe_dimensions[6].length)*TMath::Cos(-globRotationAngle)),//29055.1545;

    tmp_endpoint_z = beampipe_dimensions[6].zCenter+((0.5*beampipe_dimensions[6].length)*TMath::Cos(-globRotationAngle));
    tmp_endpoint_x = -1*beampipe_dimensions[6].xCenter+((0.5*beampipe_dimensions[6].length)*TMath::Sin(-globRotationAngle));

	//------------------------------------------------------------------------------------
    // pipe connecting taper to B2PF magnet, just past ZDC
    //------------------------------------------------------------------------------------

    //numbers here are not really correct for the full taper, just for the opening

    //length[7]            = 3000.0; // from VPC drawings
    //innerXRadius[7]      = 35.0; //
    //innerYRadius[7]      = 0.0; //NOT USED
    //outerXRadius[7]      = innerXRadius[7] + pipeThickness; //
    //outerYRadius[7]      = innerYRadius[7] + pipeThickness; //NOT USED

    //xCenter[7]           = -1*(tmp_endpoint_x+((0.5*length[7])*TMath::Sin(-globRotationAngle)));//-972.36849;
    //yCenter[7]           = 0.0;
    //zCenter[7]           = tmp_endpoint_z+((0.5*length[7])*TMath::Cos(-globRotationAngle));//29055.1545;
    //rotationAngle[7]     = globRotationAngle;

	beampipe_dimensions.push_back({

		.length            = 3000.0, // from VPC drawings
    	.innerXRadius      = 35.0, //
    	.innerYRadius      = 0.0, //NOT USED
    	.outerXRadius      = 0.0, //
    	.outerYRadius      = 0.0, //NOT USED
		.xCenter           = 0.0,//-972.36849;
    	.yCenter           = 0.0,
    	.zCenter           = 0.0,//29055.1545;
    	.rotationAngle     = globRotationAngle
	});
	
	beampipe_dimensions[7].outerXRadius      = beampipe_dimensions[7].innerXRadius + pipeThickness, //
	beampipe_dimensions[7].outerYRadius      = beampipe_dimensions[7].innerYRadius + pipeThickness, //NOT USED
	beampipe_dimensions[7].xCenter           = -1*(tmp_endpoint_x+((0.5*beampipe_dimensions[7].length)*TMath::Sin(-globRotationAngle))),//-972.36849;
	beampipe_dimensions[7].yCenter           = 0.0,
	beampipe_dimensions[7].zCenter           = tmp_endpoint_z+((0.5*beampipe_dimensions[7].length)*TMath::Cos(-globRotationAngle)),//29055.1545;

	std::cout << "All structures stored in vector container..." << std::endl;

	//------------------------------------------
	//begin building main volumes here
	//------------------------------------------

  	int pieceIdx = 0; //Larger, rectangular pipe transporting proton and neutral envelopes (neutral exit window and transfer to smaller proton line at the end)

	Box inner_subtract_outer(dd4hep::mm*beampipe_dimensions[5].outerXRadius, dd4hep::mm*beampipe_dimensions[5].outerYRadius, dd4hep::mm*beampipe_dimensions[5].length/2);
	Box inner_subtract_inner(dd4hep::mm*beampipe_dimensions[5].innerXRadius, dd4hep::mm*beampipe_dimensions[5].innerYRadius, dd4hep::mm*(beampipe_dimensions[5].length+5.0)/2);
	
	Box RP_subtract_outer(dd4hep::mm*beampipe_dimensions[1].outerXRadius, dd4hep::mm*beampipe_dimensions[1].outerYRadius, dd4hep::mm*(beampipe_dimensions[2].length+5.0)/2);
	Box RP_subtract_inner(dd4hep::mm*beampipe_dimensions[1].innerXRadius, dd4hep::mm*beampipe_dimensions[1].innerYRadius, dd4hep::mm*(beampipe_dimensions[2].length+10.0)/2);

	Box pipeAfterB1APF_outer(dd4hep::mm*beampipe_dimensions[pieceIdx].outerXRadius, dd4hep::mm*beampipe_dimensions[pieceIdx].outerYRadius, dd4hep::mm*beampipe_dimensions[pieceIdx].length/2);
	Box pipeAfterB1APF_inner(dd4hep::mm*beampipe_dimensions[pieceIdx].innerXRadius, dd4hep::mm*beampipe_dimensions[pieceIdx].innerYRadius, dd4hep::mm*(beampipe_dimensions[pieceIdx].length)/2);
	Box pipeAfterB1APF_firstEndCap(dd4hep::mm * beampipe_dimensions[pieceIdx].outerXRadius, dd4hep::mm * beampipe_dimensions[pieceIdx].outerYRadius, dd4hep::mm*5.0/2.0);
	Tube neutral_exit_window_cutout(0.0, dd4hep::mm*windowRadius, 1.0); // 1.0cm thick
	Box protonTransferWindow(dd4hep::mm * 155.0, dd4hep::mm*beampipe_dimensions[1].outerYRadius, dd4hep::mm*(5.0/2));

	SubtractionSolid tmpAfterB1APF(pipeAfterB1APF_outer, pipeAfterB1APF_inner); //This gets rid of the inner portion of the pipe, but leaves the endcaps
	SubtractionSolid tmpAfterFrontEndCap(tmpAfterB1APF, pipeAfterB1APF_firstEndCap, Position(0.0, 0.0, dd4hep::mm * (-beampipe_dimensions[pieceIdx].length)/2));
	SubtractionSolid pipeAfterProtonTransferWindow(tmpAfterFrontEndCap, protonTransferWindow, Position(dd4hep::mm * (-120.0), 0.0, dd4hep::mm * (beampipe_dimensions[pieceIdx].length)/2 ));

	SubtractionSolid pipeAfterB1APF(pipeAfterProtonTransferWindow, neutral_exit_window_cutout, Position(dd4hep::mm * 160.0 , 0.0, dd4hep::mm * 0.5*beampipe_dimensions[pieceIdx].length));

    Volume v_pipeAfterB1APF(Form("v_pipeAfterB1APF_%d", pieceIdx), pipeAfterB1APF, m_SS);
    sdet.setAttributes(det, v_pipeAfterB1APF, x_det.regionStr(), x_det.limitsStr(), vis_name);

    auto pv_pipe_0 = assembly.placeVolume(v_pipeAfterB1APF, Transform3D(RotationY(-0.025), Position(dd4hep::mm * beampipe_dimensions[pieceIdx].xCenter + 4.0, dd4hep::mm *  beampipe_dimensions[pieceIdx].yCenter, dd4hep::mm * beampipe_dimensions[pieceIdx].zCenter))); // 2353.06094)));
    pv_pipe_0.addPhysVolID("sector", 1);
    DetElement pipe_de_0(sdet, Form("sector_pipe_%d_de", pieceIdx), 1);
    pipe_de_0.setPlacement(pv_pipe_0);

	//--------------------------------------------------------------------
	
	pieceIdx = 1; //smaller rectangular pipe for the protons

	Box pipe_1_outer(dd4hep::mm*beampipe_dimensions[pieceIdx].outerXRadius, dd4hep::mm*beampipe_dimensions[pieceIdx].outerYRadius, dd4hep::mm*beampipe_dimensions[pieceIdx].length/2);
	Box pipe_1_inner(dd4hep::mm*beampipe_dimensions[pieceIdx].innerXRadius, dd4hep::mm*beampipe_dimensions[pieceIdx].innerYRadius, dd4hep::mm*(beampipe_dimensions[pieceIdx].length+5.0)/2);
	SubtractionSolid pipe_1(pipe_1_outer, pipe_1_inner);
	
    Volume v_pipe_1(Form("v_pipe_1_%d", pieceIdx), pipe_1, m_SS);
    sdet.setAttributes(det, v_pipe_1, x_det.regionStr(), x_det.limitsStr(), vis_name);

    auto pv_pipe_1 = assembly.placeVolume(v_pipe_1, Transform3D(RotationY(beampipe_dimensions[pieceIdx].rotationAngle), Position(dd4hep::mm * beampipe_dimensions[pieceIdx].xCenter, dd4hep::mm *  beampipe_dimensions[pieceIdx].yCenter, dd4hep::mm * beampipe_dimensions[pieceIdx].zCenter))); // 2353.06094)));
    pv_pipe_1.addPhysVolID("sector", 1);
    DetElement pipe_de_1(sdet, Form("sector_pipe_%d_de", pieceIdx), 1);
    pipe_de_1.setPlacement(pv_pipe_1);

	//------------------------------------------------------------------------

	pieceIdx = 3; //pipe between Roman pots station

    Box pipe_3_outer(dd4hep::mm*beampipe_dimensions[pieceIdx].outerXRadius, dd4hep::mm*beampipe_dimensions[pieceIdx].outerYRadius, dd4hep::mm*beampipe_dimensions[pieceIdx].length/2);
    Box pipe_3_inner(dd4hep::mm*beampipe_dimensions[pieceIdx].innerXRadius, dd4hep::mm*beampipe_dimensions[pieceIdx].innerYRadius, dd4hep::mm*(beampipe_dimensions[pieceIdx].length+5.0)/2);
    SubtractionSolid pipe_3(pipe_3_outer, pipe_3_inner);

    Volume v_pipe_3(Form("v_pipe_3_%d", pieceIdx), pipe_3, m_SS);
    sdet.setAttributes(det, v_pipe_3, x_det.regionStr(), x_det.limitsStr(), vis_name);

    auto pv_pipe_3 = assembly.placeVolume(v_pipe_3, Transform3D(RotationY(beampipe_dimensions[pieceIdx].rotationAngle), Position(dd4hep::mm * beampipe_dimensions[pieceIdx].xCenter, dd4hep::mm * beampipe_dimensions[pieceIdx].yCenter, dd4hep::mm * beampipe_dimensions[pieceIdx].zCenter))); // 2353.06094)));
    pv_pipe_3.addPhysVolID("sector", 1);
    DetElement pipe_de_3(sdet, Form("sector_pipe_%d_de", pieceIdx), 1);
    pipe_de_3.setPlacement(pv_pipe_3);

	//----------------------------------------------------------------

	pieceIdx = 2; //first roman pots scattering chamber
	
	Box box_rp_station_1_outer(dd4hep::mm*beampipe_dimensions[pieceIdx].outerXRadius, dd4hep::mm*beampipe_dimensions[pieceIdx].outerYRadius, dd4hep::mm*beampipe_dimensions[pieceIdx].length/2);
	Box box_rp_station_1_inner(dd4hep::mm*beampipe_dimensions[pieceIdx].innerXRadius, dd4hep::mm*beampipe_dimensions[pieceIdx].innerYRadius, dd4hep::mm*beampipe_dimensions[pieceIdx].length/2);
	SubtractionSolid tmp(box_rp_station_1_outer, box_rp_station_1_inner);
	//SubtractionSolid rpStation1(box_rp_station_1_outer, box_rp_station_1_inner);
	SubtractionSolid rpStation1(tmp, RP_subtract_outer);
	
    Volume v_rpStation1(Form("v_rpStation1_%d", pieceIdx), rpStation1, m_SS);
    sdet.setAttributes(det, v_rpStation1, x_det.regionStr(), x_det.limitsStr(), vis_name);

    auto pv_pipe_2 = assembly.placeVolume(v_rpStation1, Transform3D(RotationY(beampipe_dimensions[pieceIdx].rotationAngle), Position(dd4hep::mm * beampipe_dimensions[pieceIdx].xCenter, dd4hep::mm * beampipe_dimensions[pieceIdx].yCenter, dd4hep::mm * beampipe_dimensions[pieceIdx].zCenter))); // 2353.06094)));
    pv_pipe_2.addPhysVolID("sector", 1);
    DetElement pipe_de_2(sdet, Form("sector_pipe_%d_de", pieceIdx), 1);
    pipe_de_2.setPlacement(pv_pipe_2);

	//---------------------------------------------------------------------

	pieceIdx = 4; //Second roman pots scattering chamber

	Box box_rp_station_2_outer(dd4hep::mm*beampipe_dimensions[pieceIdx].outerXRadius, dd4hep::mm*beampipe_dimensions[pieceIdx].outerYRadius, dd4hep::mm*beampipe_dimensions[pieceIdx].length/2);
	Box box_rp_station_2_inner(dd4hep::mm*beampipe_dimensions[pieceIdx].innerXRadius, dd4hep::mm*beampipe_dimensions[pieceIdx].innerYRadius, dd4hep::mm*beampipe_dimensions[pieceIdx].length/2);
	SubtractionSolid tmp_2(box_rp_station_2_outer, box_rp_station_2_inner);
	//SubtractionSolid rpStation2(box_rp_station_2_outer, box_rp_station_2_inner);
	//SubtractionSolid tmp_22(tmp_2, pipe_3_outer);
	SubtractionSolid rpStation2(tmp_2, RP_subtract_outer);
	
    Volume v_rpStation2(Form("v_rpStation2_%d", pieceIdx), rpStation2, m_SS);
    sdet.setAttributes(det, v_rpStation2, x_det.regionStr(), x_det.limitsStr(), vis_name);

    auto pv_pipe_4 = assembly.placeVolume(v_rpStation2, Transform3D(RotationY(beampipe_dimensions[pieceIdx].rotationAngle), Position(dd4hep::mm * beampipe_dimensions[pieceIdx].xCenter, dd4hep::mm * beampipe_dimensions[pieceIdx].yCenter, dd4hep::mm * beampipe_dimensions[pieceIdx].zCenter))); // 2353.06094)));
    pv_pipe_4.addPhysVolID("sector", 1);
    DetElement pipe_de_4(sdet, Form("sector_pipe_%d_de", pieceIdx), 1);
    pipe_de_4.setPlacement(pv_pipe_4);
	
	//---------------------------------------------------------------

	pieceIdx = 5; //pipe between second RP station and the B2pf magnet

    Box pipe_5_outer(dd4hep::mm*beampipe_dimensions[pieceIdx].outerXRadius, dd4hep::mm*beampipe_dimensions[pieceIdx].outerYRadius, dd4hep::mm*beampipe_dimensions[pieceIdx].length/2);
    Box pipe_5_inner(dd4hep::mm*beampipe_dimensions[pieceIdx].innerXRadius, dd4hep::mm*beampipe_dimensions[pieceIdx].innerYRadius, dd4hep::mm*(beampipe_dimensions[pieceIdx].length+5.0)/2);
    SubtractionSolid pipe_5(pipe_5_outer, pipe_5_inner);

    Volume v_pipe_5(Form("v_pipe_5_%d", pieceIdx), pipe_5, m_SS);
    sdet.setAttributes(det, v_pipe_5, x_det.regionStr(), x_det.limitsStr(), vis_name); 
    
    auto pv_pipe_5 = assembly.placeVolume(v_pipe_5, Transform3D(RotationY(beampipe_dimensions[pieceIdx].rotationAngle), Position(dd4hep::mm * beampipe_dimensions[pieceIdx].xCenter, dd4hep::mm * beampipe_dimensions[pieceIdx].yCenter, dd4hep::mm * beampipe_dimensions[pieceIdx].zCenter))); // 2353.06094)));
    pv_pipe_5.addPhysVolID("sector", 1);
    DetElement pipe_de_5(sdet, Form("sector_pipe_%d_de", pieceIdx), 1);
    pipe_de_5.setPlacement(pv_pipe_5);

	//----------------------------------------------------------------

	pieceIdx = 6;

	Double_t trpVertices[16];
    Double_t trpVerticesInner[16];
    //(x0, y0, x1, y1, ... , x7, y7)
	//opening side - larger size
    trpVertices[0] =  -beampipe_dimensions[6].outerXRadius * dd4hep::mm;
    trpVertices[1] =  -beampipe_dimensions[6].outerYRadius * dd4hep::mm;

    trpVertices[2] = -beampipe_dimensions[6].outerXRadius * dd4hep::mm;
    trpVertices[3] = beampipe_dimensions[6].outerYRadius * dd4hep::mm;

    trpVertices[4] =  beampipe_dimensions[6].outerXRadius * dd4hep::mm;
    trpVertices[5] =  beampipe_dimensions[6].outerYRadius * dd4hep::mm;

    trpVertices[6] = beampipe_dimensions[6].outerXRadius * dd4hep::mm;
    trpVertices[7] = -beampipe_dimensions[6].outerYRadius * dd4hep::mm;
	
	//exiting side - smaller size
    
	trpVertices[8] = -beampipe_dimensions[6].outerYRadius * dd4hep::mm;
    trpVertices[9] = -beampipe_dimensions[6].outerYRadius * dd4hep::mm;

    trpVertices[10] = -beampipe_dimensions[6].outerYRadius * dd4hep::mm;
    trpVertices[11] = beampipe_dimensions[6].outerYRadius * dd4hep::mm;

    trpVertices[12] = beampipe_dimensions[6].outerYRadius * dd4hep::mm;
    trpVertices[13] = beampipe_dimensions[6].outerYRadius * dd4hep::mm;

    trpVertices[14] = beampipe_dimensions[6].outerYRadius * dd4hep::mm;
    trpVertices[15] = -beampipe_dimensions[6].outerYRadius * dd4hep::mm;

    for(int i = 0; i < 16; i++){

        if(trpVertices[i] > 0.0){trpVerticesInner[i] = trpVertices[i]-(pipeThickness * dd4hep::mm);}
        if(trpVertices[i] < 0.0){trpVerticesInner[i] = trpVertices[i]+(pipeThickness * dd4hep::mm);}

    }

	EightPointSolid taper_outer(dd4hep::mm * (0.5*beampipe_dimensions[pieceIdx].length), trpVertices);
	EightPointSolid taper_inner(dd4hep::mm * (0.5*beampipe_dimensions[pieceIdx].length), trpVerticesInner);
	
	Box taper_entrance(dd4hep::mm * beampipe_dimensions[pieceIdx].innerXRadius, dd4hep::mm * beampipe_dimensions[pieceIdx].innerYRadius, dd4hep::mm * (0.5*(pipeThickness + 5.0)));
	Box taper_exit(dd4hep::mm * beampipe_dimensions[pieceIdx].innerYRadius, dd4hep::mm * beampipe_dimensions[pieceIdx].innerYRadius, dd4hep::mm * (0.5*(pipeThickness + 5.0)));
	SubtractionSolid hollowTaper(taper_outer, taper_inner);
	SubtractionSolid taper_minus_entrance_cap(hollowTaper, taper_entrance, Position(0.0, 0.0, dd4hep::mm * (-0.5*beampipe_dimensions[pieceIdx].length)));
	SubtractionSolid finalTaper(taper_minus_entrance_cap, taper_exit, Position(0.0, 0.0, dd4hep::mm * (0.5*beampipe_dimensions[pieceIdx].length)));
	//SubtractionSolid finalTaper(taper_outer, taper_inner);

    Volume v_taper(Form("v_taper_%d", pieceIdx), finalTaper, m_SS);
    sdet.setAttributes(det, v_taper, x_det.regionStr(), x_det.limitsStr(), vis_name);

    auto pv_pipe_6 = assembly.placeVolume(v_taper, Transform3D(RotationY(beampipe_dimensions[pieceIdx].rotationAngle), Position(dd4hep::mm * beampipe_dimensions[pieceIdx].xCenter, dd4hep::mm * beampipe_dimensions[pieceIdx].yCenter, dd4hep::mm * beampipe_dimensions[pieceIdx].zCenter))); // 2353.06094)));
    pv_pipe_6.addPhysVolID("sector", 1);
    DetElement pipe_de_6(sdet, Form("sector_pipe_%d_de", pieceIdx), 1);
    pipe_de_6.setPlacement(pv_pipe_6);

	//---------------------------------------------------------------

	pieceIdx = 7; //pipe between taper and B2PF

    Tube pipe_after_taper(dd4hep::mm*beampipe_dimensions[pieceIdx].innerXRadius, dd4hep::mm*beampipe_dimensions[pieceIdx].outerXRadius, dd4hep::mm*beampipe_dimensions[pieceIdx].length/2);

    Volume v_pipe_7(Form("v_pipe_7_%d", pieceIdx), pipe_after_taper, m_SS);
    sdet.setAttributes(det, v_pipe_7, x_det.regionStr(), x_det.limitsStr(), vis_name);

    auto pv_pipe_7 = assembly.placeVolume(v_pipe_7, Transform3D(RotationY(beampipe_dimensions[pieceIdx].rotationAngle), Position(dd4hep::mm * beampipe_dimensions[pieceIdx].xCenter, dd4hep::mm * beampipe_dimensions[pieceIdx].yCenter, dd4hep::mm * beampipe_dimensions[pieceIdx].zCenter))); // 2353.06094)));
    pv_pipe_7.addPhysVolID("sector", 1);
    DetElement pipe_de_7(sdet, Form("sector_pipe_%d_de", pieceIdx), 1);
    pipe_de_7.setPlacement(pv_pipe_7);


	//--------------------------------------------------------------

	// This is the beam tube in the B0 magnet for the hadron beam
  	// doesn't use the slope information calculated before - it stands alone

	pieceIdx = 8;

  	Tube   b0_hadron_tube(b0_hadron_tube_inner_r, b0_hadron_tube_outer_r, b0_hadron_tube_length / 2.0);
  	Volume v_b0_hadron_tube("v_b0_hadron_tube", b0_hadron_tube, m_Be);
  	sdet.setAttributes(det, v_b0_hadron_tube, x_det.regionStr(), x_det.limitsStr(), vis_name);

    auto pv_pipe_8 = assembly.placeVolume(v_b0_hadron_tube, Transform3D(RotationY(-0.025), Position(-16.5, 0.0, 640.0))); // 2353.06094)));
    pv_pipe_8.addPhysVolID("sector", 1);
    DetElement pipe_de_8(sdet, Form("sector_pipe_%d_de", pieceIdx), 1);
    pipe_de_8.setPlacement(pv_pipe_6);
	
	//----------------------------------------------------------------

	pieceIdx = 9; //neutral exit window

	Box pipeAfterB1APF_LARGE(dd4hep::mm*(beampipe_dimensions[0].outerXRadius+5.0), dd4hep::mm*(beampipe_dimensions[0].outerYRadius+5.0), dd4hep::mm*(beampipe_dimensions[0].length+5.0)/2);
	Tube   neutral_exit_window(0.0, dd4hep::mm*windowRadius, 1.0); // 1.0cm thick
  	
	IntersectionSolid finalWindow(pipeAfterB1APF_outer, neutral_exit_window, Position(dd4hep::mm * 160.0 , 0.0, dd4hep::mm * 0.5*beampipe_dimensions[0].length));

	Volume v_neutral_exit_window("v_neutral_exit_window", finalWindow, m_Al);
  	sdet.setAttributes(det, v_neutral_exit_window, x_det.regionStr(), x_det.limitsStr(), "AnlRed");

    auto pv_pipe_9 = assembly.placeVolume(v_neutral_exit_window, Transform3D(RotationY(-0.025), Position( dd4hep::mm * beampipe_dimensions[0].xCenter + 4.0, 0.0, dd4hep::mm * beampipe_dimensions[0].zCenter)));
	pv_pipe_9.addPhysVolID("sector", 1);
    DetElement pipe_de_9(sdet, Form("sector_pipe_%d_de", pieceIdx), 1);
    pipe_de_9.setPlacement(pv_pipe_9);

	//-----------------------------------------------------------------
	// Build vacuum volumes here
	//-----------------------------------------------------------------


	pieceIdx = 0;

	Box vacuum_main_pipe(dd4hep::mm * beampipe_dimensions[pieceIdx].innerXRadius, dd4hep::mm*beampipe_dimensions[pieceIdx].innerYRadius, dd4hep::mm*(beampipe_dimensions[pieceIdx].length-20.0)/2);
	Box cutout_for_OMD_station(dd4hep::mm * beampipe_dimensions[pieceIdx].innerXRadius, dd4hep::mm * beampipe_dimensions[pieceIdx].innerYRadius, dd4hep::mm * 20.0);
	
	SubtractionSolid after_first_OMD_cutout(vacuum_main_pipe, cutout_for_OMD_station, Position(0.0, 0.0, dd4hep::mm * (22510.0 - beampipe_dimensions[pieceIdx].zCenter)));
	SubtractionSolid final_vacuum_main_pipe(after_first_OMD_cutout, cutout_for_OMD_station, Position(0.0, 0.0, dd4hep::mm * (24510.0 - beampipe_dimensions[pieceIdx].zCenter)));

	Volume v_vacuum_main_pipe("v_vacuum_main_pipe", final_vacuum_main_pipe, m_vac);
    sdet.setAttributes(det, v_vacuum_main_pipe, x_det.regionStr(), x_det.limitsStr(), "AnlBlue");

    auto pv_vacuum_0 = assembly.placeVolume(v_vacuum_main_pipe, Transform3D(RotationY(-0.025), Position( dd4hep::mm * beampipe_dimensions[pieceIdx].xCenter + 4.0, 0.0, dd4hep::mm * beampipe_dimensions[pieceIdx].zCenter)));
    pv_vacuum_0.addPhysVolID("sector", 1);
    DetElement vacuum_de_0(sdet, Form("sector_FF_vacuum_%d_de", pieceIdx), 1);
    vacuum_de_0.setPlacement(pv_vacuum_0);

	//------------------------------------------------------------------

	pieceIdx = 1;

    Box vacuum_pipe_before_RP_1(dd4hep::mm * beampipe_dimensions[pieceIdx].innerXRadius, dd4hep::mm*beampipe_dimensions[pieceIdx].innerYRadius, dd4hep::mm*(beampipe_dimensions[pieceIdx].length)/2);

    Volume v_pipe_before_RP_1("v_vacuum_pipe_before_RP_1", vacuum_pipe_before_RP_1, m_vac);
    sdet.setAttributes(det, v_pipe_before_RP_1, x_det.regionStr(), x_det.limitsStr(), "AnlBlue");

    auto pv_vacuum_1 = assembly.placeVolume(v_pipe_before_RP_1, Transform3D(RotationY(beampipe_dimensions[pieceIdx].rotationAngle), Position( dd4hep::mm * beampipe_dimensions[pieceIdx].xCenter, 0.0,  dd4hep::mm * beampipe_dimensions[pieceIdx].zCenter)));
    pv_vacuum_1.addPhysVolID("sector", 1);
    DetElement vacuum_de_1(sdet, Form("sector_FF_vacuum_%d_de", pieceIdx), 1);
    vacuum_de_1.setPlacement(pv_vacuum_1);

	//------------------------------------------------------------------

	pieceIdx = 2;

	//for roman pot station 1
	//to be added after Roman pot issue is solved so I can see the parts together

	//------------------------------------------------------------------

	pieceIdx = 3;

    Box vacuum_pipe_between_RP(dd4hep::mm * beampipe_dimensions[pieceIdx].innerXRadius, dd4hep::mm*beampipe_dimensions[pieceIdx].innerYRadius, dd4hep::mm*(beampipe_dimensions[pieceIdx].length)/2);

    Volume v_vacuum_pipe_between_RP("v_vacuum_pipe_between_RP", vacuum_pipe_between_RP, m_vac);
    sdet.setAttributes(det, v_vacuum_pipe_between_RP, x_det.regionStr(), x_det.limitsStr(), "AnlBlue");

    auto pv_vacuum_3 = assembly.placeVolume(v_vacuum_pipe_between_RP, Transform3D(RotationY(beampipe_dimensions[pieceIdx].rotationAngle), Position( dd4hep::mm * beampipe_dimensions[pieceIdx].xCenter, 0.0,  dd4hep::mm * beampipe_dimensions[pieceIdx].zCenter)));
    pv_vacuum_3.addPhysVolID("sector", 1);
    DetElement vacuum_de_3(sdet, Form("sector_FF_vacuum_%d_de", pieceIdx), 1);
    vacuum_de_3.setPlacement(pv_vacuum_3);

	//------------------------------------------------------------------

    pieceIdx = 4;

    //for roman pot station 2
	//to be added after Roman pot issue is solved so I can see the parts together

	//------------------------------------------------------------------

    pieceIdx = 5;

    Box vacuum_pipe_after_RP(dd4hep::mm * beampipe_dimensions[pieceIdx].innerXRadius, dd4hep::mm*beampipe_dimensions[pieceIdx].innerYRadius, dd4hep::mm*(beampipe_dimensions[pieceIdx].length)/2);

    Volume v_vacuum_pipe_after_RP("v_vacuum_pipe_after_RP", vacuum_pipe_after_RP, m_vac);
    sdet.setAttributes(det, v_vacuum_pipe_after_RP, x_det.regionStr(), x_det.limitsStr(), "AnlBlue");

    auto pv_vacuum_5 = assembly.placeVolume(v_vacuum_pipe_after_RP, Transform3D(RotationY(beampipe_dimensions[pieceIdx].rotationAngle), Position( dd4hep::mm * beampipe_dimensions[pieceIdx].xCenter, 0.0,  dd4hep::mm * beampipe_dimensions[pieceIdx].zCenter)));
    pv_vacuum_5.addPhysVolID("sector", 1);
    DetElement vacuum_de_5(sdet, Form("sector_FF_vacuum_%d_de", pieceIdx), 1);
    vacuum_de_5.setPlacement(pv_vacuum_5);


	//------------------------------------------------------------------

    pieceIdx = 6;

	EightPointSolid vacuum_taper(dd4hep::mm * (0.5*beampipe_dimensions[pieceIdx].length), trpVerticesInner);

    Volume v_vacuum_taper("v_vacuum_taper", vacuum_taper, m_vac);
    sdet.setAttributes(det, v_vacuum_taper, x_det.regionStr(), x_det.limitsStr(), "AnlBlue");

    auto pv_vacuum_6 = assembly.placeVolume(v_vacuum_taper, Transform3D(RotationY(beampipe_dimensions[pieceIdx].rotationAngle), Position( dd4hep::mm * beampipe_dimensions[pieceIdx].xCenter, 0.0,  dd4hep::mm * beampipe_dimensions[pieceIdx].zCenter)));
    pv_vacuum_6.addPhysVolID("sector", 1);
    DetElement vacuum_de_6(sdet, Form("sector_FF_vacuum_%d_de", pieceIdx), 1);
    vacuum_de_6.setPlacement(pv_vacuum_6);

	//-------------------------------------------------------------------

	pieceIdx = 7; //vacuum between taper and B2PF

    Tube vacuum_pipe_after_taper(0.0, dd4hep::mm*beampipe_dimensions[pieceIdx].innerXRadius, dd4hep::mm*beampipe_dimensions[pieceIdx].length/2);

    Volume v_vacuum_pipe_after_taper("v_vacuum_pipe_after_taper", vacuum_pipe_after_taper, m_vac);
    sdet.setAttributes(det, v_vacuum_pipe_after_taper, x_det.regionStr(), x_det.limitsStr(), "AnlBlue");

    auto pv_vacuum_7 = assembly.placeVolume(v_vacuum_pipe_after_taper, Transform3D(RotationY(beampipe_dimensions[pieceIdx].rotationAngle), Position(dd4hep::mm * beampipe_dimensions[pieceIdx].xCenter, dd4hep::mm * beampipe_dimensions[pieceIdx].yCenter, dd4hep::mm * beampipe_dimensions[pieceIdx].zCenter))); // 2353.06094)));
    pv_vacuum_7.addPhysVolID("sector", 1);
    DetElement vacuum_de_7(sdet, Form("sector_FF_vacuum_%d_de", pieceIdx), 1);
    vacuum_de_7.setPlacement(pv_vacuum_7);

	//-------------------------------------------------------------------

  	pv_assembly = det.pickMotherVolume(sdet).placeVolume(assembly); //, posAndRot);
  	pv_assembly.addPhysVolID("system", x_det.id()).addPhysVolID("barrel", 1);
  	sdet.setPlacement(pv_assembly);
  	assembly->GetShape()->ComputeBBox();
  	return sdet;
}

DECLARE_DETELEMENT(forwardBeamPipeBrazil, create_detector)
