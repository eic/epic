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

	double pipeThickness     = 5.0 * dd4hep::mm;
  
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
	double crossingAngle     = -0.025; //relevant for the neutral cone

	double b1APFEndPoint_z = 22062.3828 * dd4hep::mm; //location of proton orbit at b1APF exit -- in mm
	double b1APFEndPoint_x = 654.3372 * dd4hep::mm; //location of proton orbit at b1APF exit -- in mm

	double tmp_endpoint_z = 0.0;
	double tmp_endpoint_x = 0.0;

	//forumula -> Z = b1APFEndPoint_z+((0.5*elementLengt)*Cos(globRotationAngle))
	//forumula -> X = b1APFEndPoint_z+((0.5*elementLengt)*Sin(globRotationAngle))

	//------------------------------------------------------------------------------------
	//primary pipe after B1APF, before neutral exit window + transition to smaller pipe
	//rectangular cross-section!!!!!
	//------------------------------------------------------------------------------------

	beampipe_dimensions.push_back({

		.length            = 7615.486 * dd4hep::mm, //from VPC drawings, in mm
    	.innerXRadius      = 275.0 * dd4hep::mm,    
    	.innerYRadius      = 175.0 * dd4hep::mm,    
    	.outerXRadius      = 0.0, 
    	.outerYRadius      = 0.0, 
		.xCenter           = 0.0, 
    	.yCenter           = 0.0,
    	.zCenter           = 0.0,
    	.rotationAngle     = globRotationAngle
	});
	
	beampipe_dimensions[0].outerXRadius      = beampipe_dimensions[0].innerXRadius + pipeThickness, //
    beampipe_dimensions[0].outerYRadius      = beampipe_dimensions[0].innerYRadius + pipeThickness, //
	beampipe_dimensions[0].xCenter           = -1*(b1APFEndPoint_x+((0.5*beampipe_dimensions[0].length)*TMath::Sin(-globRotationAngle))), //-745.20328;
    beampipe_dimensions[0].yCenter           = 0.0,
    beampipe_dimensions[0].zCenter           = (b1APFEndPoint_z+((0.5*beampipe_dimensions[0].length)*TMath::Cos(-globRotationAngle))),//24060.3176;
	 
	tmp_endpoint_z = beampipe_dimensions[0].zCenter+((0.5*beampipe_dimensions[0].length)*TMath::Cos(-globRotationAngle));
	tmp_endpoint_x = -1*beampipe_dimensions[0].xCenter+((0.5*beampipe_dimensions[0].length)*TMath::Sin(-globRotationAngle));

	double windowRadius = 110.0 * dd4hep::mm;


  	//------------------------------------------------------------------------------------
  	//first small pipe section, between primary vessel and RP station 1
	//rectangular cross-section!!!!!
	//------------------------------------------------------------------------------------
  
	beampipe_dimensions.push_back({ 

		.length            = 2780.273 * dd4hep::mm, // from VPC drawings
    	.innerXRadius      = 150.0 * dd4hep::mm, //
    	.innerYRadius      = 30.0 * dd4hep::mm, //
    	.outerXRadius      = 0.0, 
    	.outerYRadius      = 0.0, 
		.xCenter           = 0.0,
    	.yCenter           = 0.0,
    	.zCenter           = 0.0,
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
  
	beampipe_dimensions.push_back({

		.length            = 200 * dd4hep::mm, // from VPC drawings
    	.innerXRadius      = 200.0 * dd4hep::mm, //
    	.innerYRadius      = 125.0 * dd4hep::mm, //
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
  
  	
	beampipe_dimensions.push_back({

		.length            = 1500.0 * dd4hep::mm, // from VPC drawings
    	.innerXRadius      = 150.0 * dd4hep::mm, //
    	.innerYRadius      = 30.0 * dd4hep::mm, //
    	.outerXRadius      = 0.0, 
    	.outerYRadius      = 0.0, 
		.xCenter           = 0.0,
    	.yCenter           = 0.0,
    	.zCenter           = 0.0,
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
  
  	
	beampipe_dimensions.push_back({

		.length            = 200 * dd4hep::mm, // from VPC drawings
    	.innerXRadius      = 200.0 * dd4hep::mm, //
    	.innerYRadius      = 125.0 * dd4hep::mm, //
    	.outerXRadius      = 0.0, 
    	.outerYRadius      = 0.0, 
		.xCenter           = 0.0,
    	.yCenter           = 0.0,
    	.zCenter           = 0.0,
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


	beampipe_dimensions.push_back({

		.length            = 100.0 * dd4hep::mm, // from VPC drawings
    	.innerXRadius      = 150.0 * dd4hep::mm, //
    	.innerYRadius      = 30.0 * dd4hep::mm, //
    	.outerXRadius      = 0.0, 
    	.outerYRadius      = 0.0, 
		.xCenter           = 0.0,
    	.yCenter           = 0.0,
    	.zCenter           = 0.0,
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

	
	beampipe_dimensions.push_back({

		.length            = 599.692 * dd4hep::mm, // from VPC drawings
    	.innerXRadius      = 150.0 * dd4hep::mm, //
    	.innerYRadius      = 30.0 * dd4hep::mm, //
    	.outerXRadius      = 0.0, 
    	.outerYRadius      = 0.0, 
		.xCenter           = 0.0,
    	.yCenter           = 0.0,
    	.zCenter           = 0.0,
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

   
	beampipe_dimensions.push_back({

		.length            = 3000.0 * dd4hep::mm, // from VPC drawings
    	.innerXRadius      = 35.0 * dd4hep::mm, //
    	.innerYRadius      = 0.0, 
    	.outerXRadius      = 0.0, 
    	.outerYRadius      = 0.0, 
		.xCenter           = 0.0,
    	.yCenter           = 0.0,
    	.zCenter           = 0.0,
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

	//This solid is to properly remove the entrance/exit portions of the RP scattering chambers (see idx 2 and 4)
	Box RP_subtract_outer(beampipe_dimensions[1].outerXRadius, beampipe_dimensions[1].outerYRadius, (beampipe_dimensions[2].length+5.0)/2);
	
	//-------------------------------------------------------------------

  	int pieceIdx = 0; //Larger, rectangular pipe transporting proton and neutral envelopes (neutral exit window and transfer to smaller proton line at the end)
	
	Box pipeAfterB1APF_outer(beampipe_dimensions[pieceIdx].outerXRadius, beampipe_dimensions[pieceIdx].outerYRadius, beampipe_dimensions[pieceIdx].length/2);
	Box pipeAfterB1APF_inner(beampipe_dimensions[pieceIdx].innerXRadius, beampipe_dimensions[pieceIdx].innerYRadius, (beampipe_dimensions[pieceIdx].length)/2);
	Box pipeAfterB1APF_firstEndCap(beampipe_dimensions[pieceIdx].outerXRadius, beampipe_dimensions[pieceIdx].outerYRadius, 5.0/2.0);
	Tube neutral_exit_window_cutout(0.0, windowRadius, 1.0); // 1.0cm thick
	//FIXME: proton transfer window is done by hand right now - not a nicer way to do it until we get the CAD drawing
	Box protonTransferWindow(155.0 * dd4hep::mm, beampipe_dimensions[1].outerYRadius, (5.0/2));

	SubtractionSolid tmpAfterB1APF(pipeAfterB1APF_outer, pipeAfterB1APF_inner); //This gets rid of the inner portion of the pipe, but leaves the endcaps
	SubtractionSolid tmpAfterFrontEndCap(tmpAfterB1APF, pipeAfterB1APF_firstEndCap, Position(0.0, 0.0, (-beampipe_dimensions[pieceIdx].length)/2));
	//FIXME: proton transfer window is done by hand right now - not a nicer way to do it until we get the CAD drawing
	SubtractionSolid pipeAfterProtonTransferWindow(tmpAfterFrontEndCap, protonTransferWindow, Position((-120.0  * dd4hep::mm), 0.0, (beampipe_dimensions[pieceIdx].length)/2 ));

	SubtractionSolid pipeAfterB1APF(pipeAfterProtonTransferWindow, neutral_exit_window_cutout, Position(160.0  * dd4hep::mm, 0.0, 0.5*beampipe_dimensions[pieceIdx].length));

    Volume v_pipeAfterB1APF(Form("v_pipeAfterB1APF_%d", pieceIdx), pipeAfterB1APF, m_SS);
    sdet.setAttributes(det, v_pipeAfterB1APF, x_det.regionStr(), x_det.limitsStr(), vis_name);

    auto pv_pipe_0 = assembly.placeVolume(v_pipeAfterB1APF, Transform3D(RotationY(crossingAngle), Position(beampipe_dimensions[pieceIdx].xCenter + 4.0,  beampipe_dimensions[pieceIdx].yCenter, beampipe_dimensions[pieceIdx].zCenter))); // 2353.06094)));
    pv_pipe_0.addPhysVolID("sector", 1);
    DetElement pipe_de_0(sdet, Form("sector_pipe_%d_de", pieceIdx), 1);
    pipe_de_0.setPlacement(pv_pipe_0);

	//--------------------------------------------------------------------
	
	double lengthDelta = 0.0; //over-length value to remove end-pieces for hollow rectangular pipes
	
	// 1 -- small pipe connecting big pipe to RP station 1
	// 2 -- roman pots scattering chamber 1
	// 3 -- small pipe connecting RP1 and RP2
	// 4 -- roman pots scattering chamber 2
	// 5 -- small pipe connecting RP2 to ZDC taper
	
	lengthDelta = 5.0; //for small beam pipes to remove endcaps
	
	for(int idx = 1; idx < 6; idx++){ //loop for the easier pieces to simplify
		
		if(idx == 2 || idx == 4){ continue;}
		
		Box outer(beampipe_dimensions[idx].outerXRadius, beampipe_dimensions[idx].outerYRadius, beampipe_dimensions[idx].length/2);
		Box inner(beampipe_dimensions[idx].innerXRadius, beampipe_dimensions[idx].innerYRadius, (beampipe_dimensions[idx].length+lengthDelta)/2);
		
		SubtractionSolid hollow_pipe(outer, inner);
		
		Volume v_hollow_pipe(Form("v_pipe_%d", idx), hollow_pipe, m_SS);
	    sdet.setAttributes(det, v_hollow_pipe, x_det.regionStr(), x_det.limitsStr(), vis_name);

	    auto pv_final = assembly.placeVolume(v_hollow_pipe, Transform3D(RotationY(beampipe_dimensions[idx].rotationAngle), Position(beampipe_dimensions[idx].xCenter,  beampipe_dimensions[idx].yCenter, beampipe_dimensions[idx].zCenter))); 
	    pv_final.addPhysVolID("sector", 1);
	    DetElement final_de(sdet, Form("sector_pipe_%d_de", idx), 1);
	    final_de.setPlacement(pv_final);
		
	}
	
	lengthDelta = 0.0; //not needed for scattering chambers
	
	for(int idx = 1; idx < 6; idx++){ //loop for the easier pieces to simplify
		
		if(idx == 1 || idx == 3 || idx == 5){ continue;}
		
		Box outer(beampipe_dimensions[idx].outerXRadius, beampipe_dimensions[idx].outerYRadius, beampipe_dimensions[idx].length/2);
		Box inner(beampipe_dimensions[idx].innerXRadius, beampipe_dimensions[idx].innerYRadius, (beampipe_dimensions[idx].length+lengthDelta)/2);
		
		SubtractionSolid tmp(outer, inner);
		SubtractionSolid hollow_pipe(tmp, RP_subtract_outer);
		
		Volume v_hollow_pipe(Form("v_pipe_%d", idx), hollow_pipe, m_SS);
	    sdet.setAttributes(det, v_hollow_pipe, x_det.regionStr(), x_det.limitsStr(), vis_name);

	    auto pv_final = assembly.placeVolume(v_hollow_pipe, Transform3D(RotationY(beampipe_dimensions[idx].rotationAngle), Position(beampipe_dimensions[idx].xCenter,  beampipe_dimensions[idx].yCenter, beampipe_dimensions[idx].zCenter))); 
	    pv_final.addPhysVolID("sector", 1);
	    DetElement final_de(sdet, Form("sector_pipe_%d_de", idx), 1);
	    final_de.setPlacement(pv_final);
		
	}
	
	
	//----------------------------------------------------------------

	pieceIdx = 6;

	Double_t trpVertices[16];
    Double_t trpVerticesInner[16];
    //(x0, y0, x1, y1, ... , x7, y7)
	//opening side - larger size
    trpVertices[0] =  -beampipe_dimensions[6].outerXRadius;
    trpVertices[1] =  -beampipe_dimensions[6].outerYRadius;

    trpVertices[2] = -beampipe_dimensions[6].outerXRadius;
    trpVertices[3] = beampipe_dimensions[6].outerYRadius;

    trpVertices[4] =  beampipe_dimensions[6].outerXRadius;
    trpVertices[5] =  beampipe_dimensions[6].outerYRadius;

    trpVertices[6] = beampipe_dimensions[6].outerXRadius;
    trpVertices[7] = -beampipe_dimensions[6].outerYRadius;
	
	//exiting side - smaller size
    
	trpVertices[8] = -beampipe_dimensions[6].outerYRadius;
    trpVertices[9] = -beampipe_dimensions[6].outerYRadius;

    trpVertices[10] = -beampipe_dimensions[6].outerYRadius;
    trpVertices[11] = beampipe_dimensions[6].outerYRadius;

    trpVertices[12] = beampipe_dimensions[6].outerYRadius;
    trpVertices[13] = beampipe_dimensions[6].outerYRadius;

    trpVertices[14] = beampipe_dimensions[6].outerYRadius;
    trpVertices[15] = -beampipe_dimensions[6].outerYRadius;

    for(int i = 0; i < 16; i++){

        if(trpVertices[i] > 0.0){trpVerticesInner[i] = trpVertices[i]-(pipeThickness);}
        if(trpVertices[i] < 0.0){trpVerticesInner[i] = trpVertices[i]+(pipeThickness);}

    }

	EightPointSolid taper_outer((0.5*beampipe_dimensions[pieceIdx].length), trpVertices);
	EightPointSolid taper_inner((0.5*beampipe_dimensions[pieceIdx].length), trpVerticesInner);
	
	Box taper_entrance(beampipe_dimensions[pieceIdx].innerXRadius, beampipe_dimensions[pieceIdx].innerYRadius, (0.5*(pipeThickness + 5.0)));
	Box taper_exit(beampipe_dimensions[pieceIdx].innerYRadius, beampipe_dimensions[pieceIdx].innerYRadius, (0.5*(pipeThickness + 5.0)));
	SubtractionSolid hollowTaper(taper_outer, taper_inner);
	SubtractionSolid taper_minus_entrance_cap(hollowTaper, taper_entrance, Position(0.0, 0.0, (-0.5*beampipe_dimensions[pieceIdx].length)));
	SubtractionSolid finalTaper(taper_minus_entrance_cap, taper_exit, Position(0.0, 0.0, (0.5*beampipe_dimensions[pieceIdx].length)));
	//SubtractionSolid finalTaper(taper_outer, taper_inner);

    Volume v_taper(Form("v_taper_%d", pieceIdx), finalTaper, m_SS);
    sdet.setAttributes(det, v_taper, x_det.regionStr(), x_det.limitsStr(), vis_name);

    auto pv_pipe_6 = assembly.placeVolume(v_taper, Transform3D(RotationY(beampipe_dimensions[pieceIdx].rotationAngle), Position(beampipe_dimensions[pieceIdx].xCenter, beampipe_dimensions[pieceIdx].yCenter, beampipe_dimensions[pieceIdx].zCenter))); // 2353.06094)));
    pv_pipe_6.addPhysVolID("sector", 1);
    DetElement pipe_de_6(sdet, Form("sector_pipe_%d_de", pieceIdx), 1);
    pipe_de_6.setPlacement(pv_pipe_6);

	//---------------------------------------------------------------

	pieceIdx = 7; //pipe between taper and B2PF

    Tube pipe_after_taper(beampipe_dimensions[pieceIdx].innerXRadius, beampipe_dimensions[pieceIdx].outerXRadius, beampipe_dimensions[pieceIdx].length/2);

    Volume v_pipe_7(Form("v_pipe_7_%d", pieceIdx), pipe_after_taper, m_SS);
    sdet.setAttributes(det, v_pipe_7, x_det.regionStr(), x_det.limitsStr(), vis_name);

    auto pv_pipe_7 = assembly.placeVolume(v_pipe_7, Transform3D(RotationY(beampipe_dimensions[pieceIdx].rotationAngle), Position(beampipe_dimensions[pieceIdx].xCenter, beampipe_dimensions[pieceIdx].yCenter, beampipe_dimensions[pieceIdx].zCenter))); // 2353.06094)));
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

    auto pv_pipe_8 = assembly.placeVolume(v_b0_hadron_tube, Transform3D(RotationY(crossingAngle), Position(-16.5, 0.0, 640.0))); // 2353.06094)));
    pv_pipe_8.addPhysVolID("sector", 1);
    DetElement pipe_de_8(sdet, Form("sector_pipe_%d_de", pieceIdx), 1);
    pipe_de_8.setPlacement(pv_pipe_6);
	
	//----------------------------------------------------------------

	pieceIdx = 9; //neutral exit window

	Box pipeAfterB1APF_LARGE((beampipe_dimensions[0].outerXRadius+5.0), (beampipe_dimensions[0].outerYRadius+5.0), (beampipe_dimensions[0].length+5.0)/2);
	Tube   neutral_exit_window(0.0, windowRadius, 1.0); // 1.0cm thick
  	
	IntersectionSolid finalWindow(pipeAfterB1APF_outer, neutral_exit_window, Position(160.0 * dd4hep::mm, 0.0, 0.5*beampipe_dimensions[0].length));

	Volume v_neutral_exit_window("v_neutral_exit_window", finalWindow, m_Al);
  	sdet.setAttributes(det, v_neutral_exit_window, x_det.regionStr(), x_det.limitsStr(), "AnlRed");

    auto pv_pipe_9 = assembly.placeVolume(v_neutral_exit_window, Transform3D(RotationY(crossingAngle), Position( beampipe_dimensions[0].xCenter + 4.0, 0.0, beampipe_dimensions[0].zCenter)));
	pv_pipe_9.addPhysVolID("sector", 1);
    DetElement pipe_de_9(sdet, Form("sector_pipe_%d_de", pieceIdx), 1);
    pipe_de_9.setPlacement(pv_pipe_9);

	//-----------------------------------------------------------------
	// Build vacuum volumes here
	//-----------------------------------------------------------------


	pieceIdx = 0;

	Box vacuum_main_pipe(beampipe_dimensions[pieceIdx].innerXRadius, beampipe_dimensions[pieceIdx].innerYRadius, (beampipe_dimensions[pieceIdx].length-2.0)/2);
	Box cutout_for_OMD_station(beampipe_dimensions[pieceIdx].innerXRadius, beampipe_dimensions[pieceIdx].innerYRadius, 2.0);
	
	SubtractionSolid after_first_OMD_cutout(vacuum_main_pipe, cutout_for_OMD_station, Position(0.0, 0.0, (2251.0 - beampipe_dimensions[pieceIdx].zCenter)));
	SubtractionSolid final_vacuum_main_pipe(after_first_OMD_cutout, cutout_for_OMD_station, Position(0.0, 0.0, (2451.0 - beampipe_dimensions[pieceIdx].zCenter)));

	Volume v_vacuum_main_pipe("v_vacuum_main_pipe", final_vacuum_main_pipe, m_vac);
    sdet.setAttributes(det, v_vacuum_main_pipe, x_det.regionStr(), x_det.limitsStr(), "AnlBlue");

    auto pv_vacuum_0 = assembly.placeVolume(v_vacuum_main_pipe, Transform3D(RotationY(crossingAngle), Position( beampipe_dimensions[pieceIdx].xCenter + 4.0, 0.0, beampipe_dimensions[pieceIdx].zCenter)));
    pv_vacuum_0.addPhysVolID("sector", 1);
    DetElement vacuum_de_0(sdet, Form("sector_FF_vacuum_%d_de", pieceIdx), 1);
    vacuum_de_0.setPlacement(pv_vacuum_0);

	//------------------------------------------------------------------

	for(int idx = 1; idx < 6; idx++){ //loop for the easier pieces to simplify
		
		if(idx == 2 || idx == 4){ continue;} //FIXME: don't fill RP chambers with vacuum yet - still an issue with RP geometry
		
		Box inner_vacuum(beampipe_dimensions[idx].innerXRadius, beampipe_dimensions[idx].innerYRadius, (beampipe_dimensions[idx].length)/2);
				
		Volume v_inner_vacuum(Form("v_vacuum_%d", idx), inner_vacuum, m_vac);
	    sdet.setAttributes(det, v_inner_vacuum, x_det.regionStr(), x_det.limitsStr(), "AnlBlue");

	    auto pv_final = assembly.placeVolume(v_inner_vacuum, Transform3D(RotationY(beampipe_dimensions[idx].rotationAngle), Position(beampipe_dimensions[idx].xCenter,  beampipe_dimensions[idx].yCenter, beampipe_dimensions[idx].zCenter))); 
	    pv_final.addPhysVolID("sector", 1);
	    DetElement final_de(sdet, Form("sector_FF_vacuum_%d_de", idx), 1);
	    final_de.setPlacement(pv_final);
		
	}

	//------------------------------------------------------------------

    pieceIdx = 6;

	EightPointSolid vacuum_taper((0.5*beampipe_dimensions[pieceIdx].length), trpVerticesInner);

    Volume v_vacuum_taper("v_vacuum_taper", vacuum_taper, m_vac);
    sdet.setAttributes(det, v_vacuum_taper, x_det.regionStr(), x_det.limitsStr(), "AnlBlue");

    auto pv_vacuum_6 = assembly.placeVolume(v_vacuum_taper, Transform3D(RotationY(beampipe_dimensions[pieceIdx].rotationAngle), Position( beampipe_dimensions[pieceIdx].xCenter, 0.0,  beampipe_dimensions[pieceIdx].zCenter)));
    pv_vacuum_6.addPhysVolID("sector", 1);
    DetElement vacuum_de_6(sdet, Form("sector_FF_vacuum_%d_de", pieceIdx), 1);
    vacuum_de_6.setPlacement(pv_vacuum_6);

	//-------------------------------------------------------------------

	pieceIdx = 7; //vacuum between taper and B2PF

    Tube vacuum_pipe_after_taper(0.0, beampipe_dimensions[pieceIdx].innerXRadius, beampipe_dimensions[pieceIdx].length/2);

    Volume v_vacuum_pipe_after_taper("v_vacuum_pipe_after_taper", vacuum_pipe_after_taper, m_vac);
    sdet.setAttributes(det, v_vacuum_pipe_after_taper, x_det.regionStr(), x_det.limitsStr(), "AnlBlue");

    auto pv_vacuum_7 = assembly.placeVolume(v_vacuum_pipe_after_taper, Transform3D(RotationY(beampipe_dimensions[pieceIdx].rotationAngle), Position(beampipe_dimensions[pieceIdx].xCenter, beampipe_dimensions[pieceIdx].yCenter, beampipe_dimensions[pieceIdx].zCenter))); // 2353.06094)));
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
