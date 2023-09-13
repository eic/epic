// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2023 Alex Jentsch


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

double getRotatedZ(double z, double x, double angle);
double getRotatedX(double z, double x, double angle);

static Ref_t create_detector(Detector& det, xml_h e, SensitiveDetector /* sens */)
{

	using namespace ROOT::Math;
	xml_det_t x_det    = e;
	string    det_name = x_det.nameStr();
	DetElement sdet(det_name, x_det.id());
	Assembly   assembly(det_name + "_assembly");
	Material   m_Vac    = det.material("Vacuum");
	string     vis_name = x_det.visStr();

	PlacedVolume pv_assembly;

	//----------------------------------------------
	// Starting point is only the magnet centers, 
	// lengths, rotations, and radii --> 
	// everything else calculated internally to 
	// make it easier to update later.
	//----------------------------------------------

	const int numGaps = 6; //number of gaps between magnets (excluding the IP to B0pf transition -- special case)
	const int numMagnets = 7; //number of actual FF magnets between IP and FF detectors
	
	bool makeIP_B0pfVacuum = true; //This is for the special gap location between IP and b0pf
	const int numDetElements = numMagnets + numGaps + 1;

	double radii_magnet[numMagnets];
	double lengths_magnet[numMagnets];
	double rotation_magnet[numMagnets];
	double x_elem_magnet[numMagnets];
	double y_elem_magnet[numMagnets];
	double z_elem_magnet[numMagnets];

	//loop to fill arrays here
	
	int idx = 0;
	for(xml_coll_t c(x_det,_U(element)); c; ++c){
		
		xml_dim_t pos       = c.child(_U(placement));
	    double    pos_x     = pos.x();
	    double    pos_y     = pos.y();
	    double    pos_z     = pos.z();
	    double    pos_theta = pos.attr<double>(_U(theta));
	    xml_dim_t dims      = c.child(_U(dimensions)); //dimensions();
	    double    dim_z     = dims.z();
	    xml_dim_t apperture = c.child(_Unicode(apperture));
	    double    app_r     = apperture.r();
		
		
		radii_magnet[idx]     = app_r; // cm
		lengths_magnet[idx]   = dim_z; //cm
		rotation_magnet[idx]  = pos_theta;  // radians
		x_elem_magnet[idx]    = pos_x*dd4hep::cm;  
		y_elem_magnet[idx]    = pos_y*dd4hep::cm;    
		z_elem_magnet[idx]    = pos_z*dd4hep::cm;
		
		idx++;
	}

	//-------------------------------------------
	//override numbers for the first element --> 
	//doesn't use the actual B0pf geometry!!!
	//-------------------------------------------
	
	radii_magnet[0]     = 2.9;     // cm
	lengths_magnet[0]   = 120.0;   // cm
	rotation_magnet[0]  = -0.025;  // radians
	x_elem_magnet[0]    = -16.5;   // cm
	y_elem_magnet[0]    = 0.0;     // cm
	z_elem_magnet[0]    = 640.0;   // cm
	
	
	
	double x_beg[numMagnets];
	double z_beg[numMagnets];
	double x_end[numMagnets];
	double z_end[numMagnets];
  
	double angle_elem_gap[numGaps];
	double z_gap[numGaps];
	double x_gap[numGaps];
	double length_gap[numGaps];

	
	DetElement *detectorElement[numDetElements];
   
	//-------------------------------------------
	//calculate entrance/exit points of magnets
	//-------------------------------------------

	for(int i = 0; i < numMagnets; i++){
	
		// need to use the common coordinate system -->
		// use x = z, and y = x to make things easier
	
		z_beg[i] = getRotatedZ(-0.5*lengths_magnet[i], 0.0, rotation_magnet[i]) + z_elem_magnet[i];
		z_end[i] = getRotatedZ( 0.5*lengths_magnet[i], 0.0, rotation_magnet[i]) + z_elem_magnet[i];
		x_beg[i] = getRotatedX(-0.5*lengths_magnet[i], 0.0, rotation_magnet[i]) + x_elem_magnet[i];
		x_end[i] = getRotatedX( 0.5*lengths_magnet[i], 0.0, rotation_magnet[i]) + x_elem_magnet[i];
			
	}

	//------------------------------------------
	// this part is a bit ugly for now - 
	// it's to make the vacuum volume between the
	// end of the IP beam pipe and the beginning of 
	// beginning of the B0pf magnet
	// 
	// -->the volume will be calculated at the end
	//-------------------------------------------
	
	double endOfCentralBeamPipe_z = 445.580*dd4hep::cm; //extracted from central_beampipe.xml, line 64
	double diameterReduce = 11.0*dd4hep::cm; //size reduction to avoid overlap with electron pipe
	double vacuumDiameterEntrance = 25.792*dd4hep::cm - diameterReduce; //extracted from central_beampipe.xml, line 64
	double vacuumDiameterExit = 17.4*dd4hep::cm; //15mrad @ entrance to magnet to not overlap electron magnet
	double crossingAngle = -0.025; //radians
	double endOfCentralBeamPipe_x = endOfCentralBeamPipe_z*crossingAngle;

	

	//-----------------------------------------------
	//calculate gap region center, length, and angle
	//-----------------------------------------------

	for(int i = 1; i < numMagnets; i++){
	
		angle_elem_gap[i-1] = (x_beg[i] - x_end[i-1])/(z_beg[i] - z_end[i-1]);
		length_gap[i-1] = TMath::Sqrt(TMath::Power(z_beg[i] - z_end[i-1], 2) + TMath::Power(x_beg[i] - x_end[i-1], 2));
		z_gap[i-1] = z_end[i-1] + 0.5*length_gap[i-1]*TMath::Cos(angle_elem_gap[i-1]);
		x_gap[i-1] = x_end[i-1] + 0.5*length_gap[i-1]*TMath::Sin(angle_elem_gap[i-1]);
			
	}   
   
	Double_t inRadius[numGaps];
	Double_t outRadius[numGaps];
	Double_t nxLow[numGaps];
	Double_t nyLow[numGaps];
	Double_t nzLow[numGaps];
	Double_t nxHigh[numGaps];
	Double_t nyHigh[numGaps];
	Double_t nzHigh[numGaps];
	Double_t phi_initial[numGaps];
	Double_t phi_final[numGaps];   

	for(int gapIdx = 0; gapIdx < numGaps; gapIdx++){

		inRadius[gapIdx]    = 0.0;
		outRadius[gapIdx]   = radii_magnet[gapIdx+1];
		phi_initial[gapIdx] = 0.0;
		phi_final[gapIdx]   = TMath::TwoPi();
		nxLow[gapIdx]       = -(length_gap[gapIdx]/2.0)*TMath::Sin(rotation_magnet[gapIdx]-angle_elem_gap[gapIdx]);
		nyLow[gapIdx]       = 0.0;
		nzLow[gapIdx]       = -(length_gap[gapIdx]/2.0)*TMath::Cos(rotation_magnet[gapIdx]-angle_elem_gap[gapIdx]);
		nxHigh[gapIdx]      = (length_gap[gapIdx]/2.0)*TMath::Sin(rotation_magnet[gapIdx+1]-angle_elem_gap[gapIdx]);
		nyHigh[gapIdx]      = 0.0;
		nzHigh[gapIdx]      = (length_gap[gapIdx]/2.0)*TMath::Cos(rotation_magnet[gapIdx+1]-angle_elem_gap[gapIdx]);

	}

	//-----------------------
	// inside magnets
	//-----------------------

	for(int pieceIdx = 0; pieceIdx < numMagnets; pieceIdx++){

		std::string piece_name      = Form("MagnetVacuum%d", pieceIdx);

		Tube magnetPiece(piece_name, 0.0, radii_magnet[pieceIdx], lengths_magnet[pieceIdx]/2);
		Volume vpiece(piece_name, magnetPiece, m_Vac);
		sdet.setAttributes(det, vpiece, x_det.regionStr(), x_det.limitsStr(), vis_name);
		
	    auto pv = assembly.placeVolume(vpiece, 
								       Transform3D(RotationY(rotation_magnet[pieceIdx]), 
								       Position(x_elem_magnet[pieceIdx], y_elem_magnet[pieceIdx], z_elem_magnet[pieceIdx])));
	    pv.addPhysVolID("sector", 1);
	    detectorElement[pieceIdx] = new DetElement(sdet, Form("sector%d_de", pieceIdx), 1);
	    detectorElement[pieceIdx]->setPlacement(pv);
	
	}
  
    //--------------------------
    //between magnets
    //--------------------------
    
	for(int pieceIdx = numMagnets; pieceIdx < numGaps + numMagnets; pieceIdx++){
    	
		int correctIdx = pieceIdx-numMagnets;
	
		std::string piece_name  = Form("GapVacuum%d", correctIdx);
	
    	CutTube gapPiece(piece_name, inRadius[correctIdx], outRadius[correctIdx], length_gap[correctIdx]/2, phi_initial[correctIdx], phi_final[correctIdx], 
										nxLow[correctIdx], nyLow[correctIdx], nzLow[correctIdx], nxHigh[correctIdx], nyHigh[correctIdx], nzHigh[correctIdx]);
				
		Volume vpiece(piece_name, gapPiece, m_Vac);		
    	sdet.setAttributes(det, vpiece, x_det.regionStr(), x_det.limitsStr(), vis_name);
    	
		auto pv = assembly.placeVolume(vpiece, 
								       Transform3D(RotationY(angle_elem_gap[correctIdx]), 
								       Position(x_gap[correctIdx], 0.0, z_gap[correctIdx])));
	    pv.addPhysVolID("sector", 1);
	    detectorElement[pieceIdx] = new DetElement(sdet, Form("sector%d_de", pieceIdx), 1);
	    detectorElement[pieceIdx]->setPlacement(pv);
  
    }
  
    //--------------------------------------------------------------
    //make and place vacuum volume to connect IP beam pipe to B0pf
    //--------------------------------------------------------------
	
  	if(makeIP_B0pfVacuum){
  
  	  	double specialGapLength = TMath::Sqrt(TMath::Power(z_beg[0] - endOfCentralBeamPipe_z, 2) + TMath::Power(x_beg[0] - endOfCentralBeamPipe_x, 2)) - 0.1;
  	  	double specialGap_z = 0.5*specialGapLength*TMath::Cos(crossingAngle) + endOfCentralBeamPipe_z;
		double specialGap_x = 0.5*specialGapLength*TMath::Sin(crossingAngle) + endOfCentralBeamPipe_x;
  
		std::string piece_name  = Form("GapVacuum%d", numGaps + numMagnets);

		Cone specialGap(piece_name, specialGapLength/2, 0.0, vacuumDiameterEntrance/2, 0.0, vacuumDiameterExit/2 );
			
		Volume specialGap_v(piece_name, specialGap, m_Vac);		
		sdet.setAttributes(det, specialGap_v, x_det.regionStr(), x_det.limitsStr(), vis_name);
	
		auto pv = assembly.placeVolume(specialGap_v, Transform3D(RotationY(crossingAngle), Position(specialGap_x, 0.0, specialGap_z)));
    	pv.addPhysVolID("sector", 1);
    	detectorElement[numGaps + numMagnets] = new DetElement(sdet, Form("sector%d_de", numGaps + numMagnets), 1);
   		detectorElement[numGaps + numMagnets]->setPlacement(pv);
  
	}
  
  	//----------------------------------------------------

  	pv_assembly = det.pickMotherVolume(sdet).placeVolume(assembly);
  	pv_assembly.addPhysVolID("system", x_det.id()).addPhysVolID("barrel", 1);
  	sdet.setPlacement(pv_assembly);
  	assembly->GetShape()->ComputeBBox();
  	return sdet;
}

double getRotatedZ(double z, double x, double angle){
	
	return z*TMath::Cos(angle) - x*TMath::Sin(angle);
}

double getRotatedX(double z, double x, double angle){
	
	return z*TMath::Sin(angle) + x*TMath::Cos(angle);
}

DECLARE_DETELEMENT(magnetElementInnerVacuum, create_detector)
