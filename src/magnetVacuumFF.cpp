// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2023 Alex Jentsch

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

double getRotatedZ(double z, double x, double angle);
double getRotatedX(double z, double x, double angle);

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

  //new code goes here!!!!!----------------------------
  
  //----------------------------------------------------------------------------
  // Starting point is only the magnet centers and OMD/RP detector locations, 
  // lengths, rotations, and radii -- everything else calculated internally to 
  // make it easier to update later.
  //----------------------------------------------------------------------------

  const int numGaps = 6;
  const int numMagnets = 7;

  TString elem_name[7] = {"b0pf", "b0apf", "q1apf", "q1bpf", "q2pf", "b1pf", "b1apf"}; //, "RP1", "RP2"};

  //-----------------------------------------------------------------------------
  // I can use the stupid unicode stuff to link this to the XML file with the 
  // magnet geometric information --> reduce chance for an error with update.
  //-----------------------------------------------------------------------------

	double radii_magnet[numMagnets];
	double lengths_magnet[numMagnets];
	double rotation_magnet[numMagnets];
	double x_elem_magnet[numMagnets];
	double y_elem_magnet[numMagnets];
	double z_elem_magnet[numMagnets];

	radii_magnet[0]     = 2.9;     // cm
	lengths_magnet[0]   = 120.0;   // 848.2683995; //290.0;    //cm
	rotation_magnet[0]  = -0.025;  // radians
	x_elem_magnet[0]    = -16.5;   // cm
	y_elem_magnet[0]    = 0.0;     // cm
	z_elem_magnet[0]    = 640.0;   // cm

	radii_magnet[1]     = 4.3;    // cm
	lengths_magnet[1]   = 60.0;   // cm
	rotation_magnet[1]  = -0.025; // radians
	x_elem_magnet[1]    = -21.0480535;
	y_elem_magnet[1]    = 0.0;
	z_elem_magnet[1]    = 819.8946015;

	radii_magnet[2]     = 5.6;     // cm
	lengths_magnet[2]   = 146.0;   // cm
	rotation_magnet[2]  = -0.0195; // radians
	x_elem_magnet[2]    = -25.4342857;
	y_elem_magnet[2]    = 0.0;
	z_elem_magnet[2]    = 962.8296939;

	radii_magnet[3]     = 7.8;    // cm
	lengths_magnet[3]   = 161.0;  // cm
	rotation_magnet[3]  = -0.015; // radians
	x_elem_magnet[3]    = -31.2840809;
	y_elem_magnet[3]    = 0.0;
	z_elem_magnet[3]    = 1156.243847;

	radii_magnet[4]     = 13.15;   // cm
	lengths_magnet[4]   = 380.0;   // cm
	rotation_magnet[4]  = -0.0148; // radians
	x_elem_magnet[4]    = -40.7362293;
	y_elem_magnet[4]    = 0.0;
	z_elem_magnet[4]    = 1466.604545;

	radii_magnet[5]     = 13.5;   // cm
	lengths_magnet[5]   = 300.0;  // cm
	rotation_magnet[5]  = -0.034; // radians
	x_elem_magnet[5]    = -50.3165042;
	y_elem_magnet[5]    = 0.0;
	z_elem_magnet[5]    = 1856.486896;

	radii_magnet[6]     = 16.8;   // cm
	lengths_magnet[6]   = 150.0;  // cm
	rotation_magnet[6]  = -0.025; // radians
	x_elem_magnet[6]    = -61.2903791;
	y_elem_magnet[6]    = 0.0;
	z_elem_magnet[6]    = 2131.298439;	

	//Off-Momentum Station 1
	//radii_magnet[7]     = 15.0;   // cm
	//lengths_magnet[7]   = 2.0;  // cm
	//rotation_magnet[7]  = -0.0454486856; // radians
	//x_elem_magnet[7]    = -100.0;
	//y_elem_magnet[7]    = 0.0;
	//z_elem_magnet[7]    = 2600.0;

	//Off-Momentum Station 2
	//radii_magnet[8]     = 15.0;   // cm
	//lengths_magnet[8]   = 2.0;  // cm
	//rotation_magnet[8]  = -0.02822; // radians
	//x_elem_magnet[8]    = -149.1239596;
	//y_elem_magnet[8]    = 0.0;
	//z_elem_magnet[8]    = 4074.293743;

	//Roman Pots Station 1
	//radii_magnet[7]     = 15.0;   // cm
	//lengths_magnet[7]   = 2.0;  // cm
	//rotation_magnet[7]  = -0.0454486856; // radians
	//x_elem_magnet[7]    = -92.3019;
	//y_elem_magnet[7]    = 0.0;
	//z_elem_magnet[7]    = 2797.0;

	//Roman Pots Station 2
	//radii_magnet[8]     = 15.0;   // cm
	//lengths_magnet[8]   = 2.0;  // cm
	//rotation_magnet[8]  = -0.0454486856; // radians
	//x_elem_magnet[8]    = -101.352;
	//y_elem_magnet[8]    = 0.0;
	//z_elem_magnet[8]    = 2996.0;


	//radii_magnet[7]     = 20.0;   // cm
	//lengths_magnet[7]   = 440.0;  // cm
	//rotation_magnet[7]  = -0.02822; // radians
	//x_elem_magnet[7]    = -149.1239596;
	//y_elem_magnet[7]    = 0.0;
	//z_elem_magnet[7]    = 4074.293743;	

	//double romanPotsStation1_z = 2797.0;
	//double romanPotsStatios1_x = -92.3019;

	//double romanPotsStation2_z = 2996.0;
	//double romanPotsStatios2_x = -1013.52;

	//double omdStation1_z = 2600.0;
	//double omdPotsStatios1_x = -100.0;

	//double omdPotsStation2_z = 2750.0;
	//double omdPotsStatios2_x = -1065.0;
 
	//double rpAndOMD_RotationAngle = -0.0454486856;

  
	double x_beg[numMagnets];
	double z_beg[numMagnets];
	double x_end[numMagnets];
	double z_end[numMagnets];
  
	double angle_elem_gap[numGaps];
	double z_gap[numGaps];
	double x_gap[numGaps];
	double length_gap[numGaps];

	//int numElements = 7;
   
	//CutTube *gapPiece[numGaps];
	//Tube *magnetPiece[numMagnets];
	//Volume  *vpiece[numMagnets + numGaps];
	DetElement *detectorElement[numMagnets + numGaps];
   
	//-------------------------------------------------------------
	//--- First step --> calculate entrance/exit points of magnets
	//-------------------------------------------------------------

	for(int i = 0; i < numMagnets; i++){
	
		// need to use the common coordinate system -->
		// use x = z, and y = x to make things easier
	
		z_beg[i] = getRotatedZ(-0.5*lengths_magnet[i], 0.0, rotation_magnet[i]) + z_elem_magnet[i];
		z_end[i] = getRotatedZ( 0.5*lengths_magnet[i], 0.0, rotation_magnet[i]) + z_elem_magnet[i];
		x_beg[i] = getRotatedX(-0.5*lengths_magnet[i], 0.0, rotation_magnet[i]) + x_elem_magnet[i];
		x_end[i] = getRotatedX( 0.5*lengths_magnet[i], 0.0, rotation_magnet[i]) + x_elem_magnet[i];
			
	}

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


  
	for(int pieceIdx = 0; pieceIdx < numMagnets; pieceIdx++){

		std::string piece_name      = Form("MagnetVacuum%d", pieceIdx);

		Tube magnetPiece(piece_name, 0.0, radii_magnet[pieceIdx], lengths_magnet[pieceIdx]/2);
	    //vpiece[pieceIdx]      = new Volume(Form("MagnetVacuum%d", pieceIdx), magnetPiece[pieceIdx], m_Vac);
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
		//vpiece[pieceIdx] = new Volume(Form("GapVacuum%d", correctIdx), gapPiece[correctIdx], m_Vac);
    	sdet.setAttributes(det, vpiece, x_det.regionStr(), x_det.limitsStr(), vis_name);
    	
		auto pv = assembly.placeVolume(vpiece, 
								       Transform3D(RotationY(angle_elem_gap[correctIdx]), 
								       Position(x_gap[correctIdx], 0.0, z_gap[correctIdx])));
	    pv.addPhysVolID("sector", 1);
	    detectorElement[pieceIdx] = new DetElement(sdet, Form("sector%d_de", pieceIdx), 1);
	    detectorElement[pieceIdx]->setPlacement(pv);
  
    }
  
  
  //----------------------------------------------------



  pv_assembly = det.pickMotherVolume(sdet).placeVolume(assembly); //, posAndRot);
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
