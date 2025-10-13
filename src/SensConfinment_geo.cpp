// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 - 2025 Dhevan Gangadharan, Simon Gardner

//==========================================================================
//
// Sensitive IR6 confinment
//
//==========================================================================

#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/Printout.h"
#include "TMath.h"
#include <XML/Helper.h>
#include "DD4hep/Shapes.h"
#include "DD4hep/Objects.h"
#include "TGeoManager.h"
#include "TVector3.h"

using namespace std;
using namespace dd4hep;

static Ref_t create_detector(Detector& det, xml_h e, SensitiveDetector /* sens */) 
{
	using namespace ROOT::Math;
	xml_det_t x_det = e;
	string det_name = x_det.nameStr();
	DetElement sdet(det_name, x_det.id());
	Assembly assembly(det_name + "_assembly");

	const double w = 1*cm;
	vector<double> r_min = {2*m,	2*m,	2*m, 	2*m,	4*m,	4*m,	2*m,	2*m,	2*m,	2*m};
	vector<double> r_max = {2*m+w,	2*m+w,	4*m+w,	4*m+w,	4*m+w,	4*m+w,	4*m+w,	4*m+w,	2*m+w,	2*m+w};
	vector<double> z_vec = {-20*m,	-5*m-w,	-5*m-w,	-5*m,	-5*m,	5*m,	5*m,	5*m+w,	5*m+w, 	30*m};

	Solid sol = Polycone(0, 2*M_PI, r_min, r_max, z_vec);
	Volume vol(det_name + "_volume", sol, det.material("Vacuum"));
	assembly.placeVolume(vol, Position(0,0,0));

	// Final placement
	auto pv_assembly =
		det.pickMotherVolume(sdet).placeVolume(assembly, Position(0.0, 0.0, 0.0));
	sdet.setPlacement(pv_assembly);
	assembly->GetShape()->ComputeBBox();

	return sdet;
}

DECLARE_DETELEMENT(SensConfinment, create_detector)
