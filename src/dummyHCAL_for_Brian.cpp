// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2023 Alex Jentsch

/****

Dummy HCAL barrel, formed as a concentric set of cylinders of Steel absorber and plastic scintillator.
For testing only.

Author: Alex Jentsch
Date: June 14th, 2023

****/

#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/OpticalSurfaces.h"
#include "DD4hep/Printout.h"
#include "DDRec/DetectorData.h"
#include "DDRec/Surface.h"
#include <XML/Helper.h>
#include <XML/Layering.h>


using namespace std;
using namespace dd4hep;

static Ref_t createDetector(Detector& desc, xml_h handle, SensitiveDetector sens)
{
	
    /// basic barrel HCAL dimensions

    //double innerRadius = 1655.0 + 640.0; //millimeter
    //double outerRadius = innerRadius + 930.0; //millimeter
    //double barrelHCALLength = 1800.0 + 1350.0 + 46.0 + 3200; //millimeter - extracted from definitions.xml -- broken up on purpose to edit later
  	//double absorberThickness = 10.0; //millimeter
    //double scintillatorThickness = 2.5; //millimeter
  	//int numLayers = 74;

	xml::DetElement detElem = handle;
	std::string     detName = detElem.nameStr();
	int             detID   = detElem.id();
	DetElement      det(detName, detID);
	sens.setType("calorimeter");
	auto      dim    = detElem.dimensions();
	auto      inner_r = dim.rmin();
	auto      outer_r = dim.rmax();
	auto      length = dim.z();
	xml_dim_t pos    = detElem.position();
	xml_dim_t rot    = detElem.rotation();

	// envelope
	Tube    envShape(inner_r, outer_r, length * 0.5);
	Volume env(detName + "_envelope", envShape, desc.material("Air"));
	env.setVisAttributes(desc.visAttributes(detElem.visStr()));

	xml_comp_t mod_x  = detElem.child(_Unicode(module));
	//auto       nbox   = mod_x.attr<int>(_Unicode(nbox));
	//auto       boxgap = mod_x.attr<int>(_Unicode(gapspace));

	xml_comp_t x_lyr = mod_x.child(_Unicode(layer));
	auto       nlyr  = x_lyr.attr<int>(_Unicode(nlayer));

	map<int, string>    v_sl_name;
	map<string, Volume> slices;
	map<string, double> sl_thickness;

	int        nsl = 0;
	double     running_thickness = 0.0;
	int layerid = 0;
	xml_coll_t ci(x_lyr, _Unicode(slice));
	for(int iLayer = 0; iLayer < nlyr; iLayer++){
		for (ci.reset(); ci; ++ci) {
			xml_comp_t x_sl    = ci;
			Material   sl_mat  = desc.material(x_sl.materialStr());
			string     sl_name = x_sl.nameStr();
			double     sl_z    = x_sl.thickness();

			Tube   sl_Shape(inner_r + running_thickness, inner_r + running_thickness + sl_z , length / 2.);
			Volume sl_Vol("slice_vol", sl_Shape, sl_mat);
			sl_Vol.setVisAttributes(desc.visAttributes(x_sl.visStr()));
			if (x_sl.isSensitive()){
		 	   sl_Vol.setSensitiveDetector(sens);
			}

			nsl++;
			running_thickness    += sl_z;
			v_sl_name[nsl]        = sl_name;
			slices[sl_name]       = sl_Vol;
			sl_thickness[sl_name] = sl_z;
			
			Position     sl_pos(0, 0, 0);
			PlacedVolume pv = env.placeVolume(slices[sl_name], sl_pos);
			if (slices[sl_name].isSensitive()){
				pv.addPhysVolID(sl_name, layerid);
				layerid++;
			}
		}
	}
	

	// detector position and rotation
	Volume       motherVol = desc.pickMotherVolume(det);
	Transform3D  tr(RotationZYX(rot.z(), rot.y(), rot.x()), Position(pos.x(), pos.y(), pos.z()));
	PlacedVolume envPV = motherVol.placeVolume(env, tr);
	envPV.addPhysVolID("system", detID);
	det.setPlacement(envPV);
	return det;
}

DECLARE_DETELEMENT(dummyBarrelHCAL, createDetector)
