// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2023 Aranya Giri

/* Date : 04/31/2023
W Scifi (EM Calorimeter) Pair Spectrometer

Scintillating fiber calorimeter with tower shape blocks
reference: https://github.com/eic/epic/blob/main/src/ScFiCalorimeter_geo.cpp
Author: Chao Peng (ANL)*/

#include "DD4hep/DetFactoryHelper.h"
#include <XML/Helper.h>
#include <algorithm>
#include <iostream>
#include <tuple>

using namespace std;
using namespace dd4hep;

// Definition of function to build the modules
static tuple<Volume, Position> build_specScifiCAL_module(const Detector& description, const xml::Component& mod_x, SensitiveDetector& sens, int mod_id);

// Driver Function
static Ref_t create_detector(Detector& description, xml_h e, SensitiveDetector sens)
{
  sens.setType("calorimeter");

  xml_det_t 	x_det 		= 	e;
  xml_comp_t    x_mod           =       x_det.child( _Unicode(module) );
  string        det_name	= 	x_det.nameStr();
  int		det_ID		=	x_det.id();

  // Create main detector element to be returned at the end
  DetElement    det( det_name, det_ID );

  // Mother volume
  Volume        motherVol = description.pickMotherVolume( det );

  // Detector assembly
  Assembly      assembly( det_name );
  assembly.setVisAttributes( description.invisible() );
 
  // Build detector components
  // loop over sectors
  for( xml_coll_t si(x_det, _Unicode(sector)); si; si++) { // sectors (top,bottom)

    xml_comp_t x_sector( si );
    int sector_id = x_sector.id();
    int mod_id = 0;

    xml_comp_t x_pos = x_sector.position();
    xml_comp_t x_rot = x_sector.rotation();

    double detSizeXY = getAttrOrDefault( x_det, _Unicode(sizeXY), 20*cm);
    double modSizeXY = getAttrOrDefault( x_mod, _Unicode(sizex), 5*mm);
    int nxyz = int( detSizeXY / modSizeXY );
    double xyzpos0 = -nxyz*modSizeXY/2.0 + modSizeXY/2.0;

    for(int iz=0; iz< nxyz; iz++){
	    
	    // Create Modules
	    auto [modVol, modSize] = build_specScifiCAL_module(description, x_mod, sens, mod_id);
	    
	    if((iz%2)==0){ //90* rotation along y-axis
		    for(int iy=0; iy< nxyz; iy++){

			    double 	mod_pos_z 	= x_pos.z() + xyzpos0 + iz*modSize.x();
			    double	mod_pos_y 	= x_pos.y() + xyzpos0 + iy*modSize.y();
			    double 	mod_pos_x 	= x_pos.x() + 0.0*cm;

			    PlacedVolume modPV = assembly.placeVolume(
					    modVol, Transform3D( RotationZYX( x_rot.x(), 90.0*degree, x_rot.z()), Position( mod_pos_x, mod_pos_y, mod_pos_z ) ) );

			    modPV.addPhysVolID( "sector", sector_id ).addPhysVolID( "module", mod_id);
			    mod_id++;
		    }//iy-close
	    }//if-close 
	    else{// 90* rotation along x-axis
	    for(int ix=0; ix< nxyz; ix++){

			    double 	mod_pos_z 	= x_pos.z() + xyzpos0 + iz*modSize.x();
			    double	mod_pos_x 	= x_pos.x() + xyzpos0 + ix*modSize.x();
			    double 	mod_pos_y 	= x_pos.y() + 0.0*cm;

			    PlacedVolume modPV = assembly.placeVolume(
					    modVol, Transform3D( RotationZYX(0.0, x_rot.y(), 90.0*degree), Position( mod_pos_x, mod_pos_y, mod_pos_z ) ) );

			    modPV.addPhysVolID( "sector", sector_id ).addPhysVolID( "module", mod_id);
			    mod_id++;
		    }//ix-loop close
	    }//else-close
	    
    }//iz-close

  } // sectors

  // Place assembly into mother volume.  Assembly is centered at origin
  PlacedVolume detPV = motherVol.placeVolume( assembly, Position(0.0, 0.0, 0.0) );
  detPV.addPhysVolID("system", det_ID);

  // Connect to system ID
  det.setPlacement(detPV);

  return det;
} //Driver class close

//--------------------------------------------------------------------
//Function for building the module
static tuple<Volume, Position> build_specScifiCAL_module( const Detector& description, const xml::Component& mod_x, SensitiveDetector& sens, int mod_id){

	//Modules
	double sx = mod_x.attr<double>(_Unicode(sizex));
	double sy = mod_x.attr<double>(_Unicode(sizey));
	double sz = mod_x.attr<double>(_Unicode(sizez));
	//double frontplatesize = mod_x.attr<double>(_Unicode(frontplateSize));

	Box    modShape( sx/2.0 , sy/2.0 , sz/2.0 );
	auto   modMat = description.material(mod_x.attr<std::string>(_Unicode(material)));

	//Scifi fibers
	auto   fiber_box  = mod_x.child(_Unicode(fiber));
	auto   fsize      = fiber_box.attr<double>(_Unicode(size));
	auto   fiberMat = description.material(fiber_box.attr<std::string>(_Unicode(material)));
	Box    fiberShape(fsize/2.0, fsize/2.0, sz/2.0);

	//make the module hollow to insert Scifi fibers
	Volume modVol("module_vol", modShape, modMat);
	modVol.setVisAttributes(description.visAttributes(mod_x.attr<std::string>(_Unicode(vis))));
	
	Volume fiberVol("fiber_vol", fiberShape, fiberMat);
	PlacedVolume detfiberPV = modVol.placeVolume(fiberVol,Transform3D( RotationZYX(0.0, 0.0, 0.0), Position(0.0, 0.0,0.0) ) );
	fiberVol.setVisAttributes(description.visAttributes(fiber_box.attr<std::string>(_Unicode(vis))));
	fiberVol.setSensitiveDetector(sens);
	detfiberPV.addPhysVolID( "fiber", mod_id);

	return make_tuple(modVol, Position{sx, sy, sz} );
}

DECLARE_DETELEMENT(LumiSpecScifiCAL, create_detector) //(det_type, driver func)
