// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Aranya Giri

// Start Date - 10/31/2022
// Homogeneous PbWO4 (EM Calorimeter) Pair Spectrometer

#include "DD4hep/DetFactoryHelper.h"
#include <XML/Helper.h>
#include <algorithm>
#include <iostream>
#include <tuple>

using namespace std;
using namespace dd4hep;

// Definition of function to build the modules
static tuple<Volume, Position> build_specHomoCAL_module(const Detector& description, const xml::Component& mod_x, SensitiveDetector& sens);

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

  // Create Modules

  auto [modVol, modSize] = build_specHomoCAL_module(description, x_mod, sens);
  double detSizeXY = getAttrOrDefault( x_det, _Unicode(sizeXY), 20 );
  int nxy = int( detSizeXY / modSize.x() );
  double xypos0 = -nxy*modSize.x()/2.0 + modSize.x()/2.0;

  // Build detector components
  // loop over sectors
  for( xml_coll_t si(x_det, _Unicode(sector)); si; si++) { // sectors (top,bottom)

    xml_comp_t x_sector( si );
    int sector_id = x_sector.id();
    int mod_id = 0;

    xml_comp_t x_pos = x_sector.position();
    xml_comp_t x_rot = x_sector.rotation();

    for(int ix=0; ix< nxy; ix++){
      for(int iy=0; iy< nxy; iy++){

        double 	mod_pos_x 	= x_pos.x() + xypos0 + ix*modSize.x();
        double	mod_pos_y 	= x_pos.y() + xypos0 + iy*modSize.y();
        double 	mod_pos_z 	= x_pos.z() + 0.0*cm;

        PlacedVolume modPV = assembly.placeVolume(
            modVol, Transform3D( RotationZYX( x_rot.x(), x_rot.y(), x_rot.z()), Position( mod_pos_x, mod_pos_y, mod_pos_z ) ) );

        modPV.addPhysVolID( "sector", sector_id ).addPhysVolID( "module", mod_id );
        mod_id++;
      }
    }


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
static tuple<Volume, Position> build_specHomoCAL_module( const Detector& description, const xml::Component& mod_x, SensitiveDetector& sens){

  double sx = mod_x.attr<double>(_Unicode(sizex));
  double sy = mod_x.attr<double>(_Unicode(sizey));
  double sz = mod_x.attr<double>(_Unicode(sizez));
  double frame_size = mod_x.attr<double>(_Unicode(frameSize));

  Box    modShape( (sx/2.0 -frame_size) , (sy/2.0 -frame_size) , sz/2.0 );
  auto   modMat = description.material(mod_x.attr<std::string>(_Unicode(material)));
  Volume modVol("module_vol", modShape, modMat);

  if (mod_x.hasAttr(_Unicode(vis))) {
    modVol.setVisAttributes(description.visAttributes(mod_x.attr<std::string>(_Unicode(vis))));
  }

  modVol.setSensitiveDetector(sens);
  return make_tuple(modVol, Position{sx, sy, sz} );
}

DECLARE_DETELEMENT(LumiSpecHomoCAL, create_detector) //(det_type, driver func)
