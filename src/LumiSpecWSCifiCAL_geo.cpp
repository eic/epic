// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Aranya Giri

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
static tuple<Volume, Position> build_specScifiCAL_module(const Detector& description, const xml::Component& mod_x, SensitiveDetector& sens);

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

  auto [modVol, modSize] = build_specScifiCAL_module(description, x_mod, sens);
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
static tuple<Volume, Position> build_specScifiCAL_module( const Detector& description, const xml::Component& mod_x, SensitiveDetector& sens){

  //--------------------Module Setup---------------------------------------------------------------------
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

  //----------------------------Sci-Fi fibers -----------------------------------------------------------
  if (mod_x.hasChild(_Unicode(fiber))) {
    auto   fiber_tube  = mod_x.child(_Unicode(fiber));
    auto   fr       = fiber_tube.attr<double>(_Unicode(radius));
    auto   fsx      = fiber_tube.attr<double>(_Unicode(spacex));
    auto   fsy      = fiber_tube.attr<double>(_Unicode(spacey));
    auto   foff     = dd4hep::getAttrOrDefault<double>(fiber_tube, _Unicode(offset), 0.5 * mm);
    auto   fiberMat = description.material(fiber_tube.attr<std::string>(_Unicode(material)));
    Tube   fiberShape(0., fr, sz / 2.);
    Volume fiberVol("fiber_vol", fiberShape, fiberMat);
    fiberVol.setVisAttributes(description.visAttributes(fiber_tube.attr<std::string>(_Unicode(vis))));
    fiberVol.setSensitiveDetector(sens);

    // Fibers are placed in a honeycomb with the radius = sqrt(3)/2. * hexagon side length
    // the parameters space x and space y are used to add additional spaces between the hexagons
    double fside  = 2. / std::sqrt(3.) * fr;
    double fdistx = 2. * fside + fsx;
    double fdisty = 2. * fr + fsy;

    // maximum numbers of the fibers, help narrow the loop range
    int nx = int(sx / (2. * fr)) + 1;
    int ny = int(sy / (2. * fr)) + 1;

    // place the fibers
    double y0      = (foff + fside);
    int    nfibers = 0;
    for (int iy = 0; iy < ny; ++iy) {
      double y = y0 + fdisty * iy;
      // about to touch the boundary
      if ((sy - y) < y0) {
        break;
      }
      double x0 = (iy % 2) ? (foff + fside) : (foff + fside + fdistx / 2.);
      for (int ix = 0; ix < nx; ++ix) {
        double x = x0 + fdistx * ix;
        // about to touch the boundary
        if ((sx - x) < x0) {
          break;
        }
        auto fiberPV = modVol.placeVolume(fiberVol, nfibers++, Position{x - sx / 2., y - sy / 2., 0});
        fiberPV.addPhysVolID("fiber_x", ix + 1).addPhysVolID("fiber_y", iy + 1);
      }
    }
    
  } //-----------------if no fibers we make the module itself sensitive----------------------------------------
  else {
    modVol.setSensitiveDetector(sens);
  }

  return make_tuple(modVol, Position{sx, sy, sz} );
}

DECLARE_DETELEMENT(LumiSpecScifiCAL, create_detector) //(det_type, driver func)
