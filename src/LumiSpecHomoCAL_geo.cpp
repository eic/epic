/*
 Written by - Aranya Giri
 University of Houston
 Start Date - 10/31/2022
 Homogeneous PbWO4 (EM Calorimeter) Pair Spectrometer
 Note : Solid, Logical & Physical Volumes are the three component of detector,
 The following abbriviations will reflect the proper usage of these 3.
*/

// Essential Header Files
#include "DD4hep/DetFactoryHelper.h"
//#include "GeometryHelpers.h"
#include <XML/Helper.h>
#include <algorithm>
#include <iostream>
#include <math.h>
#include <tuple>

//Additional Header Files

using namespace std;
using namespace dd4hep;
//using Point = ROOT::Math::XYPoint;

//Definition of function to build the modules
tuple<Volume, Position> build_specHomoCAL_module(const Detector& lccdd, const xml::Component& mod_x, SensitiveDetector& sens);

//Driver Function
static Ref_t create_detector(Detector& lccdd, xml_h e, SensitiveDetector sens){

  sens.setType("calorimeter");

  xml_det_t 	x_det 		= 	e;
  xml_dim_t 	x_dim 		= 	x_det.dimensions();
  xml_dim_t	x_pos		= 	x_det.position();
  xml_dim_t	x_rot		= 	x_det.rotation();

  string		detName		= 	x_det.nameStr();
  int		detID		=	x_det.id();

  double         	env_dx          = 	x_dim.x();
  double         	env_dy          = 	x_dim.y();
  double         	env_dz          = 	x_dim.z();

  //Create Envolope
  Box envSolid(env_dx/2.0, env_dy/2.0, env_dz/2.0);

  Volume envLogic(detName+ "_envLogic", envSolid, lccdd.material("Vacuum") );
  envLogic.setVisAttributes( x_det.visStr() );

  //Create Module

  int nxy = 8;
  //int mod_id = 1;

  auto [modVol, modSize] = build_specHomoCAL_module(lccdd, x_det.child( _Unicode(module) ) , sens);
  double xypos0 = -(nxy*modSize.x())/2.0 + modSize.x()/2.0 ;

  for(int ix=0; ix< nxy; ix++){
    for(int iy=0; iy< nxy; iy++){

      double 	mod_pos_x 	= xypos0 + ix*modSize.x();
      double	mod_pos_y 	= xypos0 + iy*modSize.y();
      double 	mod_pos_z 	= 0.0*cm;	   
      PlacedVolume modPV 		= envLogic.placeVolume(modVol, Position(mod_pos_x, mod_pos_y, mod_pos_z) );
      modPV.addPhysVolID("row", ix+1).addPhysVolID("column", iy+1);
      // mod_id++;

    }
  }

  DetElement		det(detName, detID);
  Volume       		motherVol 	= lccdd.pickMotherVolume(det);

  Transform3D  		tr( RotationZYX( x_rot.x(), x_rot.y(), x_rot.z() ), Position( x_pos.x(), x_pos.y(), x_pos.z() ) ) ;
  PlacedVolume 		envPV		= motherVol.placeVolume(envLogic, tr);

  envPV.addPhysVolID("system", detID);
  det.setPlacement(envPV);

  return det;
} //Driver class close

//--------------------------------------------------------------------
//Function for building the module
tuple<Volume, Position> build_specHomoCAL_module( const Detector& lccdd, const xml::Component& mod_x, SensitiveDetector& sens){

  double sx = mod_x.attr<double>(_Unicode(sizex));
  double sy = mod_x.attr<double>(_Unicode(sizey));
  double sz = mod_x.attr<double>(_Unicode(sizez));

  double frame_size = 0.005*cm;

  Box    modShape( (sx/2.0 -frame_size) , (sy/2.0 -frame_size) , sz/2.0 );
  auto   modMat = lccdd.material(mod_x.attr<std::string>(_Unicode(material)));
  Volume modVol("module_vol", modShape, modMat); 

  if (mod_x.hasAttr(_Unicode(vis))) {
    modVol.setVisAttributes(lccdd.visAttributes(mod_x.attr<std::string>(_Unicode(vis))));
  }

  modVol.setSensitiveDetector(sens);
  return make_tuple(modVol, Position{sx, sy, sz} );
}

DECLARE_DETELEMENT(HomoCAL, create_detector) //(det_type, driver func)
