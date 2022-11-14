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
tuple<Volume, Position> build_specHomoCAL_module(const Detector& description, const xml::Component& mod_x, SensitiveDetector& sens);

//Driver Function
static Ref_t create_detector(Detector& description, xml_h e, SensitiveDetector sens)
{
  sens.setType("calorimeter");

  xml_det_t 	x_det 		= 	e;
  //xml_dim_t 	x_dim 		= 	x_det.dimensions();
  xml_dim_t	x_pos		= 	x_det.position();
  xml_dim_t	x_rot		= 	x_det.rotation();

  string        det_name	= 	x_det.nameStr();
  int		det_ID		=	x_det.id();

  // Create main detector element to be returned at the end
  DetElement    det( det_name, det_ID );

  // Mother volume
  Volume        motherVol = description.pickMotherVolume( det );

  // Detector assembly
  Assembly      assembly( det_name );
  assembly.setVisAttributes( description.invisible() );

  //Create Modules
  int nxy = 8;
  int mod_id = 0;

  auto [modVol, modSize] = build_specHomoCAL_module(description, x_det.child( _Unicode(module) ), sens);
  double xypos0 = -(nxy*modSize.x())/2.0 + modSize.x()/2.0 ;

  for(int ix=0; ix< nxy; ix++){
    for(int iy=0; iy< nxy; iy++){

      double 	mod_pos_x 	= xypos0 + ix*modSize.x();
      double	mod_pos_y 	= xypos0 + iy*modSize.y();
      double 	mod_pos_z 	= 0.0*cm;

      PlacedVolume modPV = assembly.placeVolume(
          modVol, Position( mod_pos_x, mod_pos_y, mod_pos_z ) );

      modPV.addPhysVolID("module", mod_id + 1);
      mod_id++;
    }
  }


  Transform3D tr( RotationZYX( x_rot.x(), x_rot.y(), x_rot.z() ), Position( x_pos.x(), x_pos.y(), x_pos.z() ) ) ;
  PlacedVolume detPV = motherVol.placeVolume( assembly, tr );

  detPV.addPhysVolID("system", det_ID);

  det.setPlacement(detPV);

  return det;
} //Driver class close

//--------------------------------------------------------------------
//Function for building the module
tuple<Volume, Position> build_specHomoCAL_module( const Detector& description, const xml::Component& mod_x, SensitiveDetector& sens){

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

DECLARE_DETELEMENT(HomoSpecCAL, create_detector) //(det_type, driver func)
