// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2023 Aranya Giri, Simon Gardner

/* Date : 12/12/2023
W Scifi (EM Calorimeter) Pair Spectrometer */

#include "DD4hep/DetFactoryHelper.h"
#include <XML/Helper.h>
#include <algorithm>
#include <iostream>
#include <tuple>
#include <TVector3.h>

using namespace std;
using namespace dd4hep;

// Definition of function to build the modules
static tuple<Volume, Position> build_specScFiCAL_module(const Detector& description, const xml::Component& mod_x, SensitiveDetector& sens);

// Driver Function
static Ref_t create_detector(Detector& description, xml_h e, SensitiveDetector sens)
{
  sens.setType("calorimeter");

  xml_det_t     x_det           =       e;
  xml_comp_t    x_mod           =       x_det.child( _Unicode(module) );
  string        det_name        =       x_det.nameStr();
  int           det_ID          =       x_det.id();

  // Create main detector element to be returned at the end
  DetElement    det( det_name, det_ID );

  // Mother volume
  Volume        motherVol = description.pickMotherVolume( det );

  // Detector assembly
  Assembly      assembly( det_name );
  assembly.setVisAttributes( description.invisible() );

  double detSizeXY      = getAttrOrDefault( x_mod, _Unicode(sizeXY), 180*mm);
  double detSizeZ       = getAttrOrDefault( x_mod, _Unicode(sizeZ), 180*mm);
  int nmod_perlayer     = getAttrOrDefault( x_mod, _Unicode(nmod_perlayer), 3);
  int nlayer            = getAttrOrDefault( x_mod, _Unicode(nlayer), 20);

  auto [modVol, modSize] = build_specScFiCAL_module(description, x_mod, sens);

  // No. of module/layer, # layer in CAL
  //int nmod_perlayer   = int( detSizeXY / modSize.x() );
  //int nlayer          = int( detSizeZ / modSize.y() );

  // Position of first module, layer from center of CAL
  double mod_pos0       = -(detSizeXY/2.0) + (modSize.x()/2.0);
  double layer_pos0     = -(detSizeZ/2.0) + (modSize.y()/2.0);

  // loop over sectors(top, bottom)
  for( xml_coll_t si(x_det, _Unicode(sector)); si; si++) {

          xml_comp_t x_sector( si );
          int sector_id = x_sector.id();
          int mod_id = 0;
          xml_comp_t x_pos = x_sector.position();
          xml_comp_t x_rot = x_sector.rotation();

          for(int ilay=0; ilay< nlayer; ilay++){
                 //Modules rotation A/Odd-Even
                  if((ilay%2)==0){ // 0 : || to x-axis
                          for(int imod=0; imod< nmod_perlayer; imod++){

                                  //Build // to z-axis, then rotate
                                  double        mod_pos_z       = ( x_pos.z() - 0.0*cm ) - layer_pos0 - ilay*modSize.y();
                                  double        mod_pos_y       = x_pos.y() + mod_pos0 + imod*modSize.x();
                                  double        mod_pos_x       = x_pos.x() + 0.0*cm;

                                  PlacedVolume modPV = assembly.placeVolume(
                                                  modVol, Transform3D( RotationZYX( x_rot.z(), -90.0*degree, -90.0*degree), Position( mod_pos_x, mod_pos_y, mod_pos_z ) ) );

                                  modPV.addPhysVolID( "sector", sector_id ).addPhysVolID( "module", mod_id);
                                  mod_id++;
                          }//imod-loop close
                  }//if-close
                  else{// 1 : || to y-axis
                          for(int imod=0; imod< nmod_perlayer; imod++){

                                  double        mod_pos_z       = ( x_pos.z() - 0.0*cm ) - layer_pos0 - ilay*modSize.y();
                                  double        mod_pos_x       = x_pos.x() + mod_pos0 + imod*modSize.x();
                                  double        mod_pos_y       = x_pos.y() + 0.0*mm;

                                  PlacedVolume modPV = assembly.placeVolume(
                                                  modVol, Transform3D( RotationZYX(x_rot.z(), x_rot.y(), -90.0*degree), Position( mod_pos_x, mod_pos_y, mod_pos_z ) ) );

                                  modPV.addPhysVolID( "sector", sector_id ).addPhysVolID( "module", mod_id);
                                  mod_id++;
                          }//imod-loop close
                  }//else-close

          }//ilay-loop close
  }// sectors

  // Place assembly into mother volume.  Assembly is centered at origin
  PlacedVolume detPV = motherVol.placeVolume( assembly, Position(0.0, 0.0, 0.0) );
  detPV.addPhysVolID("system", det_ID);
  // Connect to system ID
  det.setPlacement(detPV);

  return det;
} //Driver class close

static tuple<Volume, Position> build_specScFiCAL_module( const Detector& description, const xml::Component& mod_x, SensitiveDetector& sens){

        //--------------------Module Setup---------------------------------------------------------------------
        double sx = mod_x.attr<double>(_Unicode(sizex));
        double sy = mod_x.attr<double>(_Unicode(sizey));
        double sz = mod_x.attr<double>(_Unicode(sizez));

        Position modSize(sx, sy, sz);

        Box    modShape( modSize.x()/2.0, modSize.y()/2.0, modSize.z()/2.0 );
        auto   modMat = description.material(mod_x.attr<std::string>(_Unicode(material)));
        Volume modVol("module_vol", modShape, modMat);

        modVol.setVisAttributes(description.visAttributes(mod_x.attr<std::string>(_Unicode(vis))));
        //----------------------------Scintillating fibers -----------------------------------------------------------

        //detectors will be placed in module within square shaped sub-modules inside the modules.
        /*
           | + + + + + | x (->) even
           |  + + + +  | y (up) odd
           | + + + + + | even
           |  + + + +  | odd
           fsx, fsy = space between center of two fibers or submods along x, y in fig.
           foff = space between two modules
           */

        //fibers
        auto   fiber_tube  = mod_x.child(_Unicode(fiber));
        auto   fr       = fiber_tube.attr<double>(_Unicode(radius));
        auto   fsx      = fiber_tube.attr<double>(_Unicode(spacex));
        auto   fsy      = fiber_tube.attr<double>(_Unicode(spacey));
        auto   fiberMat = description.material(fiber_tube.attr<std::string>(_Unicode(material)));
        Tube   fiberShape(0., fr, modSize.z() / 2.0);
        Volume fiberVol("fiber_vol", fiberShape, fiberMat);
        fiberVol.setVisAttributes(description.visAttributes(fiber_tube.attr<std::string>(_Unicode(vis))));
        fiberVol.setSensitiveDetector(sens);

        //double submod_sizexy = 2.0*fr; // size of square = diameter of tubes.
        int num_submodX = int (modSize.x() / (2*fr + 2.0*fsx) );
        int num_submodY = int (modSize.y() / (2*fr + 2.0*fsy) );

        double submod_xpos0 = -modSize.x()/2.0 + fr + fsx;
        double submod_ypos0 = -modSize.y()/2.0 + fr + fsy;
        int nfibers = 0;

        //Fiber Holder
        auto   fiberholder_x    = mod_x.child(_Unicode(fiberholder));
        double fh_dz            = 0.6*mm; //thickness of fiber holder

        double  fh_outerbox_y   = 2.0*fr + 2.0*fsy;
        double  fh_outerbox_x   = 2.0*fr + 2.0*fsx;
        Box fh_outerbox(fh_outerbox_x/2.0, fh_outerbox_y/2.0, fh_dz/2.0);

        double  fh_innerbox_y   = 2.0*fr;
        double  fh_innerbox_x   = 2.0*fr;
        Box fh_innerbox(fh_innerbox_x/2.0, fh_innerbox_y/2.0, fh_dz/2.0);

        SubtractionSolid fiberholder_solid(fh_outerbox, fh_innerbox, Position(0.0, 0.0,0.0));
        auto fiberholderMat = description.material(fiberholder_x.attr<std::string>(_Unicode(material)));
        Volume fiberholderVol("fiberholder_vol",fiberholder_solid, fiberholderMat);
        fiberholderVol.setVisAttributes(description.visAttributes(fiberholder_x.attr<std::string>(_Unicode(vis))));

        int nfh = 0;

        //placement of fibers and fiberholder
        for(int iy = 0; iy< num_submodY; iy++){

                for(int ix=0; ix< num_submodX; ix++){

                        double  submod_pos_x    = submod_xpos0 + ix*(2.0*fr + 2.0*fsx); //mm
                        double  submod_pos_y    = submod_ypos0 + iy*(2.0*fr + 2.0*fsy); //mm
                        double  submod_pos_z    = 0*mm ; //mm

                        //placement of fiber
                        auto fiberPV = modVol.placeVolume(fiberVol, nfibers++, Position{submod_pos_x, submod_pos_y, submod_pos_z});
                        fiberPV.addPhysVolID("fiber_x", ix+1).addPhysVolID("fiber_y", iy+1);

                        //placement of fiber holder 6.6*cm apart c-to-c
                        int num_holders = 4; // which means 4 regions
                        double fh_pos_z0 = -1*(modSize.z()/2.0) + (fh_dz/2.0);

                        for(int iz=0; iz<num_holders;iz++){
                                double fh_pos_z = fh_pos_z0 + iz*( (modSize.z() - fh_dz)/(num_holders-1) );
                                modVol.placeVolume(fiberholderVol, nfh++, Position{submod_pos_x, submod_pos_y, fh_pos_z});
                        }//iz close

                }//ix close
        }//iy close

        return make_tuple(modVol, modSize);

}// build_specScifiCAL_module function close

DECLARE_DETELEMENT(EcalLumiSpecWScFi, create_detector)
