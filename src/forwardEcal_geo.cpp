// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2025 Akio Ogawa

//==========================================================================
//  Implementation of forward calorimeter with 2025 design and 
//  the insert shape cut out
//--------------------------------------------------------------------------
//  Author: Akio Ogawa (BNL)
//==========================================================================

#include "DD4hep/DetFactoryHelper.h"
#include <XML/Helper.h>
#include <XML/Layering.h>
#include "forwardEcalMap.h"

using namespace dd4hep;

static Ref_t createDetector(Detector& desc, xml_h handle, SensitiveDetector sens) {
  forwardEcalMap *map = new forwardEcalMap();
  
  xml_det_t detElem   = handle;
  std::string detName = detElem.nameStr();
  int detID           = detElem.id();

  double blocksize = map->blockSize();         // X,Y size of block
  double blockgap  = map->spaceBetweenBlock(); // Gap between blocks 
  double nsgap  = mOffsetX[0]+mOffsetX[1];     // North-South gap
  double rmin   = 0.0;                     // Dummy variable. Set to 0 since cutting out insert
  double rmax   = 81*2.54;                 // Max radius of endcap  
  double rmaxWithGap   = rmax + nsgap/2.0; // Max radius with NS gap included
  double zmax   = map->backPlateZ();       // Face of back plate fEcal mount to = Hcal start z
  double length = 27.0;                    // Total length
  double zmin   = zmax - length;           // minimum z where detector starts
  double insert_dx[2] = {map->offsetXBeamPipe(0), map->offsetXBeamPipe(1)};
                                           // Insert x width for north and south halves
  double insert_dy    = 30.05;             // Insert y height
  double insert_dz    = 25.4;              // Insert (=Steel beam pipe protector) z depth
  double insert_thickness = 0.5*2.54;      // Insert (=Steel beam pipe protector) thickness
  double insert_x=(insert_dx[0]-insert_dx[1])/2.0;  //Insert center x

  //from compat file
  xml_dim_t dim = detElem.dimensions();
  xml_dim_t pos = detElem.position();
  if(dim.z() != length) printf("WARNING!!! forwardEcal_geo.cpp detect inconsistent Z len %f(compact) %f(map)\n",dim.z(),length);
  if(pos.z() != zmin)   printf("WARNING!!! forwardEcal_geo.cpp detect inconsistent Z pos %f(compact) %f(map)\n",pos.z(),zmin);
  printf("forwardEcal_geo : dz=%f %f zmin=%f %f\n",dim.z(),length,pos.z(),zmin);

  const double phi1[2]={-M_PI/2.0,M_PI/2.0};
  const double phi2[2]={M_PI/2.0,3.0*M_PI/2.0};
  const char* nsName[2]={"North","South"};
  const double pm[2]={1.0,-1.0}; //positive x for north, and negative for south
  
  PlacedVolume pv;
  Material air = desc.material("Air");
  Material steel = desc.material("StainlessSteel");

  // Defining envelope with full phi
  Tube envelope(rmin, rmaxWithGap, length / 2.0);

  // Removing insert shape from envelope
  Box insert((insert_dx[0] + insert_dx[1] - insert_thickness + nsgap)/2.0, (insert_dy - insert_thickness)/2.0, length / 2.);
  SubtractionSolid envelope_with_inserthole(envelope, insert, Position(insert_x,0.0,0.0));
  Volume envelopeVol(detName, envelope_with_inserthole, air);
  envelopeVol.setAttributes(desc, detElem.regionStr(), detElem.limitsStr(), detElem.visStr());

  double thickness=0.0;
  //int slice_num  = 1;
  double slice_z = -length / 2.0; // Keeps track of slices' z locations in each layer
  // Looping over each layer's slices
  for (xml_coll_t sl(detElem, _U(slice)); sl; ++sl) {
    xml_comp_t x_slice     = sl;
    double slice_thickness = x_slice.thickness();
    thickness+=slice_thickness;
    //printf("forwardEcal slice=%1d %8.4f %s \n",slice_num,slice_thickness,x_slice.materialStr().c_str());
    std::string slice_name = detName + "_" + x_slice.materialStr();
    Material slice_mat     = desc.material(x_slice.materialStr());
    slice_z += slice_thickness / 2.; // Going to slice halfway point	
    Tube slice(rmin, rmaxWithGap, slice_thickness / 2.);
	
    // Removing insert shape from each slice
    Box slice_insert((insert_dx[0] + insert_dx[1] + nsgap)/2.0, insert_dy/2.0, slice_thickness/2.0);
    SubtractionSolid slice_with_inserthole(slice, slice_insert, Position(insert_x, 0.0, 0.0));
    Volume slice_vol(slice_name, slice_with_inserthole, air); //Still air
    slice_vol.setAttributes(desc, x_slice.regionStr(), x_slice.limitsStr(), x_slice.visStr());
    pv = envelopeVol.placeVolume(slice_vol,Transform3D(RotationZYX(0, 0, 0), Position(0., 0., slice_z)));
	
    //Loop over north and south halves
    for (int ns=0; ns<2; ns++){
      Tube half(rmin, rmax, slice_thickness/2.0, phi1[ns], phi2[ns]);
      std::string half_name = slice_name + "_"  + nsName[ns];
      
      // Removing insert shape from each slice & halves   
      Box half_insert(insert_dx[ns]/2.0, insert_dy/2.0, slice_thickness/2.0);
      SubtractionSolid half_with_inserthole(half, half_insert,Position(pm[ns]*insert_dx[ns]/2.0, 0.0, 0.0));

      Material mat = slice_mat; 
      if(x_slice.isSensitive()) mat=air; //for calorimeter itself, still air to place blocks inside
      Volume half_vol(half_name, half_with_inserthole, mat);
      half_vol.setAttributes(desc, x_slice.regionStr(), x_slice.limitsStr(), x_slice.visStr());
      pv = slice_vol.placeVolume(half_vol,Transform3D(RotationZYX(0, 0, 0), Position(pm[ns]*nsgap/2.0, 0.0, 0.0)));
      pv.addPhysVolID("northsouth", ns);
      
      // For detector (sensitive) slice, placing detector blocks in col and row
      double bsize=blocksize+blockgap;
      if (x_slice.isSensitive()) {	
	//rows of blocks
	int nRowBlock = map->maxRowBlock(); //# of rows
	for(int r=0; r<nRowBlock; r++){
	  int nColBlock = map->nColBlock(ns,r);
	  double dxrow=bsize*nColBlock;
	  Box row(dxrow/2.0,bsize/2.0,slice_thickness/2.0);
	  std::string row_name = half_name +_toString(r, "_R%02d");
	  Volume row_vol(row_name, row, air);
	  row_vol.setAttributes(desc, x_slice.regionStr(), x_slice.limitsStr(), x_slice.visStr());
	  double xrow=(map->xBlock(ns,r,0) + map->xBlock(ns,r,nColBlock-1))/2.0 - pm[ns]*nsgap/2.0;
	  double yrow=map->yBlock(ns,r);
	  pv = half_vol.placeVolume(row_vol, Transform3D(RotationZYX(0, 0, 0),Position(xrow,yrow,0)));
	  pv.addPhysVolID("blockrow", r);  
	  
	  //column of blocks
	  double xcol = -pm[ns]*(dxrow/2.0 - bsize/2.0);
	  for(int c=0; c<nColBlock; c++){	    
	    Box col(bsize/2.0,bsize/2.0,slice_thickness/2.0);
	    std::string col_name = row_name +_toString(c, "C%02d");
	    Volume col_vol(col_name, col, air);
	    col_vol.setAttributes(desc, x_slice.regionStr(), x_slice.limitsStr(), x_slice.visStr());
	    pv = row_vol.placeVolume(col_vol, Transform3D(RotationZYX(0, 0, 0),Position(xcol,0,0)));
	    pv.addPhysVolID("blockcol", c);
	    //printf("r=%2d dx=%8.3f x=%8.3f y=%8.3f   c=%2d %8.3f  mapx=%8.3f\n",r,dxrow,xrow,yrow,c,xcol,map->xBlock(ns,r,c)); 
	    xcol += pm[ns]*bsize;

	    //place actual WSiFi block inside
	    Box block(blocksize/2.0,blocksize/2.0,slice_thickness/2.0);
            std::string block_name = col_name + "_WScFiBlock";
            Volume block_vol(block_name, block, slice_mat);
	    block_vol.setAttributes(desc, x_slice.regionStr(), x_slice.limitsStr(), x_slice.visStr());
            pv = col_vol.placeVolume(block_vol, Transform3D(RotationZYX(0, 0, 0),Position(0,0,0)));
	    sens.setType("calorimeter");
	    block_vol.setSensitiveDetector(sens);
	  }
	}
      }  //end if isSensitive
    }  //end loop over ns    
    slice_z += slice_thickness / 2.; // Going to end of slice
    //++slice_num;
  } //end loop over slice
  printf("forwardEcal Total thickness=%f Slice end at %f\n",thickness,slice_z);

  //Steel Beampipe Protector placed in envelope volume
  for (int ns=0; ns<2; ns++){
    Box bpp(insert_dx[ns]/2.0, insert_dy/2.0, insert_dz/2.0); 
    std::string bpp_name = detName + "_BeamPipeProtector_" + nsName[ns];
    Box bpp_hole((insert_dx[ns]-insert_thickness)/2.0, insert_dy/2.0-insert_thickness, insert_dz/2.0);
    SubtractionSolid bpp_with_hole(bpp, bpp_hole, Position(-pm[ns]*insert_thickness,0.0, 0.0));
    Volume bpp_vol(bpp_name,bpp_with_hole,steel);
    bpp_vol.setAttributes(desc, detElem.regionStr(), detElem.limitsStr(), detElem.visStr());
    pv = envelopeVol.placeVolume(bpp_vol,Transform3D(RotationZYX(0,0,0), Position(pm[ns]*(insert_dx[ns]+ nsgap)/2.0,0.0,(length-insert_dz)/2.0)));
  }
  
  DetElement det(detName, detID);
  Volume motherVol = desc.pickMotherVolume(det);

  // Placing endcap in world volume
  auto tr = Transform3D(Position(0.0, 0.0, zmin + length/2.0));

  pv = motherVol.placeVolume(envelopeVol, tr);
  pv.addPhysVolID("system", detID);
  det.setPlacement(pv);

  //printf("forwardEcal_geo Done\n");
  return det;
}
DECLARE_DETELEMENT(epic_ForwardEcal, createDetector)
