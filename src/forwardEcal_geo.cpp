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

  double blocksize = 10.0;                 // X,Y size of block  
  double blockgap  = 0.0254;               // Gap between blocks 
  double nsgap  = 0.27*2.54;               // North-South gap
  double rmin   = 0.0;                     // Dummy variable. Set to 0 since cutting out insert
  double rmax   = 81*2.54;                 // Max radius of endcap  
  double rmaxWithGap   = rmax + nsgap/2.0; // Max radius with NS gap included
  double zmax   = 362.0;                   // Face of iron plate fEcal mount to
  double length = 27.0;                    // Total length
  double zmin   = zmax - length;           // minimum z where detector starts
  double insert_dx[2] = {7.5, 20.05};      // Insert x width for north and south halves
  double insert_dy    = 30.05;             // Insert y height
  double insert_x=(insert_dx[0]-insert_dx[1])/2.0;  //Insert center x
  printf("zmin=%f\n",zmin);

  //xml_dim_t dim = detElem.dimensions();
  //double rmin   = dim.rmin(); // Dummy variable. Set to 0 since cutting out insert
  //double rmax   = dim.rmax(); // Max radius of endcap
  //double length = dim.z();    // Size along z-axis
  //xml_dim_t pos = detElem.position();
  //Getting insert dimensions
  //const xml::Component& insert_xml = detElem.child(_Unicode(insert));
  //xml_dim_t insert_dim             = insert_xml.dimensions();
  //xml_dim_t insert_local_pos       = insert_xml.position();

  PlacedVolume pv;
  Material air = desc.material("Air");

  // Defining envelope with full phi
  Tube envelope(rmin, rmaxWithGap, length / 2.0);

  // Removing insert shape from envelope
  Box insert((insert_dx[0] + insert_dx[1] + nsgap)/2.0, insert_dy/2.0, length / 2.);
  SubtractionSolid envelope_with_inserthole(envelope, insert, Position(insert_x,0.0,0.0));
  Volume envelopeVol(detName, envelope_with_inserthole, air);

  // Setting envelope attributes
  envelopeVol.setAttributes(desc, detElem.regionStr(), detElem.limitsStr(), detElem.visStr());
      
  double thickness=0.0;
  int slice_num  = 1;
  double slice_z = -length / 2.0; // Keeps track of slices' z locations in each layer
  // Looping over each layer's slices
  for (xml_coll_t sl(detElem, _U(slice)); sl; ++sl) {
    xml_comp_t x_slice     = sl;
    double slice_thickness = x_slice.thickness();
    thickness+=slice_thickness;
    printf("forwardEcal slice=%1d %8.4f %s \n",slice_num,slice_thickness,x_slice.materialStr().c_str());
    std::string slice_name = detName + "_" +_toString(slice_num, "slice%d");
    Material slice_mat     = desc.material(x_slice.materialStr());
    slice_z += slice_thickness / 2.; // Going to slice halfway point	
    Tube slice(rmin, rmaxWithGap, slice_thickness / 2.);
	
    // Removing insert shape from each slice
    Box slice_insert((insert_dx[0] + insert_dx[1] + nsgap)/2.0, insert_dy/2.0, slice_thickness/2.0);
    SubtractionSolid slice_with_inserthole(slice, slice_insert, Position(insert_x, 0.0, 0.0));
    Volume slice_vol(slice_name, slice_with_inserthole, air); //Still air
    slice_vol.setAttributes(desc, detElem.regionStr(), detElem.limitsStr(), detElem.visStr());
    pv = envelopeVol.placeVolume(slice_vol,Transform3D(RotationZYX(0, 0, 0), Position(0., 0., slice_z)));
	
    //Loop over north and south halves
    const double phi1[2]={-M_PI/2.0,M_PI/2.0};
    const double phi2[2]={M_PI/2.0,3.0*M_PI/2.0};
    const char* nsName[2]={"North","South"};
    const double pm[2]={1.0,-1.0}; //positive x for north, and negative for south
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
      pv = slice_vol.placeVolume(half_vol,Transform3D(RotationZYX(0, 0, 0), Position(pm[ns]*nsgap, 0.0, 0.0)));
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
	  double xrow=(map->xBlock(ns,r,0)+map->xBlock(ns,r,nColBlock-1))/2.0 - pm[ns]*nsgap;
	  double yrow=map->yBlock(ns,r);
	  pv = half_vol.placeVolume(row_vol, Transform3D(RotationZYX(0, 0, 0),Position(xrow,yrow,0)));
	  pv.addPhysVolID("blockrow", r);  
	  
	  //colmn of blocks
	  double xcol = -pm[ns]*(dxrow/2.0 - bsize/2.0);
	  for(int c=0; c<nColBlock; c++){	    
	    Box col(bsize,bsize,slice_thickness/2.0);
	    std::string col_name = row_name +_toString(c, "C%02d");
	    Volume col_vol(col_name, col, air);
	    col_vol.setAttributes(desc, x_slice.regionStr(), x_slice.limitsStr(), x_slice.visStr());
	    pv = row_vol.placeVolume(col_vol, Transform3D(RotationZYX(0, 0, 0),Position(xcol,0,0)));
	    pv.addPhysVolID("blockcol", c);
	    printf("r=%2d x=%8.2f y=%8.2f   c=%2d c=%8.2f\n",r,xrow,yrow,c,xcol); 
	    xcol += pm[ns]*bsize;

	    //place actual WSiFi block inside
	    Box block(blocksize,blocksize,slice_thickness/2.0);
            std::string block_name = col_name + "_WScFiBlock";
            Volume block_vol(block_name, block, slice_mat);
            pv = col_vol.placeVolume(block_vol, Transform3D(RotationZYX(0, 0, 0),Position(0,0,0)));
	    sens.setType("calorimeter");
	    block_vol.setSensitiveDetector(sens);
	  }
	}
      }
    }    
    slice_z += slice_thickness / 2.; // Going to end of slice
    ++slice_num;
  }
  printf("forwardEcal Total thickness=%f Slice end at %f\n",thickness,slice_z);

  DetElement det(detName, detID);
  Volume motherVol = desc.pickMotherVolume(det);

  // Placing endcap in world volume
  auto tr = Transform3D(Position(0.0, 0.0, zmax - length/2.0));

  pv = motherVol.placeVolume(envelopeVol, tr);
  pv.addPhysVolID("system", detID);
  det.setPlacement(pv);

  printf("forwardEcal_geo Done\n");
  return det;
}
DECLARE_DETELEMENT(epic_ForwardEcal, createDetector)
