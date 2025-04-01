// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2025 Akio Ogawa

//==========================================================================
//  Implementation of forward calorimeter with 2025 design and 
//  the insert shape cut out, with ScFi
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
  
  double blocksize = map->blockSize();         // X,Y size of block
  double blockgap  = map->spaceBetweenBlock(); // Gap between blocks 
  double nsgap  = mOffsetX[0]+mOffsetX[1];     // North-South gap
  double rmin   = 0.0;                         // Dummy variable. Set to 0 since cutting out insert
  double rmax   = 81*2.54;                     // Max radius of endcap  
  double rmaxWithGap   = rmax + nsgap/2.0;     // Max radius with NS gap included
  double zmax   = map->backPlateZ();           // Face of back plate fEcal mount to = Hcal start z
  double length = 27.0;                        // Total length
  double zmin   = zmax - length;               // minimum z where detector starts
  double insert_dx[2] = {map->offsetXBeamPipe(0), map->offsetXBeamPipe(1)};
                                               // Insert x width for north and south halves
  double insert_dy    = 30.05;                 // Insert y height
  double insert_dz    = 27.0;                  // Insert (=Al beam pipe protector) z depth
  double insert_thickness = 0.25*2.54;         // Insert (=Al beam pipe protector) thickness
  double insert_x=(insert_dx[0]-insert_dx[1])/2.0;  //Insert center x
  int nx=26;             //number of fibers in a row
  int ny=30;	         //numbers of row of fibers	
  double rFiber =0.0235; //fiber radius (PMMA outside)
  double rScfi=0.02209;  //Scintillating fiber core radius

  const double phi1[2]={-M_PI/2.0,M_PI/2.0};
  const double phi2[2]={M_PI/2.0,3.0*M_PI/2.0};
  const char* nsName[2]={"N","S"};
  const double pm[2]={1.0,-1.0}; //positive x for north, and negative for south
  
  //from compact files
  xml_det_t detElem   = handle;
  std::string detName = detElem.nameStr();
  int detID           = detElem.id();

  xml_dim_t dim = detElem.dimensions();
  xml_dim_t pos = detElem.position();
  if(dim.z() != length) printf("WARNING!!! forwardEcal_geo.cpp detect inconsistent Z len %f(compact) %f(map)\n",dim.z(),length);
  if(pos.z() != zmin)   printf("WARNING!!! forwardEcal_geo.cpp detect inconsistent Z pos %f(compact) %f(map)\n",pos.z(),zmin);
  //printf("forwardEcal_geo : dz=%f %f zmin=%f %f\n",dim.z(),length,pos.z(),zmin);

  PlacedVolume pv;
  Material air = desc.material("Air");
  Material alumi = desc.material("Aluminum5083"); //actually using 6061... does it matter?
  Material Wpowder = desc.material("WPowderplusEpoxy");
  Material PMMA    = desc.material("Plexiglass");
  Material ScFi    = desc.material("Polystyrene");
  
  // Defining envelope with full phi,with slightly increased radius for NS gap
  Tube envelope(rmin, rmaxWithGap, length / 2.0);

  // Removing insert shape from envelope
  Box insert((insert_dx[0] + insert_dx[1] - insert_thickness + nsgap)/2.0, (insert_dy - insert_thickness)/2.0, length / 2.);
  SubtractionSolid envelope_with_inserthole(envelope, insert, Position(insert_x,0.0,0.0));
  Volume envelopeVol(detName, envelope_with_inserthole, air);
  envelopeVol.setAttributes(desc, detElem.regionStr(), detElem.limitsStr(), detElem.visStr());

  //double thickness=0.0;
  //int slice_num  = 1;
  double slice_z = -length / 2.0; // Keeps track of slices' z locations in each layer
  // Looping over each layer's slices
  for (xml_coll_t sl(detElem, _U(slice)); sl; ++sl) {
    xml_comp_t x_slice     = sl;
    double slice_thickness = x_slice.thickness();
    //thickness+=slice_thickness;
    //printf("forwardEcal slice=%1d %8.4f %s \n",slice_num,slice_thickness,x_slice.materialStr().c_str());
    std::string slice_name = "fEcal" + x_slice.nameStr();
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

	//Define WSiFi block (4x4 towers)       
	Box block(blocksize/2.0,blocksize/2.0,slice_thickness/2.0);
	Volume block_vol("fEcalBlock", block, air);
	block_vol.setAttributes(desc, x_slice.regionStr(), x_slice.limitsStr(), x_slice.visStr());

	//4 rows of towers
	Box trow(blocksize/2.0,blocksize/8.0,slice_thickness/2.0);
	Volume trow_vol("fEcalTowerRow", trow, air);
	trow_vol.setAttributes(desc, x_slice.regionStr(), x_slice.limitsStr(), x_slice.visStr());
	for(int tr=0; tr<4; tr++){
	  pv = block_vol.placeVolume(trow_vol, Transform3D(RotationZYX(0, 0, 0),Position(0.0,(tr-1.5)*blocksize/4,0)));
	  pv.addPhysVolID("towery", tr);
	}
	
	//4 towers in a row - finally a W powder volume, not air
	Box tower(blocksize/8.0,blocksize/8.0,slice_thickness/2.0);
	Volume tower_vol("fEcalTower", tower, Wpowder);
	tower_vol.setAttributes(desc, x_slice.regionStr(), x_slice.limitsStr(), x_slice.visStr());
	for(int tc=0; tc<4; tc++){
	  pv = trow_vol.placeVolume(tower_vol, Transform3D(RotationZYX(0, 0, 0),Position((tc-1.5)*blocksize/4,0,0)));
	  pv.addPhysVolID("towerx", tc);
	}

	//rows of fibers
	double fiberDistanceX=blocksize/4.0/(nx+0.5); //exrea 0.5 for even/odd rows shifted by 1/2 
	Box frow(blocksize/8.0-fiberDistanceX/2.0,blocksize/8.0/ny,slice_thickness/2.0);
	Volume frow_vol("fEcalFiberRow", frow, Wpowder);
	frow_vol.setAttributes(desc, x_slice.regionStr(), x_slice.limitsStr(), x_slice.visStr());
	for(int iy=0; iy<ny; iy++){
	  double xx = 0; 
	  if(iy%2 == 1) xx+=fiberDistanceX/2.0; 
	  pv = tower_vol.placeVolume(frow_vol, Transform3D(RotationZYX(0, 0, 0),Position(xx,(iy-ny/2.0+0.5)*blocksize/4.0/ny,0)));
	  //printf("iy=%2d dy=%8.4f fiberRx2=%8.4f xx=%8.4f\n",iy,blocksize/4.0/ny,rFiber*2,xx);
	  pv.addPhysVolID("fibery", iy);
	}
	
	//columns of fibers, with 1/2 fiber distance shifted each row
	Box fcol(fiberDistanceX/2.0,blocksize/8.0/ny,slice_thickness/2.0);
	Volume fcol_vol("fEcalFiberCol", fcol, Wpowder);
	fcol_vol.setAttributes(desc, x_slice.regionStr(), x_slice.limitsStr(), x_slice.visStr());
	for(int ix=0; ix<nx; ix++){
	  double xx=(ix - nx/2.0 + 0.5)*fiberDistanceX;
	  pv = frow_vol.placeVolume(fcol_vol, Transform3D(RotationZYX(0, 0, 0),Position(xx,0,0)));
	  //printf("ix=%2d dx=%8.4f xx=%8.4f x0=%8.4f x1=%8.4f\n",ix,fiberDistanceX,xx,xx-fiberDistanceX/2,xx+fiberDistanceX/2);
	  pv.addPhysVolID("fiberx", ix);
	}

	//a fiber (with coating material, not sensitive yet) 
	Tube fiber(0,rFiber,slice_thickness/2.0);
	Volume fiber_vol("fEcalFiber", fiber, PMMA);
	fiber_vol.setAttributes(desc, x_slice.regionStr(), x_slice.limitsStr(), x_slice.visStr());
	pv = fcol_vol.placeVolume(fiber_vol, Transform3D(RotationZYX(0, 0, 0),Position(0,0,0)));		    

	//scintillating fiber core - and finally a sensitive volume
	Tube scfi(0,rScfi,slice_thickness/2.0);
	Volume scfi_vol("fEcalScFi", scfi, ScFi);
	scfi_vol.setAttributes(desc, x_slice.regionStr(), x_slice.limitsStr(), x_slice.visStr());
	pv = fiber_vol.placeVolume(scfi_vol, Transform3D(RotationZYX(0, 0, 0),Position(0,0,0)));		    
	sens.setType("calorimeter");
	scfi_vol.setSensitiveDetector(sens);

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

	    //a block inside with air gap
            pv = col_vol.placeVolume(block_vol, Transform3D(RotationZYX(0, 0, 0),Position(0,0,0)));
	  }
	}
      }  //end if isSensitive
    }  //end loop over ns    
    slice_z += slice_thickness / 2.; // Going to end of slice
    //++slice_num;
  } //end loop over slice
  //printf("forwardEcal Total thickness=%f Slice end at %f\n",thickness,slice_z);

  //Al Beampipe Protector placed in envelope volume outside slices
  for (int ns=0; ns<2; ns++){
    Box bpp(insert_dx[ns]/2.0, insert_dy/2.0, insert_dz/2.0); 
    std::string bpp_name = detName + "_BeamPipeProtector_" + nsName[ns];
    Box bpp_hole((insert_dx[ns]-insert_thickness)/2.0, insert_dy/2.0-insert_thickness, insert_dz/2.0);
    SubtractionSolid bpp_with_hole(bpp, bpp_hole, Position(-pm[ns]*insert_thickness,0.0, 0.0));
    Volume bpp_vol(bpp_name,bpp_with_hole,alumi);
    bpp_vol.setAttributes(desc, detElem.regionStr(), detElem.limitsStr(), detElem.visStr());
    pv = envelopeVol.placeVolume(bpp_vol,Transform3D(RotationZYX(0,0,0), Position(pm[ns]*(insert_dx[ns]+ nsgap)/2.0,0.0,(length-insert_dz)/2.0)));
  }
  
  // Placing in the world volume
  DetElement det(detName, detID);
  Volume motherVol = desc.pickMotherVolume(det);
  auto tr = Transform3D(Position(0.0, 0.0, zmin + length/2.0));
  pv = motherVol.placeVolume(envelopeVol, tr);
  pv.addPhysVolID("system", detID);
  det.setPlacement(pv);

  //printf("forwardEcal_geo Done\n");
  return det;
}
DECLARE_DETELEMENT(epic_ForwardEcalScFi, createDetector)
