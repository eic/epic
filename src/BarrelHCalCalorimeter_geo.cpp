//==========================================================================
//  AIDA Detector description implementation
//--------------------------------------------------------------------------
// Copyright (C) Organisation europeenne pour la Recherche nucleaire (CERN)
// All rights reserved.
//
// For the licensing terms see $DD4hepINSTALL/LICENSE.
// For the list of contributors see $DD4hepINSTALL/doc/CREDITS.
//
// Author     : M.Frank
//
//==========================================================================
//
// Specialized generic detector constructor
//
//==========================================================================
#include "DD4hep/DetFactoryHelper.h"
#include "XML/Layering.h"
#include "DD4hep/Printout.h"

using namespace std;
using namespace dd4hep;
using namespace dd4hep::detail;

static Ref_t create_detector(Detector& description, xml_h e, SensitiveDetector sens)
{

  //printout(WARNING, "BarrelHCalCalorimeter", "called create_detector ");

  xml_det_t     x_det             = e;
  int           det_id            = x_det.id();
  string        det_name          = x_det.nameStr();
  Material      air               = description.air();

  // get the solids section for this detector
  xml_comp_t    x_solids         = x_det.child("solids"); 

  DetElement       sdet(det_name, det_id);
  Volume           motherVol = description.pickMotherVolume(sdet);

  std::string   string_rmin      = dd4hep::xml::_toString(x_det.attr_value(x_det.getAttr("rmin"))); 
  std::string   string_rmax      = dd4hep::xml::_toString(x_det.attr_value(x_det.getAttr("rmax"))); 
  std::string   string_length    = dd4hep::xml::_toString(x_det.attr_value(x_det.getAttr("length"))); 

  double           rmin = (dd4hep::xml::_toDouble(string_rmin.c_str()))*dd4hep::mm;
  double           rmax = (dd4hep::xml::_toDouble(string_rmax.c_str()))*dd4hep::mm;
  double           length = (dd4hep::xml::_toDouble(string_length.c_str()))*dd4hep::mm;

  Tube             etube(rmin,rmax, length);
  Volume           envelope(det_name, etube, air);
  PlacedVolume     env_phv = motherVol.placeVolume(envelope);
  env_phv.addPhysVolID("system", det_id);
  env_phv.addPhysVolID("barrel", 0);
  sdet.setPlacement(env_phv);

  // Storage for sectors and tile assemblies
  Assembly      ChimneySector("ChimneySector"); 
  Assembly      Sector("Sector"); 
  Assembly      TileAssembly24("TileAssembly24"); 
  Assembly      TileAssembly24Chimney("TileAssembly24Chimney"); 

  xml_comp_t    det_define           = x_det.child("define");

  // Pick up the constants

  double tilePlaneRotate = 0.0; 
  double sectorRotate = 0.0; 

  for(xml_coll_t i(det_define, _Unicode(constant)); i; ++i){
    xml_comp_t  x_const = i; 
    std::string   const_name      = dd4hep::xml::_toString(x_const.attr_value(x_const.getAttr("name"))); 
    std::string   const_value     = dd4hep::xml::_toString(x_const.attr_value(x_const.getAttr("value"))); 

    if(const_name == "tilePlaneRotate")
      tilePlaneRotate = dd4hep::xml::_toDouble(const_value.c_str());
    else if(const_name == "sectorRotate")
      sectorRotate = dd4hep::xml::_toDouble(const_value.c_str());

  }

  // Loop over the defines section and pick up the tile offsets

  std::vector<double> xposOuter; 
  std::vector<double> yposOuter; 

  std::vector<double> xposTileS; 
  std::vector<double> yposTileS; 
  std::vector<double> zposTileS; 

  std::vector<double> xposTileN; 
  std::vector<double> yposTileN; 
  std::vector<double> zposTileN; 

  std::vector<double> xposChimneyTileS; 
  std::vector<double> yposChimneyTileS; 
  std::vector<double> zposChimneyTileS; 


  for(xml_coll_t i(det_define, _Unicode(matrix)); i; ++i){
    xml_comp_t  x_mtrx = i; 
    std::string   mtrx_name       = dd4hep::xml::_toString(x_mtrx.attr_value(x_mtrx.getAttr("name"))); 
    std::string   mtrx_values     = dd4hep::xml::_toString(x_mtrx.attr_value(x_mtrx.getAttr("values"))); 
    
    std::vector<double> *aptr = NULL; 

    if(mtrx_name == "xposOuter") 
      aptr = &xposOuter; 
    else if(mtrx_name == "yposOuter") 
      aptr = &yposOuter;
    else if(mtrx_name == "xposTileS") 
      aptr = &xposTileS;
    else if(mtrx_name == "yposTileS") 
      aptr = &yposTileS;
    else if(mtrx_name == "zposTileS") 
      aptr = &zposTileS;
    else if(mtrx_name == "xposTileN") 
      aptr = &xposTileN;
    else if(mtrx_name == "yposTileN") 
      aptr = &yposTileN;
    else if(mtrx_name == "zposTileN") 
      aptr = &zposTileN;
    else if(mtrx_name == "xposChimneyTileS") 
      aptr = &xposChimneyTileS;
    else if(mtrx_name == "yposChimneyTileS") 
      aptr = &yposChimneyTileS;
    else if(mtrx_name == "zposChimneyTileS") 
      aptr = &zposChimneyTileS;
    else{
      printout(WARNING, "BarrelHCalCalorimeter", "unknown <matrix> data!");
      continue;
    }
      
    std::string delimiter = " "; 
    size_t pos = 0;
    std::string token;
    while ((pos = mtrx_values.find(delimiter)) != std::string::npos) {
      token = mtrx_values.substr(0, pos);
      aptr->push_back(dd4hep::xml::_toDouble(token.c_str())); 
      mtrx_values.erase(0, pos + delimiter.length());
    }
    aptr->push_back(dd4hep::xml::_toDouble(mtrx_values.c_str())); 

  }

  // Loop over the solids, create them and add them to the detector volume

  for(xml_coll_t k(x_solids, _Unicode(solid)); k; ++k){
    
    xml_comp_t    x_solid = k; 

    // get the sector solid definitions
    xml_comp_t    define           = x_solid.child("define");
    xml_comp_t    tessellated      = x_solid.child("tessellated"); 
    std::string   solid_name       = dd4hep::xml::_toString(x_solid.attr_value(x_solid.getAttr("name"))); 
    std::string   solidMatString   = dd4hep::xml::_toString(x_solid.attr_value(x_solid.getAttr("material"))); 
    Material      solid_material   = description.material(solidMatString); 

    double offset_x = dd4hep::xml::_toDouble(x_solid.attr_value(x_solid.getAttr("x")))*dd4hep::mm; 
    double offset_y = dd4hep::xml::_toDouble(x_solid.attr_value(x_solid.getAttr("y")))*dd4hep::mm; 
    double offset_z = dd4hep::xml::_toDouble(x_solid.attr_value(x_solid.getAttr("z")))*dd4hep::mm; 

    // Get the vertices
    std::vector<Tessellated::Vertex_t> vertices; 
    for(xml_coll_t j(define, _Unicode(position)); j; ++j){
      xml_comp_t pos = j;

      // create the vertex point
      // Note that we trivially rotate everything about y by taking the negative of the x,z coordinate
      Tessellated::Vertex_t thisPoint(-(dd4hep::xml::_toDouble(pos.attr_value(pos.getAttr("x")))*dd4hep::mm - offset_x),
      				      dd4hep::xml::_toDouble(pos.attr_value(pos.getAttr("y")))*dd4hep::mm - offset_y,
      				      -(dd4hep::xml::_toDouble(pos.attr_value(pos.getAttr("z")))*dd4hep::mm - offset_z)); 
    
      vertices.push_back(thisPoint);

    }

    TessellatedSolid solid(solid_name.c_str(),vertices);

    for(xml_coll_t i(tessellated, _Unicode(triangular)); i; ++i){
      xml_comp_t triang = i; 

      int vtx1 = -1; 
      int vtx2 = -1; 
      int vtx3 = -1; 

      std::string facetName1 = dd4hep::xml::_toString(triang.attr_value(triang.getAttr("vertex1"))); 
      std::string facetName2 = dd4hep::xml::_toString(triang.attr_value(triang.getAttr("vertex2"))); 
      std::string facetName3 = dd4hep::xml::_toString(triang.attr_value(triang.getAttr("vertex3"))); 

      // Search the define collection to match things up
      int idx = 0; 
      for(xml_coll_t j(define, _Unicode(position)); j; ++j){
	xml_comp_t pos = j;
	std::string posName = dd4hep::xml::_toString(pos.attr_value(pos.getAttr("name"))); 

	if( posName == facetName1 ) vtx1 = idx; 
	if( posName == facetName2 ) vtx2 = idx; 
	if( posName == facetName3 ) vtx3 = idx; 
      
	if( (vtx1>=0) && (vtx2>=0) && (vtx3>=0) ) break; 

	idx++; 

      }

      // Add the facet to the solid

      if( (vtx1>=0) && (vtx2>=0) && (vtx3>=0) && (vtx1!=vtx2) && (vtx1!=vtx3) && (vtx2!=vtx3) )
	solid->AddFacet(vtx1,vtx2,vtx3); 
      else
	printout(WARNING, "BarrelHCalCalorimeter", "bad facet! %d %d %d", vtx1, vtx2, vtx3);

    }

    // Complete the shape
    solid->CloseShape(true,true,false); 
    Volume           solidVolume(solid_name, solid, solid_material);
    solidVolume.setVisAttributes(description, x_det.visStr());

    //printout(WARNING, "BarrelHCalCalorimeter", "tesselated solid name %s", solid_name.c_str());

    if(solid_name == "HCAL_Chimney_Sector_Half_Plate"){

      ChimneySector.placeVolume(solidVolume, 0, Transform3D(RotationZ(sectorRotate*dd4hep::deg), Translation3D(0, 0, 0) ));
 
      //ChimneySector.placeVolume(solidVolume, 1, Transform3D(RotationZ((9.5*2*M_PI / 320) + sectorRotate*dd4hep::deg), Translation3D(0, 0, 0) ));
 
    }
    else if(solid_name == "HCAL_Chimney_Sector_Plate"){

      for(int i=5; i<6; i++)
        ChimneySector.placeVolume(solidVolume, i, Transform3D(RotationZ((i*2*M_PI / 320) + sectorRotate*dd4hep::deg), Translation3D(0, 0, 0) ));

    }
    else if(solid_name == "HCAL_Sector_Half_Plate"){

      Sector.placeVolume(solidVolume, 0, Transform3D(RotationZ(sectorRotate*dd4hep::deg), Translation3D(0, 0, 0) ));
 
      Sector.placeVolume(solidVolume, 1, Transform3D(RotationZ((9.5*2*M_PI / 320) + sectorRotate*dd4hep::deg), Translation3D(0, 0, 0) ));

    }
    else if(solid_name == "HCAL_Sector_Plate"){

      for(int i=0; i<9; i++)
	Sector.placeVolume(solidVolume, i, Transform3D(RotationZ((i*2*M_PI / 320) + sectorRotate*dd4hep::deg), Translation3D(0, 0, 0) ));

    }
    else{

      // If it's not sectors then it's a tile - for these we build an assembly to get the full array of tiles
      // Offsets and rotation are to properly orient the tiles in the assembly

      if(solid_name.size()>0){

	std::string type = solid_name.substr(0,solid_name.size()-2); 

	if( type=="OuterHCalTile" || type=="OuterHCalChimneyTile" ){

	  std::string stnum = solid_name.substr(solid_name.size()-2,2);
	  int tnum = dd4hep::xml::_toInt(stnum.c_str())-1; 
	  
	  if(type=="OuterHCalTile"){

	    TileAssembly24.placeVolume(solidVolume,0,Transform3D(RotationY(0.0), 
								 Translation3D(xposTileS[tnum]*dd4hep::mm, yposTileS[tnum]*dd4hep::mm, zposTileS[tnum]*dd4hep::mm) ));

	    TileAssembly24.placeVolume(solidVolume,1,Transform3D(RotationY(180.0*dd4hep::deg),
	    							 Translation3D(xposTileN[tnum]*dd4hep::mm, yposTileN[tnum]*dd4hep::mm, zposTileN[tnum]*dd4hep::mm)  )); 

	    if(tnum<=7){
	      TileAssembly24Chimney.placeVolume(solidVolume,2,Transform3D(RotationY(0.0), 
									  Translation3D(xposTileS[tnum]*dd4hep::mm, yposTileS[tnum]*dd4hep::mm, zposTileS[tnum]*dd4hep::mm) ));
	      TileAssembly24Chimney.placeVolume(solidVolume,3,Transform3D(RotationY(180.0*dd4hep::deg), 
	      							  Translation3D(xposTileN[tnum]*dd4hep::mm, yposTileN[tnum]*dd4hep::mm, zposTileN[tnum]*dd4hep::mm) ));
	    }
	    else{
	      TileAssembly24Chimney.placeVolume(solidVolume,3,Transform3D(RotationY(180.0*dd4hep::deg), 
									  Translation3D(xposTileN[tnum]*dd4hep::mm, yposTileN[tnum]*dd4hep::mm, zposTileN[tnum]*dd4hep::mm) ));
	    }
	      

	  }

	  if( (tnum>7) && (type=="OuterHCalChimneyTile") )
	    TileAssembly24Chimney.placeVolume(solidVolume,0,Transform3D(RotationY(0.0),
	  							      Translation3D(xposChimneyTileS[tnum-8]*dd4hep::mm, yposChimneyTileS[tnum-8]*dd4hep::mm, zposChimneyTileS[tnum-8]*dd4hep::mm))); 
	  

	}
	else	 
          printout(WARNING, "BarrelHCalCalorimeter", "invalid solid_name, not a tile type?");

      }
      else
        printout(WARNING, "BarrelHCalCalorimeter", "solid_name.size() invalid! ");


    }

  }

  // Rotate the tile assemblies

  Assembly TileAssembly24Rotated("TileAssembly24Rotated"); 
  TileAssembly24Rotated.placeVolume(TileAssembly24,0,Transform3D(RotationX(tilePlaneRotate*dd4hep::deg), Translation3D(0.0,0.0,0.0) ));

  Assembly TileAssembly24ChimneyRotated("TileAssembly24ChimneyRotated"); 
  TileAssembly24ChimneyRotated.placeVolume(TileAssembly24Chimney,0,Transform3D(RotationX(tilePlaneRotate*dd4hep::deg), Translation3D(0.0,0.0,0.0) ));

  // Place the sectors into the mother volume

  sens.setType("calorimeter");
  DetElement sector_det("sector0", det_id);

  // Chimney sectors
  for(int i=-1; i<0; i++){
    PlacedVolume     sect_phv = envelope.placeVolume(ChimneySector, i+1, Transform3D(RotationZ((i*2*M_PI / 32)), Translation3D(0, 0, 0) ));
    sect_phv.addPhysVolID("system", det_id);
    sect_phv.addPhysVolID("barrel", 0);
    sect_phv.addPhysVolID("module", i+1);
    DetElement sd = sector_det.clone(_toString(i+1, "sector%d")); 
    sd.setPlacement(sect_phv);
    sdet.add(sd);
  }

  /*
  // Normal sectors
  for(int i=4; i<33; i++){
    PlacedVolume     env_phv = motherVol.placeVolume(Sector, i-1, RotationZYX((i*2*M_PI/32)+ sectorRotate*dd4hep::deg,0,0) );
    env_phv.addPhysVolID("system", det_id);
    env_phv.addPhysVolID("barrel", 0);
    env_phv.addPhysVolID("sector", i-1);
    sdet.setPlacement(env_phv);
  }
  */

  // Place the tile assemblies

  /*
  for(int i=0; i<18; i++){

    PlacedVolume     tile_phv = envelope.placeVolume(TileAssembly24ChimneyRotated, i, 
						       RotationZ((i-0.5)*(360.0/320.0)*dd4hep::deg)*
						       Transform3D(RotationY(90.0*dd4hep::deg), Translation3D(xposOuter[0]*dd4hep::mm, yposOuter[0]*dd4hep::mm, 0.0)) );
    tile_phv.addPhysVolID("system", det_id);
    tile_phv.addPhysVolID("barrel", 0);
    tile_phv.addPhysVolID("tilerow", i+1);
    DetElement sd = sector_det.clone(_toString(i+1, "tilerow%d")); 
    sd.setPlacement(tile_phv);
    sdet.add(sd);

  }
  */
  
  /*

  for(int i=18; i<308; i++){

    PlacedVolume     env_phv_c = motherVol.placeVolume(TileAssembly24Rotated, i, 
						       RotationZ((i-0.5)*(360.0/320.0)*dd4hep::deg)*
						       Transform3D(RotationY(90.0*dd4hep::deg), Translation3D(xposOuter[0]*dd4hep::mm, yposOuter[0]*dd4hep::mm, 0.0)) );
    env_phv_c.addPhysVolID("system", det_id);
    env_phv_c.addPhysVolID("barrel", 0);
    env_phv_c.addPhysVolID("row", i);
    sdet.setPlacement(env_phv_c);

  }

  for(int i=309; i<320; i++){

    PlacedVolume     env_phv_c = motherVol.placeVolume(TileAssembly24ChimneyRotated, i, 
						       RotationZ((i-0.5)*(360.0/320.0)*dd4hep::deg)*
						       Transform3D(RotationY(90.0*dd4hep::deg), Translation3D(xposOuter[0]*dd4hep::mm, yposOuter[0]*dd4hep::mm, 0.0)) );
    env_phv_c.addPhysVolID("system", det_id);
    env_phv_c.addPhysVolID("barrel", 0);
    env_phv_c.addPhysVolID("row", i);
    sdet.setPlacement(env_phv_c);

  }

  */

  envelope.setAttributes(description, x_det.regionStr(), x_det.limitsStr(), x_det.visStr());
  return sdet;

}

#ifdef EPIC_ECCE_LEGACY_COMPAT
DECLARE_DETELEMENT(ecce_oHcalBarrel, create_detector)
#endif
DECLARE_DETELEMENT(epic_oHcalBarrel, create_detector)
#ifdef EPIC_ECCE_LEGACY_COMPAT
DECLARE_DETELEMENT(ecce_HcalBarrel, create_detector)
#endif
DECLARE_DETELEMENT(epic_HcalBarrel, create_detector)

