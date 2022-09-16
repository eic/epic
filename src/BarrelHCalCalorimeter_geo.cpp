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

  // get the solids section for this detector
  xml_comp_t    x_solids         = x_det.child("solids"); 

  DetElement       sdet(det_name, det_id);
  Volume           motherVol = description.pickMotherVolume(sdet);

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
      // Note that we trivially rotate everything by swapping the z coordinate
      Tessellated::Vertex_t thisPoint(dd4hep::xml::_toDouble(pos.attr_value(pos.getAttr("x")))*dd4hep::mm - offset_x,
				      dd4hep::xml::_toDouble(pos.attr_value(pos.getAttr("y")))*dd4hep::mm - offset_y,
				      -(dd4hep::xml::_toDouble(pos.attr_value(pos.getAttr("z")))*dd4hep::mm - offset_z)); 
    
      vertices.push_back(thisPoint);

    }

    TessellatedSolid solid(vertices);

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
    //printout(WARNING, "BarrelHCalCalorimeter", "calling sector->CloseShape() ");
    solid->CloseShape(true,true,false); 
    Volume           solidVolume(solid_name, solid, solid_material);
 
    if(solid_name == "HCAL_Chimney_Sector"){
      for(int i=-1; i<2; i++){
	PlacedVolume     env_phv = motherVol.placeVolume(solidVolume, Transform3D(RotationZ((i*2*M_PI / 32)*dd4hep::rad), Translation3D(0, 0, 0) ));
	env_phv.addPhysVolID("system", det_id);
	env_phv.addPhysVolID("barrel", 0);
	env_phv.addPhysVolID("sector", i+1);
	sdet.setPlacement(env_phv);
      }
    }
    else if(solid_name == "HCAL_Sector"){
      for(int i=4; i<33; i++){
	PlacedVolume     env_phv = motherVol.placeVolume(solidVolume, RotationZYX(i*2*M_PI/32,0,0) );
	env_phv.addPhysVolID("system", det_id);
	env_phv.addPhysVolID("barrel", 0);
	env_phv.addPhysVolID("sector", i-1);
	sdet.setPlacement(env_phv);
      }
    }

  }
  
  sens.setType("calorimeter");

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

