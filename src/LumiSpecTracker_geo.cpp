// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Dhevan Gangadharan

#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/Printout.h"
#include <XML/Helper.h>

using namespace std;
using namespace dd4hep;

static Ref_t create_detector(Detector& description, xml_h e, SensitiveDetector sens)
{
  sens.setType("tracker");

  xml_det_t     x_det           = e;
  string        det_name        = x_det.nameStr();
  int           det_ID          = x_det.id();
  string        vis_Si          = getAttrOrDefault<string>(x_det, _Unicode(vis), "TrackerVis");
  string        vis_Cu          = getAttrOrDefault<string>(x_det, _Unicode(visCu), "TrackerServiceVis");
  string        Si_name         = getAttrOrDefault<string>(x_det, _Unicode(materialSi), "SiliconOxide");
  string        Cu_name         = getAttrOrDefault<string>(x_det, _Unicode(materialCu), "Copper");
  double        Si_DZ           = getAttrOrDefault<double>(x_det, _Unicode(thicknessSi), 0.03/2.0);
  double        Cu_DZ           = getAttrOrDefault<double>(x_det, _Unicode(thicknessCu), 0.014/2.0);

  Material      m_Si     = description.material( Si_name );
  Material      m_Cu     = description.material( Cu_name );

  // Create main detector element to be returned at the end
  DetElement    det( det_name, det_ID );

  // Mother volume
  Volume        motherVol = description.pickMotherVolume( det );

  // Detector assembly
  Assembly      assembly( det_name );
  assembly.setVisAttributes( description.invisible() );

  // Build detector components
  // loop over modules
  for( xml_coll_t mi(x_det, _Unicode(module)); mi; mi++) { // modules

    xml_comp_t x_mod( mi );
    int module_id = x_mod.id();

    // loop over sectors within each module
    for( xml_coll_t si( mi, _Unicode(sector)); si; si++) { // sectors

      xml_comp_t x_sector( si );
      int sector_id = x_sector.id();
      string name = getAttrOrDefault<string>(x_sector, _Unicode(name), "");

      double posX = x_sector.position().x();
      double posY = x_sector.position().y();
      double posZ = x_sector.position().z();
      double sizeX = x_sector.dimensions().x();
      double sizeY = x_sector.dimensions().y();

      // Silicon sensor
      Box box_Si( sizeX, sizeY, Si_DZ );
      Volume vol_Si( det_name + "_" + name, box_Si, m_Si );
      vol_Si.setVisAttributes( description.visAttributes( vis_Si ) );
      vol_Si.setSensitiveDetector( sens );

      // Cu layer to approximate ASICs/cooling
      Box box_Cu( sizeX, sizeY, Cu_DZ );
      Volume vol_Cu( det_name + "_" + name + "_Cu", box_Cu, m_Cu );
      vol_Cu.setVisAttributes( description.visAttributes( vis_Cu ) );

      // place into assembly
      PlacedVolume pv = assembly.placeVolume(
          vol_Si, Transform3D( RotationZYX(0.0,0.0,0.0), Position( posX, posY, posZ ) ) );
      assembly.placeVolume(
          vol_Cu, Transform3D( RotationZYX(0.0,0.0,0.0), Position( posX, posY, posZ - (Si_DZ+Cu_DZ) ) ) );

      // Connect sector and module IDs
      pv.addPhysVolID("sector", sector_id).addPhysVolID("module", module_id);

    } // sectors
  } // modules

  // Place assembly into mother volume.  Assembly is centered at origin
  PlacedVolume detPV = motherVol.placeVolume( assembly, Position(0.0, 0.0, 0.0) );

  // Connect system ID
  detPV.addPhysVolID( "system", det_ID );

  det.setPlacement( detPV );

  return det;
}

DECLARE_DETELEMENT(LumiSpecTracker, create_detector)
