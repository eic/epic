#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/Printout.h"
#include "TMath.h"
#include <XML/Helper.h>

using namespace std;
using namespace dd4hep;

static Ref_t create_detector(Detector& description, xml_h e, SensitiveDetector sens)
{
  sens.setType("tracker");

  xml_det_t     x_det           = e;
  string        det_name        = x_det.nameStr();
  int           det_ID          = x_det.id();

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

      double posX = x_sector.position().x();
      double posY = x_sector.position().y();
      double posZ = x_sector.position().z();
      double sizeX = x_sector.dimensions().x();
      double sizeY = x_sector.dimensions().y();
      double sizeZ = x_sector.dimensions().z();

      Box box( sizeX, sizeY, sizeZ );
      Volume vol( det_name + "_vol", box, description.material( "Silicon" ) );
      vol.setVisAttributes( description.visAttributes( x_det.visStr() ) );
      vol.setSensitiveDetector( sens );

      // place into assembly
      PlacedVolume pv = assembly.placeVolume(
          vol, Transform3D( RotationZYX(0.0,0.0,0.0), Position( posX, posY, posZ ) ) );

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
