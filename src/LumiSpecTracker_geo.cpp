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
  // loop over layers
  for( xml_coll_t li(x_det, _U(layer)); li; li++) { // layers

    xml_comp_t x_layer( li );
    int layer_id = x_layer.id();

    // loop over sensors within each layer
    for( xml_coll_t sensor( li, _U(sensor)); sensor; sensor++) { // sensors

      xml_comp_t x_sensor( sensor );
      int sensor_id = x_sensor.id();

      double posX = x_sensor.position().x();
      double posY = x_sensor.position().y();
      double posZ = x_sensor.position().z();
      double sizeX = x_sensor.dimensions().x();
      double sizeY = x_sensor.dimensions().y();
      double sizeZ = x_sensor.dimensions().z();

      Box box( sizeX, sizeY, sizeZ );
      Volume vol( det_name + "_vol", box, description.material( "Silicon" ) );
      vol.setVisAttributes( description.visAttributes( x_det.visStr() ) );
      vol.setSensitiveDetector( sens );

      // place sensor into assembly
      PlacedVolume sensor_pv = assembly.placeVolume(
          vol, Transform3D( RotationZYX(0.0,0.0,0.0), Position( posX, posY, posZ ) ) );

      // Connect layer and sensor IDs
      sensor_pv.addPhysVolID("module", layer_id).addPhysVolID("layer", sensor_id);

    } // sensors
  } // layers

  // Place assembly into mother volume.  Assembly is centered at origin
  PlacedVolume detPV = motherVol.placeVolume( assembly, Position(0.0, 0.0, 0.0) );

  // Connect system ID
  detPV.addPhysVolID( "system", det_ID );

  det.setPlacement( detPV );

  return det;
}

DECLARE_DETELEMENT(LumiSpecTracker, create_detector)
