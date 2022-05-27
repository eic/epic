
#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/OpticalSurfaces.h"
#include "DD4hep/Printout.h"
#include "DDRec/DetectorData.h"
#include "DDRec/Surface.h"
#include <XML/Helper.h>

//////////////////////////////////////////////////
// Low Q2 Tagger
//////////////////////////////////////////////////

using namespace std;
using namespace dd4hep;

static Ref_t createDetector(Detector& desc, xml_h e, SensitiveDetector sens)
{
  xml_det_t x_det   = e;
  string    detName = x_det.nameStr();
  int       detID   = x_det.id();

  xml_dim_t dim       = x_det.dimensions();
  double    Width     = dim.x();
  double    Height    = dim.y();
  double    Thickness = dim.z() / 2;

  xml_dim_t pos = x_det.position();
  // xml_dim_t  rot        = x_det.rotation();

  Material Vacuum   = desc.material("Vacuum");
  Material Silicon  = desc.material("Silicon");
  Material Absorber = desc.material("TungstenDens24");

  sens.setType("calorimeter");

  // Create Global Volume
  Box    Tagger_Box(Width, Height, Thickness);
  Volume detVol("Tagger_Box", Tagger_Box, Vacuum);
  detVol.setVisAttributes(desc.visAttributes(x_det.visStr()));

  // Calorimeter
  double abso_z = 8.5 * mm; // for 20 layers
  double sens_z = 0.3 * mm;
  // double calo_start = 30*mm;

  for (int i = 0; i < 20; i++) {

    // absorber
    Box    Abso_Box(Width, Height, abso_z / 2);
    Volume absoVol("AbsorberVolume", Abso_Box, Absorber);
    absoVol.setVisAttributes(desc.visAttributes("BlueVis"));

    detVol.placeVolume(absoVol, Position(0, 0, Thickness - (i) * (abso_z + sens_z) - abso_z / 2));

    // sensitive layer
    Box    Cal_Box(Width, Height, sens_z / 2);
    Volume calVol("SensVolume", Cal_Box, Silicon);
    calVol.setSensitiveDetector(sens);
    calVol.setVisAttributes(desc.visAttributes("RedVis"));

    PlacedVolume pv_mod =
        detVol.placeVolume(calVol, Position(0, 0, Thickness - (i) * (abso_z + sens_z) - abso_z - sens_z / 2));
    pv_mod.addPhysVolID("layer", i + 3); // leave room for tracking IDs
  }

  // mother volume for calorimeter
  std::string   mother_nam = dd4hep::getAttrOrDefault(x_det, _Unicode(place_into), "");
  VolumeManager man        = VolumeManager::getVolumeManager(desc);
  DetElement    mdet       = man.detector().child(mother_nam);

  // placement in envelope
  Transform3D  tr(RotationZYX(0, 0, 0), Position(pos.x(), pos.y(), pos.z()));
  PlacedVolume detPV = mdet.volume().placeVolume(detVol, tr);
  detPV.addPhysVolID("system", detID);
  DetElement det(detName, detID);
  det.setPlacement(detPV);

  return det;
}

DECLARE_DETELEMENT(TaggerCalWSi, createDetector)
