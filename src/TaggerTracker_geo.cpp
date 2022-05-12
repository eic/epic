#include "DDRec/Surface.h"
#include "DDRec/DetectorData.h"
#include "DD4hep/OpticalSurfaces.h"
#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/Printout.h"
#include <XML/Helper.h>

//////////////////////////////////////////////////
// Low Q2 Tagger, tracking layers
//////////////////////////////////////////////////

using namespace std;
using namespace dd4hep;

static Ref_t createDetector(Detector& desc, xml_h e, SensitiveDetector sens)
{
  xml_det_t  x_det      = e;
  string     detName    = x_det.nameStr();
  int        detID      = x_det.id();

  xml_dim_t  dim        = x_det.dimensions();
  double     Width      = dim.x();
  double     Height     = dim.y();
  double     Thickness  = dim.z()/2;
  
  xml_dim_t  pos        = x_det.position();

  Material   Vacuum     = desc.material("Vacuum");
  Material   Silicon    = desc.material("Silicon");

  sens.setType("tracker");

  // Create Global Volume 
  Box Tagger_Box(Width, Height, Thickness);
  Volume detVol("Tagger_Box", Tagger_Box, Vacuum);
  detVol.setVisAttributes(desc.visAttributes(x_det.visStr()));

  //Add Hodoscope layers
  int N_layers      = 0;
  for(xml_coll_t lay( x_det, _Unicode(trackLayer) ); lay; ++lay, ++N_layers)  {
    
    string layerType      = dd4hep::getAttrOrDefault(lay,  _Unicode(type)             , "timepix"    );
    string layerVis       = dd4hep::getAttrOrDefault(lay,  _Unicode(vis)              , "GreenVis" );
    double layerZ         = dd4hep::getAttrOrDefault(lay,  _Unicode(z)                , 0            );
    double layerThickness = dd4hep::getAttrOrDefault(lay,  _Unicode(sensor_thickness) , 200*um       );

    Box Layer_Box(Width, Height, layerThickness/2);
    Volume layVol("LayerVolume", Layer_Box, Silicon);
    layVol.setSensitiveDetector(sens);
    layVol.setVisAttributes(desc.visAttributes(layerVis));

    //string module_name = detName + _toString(N_layers,"_TrackLayer_%d"); 
   
    PlacedVolume pv_mod = detVol.placeVolume(layVol, Position(0,0,Thickness-layerZ-layerThickness));
    pv_mod.addPhysVolID("layer",N_layers);
    
  }

  //mother volume for the tracker
  std::string mother_nam = dd4hep::getAttrOrDefault(x_det, _Unicode(place_into), "");
  VolumeManager man = VolumeManager::getVolumeManager(desc);
  DetElement mdet = man.detector().child(mother_nam);

  //placement in mother volume
  Transform3D tr(RotationZYX(0, 0, 0), Position(pos.x(), pos.y(), pos.z()));
  PlacedVolume detPV = mdet.volume().placeVolume(detVol, tr);
  detPV.addPhysVolID("system", detID);
  DetElement det(detName, detID);
  det.setPlacement(detPV);

  return det;
}

DECLARE_DETELEMENT(TaggerTracker, createDetector)
















