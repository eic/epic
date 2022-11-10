// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Simon Gardner, Wouter Deconinck

#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/OpticalSurfaces.h"
#include "DD4hep/Printout.h"
#include "DDRec/DetectorData.h"
#include "DDRec/Surface.h"
#include <XML/Helper.h>

//////////////////////////////////////////////////
// Vacuum drift volume in far backwards region
//////////////////////////////////////////////////

using namespace std;
using namespace dd4hep;

static Ref_t create_detector(Detector& desc, xml_h e, SensitiveDetector sens)
{
  xml_det_t x_det   = e;
  string    detName = x_det.nameStr();
  int       detID   = x_det.id();

  sens.setType("tracker");

  std::string   mother_name  = dd4hep::getAttrOrDefault<std::string>(x_det, _Unicode(place_into), "");
  std::string   mother_name2 = dd4hep::getAttrOrDefault<std::string>(x_det, _Unicode(place_into2), "");
  VolumeManager man          = VolumeManager::getVolumeManager(desc);

  // Add check tha mother volumes exist
  DetElement mdet1    = man.detector().child(mother_name);
  int        motherid = mdet1.id();

  DetElement mdet = man.detector().child(mother_name).child(mother_name2);

  DetElement det(detName, detID);

  // Materials
  // Material Vacuum  = desc.material("Vacuum");
  Material Air     = desc.material("Air");
  Material Silicon = desc.material("Silicon");
  // Material Steel   = desc.material("StainlessSteel");

  double tag_h   = dd4hep::getAttrOrDefault(x_det, _Unicode(height), 100 * mm);
  double tag_w   = dd4hep::getAttrOrDefault(x_det, _Unicode(width), 100 * mm);
  double tagboxL = dd4hep::getAttrOrDefault(x_det, _Unicode(length), 100 * mm);

  Assembly taggerAssembly("tagAssembly");

  Volume Tagger_Air;

  double airThickness    = 0;
  double vacuumThickness = tagboxL;

  // Add window layer and air-vacuum boxes
  for (xml_coll_t lay(x_det, _Unicode(windowLayer)); lay; ++lay) {

    string layerType      = dd4hep::getAttrOrDefault<std::string>(lay, _Unicode(type), "window");
    string layerVis       = dd4hep::getAttrOrDefault<std::string>(lay, _Unicode(vis), "FFTrackerShieldingVis");
    double layerZ         = dd4hep::getAttrOrDefault<double>(lay, _Unicode(z), 0 * mm);
    double layerThickness = dd4hep::getAttrOrDefault<double>(lay, _Unicode(sensor_thickness), 1 * mm);
    string layerMaterial  = dd4hep::getAttrOrDefault<std::string>(lay, _Unicode(material), "Copper");

    Material WindowMaterial = desc.material(layerMaterial);

    airThickness    = tagboxL - layerZ;
    vacuumThickness = tagboxL - airThickness;

    Box Box_Air(tag_w, tag_h, airThickness / 2);
    Tagger_Air = Volume("AirVolume", Box_Air, Air);
    Tagger_Air.setVisAttributes(desc.visAttributes("BackwardsAir"));
    // mdet.volume().placeVolume(Tagger_Air, Position(0, 0, -tagboxL/2 + airThickness/2));

    Box    Window_Box(tag_w, tag_h, layerThickness / 2);
    Volume layVol("WindowVolume", Window_Box, WindowMaterial);
    layVol.setVisAttributes(desc.visAttributes(layerVis));

    Tagger_Air.placeVolume(layVol, Position(0, 0, airThickness / 2 - layerThickness / 2));

    taggerAssembly.placeVolume(Tagger_Air, Position(0, 0, tagboxL - airThickness / 2));
  }

  // Add Hodoscope layers
  int N_layers = 0;
  for (xml_coll_t lay(x_det, _Unicode(trackLayer)); lay; ++lay, ++N_layers) {

    int    layerID        = dd4hep::getAttrOrDefault<int>(lay, _Unicode(id), 0);
    string layerType      = dd4hep::getAttrOrDefault<std::string>(lay, _Unicode(type), "timepix");
    string layerVis       = dd4hep::getAttrOrDefault<std::string>(lay, _Unicode(vis), "FFTrackerLayerVis");
    double layerZ         = dd4hep::getAttrOrDefault<double>(lay, _Unicode(z), 0 * mm);
    double layerThickness = dd4hep::getAttrOrDefault<double>(lay, _Unicode(sensor_thickness), 200 * um);

    Volume mother          = taggerAssembly;
    double MotherThickness = tagboxL;

    if (layerZ > vacuumThickness) {
      mother = Tagger_Air;
      layerZ -= vacuumThickness;
      MotherThickness = airThickness / 2;
    }

    Box    Layer_Box(tag_w, tag_h, layerThickness / 2);
    Volume layVol("TrackerVolume", Layer_Box, Silicon);
    layVol.setSensitiveDetector(sens);
    layVol.setVisAttributes(desc.visAttributes(layerVis));

    // string module_name = detName + _toString(N_layers,"_TrackLayer_%d");

    PlacedVolume pv_layer = mother.placeVolume(layVol, Position(0, 0, MotherThickness - layerZ + layerThickness / 2));
    pv_layer.addPhysVolID("layer", layerID);
  }

  //     // Add Calorimeter layers
  //     for (xml_coll_t lay(mod, _Unicode(calorimeter)); lay; ++lay, ++N_layers) {

  //       int    layerID        = dd4hep::getAttrOrDefault(lay, _Unicode(id),   0);
  //       string layerType      = dd4hep::getAttrOrDefault(lay, _Unicode(type), "TaggerCalPbWO4");
  //       string layerVis       = dd4hep::getAttrOrDefault(lay, _Unicode(vis),  "RedVis");
  //       double layerThickness = dd4hep::getAttrOrDefault(lay, _Unicode(thickness), 180 * mm);

  //       Box    Layer_Box(tag_w,tag_h, layerThickness / 2);
  //       Volume layVol("CalorimeterVolume", Layer_Box, desc.material("PbWO4"));
  //       layVol.setSensitiveDetector(sens);
  //       layVol.setVisAttributes(desc.visAttributes(layerVis));

  //       // string module_name = detName + _toString(N_layers,"_TrackLayer_%d");

  //       PlacedVolume pv_mod = Tagger_Air.placeVolume(layVol, Position(-wall/2, 0, -airThickness/2 +
  //       layerThickness/2)); pv_mod.addPhysVolID("layer", layerID+100);
  //     }
  // return mdet;

  PlacedVolume pv_mod = mdet.volume().placeVolume(taggerAssembly);

  pv_mod.addPhysVolID("module", detID);
  pv_mod.addPhysVolID("sensor", 1);
  pv_mod.addPhysVolID("system", motherid);
  det.setPlacement(pv_mod);
  desc.declareParent(detName, mdet);

  return det;
}

DECLARE_DETELEMENT(BackwardsTracker, create_detector)
