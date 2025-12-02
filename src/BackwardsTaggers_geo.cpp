// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022-2025 Simon Gardner

#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/OpticalSurfaces.h"
#include "DD4hep/Printout.h"
#include "DDRec/DetectorData.h"
#include "DDRec/Surface.h"
#include "XML/Utilities.h"

//////////////////////////////////////////////////
// Low Q2 tagger trackers
//////////////////////////////////////////////////

using namespace std;
using namespace dd4hep;
using namespace dd4hep::rec;

// Helper function to make the tagger tracker detectors
static void Make_Tagger(Detector& desc, xml_coll_t& mod, Assembly& env, DetElement modElement,
                        SensitiveDetector& sens);

static Ref_t create_detector(Detector& desc, xml_h e, SensitiveDetector sens) {

  xml_det_t x_det = e;
  string detName  = x_det.nameStr();
  int detID       = x_det.id();

  DetElement det(detName, detID);

  string vis_name = dd4hep::getAttrOrDefault<std::string>(x_det, _Unicode(vis), "BackwardsBox");

  // apply any detector type flags set in XML
  dd4hep::xml::setDetectorTypeFlag(x_det, det);

  // Create an Assembly to hold all modules
  Assembly taggerAssembly("BackwardsTaggerAssembly");

  // Add Tagger stations (modules) to the Assembly
  for (xml_coll_t mod(x_det, _Unicode(module)); mod; ++mod) {
    int moduleID      = dd4hep::getAttrOrDefault<int>(mod, _Unicode(id), 0);
    string moduleName = dd4hep::getAttrOrDefault<std::string>(mod, _Unicode(name), "Tagger0");

    // Size of the actual tagger box, replicated in BackwardsTagger
    xml_dim_t pos = mod.child(_U(position));
    xml_dim_t rot = mod.child(_U(rotation));

    // Create a DetElement for the module
    DetElement moddet(det, moduleName, moduleID);

    // Create an Assembly for the module contents
    Assembly moduleAssembly(moduleName + "_assembly");

    // Fill the module assembly with tagger layers
    Make_Tagger(desc, mod, moduleAssembly, moddet, sens);

    // Place the module assembly into the taggerAssembly at the specified position/rotation
    PlacedVolume pv_mod = taggerAssembly.placeVolume(
        moduleAssembly, Transform3D(RotationY(rot.y()), Position(pos.x(), pos.y(), pos.z())));
    pv_mod.addPhysVolID("module", moduleID);
    moddet.setPlacement(pv_mod);
  }

  // Place the taggerAssembly into the mother volume
  Volume motherVol       = desc.pickMotherVolume(det);
  PlacedVolume pv_tagger = motherVol.placeVolume(taggerAssembly, Position(0, 0, 0));
  pv_tagger.addPhysVolID("system", detID);
  det.setPlacement(pv_tagger);

  return det;
}

static void Make_Tagger(Detector& desc, xml_coll_t& mod, Assembly& env, DetElement modElement,
                        SensitiveDetector& sens) {

  sens.setType("tracker");

  Material Silicon = desc.material("Silicon");

  xml_dim_t moddim = mod.child(_Unicode(dimensions));
  double tag_w     = moddim.x() / 2;
  double tag_h     = moddim.y() / 2;

  // Add Hodoscope layers
  int N_layers = 0;
  for (xml_coll_t lay(mod, _Unicode(trackLayer)); lay; ++lay, ++N_layers) {

    int layerID      = dd4hep::getAttrOrDefault<int>(lay, _Unicode(id), 0);
    string layerType = dd4hep::getAttrOrDefault<std::string>(lay, _Unicode(type), "timepix");
    string layerVis =
        dd4hep::getAttrOrDefault<std::string>(lay, _Unicode(vis), "FFTrackerLayerVis");
    double layerRot = dd4hep::getAttrOrDefault<double>(lay, _Unicode(angle), 0);
    double layerZ   = dd4hep::getAttrOrDefault<double>(lay, _Unicode(z), 0 * mm);
    double layerThickness =
        dd4hep::getAttrOrDefault<double>(lay, _Unicode(sensor_thickness), 200 * um);

    RotationY rotate(layerRot);

    Box Layer_Box(tag_w, tag_h, layerThickness / 2);
    Volume layVol("Tagger_tracker_layer", Layer_Box, Silicon);
    layVol.setSensitiveDetector(sens);
    layVol.setVisAttributes(desc.visAttributes(layerVis));

    PlacedVolume pv_layer =
        env.placeVolume(layVol, Transform3D(rotate, Position(0, 0, layerZ - layerThickness / 2)));
    pv_layer.addPhysVolID("layer", layerID);

    DetElement laydet(modElement, "layerName" + std::to_string(layerID), layerID);
    laydet.setPlacement(pv_layer);
  }
}

DECLARE_DETELEMENT(BackwardsTagger, create_detector)
