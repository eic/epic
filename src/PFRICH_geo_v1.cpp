// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Christopher Dilks, Sylvester Joosten

//----------------------------------
//  pfRICH: Proximity Focusing RICH
//  Author: C. Dilks
//----------------------------------

#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/OpticalSurfaces.h"
#include "DD4hep/Printout.h"
#include "DDRec/DetectorData.h"
#include "DDRec/Surface.h"

#include <XML/Helper.h>
#include "XML/Layering.h"

#include "TVector3.h"

// ROOT includes
#include "TGDMLParse.h"
#include "TGDMLWrite.h"
#include "TGeoElement.h"
#include "TGeoManager.h"
#include "TInterpreter.h"
#include "TUri.h"

using namespace std;
using namespace dd4hep;
using namespace dd4hep::rec;

using namespace dd4hep::detail;

static Ref_t createDetector(Detector& description, xml_h e, SensitiveDetector sens) {

///*--------------------------------------------------*/ 
///*--------------------------------------------------*/ 
///*--------------------------------------------------*/ 
 //  xml::DetElement       detElem = e;
 //  int                   detID   = detElem.id();
 ////  xml::Component        dims    = detElem.dimensions();
 
 //  DetElement            det(detName, detID);
 //  sens.setType("tracker");
 //
 
   xml_det_t x_det    = e;
   int       det_id   = x_det.id();
   int       detID   = det_id;
 
   string det_name = x_det.nameStr();
 //  Material  air      = description.air();
 
   DetElement sdet(det_name, det_id);
 
 //  std::vector<double> rmins = {rmin2, rmin2, rmin1, rmin1, rmin2, rmin2};
 //  std::vector<double> rmaxs = {rmax, rmax, rmax, rmax, rmax, rmax};
 //  std::vector<double> zs    = {-length2 / 2., -length1 / 2., -length1 / 2., length1 / 2., length1 / 2., length2 / 2.};
 
   std::vector<double> rmins = {100,  100};
   std::vector<double> rmaxs = {1300, 1300};
   std::vector<double> zs    = {-5000, 5000};
 
 
 /*--------------------------------------------------*/ 
 /*--------------------------------------------------*/ 
   Volume     motherVol = description.pickMotherVolume(sdet);
 //  Material  air      = description.air();
 //
 //  Polycone ptube(0.0, 2.0 * M_PI, rmins, rmaxs, zs);
 //  // Tube ptube(rmin1,rmax,length2);
 //  Volume envelope(det_name, ptube, air);
 //
 //  PlacedVolume env_phv = motherVol.placeVolume(envelope);
 //  env_phv.addPhysVolID("system", det_id);
 //  sdet.setPlacement(env_phv);
 
 
 
 
   sens.setType("tracker");
   description.invisible(); 
 
   xml::DetElement       detElem = e;
   std::string           detName = detElem.nameStr();
   xml::Component        dims    = detElem.dimensions();
 
 //  double vessel_radius   =  dims.attr<double>(_Unicode(rmax));;
 
   ///*--------------------------------------------------*/
   /// Porting in GDML volume
  
   xml_dim_t  x_gdml(x_det.child(_U(gdmlFile)));
   xml_dim_t  x_par(x_det.child(_U(parent)));
 
  string     name         = x_det.nameStr();
   string     gdml         = x_gdml.attr<string>(_U(ref));
   string     par_nam        = x_par.nameStr();
   DetElement det_parent   = description.detector(par_nam);
 
   TGDMLParse parser;
   if (!gdml.empty() && gdml[0] == '/') {
     TUri uri(gdml.c_str());
     gdml = uri.GetRelativePart();
   } else {
     string path = xml::DocumentHandler::system_path(e, gdml);
     TUri   uri(path.c_str());
     gdml = uri.GetRelativePart();
   }
 
   if (!det_parent.isValid()) {
     except(name, "+++ Cannot access detector parent: %s", par_nam.c_str());
   }
 
 //  DetElement sdet(name, id);
   Volume     volume = parser.GDMLReadFile(gdml.c_str());
   if (!volume.isValid()) {
     except("ROOTGDMLParse", "+++ Failed to parse GDML file:%s", gdml.c_str());
   }
   volume.import(); // We require the extensions in dd4hep.
   printout(INFO, "ROOTGDMLParse", "+++ Attach GDML volume %s", volume.name());
   Volume       mother = det_parent.volume();
   PlacedVolume pv;
 
   string     gdml_physvol = dd4hep::getAttrOrDefault<string>(x_gdml, _Unicode(physvol), "");
 
   if (!gdml_physvol.empty()) {
     PlacedVolume node = volume->FindNode(gdml_physvol.c_str());
     if (!node.isValid()) {
       printout(ERROR, "ROOTGDMLParse", "+++ Invalid gdml placed volume %s", gdml_physvol.c_str());
       printout(ERROR, "ROOTGDMLParse", "+++ Valid top-level nodes are:");
       volume->PrintNodes();
       except("ROOTGDMLParse", "+++ Failed to parse GDML file:%s for node:%s", gdml.c_str(), gdml_physvol.c_str());
     }
     volume = node.volume();
   }
 
 
   cout << volume.name() << endl;
 //  exit(0);
 
   int        id    = x_det.hasAttr(_U(id)) ? x_det.id() : 0;
   xml_dim_t  x_pos(x_det.child(_U(position), false));
   xml_dim_t  x_rot(x_det.child(_U(rotation), false));
 
 //  xml_dim_t  x_rot = detElem.rotation();
 
 ///*--------------------------------------------------*/ 
 //  if (x_pos && x_rot) {
 //    Rotation3D  rot(RotationZYX(x_rot.z(), x_rot.y(), x_rot.x()));
 //    Transform3D transform(rot, Position(x_pos.x(), x_pos.y(), x_pos.z()));
 //    pv = mother.placeVolume(volume, transform);
 //  } else if (x_rot) {
 //    Rotation3D  rot(RotationZYX(x_rot.z(), x_rot.y(), x_rot.x()));
 //    Transform3D transform(rot, Position(0, 0, 0));
 //    pv = mother.placeVolume(volume, transform);
 //  } else if (x_pos) {
 //    pv = mother.placeVolume(volume, Position(x_pos.x(), x_pos.y(), x_pos.z()));
 //  } else {
 //    pv = mother.placeVolume(volume);
 //  }
 
 
     Rotation3D  rot(RotationZYX(0, M_PI, 0));
     Transform3D transform(rot, Position(0, 0, -149));
     pv = mother.placeVolume(volume, transform);
 
 
 
 
   volume.setVisAttributes(description, x_det.visStr());
   volume.setLimitSet(description, x_det.limitsStr());
   volume.setRegion(description, x_det.regionStr());
   if (id != 0) {
     pv.addPhysVolID("system", id);
   }
   sdet.setPlacement(pv);
 
   // BUILD SENSORS ///////////////////////
 
   // solid and volume: single sensor module
 
   OpticalSurfaceManager surfMgr = description.surfaceManager();
 
//   auto sensorElem = detElem.child(_Unicode(sensors)).child(_Unicode(module));
//   auto sensorMat = description.material("VacuumOptical");
//   auto sensorVis = description.visAttributes(sensorElem.attr<std::string>(_Unicode(vis)));
 
  // - sensor module
  auto   sensorElem      = detElem.child(_Unicode(sensors)).child(_Unicode(module));
  auto   sensorMat       = description.material(sensorElem.attr<std::string>(_Unicode(material)));
  auto   sensorVis       = description.visAttributes(sensorElem.attr<std::string>(_Unicode(vis)));
  auto   sensorSurf      = surfMgr.opticalSurface(sensorElem.attr<std::string>(_Unicode(surface)));
  double sensorSide      = sensorElem.attr<double>(_Unicode(side));
  double sensorGap       = sensorElem.attr<double>(_Unicode(gap));
  double sensorThickness = sensorElem.attr<double>(_Unicode(thickness));
  auto   readoutName     = detElem.attr<std::string>(_Unicode(readout));



   double vesselZmin      = dims.attr<double>(_Unicode(zmin));
   double vesselRmin0     = dims.attr<double>(_Unicode(rmin0));
   double vesselRmin1     = dims.attr<double>(_Unicode(rmin1));
   double vesselRmax0     = dims.attr<double>(_Unicode(rmax0));
   double vesselRmax1     = dims.attr<double>(_Unicode(rmax1));
 
   int    imod           = 0;           // module number
   double tBoxMax        = vesselRmax1; // sensors will be tiled in tBox, within annular limits
 
 //  auto   sensorElem      = detElem.child(_Unicode(sensors)).child(_Unicode(module));
//   double sensorGap      = sensorElem.attr<double>(_Unicode(gap));

   auto   sensorPlaneElem = detElem.child(_Unicode(sensors)).child(_Unicode(plane));
   double sensorPlaneRmin = sensorPlaneElem.attr<double>(_Unicode(rmin));
   double sensorPlaneRmax = sensorPlaneElem.attr<double>(_Unicode(rmax));
 
//   double sensorSide      = sensorElem.attr<double>(_Unicode(side));
//   double sensorThickness = sensorElem.attr<double>(_Unicode(thickness));
 
   auto   gasvolMat       = description.material("C4F10_PFRICH");
   auto   gasvolVis       = description.visAttributes("DRICH_gas_vis");
   auto   vesselVis       = description.visAttributes("DRICH_gas_vis");
  
   double windowThickness = dims.attr<double>(_Unicode(window_thickness));
   double wallThickness = dims.attr<double>(_Unicode(wall_thickness));
 
   double proximityGap    = dims.attr<double>(_Unicode(proximity_gap));
 
   long debug_optics_mode = description.constantAsLong("PFRICH_debug_optics");
 
   // if debugging optics, override some settings
   bool debug_optics = debug_optics_mode > 0;
 
   // /*--------------------------------------------------*/ 
   // /*--------------------------------------------------*/ 
   // /*--------------------------------------------------*/ 
 
 
   // sensor plane positioning: we want `proximityGap` to be the distance between the
   // aerogel backplane (i.e., aerogel/filter boundary) and the sensor active surface (e.g, photocathode)
 
   auto   radiatorElem       = detElem.child(_Unicode(radiator));
   double radiatorFrontplane = radiatorElem.attr<double>(_Unicode(frontplane));
 
   auto   aerogelElem      = radiatorElem.child(_Unicode(aerogel));
   double aerogelThickness = aerogelElem.attr<double>(_Unicode(thickness));
 
   double radiatorRmin     = radiatorElem.attr<double>(_Unicode(rmin));
   double radiatorRmax     = radiatorElem.attr<double>(_Unicode(rmax));
 
   double radiatorPitch    = radiatorElem.attr<double>(_Unicode(pitch));
 
   double airgapThickness = 0.1;
   double filterThickness = 1;
 
   auto   aerogelMat       = description.material("C4F10_PFRICH");
   auto   filterMat       = description.material("C4F10_PFRICH");
 
 
   double vesselLength    = dims.attr<double>(_Unicode(length));
   auto originFront = Position(0., 0., vesselLength / 2.0);
   double sensorZpos     = radiatorFrontplane - aerogelThickness - proximityGap - 0.5 * sensorThickness;
   auto   sensorPlanePos = Position(0., 0., sensorZpos) + originFront; // reference position
 
   ///*--------------------------------------------------*/ 
   // readout coder <-> unique sensor ID
   /* - `sensorIDfields` is a list of readout fields used to specify a unique sensor ID
    * - `cellMask` is defined such that a hit's `cellID & cellMask` is the corresponding sensor's unique ID
    * - this redundant generalization is for future flexibility, and consistency with dRICH
    */
 
//   auto   readoutName     = "PFRICHHits";
 
   std::vector<std::string> sensorIDfields = {"module"};
   const auto&              readoutCoder   = *description.readout(readoutName).idSpec().decoder();
   // determine `cellMask` based on `sensorIDfields`
   uint64_t cellMask = 0;
   for (const auto& idField : sensorIDfields)
     cellMask |= readoutCoder[idField].mask();
   description.add(Constant("PFRICH_cell_mask", std::to_string(cellMask)));
   // create a unique sensor ID from a sensor's PlacedVolume::volIDs
   auto encodeSensorID = [&readoutCoder](auto ids) {
     uint64_t enc = 0;
     for (const auto& [idField, idValue] : ids)
       enc |= uint64_t(idValue) << readoutCoder[idField].offset();
     return enc;
   };
 
 //  // miscellaneous
 //  double vesselRmin0     = 63;
 //  double vesselRmin1     = 63;
 //  double vesselRmax0     = 65;
 //  double vesselRmax1     = 65;
 
 
 //  double radiatorFrontplane = -150;
 //  double aerogelThickness = 10;
 ////  double vesselLength    = dims.attr<double>(_Unicode(length));
 ////  double vesselRmax1     = 130;
 //  double proximityGap    = 30;
 
 
 
 //  ///*--------------------------------------------------*/ 
 //  // tank solids
 
   double boreDelta = vesselRmin1 - vesselRmin0;
   Cone   vesselTank(vesselLength / 2.0, vesselRmin1, vesselRmax1, vesselRmin0, vesselRmax0);
   Cone   gasvolTank(vesselLength / 2.0 - windowThickness, vesselRmin1 + wallThickness, vesselRmax1 - wallThickness,
                     vesselRmin0 + wallThickness, vesselRmax0 - wallThickness);
 //
 //  //  extra solids for `debug_optics` only
 //  Box vesselBox(1001, 1001, 1001);
   Box gasvolBox(1000, 1000, 1000);
 //
 //  // choose vessel and gasvol solids (depending on `debug_optics_mode` (0=disabled))
 //  Solid vesselSolid, gasvolSolid;
 //  switch (debug_optics_mode) {
 //  case 0:
 //    vesselSolid = vesselTank;
 //    gasvolSolid = gasvolTank;
 //    break; // `!debug_optics`
 //  case 1:
 //    vesselSolid = vesselBox;
 //    gasvolSolid = gasvolBox;
 //    break;
 //  case 2:
 //    vesselSolid = vesselBox;
 //    gasvolSolid = gasvolTank;
 //    break;
 //  };
 //
 
   Solid gasvolSolid;
   gasvolSolid = gasvolTank;
 
   Solid vesselSolid;
   vesselSolid = vesselTank;
 
 
   auto   vesselMat = description.material("VacuumOptical");
 
   Volume vesselVol(detName, vesselSolid, vesselMat);
   Volume gasvolVol(detName + "_gas", gasvolSolid, gasvolMat);
   vesselVol.setVisAttributes(vesselVis);
   gasvolVol.setVisAttributes(gasvolVis);
 
   // place mother volume (vessel)
 //  Volume       motherVol = description.pickMotherVolume(sdet);
 
   // place gas volume
   PlacedVolume gasvolPV = vesselVol.placeVolume(gasvolVol, Position(0, 0, 0));
   DetElement   gasvolDE(sdet, "gasvol_de", 0);
   gasvolDE.setPlacement(gasvolPV);
 
   auto vesselPos = Position(0, 0, vesselZmin) - originFront;
 
 //  PlacedVolume vesselPV  = motherVol.placeVolume(vesselVol, vesselPos);
 //  vesselPV.addPhysVolID("system", detID);
 //  sdet.setPlacement(vesselPV);
 
 
   // BUILD RADIATOR //////////////////////////////////////
 
   // solid and volume: create aerogel and filter
   Cone aerogelSolid(aerogelThickness / 2, radiatorRmin + boreDelta * aerogelThickness / vesselLength, /* at backplane */
                     radiatorRmax, radiatorRmin, /* at frontplane */
                     radiatorRmax);
   Cone filterSolid(filterThickness / 2,
                    radiatorRmin + boreDelta * (aerogelThickness + airgapThickness + filterThickness) /
                                       vesselLength,                                                /* at backplane */
                    radiatorRmax,
                    radiatorRmin + boreDelta * (aerogelThickness + airgapThickness) / vesselLength, /* at frontplane */
                    radiatorRmax);
   Volume aerogelVol(detName + "_aerogel", aerogelSolid, aerogelMat);
   Volume filterVol(detName + "_filter", filterSolid, filterMat);
 
 //  aerogelVol.setVisAttributes(aerogelVis);
 //  filterVol.setVisAttributes(filterVis);
 
 
   ///*--------------------------------------------------*/ 
   ///*--------------------------------------------------*/ 
   ///*--------------------------------------------------*/ 
   /// Radiator
 
   // aerogel placement and surface properties
   // TODO [low-priority]: define skin properties for aerogel and filter
   // FIXME: radiatorPitch might not be working correctly (not yet used)
   auto radiatorPos      = Position(0., 0., radiatorFrontplane - 0.5 * aerogelThickness) + originFront;
   auto aerogelPlacement = Translation3D(radiatorPos.x(), radiatorPos.y(), radiatorPos.z()-43) * // re-center to originFront
                           RotationY(radiatorPitch); // change polar angle to specified pitch
 //  auto       aerogelPV = gasvolVol.placeVolume(aerogelVol, aerogelPlacement);
 
   auto       aerogelPV = volume.placeVolume(aerogelVol, aerogelPlacement);
 
   DetElement aerogelDE(sdet, "aerogel_de", 0);
   aerogelDE.setPlacement(aerogelPV);
   // SkinSurface aerogelSkin(desc, aerogelDE, "mirror_optical_surface", aerogelSurf, aerogelVol);
   // aerogelSkin.isValid();
 
   // filter placement and surface properties
   PlacedVolume filterPV;
   if (!debug_optics) {
     auto filterPlacement =
         Translation3D(0., 0., -airgapThickness) *                          // add an airgap
         Translation3D(radiatorPos.x(), radiatorPos.y(), radiatorPos.z()-39) * // re-center to originFront
         RotationY(radiatorPitch) *                                         // change polar angle
         Translation3D(0., 0., -(aerogelThickness + filterThickness) / 2.); // move to aerogel backplane
 //    filterPV = gasvolVol.placeVolume(filterVol, filterPlacement);
     filterPV = volume.placeVolume(filterVol, filterPlacement);
 
 //    DetElement filterDE(sdet, "filter_de", 0);
 //    filterDE.setPlacement(filterPV);
 //
 //    // SkinSurface filterSkin(desc, filterDE, "mirror_optical_surface", filterSurf, filterVol);
 //    // filterSkin.isValid();
   };
 
   // radiator z-positions (w.r.t. IP)
   double aerogelZpos = vesselPos.z() + aerogelPV.position().z();
   double filterZpos  = vesselPos.z() + filterPV.position().z();
   description.add(Constant("PFRICH_aerogel_zpos", std::to_string(aerogelZpos)));
   description.add(Constant("PFRICH_filter_zpos", std::to_string(filterZpos)));
 
   // radiator material names
   description.add(Constant("PFRICH_aerogel_material", aerogelMat.ptr()->GetName(), "string"));
   description.add(Constant("PFRICH_filter_material", filterMat.ptr()->GetName(), "string"));
   description.add(Constant("PFRICH_gasvol_material", gasvolMat.ptr()->GetName(), "string"));
 
 
 // /*--------------------------------------------------*/ 
 
 ////
 //
 //
   // SENSOR MODULE LOOP ------------------------
   /* cartesian tiling loop
    * - start at (x=0,y=0), to center the grid
    * - loop over positive-x positions; for each, place the corresponding negative-x sensor too
    * - nested similar loop over y positions
    */
 
   Box    sensorSolid(sensorSide / 2., sensorSide / 2., sensorThickness / 2.);
   Volume sensorVol(detName + "_sensor", sensorSolid, sensorMat);
   sensorVol.setVisAttributes(sensorVis);
 
   // sensitivity
   if (!debug_optics)
     sensorVol.setSensitiveDetector(sens);
 
   double sx, sy;
   for (double usx = 0; usx <= tBoxMax; usx += sensorSide + sensorGap) {
     for (int sgnx = 1; sgnx >= (usx > 0 ? -1 : 1); sgnx -= 2) {
       for (double usy = 0; usy <= tBoxMax; usy += sensorSide + sensorGap) {
         for (int sgny = 1; sgny >= (usy > 0 ? -1 : 1); sgny -= 2) {
 
           // sensor (x,y) center
           sx = sgnx * usx;
           sy = sgny * usy;
 
           // annular cut
           if (std::hypot(sx, sy) < sensorPlaneRmin || std::hypot(sx, sy) > sensorPlaneRmax)
             continue;
 
           // placement (note: transformations are in reverse order)
           auto sensorPlacement = Transform3D(
               Translation3D(sensorPlanePos.x(), sensorPlanePos.y(), sensorPlanePos.z() + 31)  * // move to reference position
               Translation3D(sx, sy, 0.)                                                   // move to grid position
           );
 
 //		 cout << sensorPlanePos.x() << "  " << sensorPlanePos.y() << "  "<< sensorPlanePos.z() << endl;
 //          auto sensorPV = gasvolVol.placeVolume(sensorVol, sensorPlacement);
           auto sensorPV = volume.placeVolume(sensorVol, sensorPlacement);
 
           // generate LUT for module number -> sensor position, for readout mapping tests
           // printf("%d %f %f\n",imod,sensorPV.position().x(),sensorPV.position().y());
 
           // properties
           sensorPV.addPhysVolID("module", imod); // NOTE: must be consistent with `sensorIDfields`
           auto       imodEnc = encodeSensorID(sensorPV.volIDs());
           DetElement sensorDE(sdet, "sensor_de_" + std::to_string(imod), imodEnc);
           sensorDE.setPlacement(sensorPV);
           if (!debug_optics) {
             SkinSurface sensorSkin(description, sensorDE, "sensor_optical_surface_" + std::to_string(imod), sensorSurf,
                                    sensorVol);
             sensorSkin.isValid();
           };
 
           // increment sensor module number
           imod++;
         };
       };
     };
   };
   // END SENSOR MODULE LOOP ------------------------
 
 
 
 
 
 
 ///*--------------------------------------------------*/ 
 ///*--------------------------------------------------*/ 
 
   cout << "aaaaaaaaaaaaaaaaaaaaaaaa " << "    " << vesselRmax0 << endl;
 
   return sdet;
 







}

// clang-format off
DECLARE_DETELEMENT(epic_PFRICH_v1, createDetector)
