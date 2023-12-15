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
   Material  air      = description.air();
 
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

   cout << volume.name() << endl;


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

       cout << "aaaa" << endl;
     }
     volume = node.volume();

     cout << "bbbb" << endl;

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
 
  auto mirrorElem      = detElem.child(_Unicode(mirror));
  auto mirrorMat       = description.material(mirrorElem.attr<std::string>(_Unicode(material)));
  auto mirrorVis       = description.visAttributes(mirrorElem.attr<std::string>(_Unicode(vis)));
  auto mirrorSurf      = surfMgr.opticalSurface(mirrorElem.attr<std::string>(_Unicode(surface)));

   Cone   mirror_cone(vesselLength / 2.0, vesselRmax1-7, vesselRmax1-7+0.3, vesselRmax1-13, vesselRmax1-13+0.3);



 
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
 

   Solid mirrorSolid;
   mirrorSolid = mirror_cone;
 
   auto   vesselMat = description.material("VacuumOptical");
 
   Volume vesselVol(detName, vesselSolid, vesselMat);
   Volume gasvolVol(detName + "_gas", gasvolSolid, gasvolMat);
   vesselVol.setVisAttributes(vesselVis);
   gasvolVol.setVisAttributes(gasvolVis);

   Volume mirrorVol(detName, mirrorSolid, mirrorMat);
   mirrorVol.setVisAttributes(mirrorVis);

   // place gas volume
   PlacedVolume gasvolPV = vesselVol.placeVolume(gasvolVol, Position(0, 0, 0));
   DetElement   gasvolDE(sdet, "gasvol_de", 0);
   gasvolDE.setPlacement(gasvolPV);

   // place mother volume (vessel)
 //  Volume       motherVol = description.pickMotherVolume(sdet);
 
   // place mirror volume
   PlacedVolume mirrorPV = volume.placeVolume(mirrorVol, Position(0, 0, 0));
   DetElement   mirrorDE(sdet, "mirror_de", 0);
   mirrorDE.setPlacement(mirrorPV);

   SkinSurface mirrorSkin(description, mirrorDE, "mirror_optical_surface_cone", mirrorSurf, mirrorVol);
   mirrorSkin.isValid();

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

   float _FLANGE_EPIPE_DIAMETER_ = 10.53;   // in cm
   float _FLANGE_HPIPE_DIAMETER_ = 4.47;    // in cm
   float _FLANGE_HPIPE_OFFSET_ = 6.76;    // in cm
   float clearance = 0.5;                   // in cm
   float length = 45;                       // in cm


/////*--------------------------------------------------*/ 
//  /// Inner mirror cone
//  //
//  // FIXME: do I really care about re-using the same names for these shapes?;
//  // A wedge bridging two cylinders;

  Tube   eflange(0.0, _FLANGE_EPIPE_DIAMETER_/2 + clearance, length/2);
  Tube   hflange(0.0, _FLANGE_HPIPE_DIAMETER_/2 + clearance, length/2);

  double r0 = _FLANGE_EPIPE_DIAMETER_/2 + clearance; 
  double r1 = _FLANGE_HPIPE_DIAMETER_/2 + clearance;
  double L = _FLANGE_HPIPE_OFFSET_;
  double a = r0*L/(r0-r1);
  double b = r0*r0/a;
  double c = r1*(a-b)/r0;

  // GEANT variables to define G4Trap;
  double pDz = length/2, pTheta = 0.0, pPhi = 0.0, pDy1 = (a - b - c)/2, pDy2 = pDy1; 
  double pDx1 = sqrt(r0*r0 - b*b), pDx2 = pDx1*r1/r0, pDx3 = pDx1, pDx4 = pDx2, pAlp1 = 0.0, pAlp2 = 0.0;

  Trap wedge(pDz, pTheta, pPhi, pDy1, pDx1, pDx2, pAlp1, pDy2, pDx3, pDx4, pAlp2);
 
  // flange_shape = new G4UnionSolid("Tmp", eflange, hflange, 0, G4ThreeVector(-_FLANGE_HPIPE_OFFSET_, 0.0, 0.0));
//  Rotation3D  rot(RotationZYX(0, M_PI, 0));
//  Transform3D transform(rot, Position(0, 0, -149));

  UnionSolid flange_shape(eflange, hflange, Position(-_FLANGE_HPIPE_OFFSET_, 0.0, 0.0));

  Rotation3D rZ(RotationZYX(M_PI/2.0, 0.0, 0.0));

  Transform3D transform_flange(rZ, Position(-b - pDy1, 0.0, 0.0));

  UnionSolid flange_final_shape(flange_shape, wedge, transform_flange);

//  Volume flangeVol(detName + "_flange", flange_shape, mirrorMat);
  Volume flangeVol(detName + "_flange", flange_final_shape, mirrorMat);
  flangeVol.setVisAttributes(mirrorVis);

   PlacedVolume flangePV = volume.placeVolume(flangeVol, Position(0, 0, 0));
   DetElement   flangeDE(sdet, "flange_de", 0);
   flangeDE.setPlacement(flangePV);

  ///*--------------------------------------------------*/ 
  /// SENSOR MODULE LOOP ------------------------
  /* cartesian tiling loop
   * - start at (x=0,y=0), to center the grid
   * - loop over positive-x positions; for each, place the corresponding negative-x sensor too
   * - nested similar loop over y positions
   */
 
   Box    sensorSolid(sensorSide / 2., sensorSide / 2., sensorThickness / 2.);
   Volume sensorVol(detName + "_sensor", sensorSolid, sensorMat);
   sensorVol.setVisAttributes(sensorVis);
 

  ///*--------------------------------------------------*/ 
  ///*--------------------------------------------------*/ 
  /// Detailed sensor description 
 
  double _FIDUCIAL_VOLUME_LENGTH_ = 49.1;         // cm
  double _SENSOR_AREA_LENGTH_ = 5;                // cm
  double _HRPPD_CENTRAL_ROW_OFFSET_ = 4.0;        // cm
  double _HRPPD_WINDOW_THICKNESS_ = 0.38;         // cm
  double _HRPPD_CONTAINER_VOLUME_HEIGHT_ = 3.2;   // cm
  double _HRPPD_INSTALLATION_GAP_ = 2.5;          // cm

  double _HRPPD_TILE_SIZE_ = 12.0;                // cm
  double _HRPPD_OPEN_AREA_SIZE_ = 11.4;           // cm
  double _HRPPD_ACTIVE_AREA_SIZE_  = 10.8;        // cm
  double _HRPPD_CERAMIC_BODY_THICKNESS_ = 0.9;    // cm
  double _HRPPD_BASEPLATE_THICKNESS_ = 0.3;       // cm
  double _HRPPD_PLATING_LAYER_THICKNESS_ = 0.006; // cm
  double _EFFECTIVE_MCP_THICKNESS_ = 2*0.06*0.3;  // cm

  double _READOUT_PCB_THICKNESS_ = 0.2;
  double _READOUT_PCB_SIZE_ = _HRPPD_OPEN_AREA_SIZE_ - 0.2;

  double _ASIC_SIZE_XY_ = 1.6;
  double _ASIC_THICKNESS_ = 0.1;

  double azOffset = _FIDUCIAL_VOLUME_LENGTH_/2 - _SENSOR_AREA_LENGTH_;
  double xysize = _HRPPD_TILE_SIZE_, wndthick = _HRPPD_WINDOW_THICKNESS_, zwnd = azOffset + wndthick/2;
  
  // HRPPD assembly container volume; 
  double hrppd_container_volume_thickness = _HRPPD_CONTAINER_VOLUME_HEIGHT_;
  double zcont = azOffset + hrppd_container_volume_thickness/2;
  double zcables = azOffset + hrppd_container_volume_thickness + 1.0;///2;

  // XY-size is calculated based on the ASIC size and location;
  double _COLD_PLATE_THICKNESS_  = 0.15;
  
  // Water cooling pipe;
  double  _COOLING_PIPE_INNER_DIAMETER_ = 0.55372;  // cm
  double  _COOLING_PIPE_WALL_THICKNESS_ = 0.04064;  // cm




  ///*--------------------------------------------------*/ 
  // HRPPD

  Box hrppd_Solid(xysize/2, xysize/2, hrppd_container_volume_thickness/2);

  Volume hrppdVol_air(detName + "_air_hrppd", hrppd_Solid, air);
  Volume hrppdVol(detName + "_hrppd", hrppd_Solid, sensorMat);

  hrppdVol_air.setVisAttributes(gasvolVis);

//  VisAttr attr("visible");
//  attr.setColor(1.0, 0.5, 0.5, 0.5);
//  attr.setLineStyle(VisAttr::SOLID);
//  attr.setDrawingStyle(VisAttr::SOLID);
//  attr.setVisible(true);
//  attr.setShowDaughters(true);

//  hrppdVol.setVisAttributes(attr);

//  PlacedVolume hrppdairPV = volume.placeVolume(hrppdVol_air, Position(0, 0, 30));
//  DetElement   hrppdairDE(sdet, "hrppdair_de", 0);
//  hrppdairDE.setPlacement(hrppdairPV);

//  PlacedVolume hrppdPV = hrppdVol_air.placeVolume(hrppdVol, Position(0, 0, 0));
//  DetElement   hrppdDE(sdet, "hrppd_de", 0);

//  hrppdVol.setVisAttributes(sensorVis);
//  hrppdDE.setPlacement(hrppdairPV);

  ///*--------------------------------------------------*/ 
  // Quartz Window
  Box wnd_Solid(xysize/2, xysize/2, wndthick/2);

  Volume wndVol(detName + "_wnd", wnd_Solid, gasvolMat);
  wndVol.setVisAttributes(gasvolVis);
//  wndVol.setVisAttributes(attr);

  double accu = -hrppd_container_volume_thickness/2;

  PlacedVolume wndPV = hrppdVol_air.placeVolume(wndVol, Position(0, 0, accu + wndthick/2));
//  DetElement   wndDE(hrppdDE, "wnd_de", 0);
//  wndDE.setPlacement(wndPV);

  double pitch = xysize + _HRPPD_INSTALLATION_GAP_, xyactive = _HRPPD_ACTIVE_AREA_SIZE_;
  double xyopen = _HRPPD_OPEN_AREA_SIZE_;
  double certhick = _HRPPD_CERAMIC_BODY_THICKNESS_;//, zcer = azOffset + wndthick + certhick/2;

//  G4Box *cerbox  = new G4Box("CeramicBox", xysize/2, xysize/2, certhick/2);
//  G4Box *cut_box  = new G4Box("CeramicCut", xyopen/2, xyopen/2, certhick/2);
//  auto ceramic = new G4SubtractionSolid("CeramicBody", cerbox, cut_box, 0, 
//					G4ThreeVector(0,0, -_HRPPD_BASEPLATE_THICKNESS_));

  accu += wndthick;


  ///*--------------------------------------------------*/ 
  // Ceramic body
  Box cerbox(xysize/2, xysize/2, certhick/2);
  Box cut_box(xyopen/2, xyopen/2, certhick/2);

  SubtractionSolid ceramic (cerbox, cut_box, Position(0,0, -_HRPPD_BASEPLATE_THICKNESS_));

  Volume ceramicVol(detName + "_ceramic", ceramic, air);
  ceramicVol.setVisAttributes(gasvolVis);

  PlacedVolume ceramicPV = hrppdVol_air.placeVolume(ceramicVol, Position(0.0, 0.0, accu + certhick/2));
  DetElement   ceramicDE(sdet, "ceramic_de", 0);
  ceramicDE.setPlacement(ceramicPV);

  ///*--------------------------------------------------*/ 
  // Plating body

  Box plating_solid(xyopen/2, xyopen/2, _HRPPD_PLATING_LAYER_THICKNESS_/2);
  Volume platingVol( detName + "_plating", plating_solid, air);

  platingVol.setVisAttributes(gasvolVis);
  PlacedVolume platingPV = hrppdVol_air.placeVolume(platingVol, Position(0.0, 0.0, accu + certhick/2));
  DetElement   platingDE(sdet, "plating_de", 0);
  platingDE.setPlacement(platingPV);

  ///*--------------------------------------------------*/ 
  // MCP body

  Box mcp_solid( xyopen/2, xyopen/2, _EFFECTIVE_MCP_THICKNESS_/2);
  Volume mcpVol( detName + "_mcp", mcp_solid, air);

  mcpVol.setVisAttributes(gasvolVis);
  PlacedVolume mcpPV = hrppdVol_air.placeVolume(mcpVol, Position(0.0, 0.0, accu + certhick/2 + 
					 _HRPPD_PLATING_LAYER_THICKNESS_/2 + _EFFECTIVE_MCP_THICKNESS_/2));
  DetElement   mcpDE(sdet, "mcp_de", 0);
  mcpDE.setPlacement(mcpPV);



  ///*--------------------------------------------------*/ 
  // 

  double pdthick = 0.001, zpdc = azOffset + wndthick + pdthick/2;

  Box pdbox_solid(xyactive/2, xyactive/2, pdthick/2);
  Volume pdboxVol( detName + "_pd", pdbox_solid, air);

  pdboxVol.setVisAttributes(gasvolVis);
  PlacedVolume pdboxPV = hrppdVol_air.placeVolume(pdboxVol, Position(0.0, 0.0, accu + pdthick + pdthick/2));

  DetElement   pdboxDE(sdet, "pdbox_de", 0);
  pdboxDE.setPlacement(pdboxPV);

  ///*--------------------------------------------------*/ 

  Box qdbox_solid(xyactive/2, xyactive/2, pdthick/2);
  Volume qdboxVol( detName + "_qd", qdbox_solid, air);

  qdboxVol.setVisAttributes(gasvolVis);
  PlacedVolume qdboxPV = hrppdVol_air.placeVolume(qdboxVol, Position(0.0, 0.0, accu + pdthick/2));

  DetElement   qdboxDE(sdet, "qdbox_de", 0);
  pdboxDE.setPlacement(qdboxPV);

  accu += certhick + 1*mm;


  ///*--------------------------------------------------*/ 
  // PCB Board

  Box pcb_solid (_READOUT_PCB_SIZE_/2, _READOUT_PCB_SIZE_/2, _READOUT_PCB_THICKNESS_/2);
  Volume pcbVol( detName + "_pcb", pcb_solid, air);

  pcbVol.setVisAttributes(gasvolVis);
  PlacedVolume pcbPV = hrppdVol_air.placeVolume(pcbVol, Position(0.0, 0.0, accu + _READOUT_PCB_THICKNESS_ /2));

  DetElement   pcbDE(sdet, "pcb_de", 0);
  pcbDE.setPlacement(pcbPV);
 
  accu += _READOUT_PCB_THICKNESS_ + 0.001;


  ///*--------------------------------------------------*/ 
  // PCB Board

  Box asic_solid( _ASIC_SIZE_XY_/2, _ASIC_SIZE_XY_/2, _ASIC_THICKNESS_/2);
  Volume asicVol( detName + "_asic", asic_solid, mirrorMat);
  asicVol.setVisAttributes(mirrorVis);

  double asic_pitch = _READOUT_PCB_SIZE_/2;

  for(unsigned ix=0; ix<2; ix++) {
	double xOffset = asic_pitch*(ix - (2-1)/2.);
	
	for(unsigned iy=0; iy<2; iy++) {
	  double yOffset = asic_pitch*(iy - (2-1)/2.);
	  
//	  new G4PVPlacement(0, G4ThreeVector(xOffset, yOffset, accu + _ASIC_THICKNESS_/2), asic_log, "ASIC", 
//			    hrppd_log, false, ix*2 + iy);

       auto sensorPV = hrppdVol_air.placeVolume(asicVol,  Position(xOffset, yOffset, accu + _ASIC_THICKNESS_/2));

	} //for iy
  } //for ix

  accu += _ASIC_THICKNESS_ + 0.01*mm;



//  ///*--------------------------------------------------*/ 
// Cold plates and pipes;

  double pcb_pitch = _READOUT_PCB_SIZE_/2;

  Box cplate_solid(_ASIC_SIZE_XY_/2, _ASIC_SIZE_XY_/2, _COLD_PLATE_THICKNESS_/2);
  Volume cplateVol(detName + "_cplate", cplate_solid, mirrorMat);

  // FIXME: yes, have to make the pipes 1.5mm shorter to fit into the HRPPD container volume;
  // must be a minor simplification I guess;
  double cooling_length = _HRPPD_TILE_SIZE_ /*+ _HRPPD_INSTALLATION_GAP_*/, iradius = _COOLING_PIPE_INNER_DIAMETER_/2; 
  double oradius = iradius + _COOLING_PIPE_WALL_THICKNESS_;
  
  Tube pipe_solid(0.0, oradius, cooling_length/2, 0*degree, 360*degree);
  Volume pipeVol(detName + "_pipe", pipe_solid, mirrorMat);

  Tube Water_solid(0.0, oradius, cooling_length/2, 0*degree, 360*degree);
  Volume WaterVol(detName + "_pipe", Water_solid, mirrorMat);

  auto WaterPV = pipeVol.placeVolume(WaterVol,  Position(0.0, 0.0, 0.0));

  Rotation3D rY(RotationZYX(0.0, M_PI/2.0, 0.0));

      for(unsigned bt=0; bt<2; bt++) {
	    double yOffset = pcb_pitch*(bt - (2-1)/2.);

   auto water1PV = hrppdVol_air.placeVolume(pipeVol, Transform3D(rY, Position(0.0, yOffset, accu + _COLD_PLATE_THICKNESS_/2)));
	  
	
//	new G4PVPlacement(rY, G4ThreeVector(0.0, yOffset, accu + _COLD_PLATE_THICKNESS_ + oradius), 
//			  pipe_log, "CoolingPipe", hrppd_log, false, 0);
	
	for(unsigned lr=0; lr<2; lr++) {
	    double xOffset = pcb_pitch*(lr - (2-1)/2.);


        cout << xOffset << endl;

      auto cplatePV = hrppdVol_air.placeVolume(cplateVol, Transform3D(rY, Position(xOffset, yOffset, accu + _COLD_PLATE_THICKNESS_/2)));
	  
//	  new G4PVPlacement(0, G4ThreeVector(xOffset, yOffset, accu + _COLD_PLATE_THICKNESS_/2), cplate_log, 
//			    "ColdPlate", hrppd_log, false, 0);



	} //for lr
      } //for bt
  
  accu += _COLD_PLATE_THICKNESS_ + 2*oradius + 0.01*mm;



//  ///*--------------------------------------------------*/ 
//
//   // sensitivity
//   if (!debug_optics)
//     sensorVol.setSensitiveDetector(sens);
// 
//
//   sensorSide =  _HRPPD_TILE_SIZE_;
//   sensorGap = _HRPPD_INSTALLATION_GAP_;
//
//   double sx, sy;
//   for (double usx = 0; usx <= tBoxMax; usx += sensorSide + sensorGap) {
//     for (int sgnx = 1; sgnx >= (usx > 0 ? -1 : 1); sgnx -= 2) {
//       for (double usy = 0; usy <= tBoxMax; usy += sensorSide + sensorGap) {
//         for (int sgny = 1; sgny >= (usy > 0 ? -1 : 1); sgny -= 2) {
// 
//           // sensor (x,y) center
//           sx = sgnx * usx - sensorSide/2. - sensorGap/2.;
//           sy = sgny * usy;
// 
//           // annular cut
//           if (std::hypot(sx, sy) < sensorPlaneRmin || std::hypot(sx, sy) > sensorPlaneRmax)
//             continue;
// 
//           // placement (note: transformations are in reverse order)
//           auto sensorPlacement = Transform3D(
//               Translation3D(sensorPlanePos.x(), sensorPlanePos.y(), sensorPlanePos.z() + 34)  * // move to reference position
//               Translation3D(sx, sy, 0.)                                                   // move to grid position
//           );
// 
// //		 cout << sensorPlanePos.x() << "  " << sensorPlanePos.y() << "  "<< sensorPlanePos.z() << endl;
// //          auto sensorPV = gasvolVol.placeVolume(sensorVol, sensorPlacement);
////           auto sensorPV = volume.placeVolume(sensorVol, sensorPlacement);
//           auto sensorPV = volume.placeVolume(hrppdVol_air, sensorPlacement);
// 
//           // generate LUT for module number -> sensor position, for readout mapping tests
//           // printf("%d %f %f\n",imod,sensorPV.position().x(),sensorPV.position().y());
// 
//           // properties
//           sensorPV.addPhysVolID("module", imod); // NOTE: must be consistent with `sensorIDfields`
//           auto       imodEnc = encodeSensorID(sensorPV.volIDs());
//           DetElement sensorDE(sdet, "sensor_de_" + std::to_string(imod), imodEnc);
//           sensorDE.setPlacement(sensorPV);
//           if (!debug_optics) {
//             SkinSurface sensorSkin(description, sensorDE, "sensor_optical_surface_" + std::to_string(imod), sensorSurf,
//                                    sensorVol);
//             sensorSkin.isValid();
//           };
// 
//           // increment sensor module number
//           imod++;
//         };
//       };
//     };
//   };
//   // END SENSOR MODULE LOOP ------------------------
 

   ///*--------------------------------------------------*/ 
   /// Loading the coordinates
   ///

    unsigned const hdim = 9;
    const unsigned flags[hdim][hdim] = {
      // NB: WYSIWIG fashion; well, it is top/ bottom and left/right symmetric;
      {0, 0, 1, 1, 1, 1, 1, 0, 0},
      {0, 1, 1, 1, 1, 1, 1, 1, 0},
      {1, 1, 1, 1, 1, 1, 1, 1, 1},
      {1, 1, 1, 1, 2, 1, 1, 1, 1},
      {3, 3, 3, 4, 0, 2, 1, 1, 1},
      {1, 1, 1, 1, 2, 1, 1, 1, 1},
      {1, 1, 1, 1, 1, 1, 1, 1, 1},
      {0, 1, 1, 1, 1, 1, 1, 1, 0},
      {0, 0, 1, 1, 1, 1, 1, 0, 0}
    };

    std::vector<std::pair<TVector2, bool>> coord;

    for(unsigned ix=0; ix<hdim; ix++) {
      double xOffset = (_HRPPD_TILE_SIZE_ + _HRPPD_INSTALLATION_GAP_)*(ix - (hdim-1)/2.); 
      
      for(unsigned iy=0; iy<hdim; iy++) {
	       double yOffset = (_HRPPD_TILE_SIZE_ + _HRPPD_INSTALLATION_GAP_)*(iy - (hdim-1)/2.); 
           unsigned flag = flags[hdim-iy-1][ix];

           if (!flag) continue;

  	       double qxOffset = xOffset + (flag >= 3 ? -_HRPPD_CENTRAL_ROW_OFFSET_ : 0.0);
	       coord.push_back(std::make_pair(TVector2(qxOffset, yOffset), flag%2));
      } //for iy
    } //for ix

//    printf("%lu sensors total\n", coord.size());

   ///*--------------------------------------------------*/ 
   /// Set sensors into the coordinates
   ///
    for(auto xyptr: coord) {
	   auto &xy = xyptr.first;

           double sx = xy.X();
           double sy = xy.Y();
 
           // placement (note: transformations are in reverse order)
           auto sensorPlacement = Transform3D(
               Translation3D(sensorPlanePos.x(), sensorPlanePos.y(), sensorPlanePos.z() + 34)  * // move to reference position
               Translation3D(sx, sy, 0.)                                                   // move to grid position
           );
 
           auto sensorPV = volume.placeVolume(hrppdVol_air, sensorPlacement);
 
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


      }
 
 ///*--------------------------------------------------*/ 
 ///*--------------------------------------------------*/ 
 
   cout << "aaaaaaaaaaaaaaaaaaaaaaaa " << "    " << vesselRmax0 << endl;
 
   return sdet;

}

// clang-format off
DECLARE_DETELEMENT(epic_PFRICH_v1, createDetector)
