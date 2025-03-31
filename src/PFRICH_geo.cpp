// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2023  Wenliang (Bill) Li, Alexander Kiselev, Karthik Suresh

//----------------------------------
// pfRICH: Proximity Focusing RICH
// Author: Wenliang (Bill) Li
//
// - Design Adapted from standalone Geant4 description by
//   Alexander Kiselev and Chandradoy Chatterjee
//----------------------------------

#define _WITH_OPTICS_

#include "DD4hep/DetFactoryHelper.h"
//#include "DD4hep/OpticalSurfaces.h"
//#include "DD4hep/Printout.h"
//#include "DDRec/DetectorData.h"
//#include "DDRec/Surface.h"

//#include <XML/Helper.h>
//#include "XML/Layering.h"
//#include "TVector3.h"
#include "TFile.h"
#include "TVector2.h"

//#include "TGeoElement.h"
//#include "TGeoManager.h"
//#include "TInterpreter.h"
//#include "TUri.h"

//using namespace std;
using namespace dd4hep;
//using namespace dd4hep::rec;
//using namespace dd4hep::detail;

#ifdef _WITH_OPTICS_
#include "IRT/CherenkovDetectorCollection.h"
#include "IRT/CylindricalSurface.h"
#endif

//#define _HRPPD_PITCH_                          (3.25*mm)
//#define _HRPPD_ACTIVE_AREA_SIZE_      (32*_HRPPD_PITCH_)

#define _SENSOR_PLANE_GEOMETRIC_EFFICIENCY_       (1.00)
#define _SAFETY_FACTOR_                           (0.70)

// -------------------------------------------------------------------------------------

static UnionSolid FlangeCut(Detector& description, double length, double clearance)
{
  // FIXME: not the most efficient way to recalculate them every time new,
  // but code is more readable, and overhead is negligible anyway;
  auto _FLANGE_EPIPE_DIAMETER_ = description.constant<double>("FLANGE_EPIPE_DIAMETER");
  auto _FLANGE_HPIPE_DIAMETER_ = description.constant<double>("FLANGE_HPIPE_DIAMETER");
  auto _FLANGE_HPIPE_OFFSET_   = description.constant<double>("FLANGE_HPIPE_OFFSET");
  
  // A wedge bridging two cylinders;
  Tube eflange(0.0, _FLANGE_EPIPE_DIAMETER_ / 2 + clearance, length/2);//25);
  Tube hflange(0.0, _FLANGE_HPIPE_DIAMETER_ / 2 + clearance, length/2);//25);

  double r0 = _FLANGE_EPIPE_DIAMETER_ / 2 + clearance;
  double r1 = _FLANGE_HPIPE_DIAMETER_ / 2 + clearance;
  double L  = _FLANGE_HPIPE_OFFSET_;
  double a  = r0 * L / (r0 - r1);
  double b  = r0 * r0 / a;
  double c  = r1 * (a - b) / r0;

  // GEANT variables to define G4Trap;
  double pDz = /*25*/length/2, pTheta = 0.0, pPhi = 0.0, pDy1 = (a - b - c) / 2, pDy2 = pDy1;
  double pDx1 = sqrt(r0 * r0 - b * b), pDx2 = pDx1 * r1 / r0, pDx3 = pDx1, pDx4 = pDx2, pAlp1 = 0.0,
         pAlp2 = 0.0;

  Trap wedge(pDz, pTheta, pPhi, pDy1, pDx1, pDx2, pAlp1, pDy2, pDx3, pDx4, pAlp2);

  UnionSolid flange_shape(eflange, hflange, Position(-_FLANGE_HPIPE_OFFSET_, 0.0, 0.0));
  Rotation3D rZ(RotationZYX(M_PI / 2.0, 0.0, 0.0));
  Transform3D transform_flange(rZ, Position(-b - pDy1, 0.0, 0.0));
  
  return UnionSolid(flange_shape, wedge, transform_flange);
} // FlangeCut()

// -------------------------------------------------------------------------------------

static Ref_t createDetector(Detector& description, xml_h e, SensitiveDetector sens)
{
  xml::DetElement detElem = e;
  int det_id      = detElem.id();

  std::string det_name = detElem.nameStr();
  Material air    = description.air();

  DetElement sdet(det_name, det_id);

  sens.setType("tracker");
  description.invisible();

  std::string detName     = detElem.nameStr();
  //+xml::Component dims     = detElem.dimensions();

  int id = detElem.hasAttr(_U(id)) ? /*x_det*/detElem.id() : 0;

  // Start optical configuration if needed;
#ifdef _WITH_OPTICS_
  // Output file with optics description;
  auto fname = detElem.attr<std::string>(_Unicode(optics));
  
  auto output_file = TFile::Open(fname.c_str(), "RECREATE");
  auto geometry = new CherenkovDetectorCollection();
  auto cdet = geometry->AddNewDetector(detName.c_str());
#endif

  OpticalSurfaceManager surfMgr = description.surfaceManager();

  // - sensor module
  auto sensorElem        = detElem.child(_Unicode(sensors)).child(_Unicode(module));
#if _TODAY_
  auto sensorMat         = description.material(sensorElem.attr<std::string>(_Unicode(material)));
  auto sensorVis         = description.visAttributes(sensorElem.attr<std::string>(_Unicode(vis)));
  auto sensorSurf        = surfMgr.opticalSurface(sensorElem.attr<std::string>(_Unicode(surface)));
  double sensorSide      = sensorElem.attr<double>(_Unicode(side));
  double sensorThickness = sensorElem.attr<double>(_Unicode(thickness));
#endif
  auto readoutName       = detElem.attr<std::string>(_Unicode(readout));
  auto HRPPD_WindowMat = description.material(sensorElem.attr<std::string>(_Unicode(windowmat)));
  auto HRPPD_pcMat     = description.material(sensorElem.attr<std::string>(_Unicode(photocathode)));
  auto HRPPD_MPDMat    = description.material(sensorElem.attr<std::string>(_Unicode(mpdmat)));
  auto HRPPD_PCBMat    = description.material(sensorElem.attr<std::string>(_Unicode(pcbmat)));
  auto HRPPD_ASICMat   = description.material(sensorElem.attr<std::string>(_Unicode(asicmat)));

#if _TODAY_
  double vesselRmin0 = dims.attr<double>(_Unicode(rmin0));
  double vesselRmin1 = dims.attr<double>(_Unicode(rmin1));
  double vesselRmax0 = dims.attr<double>(_Unicode(rmax0));
  double vesselRmax1 = dims.attr<double>(_Unicode(rmax1));

  int imod = 0; // module number

  auto gasvolMat = description.material(detElem.attr<std::string>(_Unicode(gas)));
#endif
  auto gasvolVis = description.visAttributes("DRICH_gas_vis");
  auto sensorVis = description.visAttributes("PFRICH_sensor_vis");//DRICH_gas_vis");
#if _TODAY_
  auto vesselVis = description.visAttributes(detElem.attr<std::string>(_Unicode(vis_vessel)));

  double windowThickness = dims.attr<double>(_Unicode(window_thickness));
  double wallThickness   = dims.attr<double>(_Unicode(wall_thickness));

  double proximityGap = dims.attr<double>(_Unicode(proximity_gap));

  long debug_optics_mode = description.constantAsLong("PFRICH_debug_optics");

  bool debug_optics = debug_optics_mode > 0;

  double radiatorRmin = radiatorElem.attr<double>(_Unicode(rmin));
  double radiatorRmax = radiatorElem.attr<double>(_Unicode(rmax));

  auto filterElem = radiatorElem.child(_Unicode(filter));

  double airgapThickness = 0.1;
  double filterThickness = 1;
#endif

  auto radiatorElem         = detElem.child(_Unicode(radiator));
  //?double radiatorFrontplane = radiatorElem.attr<double>(_Unicode(frontplane));
  auto aerogelElem        = radiatorElem.child(_Unicode(aerogel));
  auto aerogelThickness = aerogelElem.attr<double>(_Unicode(thickness));
  auto aerogelMat = description.material(aerogelElem.attr<std::string>(_Unicode(material)));
  
#if _TODAY_
  auto filterMat  = description.material(filterElem.attr<std::string>(_Unicode(material)));

  double vesselLength = dims.attr<double>(_Unicode(length));
  auto originFront    = Position(0., 0., vesselLength / 2.0);
  double sensorZpos = radiatorFrontplane - aerogelThickness - proximityGap - 0.5 * sensorThickness;
  auto sensorPlanePos = Position(0., 0., sensorZpos) + originFront; // reference position
#endif
  
  // readout coder <-> unique sensor ID
  /* - `sensorIDfields` is a list of readout fields used to specify a unique sensor ID
     * - `cellMask` is defined such that a hit's `cellID & cellMask` is the corresponding sensor's unique ID
     * - this redundant generalization is for future flexibility, and consistency with dRICH
     */
  std::vector<std::string> sensorIDfields = {"hrppd"};//module"};
  const auto& readoutCoder                = *description.readout(readoutName).idSpec().decoder();
  // determine `cellMask` based on `sensorIDfields`
  uint64_t cellMask = 0;
  for (const auto& idField : sensorIDfields)
    cellMask |= readoutCoder[idField].mask();
  description.add(Constant("PFRICH_cell_mask", std::to_string(cellMask)));
#ifdef _WITH_OPTICS_
  // Do not mind to store it twice;
  cdet->SetReadoutCellMask(cellMask);
#endif
  // create a unique sensor ID from a sensor's PlacedVolume::volIDs
  auto encodeSensorID = [&readoutCoder](auto ids) {
    uint64_t enc = 0;
    for (const auto& [idField, idValue] : ids)
      enc |= uint64_t(idValue) << readoutCoder[idField].offset();
    return enc;
  };

#ifdef _WITH_OPTICS_
  const bool flip = true;
  const double sign = flip ? -1.0 : 1.0;
#endif
  
  auto mirrorElem = detElem.child(_Unicode(mirror));
  auto mirrorMat  = description.material(mirrorElem.attr<std::string>(_Unicode(material)));
  auto mirrorVis  = description.visAttributes(mirrorElem.attr<std::string>(_Unicode(vis)));
  auto asicVis  = description.visAttributes("PFRICH_service_vis");//mirrorElem.attr<std::string>(_Unicode(vis)));
  // auto vesselMat = description.material(detElem.attr<std::string>(_Unicode(material)));
  auto vesselGas = description.material(detElem.attr<std::string>(_Unicode(gas)));

  //
  // Fiducial volume; will be split into three subvolumes along the beam line: front wall, gas volume
  // and sensor compartment; FIXME: called 'air' but presently filled with nitrogen;
  //
  auto _FIDUCIAL_VOLUME_LENGTH_ = description.constant<double>("FIDUCIAL_VOLUME_LENGTH");
  auto _VESSEL_OUTER_RADIUS_    = description.constant<double>("VESSEL_OUTER_RADIUS");
  
  Tube pfRICH_air_volume(0.0, _VESSEL_OUTER_RADIUS_, _FIDUCIAL_VOLUME_LENGTH_/2);

  auto _FLANGE_CLEARANCE_   = description.constant<double>("FLANGE_CLEARANCE");
  SubtractionSolid pfRICH_volume_shape(pfRICH_air_volume,
				       FlangeCut(description, _FIDUCIAL_VOLUME_LENGTH_ + 1*mm, _FLANGE_CLEARANCE_));
  Volume pfRICH_volume(detName + "_Vol", pfRICH_volume_shape, vesselGas); 
  
  Volume mother = description.pickMotherVolume(sdet);
  auto _FIDUCIAL_VOLUME_OFFSET_ = description.constant<double>("FIDUCIAL_VOLUME_OFFSET");
  Transform3D transform(RotationZYX(0, M_PI, 0), Position(0, 0, _FIDUCIAL_VOLUME_OFFSET_));
  PlacedVolume pv = mother.placeVolume(pfRICH_volume, transform);
  if (id != 0) pv.addPhysVolID("system", id);
  sdet.setPlacement(pv);

  double fvOffset = _FIDUCIAL_VOLUME_OFFSET_;
  
  //
  // Gas volume;
  //
  auto _VESSEL_FRONT_SIDE_THICKNESS_ = description.constant<double>("VESSEL_FRONT_SIDE_THICKNESS");
  auto _SENSOR_AREA_LENGTH_          = description.constant<double>("SENSOR_AREA_LENGTH");
  auto _VESSEL_OUTER_WALL_THICKNESS_ = description.constant<double>("VESSEL_OUTER_WALL_THICKNESS");
  
  double gas_volume_length =  _FIDUCIAL_VOLUME_LENGTH_ - _VESSEL_FRONT_SIDE_THICKNESS_ - _SENSOR_AREA_LENGTH_;
  double gas_volume_radius = _VESSEL_OUTER_RADIUS_ - _VESSEL_OUTER_WALL_THICKNESS_;
  double gas_volume_offset = -(_SENSOR_AREA_LENGTH_ - _VESSEL_FRONT_SIDE_THICKNESS_)/2, gvOffset = gas_volume_offset; 
  Tube gasTube(0.0, gas_volume_radius, gas_volume_length/2);
  SubtractionSolid gasSolid(gasTube, FlangeCut(description, gas_volume_length + 1*mm, _FLANGE_CLEARANCE_));
  Volume gasVolume(detName + "_GasVol", gasSolid, vesselGas);
  /*PlacedVolume gasPV =*/ pfRICH_volume.placeVolume(gasVolume, Position(0, 0, gas_volume_offset));
#ifdef _WITH_OPTICS_
  {
    // FIXME: Z-location does not really matter here, right?;
    auto boundary = new FlatSurface(TVector3(0,0,0), sign*TVector3(1,0,0), TVector3(0,-1,0));
    
    auto radiator = geometry->SetContainerVolume(cdet, "GasVolume", 0, (G4LogicalVolume*)(0x0), 0, boundary);
    radiator->SetAlternativeMaterialName("N2cherenkov");//gasvolMaterialName.c_str());
  }
#endif

  //
  // Vessel walls -> cut'n'paste from pfRICH standalone code;;
  //

  auto _BUILDING_BLOCK_CLEARANCE_ = description.constant<double>("BUILDING_BLOCK_CLEARANCE");
  auto _VESSEL_INNER_WALL_THICKNESS_ = description.constant<double>("VESSEL_INNER_WALL_THICKNESS");
  
  // To be used in boolean operations in several places;
  auto flange = FlangeCut(description, gas_volume_length + 1*mm, _FLANGE_CLEARANCE_ + _VESSEL_INNER_WALL_THICKNESS_ + 
			  _BUILDING_BLOCK_CLEARANCE_);
  
  double gzOffset = -gas_volume_length / 2 + _BUILDING_BLOCK_CLEARANCE_;// + agthick / 2;
    
  auto _FLANGE_EPIPE_DIAMETER_ = description.constant<double>("FLANGE_EPIPE_DIAMETER");
  float m_r0min = _FLANGE_EPIPE_DIAMETER_ / 2 + _FLANGE_CLEARANCE_ + _VESSEL_INNER_WALL_THICKNESS_ +
    _BUILDING_BLOCK_CLEARANCE_;
  float m_r0max = gas_volume_radius - _BUILDING_BLOCK_CLEARANCE_;
    
  //
  // Aerogel;
  //
  {
    const int _AEROGEL_BAND_COUNT_ = 3;
    auto _AEROGEL_INNER_WALL_THICKNESS_ =
      description.constant<double>("AEROGEL_INNER_WALL_THICKNESS");
    //float _VESSEL_OUTER_RADIUS_         = description.constant<double>("VESSEL_OUTER_RADIUS");
    //float _FLANGE_CLEARANCE_         = description.constant<double>("FLANGE_CLEARANCE");
    //float _AEROGEL_BAND_COUNT_ = aerogel_band_count;
    auto _AEROGEL_SEPARATOR_WALL_THICKNESS_ =
      description.constant<double>("AEROGEL_SEPARATOR_WALL_THICKNESS");
    auto _AEROGEL_OUTER_WALL_THICKNESS_ =
      description.constant<double>("AEROGEL_OUTER_WALL_THICKNESS");
    
    const unsigned adim[_AEROGEL_BAND_COUNT_] = {9, 14, 20};
    double rheight =
      (m_r0max - m_r0min - (_AEROGEL_BAND_COUNT_ - 1) * _AEROGEL_SEPARATOR_WALL_THICKNESS_ -
       _AEROGEL_INNER_WALL_THICKNESS_ - _AEROGEL_OUTER_WALL_THICKNESS_) /
      _AEROGEL_BAND_COUNT_;
    
    double agthick = aerogelThickness;//description.constant<double>("PFRICH_aerogel_thickness");//2.5; // cm
    
    std::string aerogel_name = "a1040";
    int kkcounter = 0;

    for (unsigned ir = 0; ir < _AEROGEL_BAND_COUNT_; ir++) {
      int counter       = ir ? -1 : 0;
      double apitch     = 360 * degree / adim[ir];
      double aerogel_r0 = m_r0min + _AEROGEL_INNER_WALL_THICKNESS_ +
	ir * (_AEROGEL_SEPARATOR_WALL_THICKNESS_ + rheight);
      double aerogel_r1 = aerogel_r0 + rheight;
      double rm         = (aerogel_r0 + aerogel_r1) / 2;
      
      // Calculate angular space occupied by the spacers and by the tiles; no gas gaps for now;
      // assume that a wegde shape is good enough (GEANT visualization does not like boolean objects),
      // rather than creating constant thicjkess azimuthal spacers; just assume that spacer thickness is
      // _AEROGEL_FRAME_WALL_THICKNESS_ at r=rm;
      double l0   = 2 * M_PI * rm / adim[ir];
      double l1   = _AEROGEL_SEPARATOR_WALL_THICKNESS_;
      double lsum = l0 + l1;
      
      // FIXME: names overlap in several places!;
      double wd0      = (l0 / lsum) * (360 * degree / adim[ir]);
      //+double wd1      = (l1 / lsum) * (360 * degree / adim[ir]);
      TString ag_name = "Tmp";//+, sp_name = "Tmp";
      
      if (ir) ag_name.Form("%s-%d-00", aerogel_name.c_str(), ir);
      //+if (ir) sp_name.Form("A-Spacer--%d-00", ir);
      
      Tube agtube(aerogel_r0, aerogel_r1, agthick / 2, 0 * degree, wd0);
      //+Tube sptube(aerogel_r0, aerogel_r1, agthick / 2, wd0, wd0 + wd1);
      
      for (unsigned ia = 0; ia < adim[ir]; ia++) {
	
	Rotation3D r_aerogel_Z(RotationZYX(ia * apitch, 0.0, 0.0));
	Rotation3D r_aerogel_Zinv(RotationZYX(-1. * ia * apitch, 0.0, 0.0));
	
	if (ir) {
	  ag_name.Form("%s-%d-%02d", "aerogel", ir, ia);
	  
	  Volume agtubeVol(ag_name.Data(), agtube, aerogelMat);//gasvolMat);
	  auto aerogelTilePlacement = Transform3D(r_aerogel_Z, Position(0.0, 0.0, gzOffset + agthick/2));//-m_gzOffset));
	  auto aerogelTilePV        = /*pfRICH_volume*/gasVolume.placeVolume(agtubeVol, aerogelTilePlacement);
	  DetElement aerogelDE(sdet, "aerogel_de_" + std::to_string(kkcounter), 0);
	  aerogelDE.setPlacement(aerogelTilePV);
#if _LATER_
	  Volume sptubeVol(detName + "_sptube", sptube, mirrorMat);
	  auto sptubePlacement = Transform3D(r_aerogel_Z, Position(0.0, 0.0, gzOffset + agthick/2));//-m_gzOffset));
	  auto sptubePV        = . pfRICH_volume.placeVolume(sptubeVol, sptubePlacement);
	  DetElement sptubeDE(sdet, "sptube_de_" + std::to_string(kkcounter), 0);
	  sptubeDE.setPlacement(sptubePV);
#endif
	} else {
	  
	  ag_name.Form("%s-%d-%02d", "aerogel_inner", ir, ia);
	  
	  Tube agtube_inner(aerogel_r0, aerogel_r1, agthick / 2, 0 * degree + ia * apitch,
			    wd0 + ia * apitch);
	  SubtractionSolid agsub(agtube_inner, flange);//_final_shape);
	  Volume agsubtubeVol(ag_name.Data(), agsub, aerogelMat);//gasvolMat);
	  auto aerogelTilePlacement =
            Transform3D(RotationZYX(0.0, 0.0, 0.0), Position(0.0, 0.0, gzOffset + agthick/2));//-m_gzOffset));
	  auto agsubTilePV = /*pfRICH_volume*/gasVolume.placeVolume(agsubtubeVol, aerogelTilePlacement);
	  
	  DetElement aerogelDE(sdet, "agsubTile_de_" + std::to_string(counter), 0);
	  aerogelDE.setPlacement(agsubTilePV);
#if _LATER_
	  sp_name.Form("%s-%d-%02d", "sp_inner", ir, ia);
	  Tube sptube_inner(aerogel_r0, aerogel_r1, agthick / 2, wd0 + ia * apitch,
			    wd0 + wd1 + ia * apitch);
	  
	  SubtractionSolid spsub(sptube_inner, flange_final_shape);
	  Volume spsubtubeVol(sp_name.Data(), spsub, mirrorMat);
	  auto spTilePlacement =
            Transform3D(RotationZYX(0.0, 0.0, 0.0), Position(0.0, 0.0, -m_gzOffset));
	  auto spsubTilePV = . pfRICH_volume.placeVolume(spsubtubeVol, spTilePlacement);
	  
	  DetElement sptubeTileDE(sdet, "sptubeTile_de_" + std::to_string(counter), 0);
	  sptubeTileDE.setPlacement(spsubTilePV);
	  
	  sp_name.Form("A-Spacer--%d-%02d", ir, ia);
#endif
	} //if
	
	counter++;
	kkcounter++;
	
      } //for ia
    } // for ir
  
    // Placing radial spacer
#if _LATER_
    double sp_accu   = m_r0min;
    int tube_counter = 0;
    
    for (unsigned ir = 0; ir < _AEROGEL_BAND_COUNT_ + 1; ir++) {
      double thickness = ir ? (ir == _AEROGEL_BAND_COUNT_ ? _AEROGEL_OUTER_WALL_THICKNESS_
			       : _AEROGEL_SEPARATOR_WALL_THICKNESS_)
	: _AEROGEL_INNER_WALL_THICKNESS_;
      double sp_r0     = sp_accu;
      double sp_r1     = sp_r0 + thickness;
      
      TString sp_name = "Tmp";
      if (ir)
	sp_name.Form("R-Spacer--%d-00", ir);
      
      Tube sptube(sp_r0, sp_r1, agthick / 2, 0 * degree, 360 * degree);
      Volume sptubeVol(detName + "_radial_sptube", sptube, sensorMat);
      auto sptubePlacement = Transform3D(RotationZYX(0.0, 0.0, 0.0), Position(0.0, 0.0, -m_gzOffset));
      
      if (ir) {
	
	auto sptubePV = . pfRICH_volume.placeVolume(sptubeVol, sptubePlacement);
	
	DetElement sptubeDE(sdet, "sptube_de_" + std::to_string(tube_counter), 0);
	sptubeDE.setPlacement(sptubePV);
	
      }
      
      else {
	
	SubtractionSolid spsub(sptube, flange_final_shape);
	Volume agsubtubeVol(detName + "_radial_sptube_inner", spsub, gasvolMat);
	
	auto sptubePV = . pfRICH_volume.placeVolume(agsubtubeVol, sptubePlacement);
	
	DetElement sptubeDE(sdet, "sptube_de_" + std::to_string(tube_counter), 0);
	sptubeDE.setPlacement(sptubePV);
	
      } //if
      
      sp_accu += thickness + rheight;
      
      tube_counter++;
      
    } //for ir
#endif
    
#ifdef _WITH_OPTICS_
    {
      TVector3 nx(1*sign,0,0), ny(0,-1,0);
      
      auto surface = new FlatSurface(sign*(1/mm)*TVector3(0,0,fvOffset + gvOffset + gzOffset), nx, ny);
      
      auto radiator = geometry->AddFlatRadiator(cdet, "Aerogel", CherenkovDetector::Upstream, 
						0, (G4LogicalVolume*)(0x1), 0, surface, agthick/mm);//agThick/mm);
      radiator->SetAlternativeMaterialName("Aerogel_QRICH");//aerogelMaterialName.c_str());
      // FIXME: what is it good for in ePIC IRT 2.0 implementation?;
      geometry->AddRadiatorLogicalVolume(radiator, (G4LogicalVolume*)(0x1));
    }
#endif
	
    gzOffset += agthick;
  }

  
  //
  // Acrylic filter;
  //
  // FIXME: XML;
  double _ACRYLIC_THICKNESS_ = 3*mm;
  {
    double acthick = _ACRYLIC_THICKNESS_;// =*/ 3*mm;
    auto acrylicMat = description.material("Acrylic_QRICH");
    // m_gzOffset += acthick/2;
    
    Tube ac_tube(m_r0min + 3, m_r0max - 1, acthick / 2, 0 * degree, 360 * degree);
    SubtractionSolid ac_shape(ac_tube, flange);//_final_shape);
    Volume acVol(detName + "_ac", ac_shape, acrylicMat);
    
    //PlacedVolume ac_PV =
    gasVolume.placeVolume(acVol, Position(0, 0, gzOffset + acthick/2));//-21.3));
    
    //DetElement acDE(sdet, "ac_de", 0);
    //acDE.setPlacement(ac_PV);
    
#ifdef _WITH_OPTICS_
    {
      TVector3 nx(1*sign,0,0), ny(0,-1,0);
      
      auto surface = new FlatSurface(sign*(1/mm)*TVector3(0,0,fvOffset + gvOffset + gzOffset), nx, ny);
      
      auto radiator = geometry->AddFlatRadiator(cdet, "Acrylic", CherenkovDetector::Upstream, 
						0, (G4LogicalVolume*)(0x2), 0, surface, acthick/mm);//acThick/mm);
      radiator->SetAlternativeMaterialName("Acrylic_QRICH");//acrylicMaterialName.c_str());
    }
#endif
  }

  //
  // Mirrors;
  //
  {
    auto mirrorSurf          = surfMgr.opticalSurface("MirrorSurface_QRICH");//detElem.attr<std::string>(_Unicode(mirror)));
  
    auto _CONICAL_MIRROR_INNER_RADIUS_ = description.constant<double>("CONICAL_MIRROR_INNER_RADIUS");
    auto _CONICAL_MIRROR_OUTER_RADIUS_ = description.constant<double>("CONICAL_MIRROR_OUTER_RADIUS");
    auto _INNER_MIRROR_THICKNESS_      = description.constant<double>("INNER_MIRROR_THICKNESS");
    auto _OUTER_MIRROR_THICKNESS_      = description.constant<double>("OUTER_MIRROR_THICKNESS");
  
    double mlen = gas_volume_length - 4*_BUILDING_BLOCK_CLEARANCE_ - aerogelThickness - _ACRYLIC_THICKNESS_;
    //mlen -= _BUILDING_BLOCK_CLEARANCE_;// + _HRPPD_SUPPORT_GRID_BAR_HEIGHT_;
    //mlen -= _BUILDING_BLOCK_CLEARANCE_;// + _HRPPD_SUPPORT_GRID_BAR_HEIGHT_;
    double mzoffset = (aerogelThickness + _ACRYLIC_THICKNESS_ + 2*_BUILDING_BLOCK_CLEARANCE_ )/2;
    
    double mirror_r0[2] = {m_r0min, m_r0max - _BUILDING_BLOCK_CLEARANCE_};
    double mirror_r1[2] = {_CONICAL_MIRROR_INNER_RADIUS_, _CONICAL_MIRROR_OUTER_RADIUS_};
    
    for (unsigned im = 0; im < 2; im++) {
      
      double mirror_thickness = im ? _OUTER_MIRROR_THICKNESS_ : _INNER_MIRROR_THICKNESS_;
      
      if (im) {
	
	Cone mirror_outer_cone_shape(mlen / 2.0, mirror_r0[im], mirror_r0[im] + mirror_thickness,
				     mirror_r1[im], mirror_r1[im] + mirror_thickness);
	
	Volume outer_mirrorVol(detName + "_outer_mirror", mirror_outer_cone_shape, mirrorMat);
	outer_mirrorVol.setVisAttributes(mirrorVis);
	
	PlacedVolume mirror_outerPV = gasVolume.placeVolume(outer_mirrorVol, Position(0, 0, mzoffset));
	DetElement mirror_outerDE(sdet, "_outer_mirror_de", 0);
	mirror_outerDE.setPlacement(mirror_outerPV);
	SkinSurface mirrorSkin(description, mirror_outerDE, "outer_mirror_optical_surface_"/* + secName*/, mirrorSurf,
			       outer_mirrorVol);
	mirrorSkin.isValid();
	
      } else {
	
	Cone mirror_inner_cone_shape(mlen / 2., mirror_r0[im], mirror_r0[im] + mirror_thickness,
				     mirror_r1[im], mirror_r1[im] + mirror_thickness);
	
	SubtractionSolid mirror_inner_sub(mirror_inner_cone_shape, flange);//_final_shape);
	
	Volume inner_mirrorVol(detName + "_inner_mirror", mirror_inner_sub, mirrorMat);
	inner_mirrorVol.setVisAttributes(mirrorVis);
	
	PlacedVolume mirror_innerPV = gasVolume.placeVolume(inner_mirrorVol, Position(0, 0, mzoffset));
	DetElement mirror_innerDE(sdet, "_inner_mirror_de", 0);
	mirror_innerDE.setPlacement(mirror_innerPV);
	SkinSurface mirrorSkin(description, mirror_innerDE, "inner_mirror_optical_surface_"/* + secName*/, mirrorSurf,
			       inner_mirrorVol);
	mirrorSkin.isValid();
      }
    } //for im
  }
      
  //
  // HRPPDs;
  //
  {
    auto _HRPPD_CENTRAL_ROW_OFFSET_ = description.constant<double>("HRPPD_CENTRAL_ROW_OFFSET");
    auto _HRPPD_WINDOW_THICKNESS_   = description.constant<double>("HRPPD_WINDOW_THICKNESS");
    auto _HRPPD_CONTAINER_VOLUME_HEIGHT_ =
      description.constant<double>("HRPPD_CONTAINER_VOLUME_HEIGHT");
    auto _HRPPD_INSTALLATION_GAP_ = description.constant<double>("HRPPD_INSTALLATION_GAP");
    
    //auto _HRPPD_SUPPORT_GRID_BAR_HEIGHT_ =
    //description.constant<double>("HRPPD_SUPPORT_GRID_BAR_HEIGHT");
    
    auto _HRPPD_TILE_SIZE_        = description.constant<double>("HRPPD_TILE_SIZE");
    auto _HRPPD_OPEN_AREA_SIZE_   = description.constant<double>("HRPPD_OPEN_AREA_SIZE");
    auto _HRPPD_ACTIVE_AREA_SIZE_ = description.constant<double>("HRPPD_ACTIVE_AREA_SIZE");
    auto _HRPPD_CERAMIC_BODY_THICKNESS_ =
      description.constant<double>("HRPPD_CERAMIC_BODY_THICKNESS");
    auto _HRPPD_BASEPLATE_THICKNESS_ = description.constant<double>("HRPPD_BASEPLATE_THICKNESS");
    auto _HRPPD_PLATING_LAYER_THICKNESS_ =
      description.constant<double>("HRPPD_PLATING_LAYER_THICKNESS");
    auto _EFFECTIVE_MCP_THICKNESS_ = description.constant<double>("EFFECTIVE_MCP_THICKNESS");
        
    double xysize = _HRPPD_TILE_SIZE_, wndthick = _HRPPD_WINDOW_THICKNESS_;

#ifdef _WITH_OPTICS_
    // [0,0]: have neither access to G4VSolid nor to G4Material; IRT code does not care; fine;
    auto pd = new CherenkovPhotonDetector(0, 0);
    
    // FIXME: '0' stands for the unknown (and irrelevant) G4LogicalVolume;
    geometry->AddPhotonDetector(cdet, 0, pd);
    
    // FIXME: really needed in EICrecon?;
    //pd->SetCopyIdentifierLevel(1);
    //?pd->DefineLogicalVolume();
    // Cannot access GEANT shapes in the reconstruction code -> store this value;
    pd->SetActiveAreaSize(_HRPPD_ACTIVE_AREA_SIZE_/mm);
    
    pd->SetGeometricEfficiency(_SENSOR_PLANE_GEOMETRIC_EFFICIENCY_ * _SAFETY_FACTOR_);
#endif

    
    // HRPPD assembly container volume;
    double hrppd_container_volume_thickness = _HRPPD_CONTAINER_VOLUME_HEIGHT_;

    // For now assume it is a unique surface, same for all HRPPDs;
#ifdef _WITH_OPTICS_
    {	
      TVector3 nx(1*sign,0,0), ny(0,-1,0);
      
      auto surface = 
	new FlatSurface(sign*(1/mm)*TVector3(0,0,/*fvzOffset + wzOffset + _HRPPD_WINDOW_THICKNESS_/2*/-1700.*mm + 2.5*mm), nx, ny);
      
      auto radiator = geometry->AddFlatRadiator(cdet, "QuartzWindow", CherenkovDetector::Downstream, 
						0, (G4LogicalVolume*)(0x3), 0, /*wndVol, m_FusedSilica,*/ surface,
						_HRPPD_WINDOW_THICKNESS_/mm);
      radiator->SetAlternativeMaterialName("AirOptical");//windowMaterialName.c_str());
    }	
#endif
    
    // Create a vector of HRPPD XY-coordinates in a separate loop (it will be used place HRPPDs
    // but also create cutouts in the sensor plane later, etc);
    std::vector<std::pair<TVector2, bool>> coord;
    {
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
	{0, 0, 1, 1, 1, 1, 1, 0, 0}};
      
      for (unsigned ix = 0; ix < hdim; ix++) {
	double xOffset = (_HRPPD_TILE_SIZE_ + _HRPPD_INSTALLATION_GAP_) * (ix - (hdim - 1) / 2.);
	
	for (unsigned iy = 0; iy < hdim; iy++) {
	  double yOffset = (_HRPPD_TILE_SIZE_ + _HRPPD_INSTALLATION_GAP_) * (iy - (hdim - 1) / 2.);
	  unsigned flag  = flags[hdim - iy - 1][ix];
	  
	  if (!flag) continue;
	  
	  double qxOffset = xOffset + (flag >= 3 ? -_HRPPD_CENTRAL_ROW_OFFSET_ : 0.0);
	  coord.push_back(std::make_pair(TVector2(qxOffset, yOffset), flag % 2));
	} //for iy
      } //for ix
    }

    unsigned imod = 0;

    //
    // It looks like each HRPPD should be a new object rather than a copy of the same volume,
    // because otherwise one cannot make photocathodes sensitive (this is done on a PlacedVolume
    // level);
    //
    for (auto xyptr : coord) {
      auto& xy = xyptr.first;

      double pdthick = 0.01*mm;

      uint64_t sensorID = 0x0;
	
      // HRPPD container volume;
      //
      Box hrppd_Solid(xysize / 2, xysize / 2, hrppd_container_volume_thickness / 2);
      TString hrppdName; hrppdName.Form("%s-hrppd-%02d", detName.c_str(), imod);
      // FIXME: may want to use AirOptical here, but then return a dummy absorber
      // layer behind the actual phtocathode;
      Volume hrppdVol_air(hrppdName.Data(), hrppd_Solid, air);
      hrppdVol_air.setVisAttributes(gasvolVis);

      // A running variable to pack layers one after the other one;
      double accu = -hrppd_container_volume_thickness / 2;
      
      // Quartz Window;
      //
      Box wnd_Solid(xysize / 2, xysize / 2, wndthick / 2);
      TString wndName; wndName.Form("%s-window-%02d", detName.c_str(), imod);
      Volume wndVol(wndName.Data(), wnd_Solid, HRPPD_WindowMat);
      wndVol.setVisAttributes(gasvolVis);     
      //PlacedVolume wndPV =
      hrppdVol_air.placeVolume(wndVol, Position(0, 0, accu + wndthick / 2));
      //DetElement wndDE(hrppdDE, "wnd_de", 0);
      //wndDE.setPlacement(wndPV);
      
      //  double pitch = xysize + _HRPPD_INSTALLATION_GAP_;
      double xyactive = _HRPPD_ACTIVE_AREA_SIZE_;
      
      accu += wndthick;

      // Photocathode layer (sensitive volume);
      //
      auto pcBox = Box(xyactive / 2, xyactive / 2, pdthick / 2);
      TString pcName; pcName.Form("%s-photocathode-%02d", detName.c_str(), imod);
      Volume pcVol(pcName.Data(), pcBox, HRPPD_pcMat);
      pcVol.setSensitiveDetector(sens);
      
      pcVol.setVisAttributes(sensorVis);
      PlacedVolume pcPV = hrppdVol_air.placeVolume(pcVol, Position(0.0, 0.0, accu + pdthick / 2));
      {
	pcPV.addPhysVolID("hrppd", imod);
	    
	// sensor DetElement
	/*auto*/ sensorID = encodeSensorID(pcPV.volIDs());
	//printf("@S@: %lu\n", sensorID);
	TString deName; deName.Form("%s-sensor-%02d", detName.c_str(), imod);
	DetElement pcDE(sdet, deName.Data(), sensorID);
	pcDE.setPlacement(pcPV);
      }

      //accu += pdthick;
      
      double xyopen   = _HRPPD_OPEN_AREA_SIZE_;
      double certhick = _HRPPD_CERAMIC_BODY_THICKNESS_; //, zcer = azOffset + wndthick + certhick/2;
      
      // Ceramic body (sidewall and anode);
      //
      Box cerbox(xysize / 2, xysize / 2, certhick / 2);
      Box cut_box(xyopen / 2, xyopen / 2, certhick / 2);
      
      SubtractionSolid ceramic(cerbox, cut_box, Position(0, 0, -_HRPPD_BASEPLATE_THICKNESS_));
      
      TString cerName; cerName.Form("%s-ceramic-%02d", detName.c_str(), imod);
      Volume ceramicVol(cerName.Data(), ceramic, HRPPD_MPDMat);
      //ceramicVol.setVisAttributes(gasvolVis);
      
      //PlacedVolume ceramicPV =
      hrppdVol_air.placeVolume(ceramicVol, Position(0.0, 0.0, accu + certhick / 2));
      //DetElement ceramicDE(sdet, "ceramic_de", 0);
      //ceramicDE.setPlacement(ceramicPV);
    
      // Effective anode plating layer
      //
      Box plating_solid(xyopen / 2, xyopen / 2, _HRPPD_PLATING_LAYER_THICKNESS_ / 2);
      TString pltName; pltName.Form("%s-plating-%02d", detName.c_str(), imod);
      Volume platingVol(pltName.Data(), plating_solid, HRPPD_MPDMat);
      //platingVol.setVisAttributes(gasvolVis);
      //PlacedVolume platingPV =
      // Place somewhere in the middle of the ceramic body gap;
      hrppdVol_air.placeVolume(platingVol, Position(0.0, 0.0, accu + certhick / 2));
      //DetElement platingDE(sdet, "plating_de", 0);
      //platingDE.setPlacement(platingPV);
      
      // Effective MCP layer
      //
      Box mcp_solid(xyopen / 2, xyopen / 2, _EFFECTIVE_MCP_THICKNESS_ / 2);
      TString mcpName; mcpName.Form("%s-mcp-%02d", detName.c_str(), imod);
      Volume mcpVol(mcpName.Data(), mcp_solid, HRPPD_MPDMat);
      //mcpVol.setVisAttributes(gasvolVis);
      //PlacedVolume mcpPV =
      hrppdVol_air.placeVolume(
			       mcpVol, Position(0.0, 0.0,
						accu + certhick / 2 + _HRPPD_PLATING_LAYER_THICKNESS_ +
						_EFFECTIVE_MCP_THICKNESS_ / 2));
      //DetElement mcpDE(sdet, "mcp_de", 0);
      //mcpDE.setPlacement(mcpPV);
    
      accu += certhick;// + 1 * mm;

      {
	auto _READOUT_PCB_THICKNESS_ = description.constant<double>("READOUT_PCB_THICKNESS");
	auto _READOUT_PCB_SIZE_      = description.constant<double>("READOUT_PCB_SIZE");
	auto _ASIC_SIZE_XY_          = description.constant<double>("ASIC_SIZE_XY");
	auto _ASIC_THICKNESS_        = description.constant<double>("ASIC_THICKNESS");
	
	// PCB Board
	//
	Box pcb_solid(_READOUT_PCB_SIZE_ / 2, _READOUT_PCB_SIZE_ / 2, _READOUT_PCB_THICKNESS_ / 2);
	TString pcbName; pcbName.Form("%s-pcb-%02d", detName.c_str(), imod);
	Volume pcbVol(pcbName.Data(), pcb_solid, HRPPD_PCBMat);
	//pcbVol.setVisAttributes(gasvolVis);
	//PlacedVolume pcbPV =
	hrppdVol_air.placeVolume(pcbVol, Position(0.0, 0.0, accu + _READOUT_PCB_THICKNESS_ / 2));
	//DetElement pcbDE(sdet, "pcb_de", 0);
	//pcbDE.setPlacement(pcbPV);
	
	accu += _READOUT_PCB_THICKNESS_;// + 0.001;
	
	// ASICs
	//
	Box asic_solid(_ASIC_SIZE_XY_ / 2, _ASIC_SIZE_XY_ / 2, _ASIC_THICKNESS_ / 2);
	TString asicName; asicName.Form("%s-asic-%02d", detName.c_str(), imod);
	Volume asicVol(asicName.Data(), asic_solid, HRPPD_ASICMat);
	asicVol.setVisAttributes(asicVis);
	
	double asic_pitch = _READOUT_PCB_SIZE_ / 2;
	
	//imod = 0;
	
	for (unsigned ix = 0; ix < 2; ix++) {
	  double xOffset = asic_pitch * (ix - (2 - 1) / 2.);
	  
	  for (unsigned iy = 0; iy < 2; iy++) {
	    double yOffset = asic_pitch * (iy - (2 - 1) / 2.);
	    
	    //auto asicPV =
	    hrppdVol_air.placeVolume(asicVol, Position(xOffset, yOffset, accu + _ASIC_THICKNESS_ / 2));
	    
	    //DetElement asicDE(sdet, "asic_de_" + std::to_string(imod), 0);
	    //asicDE.setPlacement(asicPV);
	    
	    //imod++;
	  } //for iy
	} //for ix
	
	accu += _ASIC_THICKNESS_ + 0.01 * mm;
      }

      // Eventually place the whole HRPPD container volume;
      double dz = _FIDUCIAL_VOLUME_LENGTH_/2 - _SENSOR_AREA_LENGTH_ + _HRPPD_CONTAINER_VOLUME_HEIGHT_/2;
      /*auto sensorPV =*/ pfRICH_volume.placeVolume(hrppdVol_air, Position(xy.X(), xy.Y(), dz));
      
#if _LATER_
      if (!debug_optics) {
	SkinSurface sensorSkin(description, sensorDE,
			       "sensor_optical_surface_" + std::to_string(imod), sensorSurf,
			       sensorVol);
	sensorSkin.isValid();
      };
#endif
      
#ifdef _WITH_OPTICS_
      {
	// Photocathode surface;
	double xOffset = xy.X(), yOffset = xy.Y();
	auto surface = new FlatSurface((1/mm)*TVector3(sign*xOffset, yOffset, -1700.*mm),//sign*(fvzOffset + zpdc)),
				       TVector3(1*sign,0,0), TVector3(0,-1,0));
	
	{
	  // '0': QRICH has no division in sectors (unlike e.g. dRICH);
	  unsigned sector = 0;
	  
	  // Just two configurations: with and without a reflection on a conical mirror;
	  //#ifdef _WITH_MIRROR_
	  //for(unsigned iq=0; iq<2; iq++) {
	  //#else
	  for(unsigned iq=0; iq<1; iq++) {
	    //#endif
	    auto irt = pd->AllocateIRT(sector, sensorID);
	      
	    // Aerogel and acrylic;
	    if (cdet->m_OpticalBoundaries[CherenkovDetector::Upstream].find(sector) != 
		cdet->m_OpticalBoundaries[CherenkovDetector::Upstream].end())
	      for(auto boundary: cdet->m_OpticalBoundaries[CherenkovDetector::Upstream][sector])
		irt->AddOpticalBoundary(boundary);
	    
	    // Optional mirror reflection;
	    //#ifdef _WITH_MIRROR_
	    //if (iq) irt->AddOpticalBoundary(mboundary);
	    //#endif
	      
	    // Fused silica windows;
	    if (cdet->m_OpticalBoundaries[CherenkovDetector::Downstream].find(sector) != 
		cdet->m_OpticalBoundaries[CherenkovDetector::Downstream].end())
	      for(auto boundary: cdet->m_OpticalBoundaries[CherenkovDetector::Downstream][sector])
		irt->AddOpticalBoundary(boundary);
	    
	    // Terminate the optical path;
	    pd->AddItselfToOpticalBoundaries(irt, surface);
	  } //for iq
	}
      } 
#endif
	  
      imod++;
    } //for coord
  } 
  
  // Write out optical configuration;
#ifdef _WITH_OPTICS_
  geometry->Write();
  output_file->Close();
#endif
  
  return sdet;
} // createDetector()

// -------------------------------------------------------------------------------------

// clang-format off
DECLARE_DETELEMENT(epic_PFRICH, createDetector)
