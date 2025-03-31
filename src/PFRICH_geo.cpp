// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2023  Wenliang (Bill) Li, Alexander Kiselev, Karthik Suresh

//----------------------------------
// pfRICH: Proximity Focusing RICH
// Author: Wenliang (Bill) Li
//
// - Design Adapted from standalone Geant4 description by
//   Alexander Kiselev and Chandradoy Chatterjee
//----------------------------------

#include "DD4hep/DetFactoryHelper.h"
//#include "DD4hep/OpticalSurfaces.h"
//#include "DD4hep/Printout.h"
//#include "DDRec/DetectorData.h"
//#include "DDRec/Surface.h"

//#include <XML/Helper.h>
//#include "XML/Layering.h"
//#include "TVector3.h"
#include "TVector2.h"

//#include "TGeoElement.h"
//#include "TGeoManager.h"
//#include "TInterpreter.h"
//#include "TUri.h"

//using namespace std;
using namespace dd4hep;
//using namespace dd4hep::rec;
//using namespace dd4hep::detail;

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

  // BUILD SENSORS ///////////////////////
  // solid and volume: single sensor module

  //+OpticalSurfaceManager surfMgr = description.surfaceManager();

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
#if _TODAY_
  auto HRPPD_PCBMat    = description.material(sensorElem.attr<std::string>(_Unicode(pcbmat)));
  auto HRPPD_ASICMat   = description.material(sensorElem.attr<std::string>(_Unicode(asicmat)));

  double vesselRmin0 = dims.attr<double>(_Unicode(rmin0));
  double vesselRmin1 = dims.attr<double>(_Unicode(rmin1));
  double vesselRmax0 = dims.attr<double>(_Unicode(rmax0));
  double vesselRmax1 = dims.attr<double>(_Unicode(rmax1));

  int imod = 0; // module number

  auto gasvolMat = description.material(detElem.attr<std::string>(_Unicode(gas)));
#endif
  auto gasvolVis = description.visAttributes("DRICH_gas_vis");
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
  // create a unique sensor ID from a sensor's PlacedVolume::volIDs
  auto encodeSensorID = [&readoutCoder](auto ids) {
    uint64_t enc = 0;
    for (const auto& [idField, idValue] : ids)
      enc |= uint64_t(idValue) << readoutCoder[idField].offset();
    return enc;
  };

#if _TODAY_
  auto mirrorElem = detElem.child(_Unicode(mirror));
  auto mirrorMat  = description.material(mirrorElem.attr<std::string>(_Unicode(material)));
  auto mirrorVis  = description.visAttributes(mirrorElem.attr<std::string>(_Unicode(vis)));
#endif

  /*--------------------------------------------------*/
  // Vessel
  // auto vesselMat = description.material(detElem.attr<std::string>(_Unicode(material)));
  auto vesselGas = description.material(detElem.attr<std::string>(_Unicode(gas)));

  /*--------------------------------------------------*/

  // Mirrors
#if _TODAY_
  float _CONICAL_MIRROR_INNER_RADIUS_ = description.constant<double>("CONICAL_MIRROR_INNER_RADIUS");
  float _CONICAL_MIRROR_OUTER_RADIUS_ = description.constant<double>("CONICAL_MIRROR_OUTER_RADIUS");
  float _INNER_MIRROR_THICKNESS_      = description.constant<double>("INNER_MIRROR_THICKNESS");
  float _OUTER_MIRROR_THICKNESS_      = description.constant<double>("OUTER_MIRROR_THICKNESS");


  // Aerogel
#endif


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

  //
  // Gas volume;
  //
  auto _VESSEL_FRONT_SIDE_THICKNESS_ = description.constant<double>("VESSEL_FRONT_SIDE_THICKNESS");
  auto _SENSOR_AREA_LENGTH_          = description.constant<double>("SENSOR_AREA_LENGTH");
  auto _VESSEL_OUTER_WALL_THICKNESS_ = description.constant<double>("VESSEL_OUTER_WALL_THICKNESS");
  
  double gas_volume_length =  _FIDUCIAL_VOLUME_LENGTH_ - _VESSEL_FRONT_SIDE_THICKNESS_ - _SENSOR_AREA_LENGTH_;
  double gas_volume_radius = _VESSEL_OUTER_RADIUS_ - _VESSEL_OUTER_WALL_THICKNESS_;
  double gas_volume_offset = -(_SENSOR_AREA_LENGTH_ - _VESSEL_FRONT_SIDE_THICKNESS_)/2;
  Tube gasTube(0.0, gas_volume_radius, gas_volume_length/2);
  SubtractionSolid gasSolid(gasTube, FlangeCut(description, gas_volume_length + 1*mm, _FLANGE_CLEARANCE_));
  Volume gasVolume(detName + "_GasVol", gasSolid, vesselGas);
  /*PlacedVolume gasPV =*/ pfRICH_volume.placeVolume(gasVolume, Position(0, 0, gas_volume_offset));

  //
  // Vessel walls -> cut'n'paste from pfRICH standalone code;;
  //

  auto _BUILDING_BLOCK_CLEARANCE_ = description.constant<double>("BUILDING_BLOCK_CLEARANCE");
  auto _VESSEL_INNER_WALL_THICKNESS_ = description.constant<double>("VESSEL_INNER_WALL_THICKNESS");
  
  // To be used in boolean operations in several places;
  auto flange = FlangeCut(description, gas_volume_length + 1*mm, _FLANGE_CLEARANCE_ + _VESSEL_INNER_WALL_THICKNESS_ + 
			  _BUILDING_BLOCK_CLEARANCE_);
  
  //
  // Aerogel;
  //
  {
    const int _AEROGEL_BAND_COUNT_ = 3;
    auto _FLANGE_EPIPE_DIAMETER_ = description.constant<double>("FLANGE_EPIPE_DIAMETER");
    auto _AEROGEL_INNER_WALL_THICKNESS_ =
      description.constant<double>("AEROGEL_INNER_WALL_THICKNESS");
    //float _VESSEL_OUTER_RADIUS_         = description.constant<double>("VESSEL_OUTER_RADIUS");
    //float _FLANGE_CLEARANCE_         = description.constant<double>("FLANGE_CLEARANCE");
    //float _AEROGEL_BAND_COUNT_ = aerogel_band_count;
    auto _AEROGEL_SEPARATOR_WALL_THICKNESS_ =
      description.constant<double>("AEROGEL_SEPARATOR_WALL_THICKNESS");
    auto _AEROGEL_OUTER_WALL_THICKNESS_ =
      description.constant<double>("AEROGEL_OUTER_WALL_THICKNESS");
    float m_r0min = _FLANGE_EPIPE_DIAMETER_ / 2 + _FLANGE_CLEARANCE_ + _VESSEL_INNER_WALL_THICKNESS_ +
      _BUILDING_BLOCK_CLEARANCE_;
    float m_r0max = gas_volume_radius - _BUILDING_BLOCK_CLEARANCE_;
    
    const unsigned adim[_AEROGEL_BAND_COUNT_] = {9, 14, 20};
    double rheight =
      (m_r0max - m_r0min - (_AEROGEL_BAND_COUNT_ - 1) * _AEROGEL_SEPARATOR_WALL_THICKNESS_ -
       _AEROGEL_INNER_WALL_THICKNESS_ - _AEROGEL_OUTER_WALL_THICKNESS_) /
      _AEROGEL_BAND_COUNT_;
    
    double agthick = aerogelThickness;//description.constant<double>("PFRICH_aerogel_thickness");//2.5; // cm
    double gzOffset = -gas_volume_length / 2 + _BUILDING_BLOCK_CLEARANCE_;// + agthick / 2;
    
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
  }

  //
  // Acrylic filter;
  //
  //+double _ACRYLIC_THICKNESS_ = 0.3;
#if _TODAY_
  {
    double acthick = _ACRYLIC_THICKNESS_;
    // m_gzOffset += acthick/2;
    
    Tube ac_tube(m_r0min + 3, m_r0max - 1, acthick / 2, 0 * degree, 360 * degree);
    SubtractionSolid ac_shape(ac_tube, flange_final_shape);
    
    Volume acVol(detName + "_ac", ac_shape, gasvolMat);
    
    PlacedVolume ac_PV = pfRICH_volume.placeVolume(acVol, Position(0, 0, -21.3));
    
    DetElement acDE(sdet, "ac_de", 0);
    acDE.setPlacement(ac_PV);
  }
#endif

  //
  // Mirrors;
  //
  /// Mirror construction
#if _TODAY_
  double mlen = m_gas_volume_length - _BUILDING_BLOCK_CLEARANCE_;

  mlen -= _BUILDING_BLOCK_CLEARANCE_ + _HRPPD_SUPPORT_GRID_BAR_HEIGHT_;

  double mirror_r0[2] = {m_r0min, m_r0max};
  double mirror_r1[2] = {_CONICAL_MIRROR_INNER_RADIUS_, _CONICAL_MIRROR_OUTER_RADIUS_};

  for (unsigned im = 0; im < 2; im++) {

    double mirror_thickness = im ? _OUTER_MIRROR_THICKNESS_ : _INNER_MIRROR_THICKNESS_;

    if (im) {

      Cone mirror_outer_cone_shape(mlen / 2.0, mirror_r0[im], mirror_r0[im] + mirror_thickness,
                                   mirror_r1[im], mirror_r1[im] + mirror_thickness);

      Volume outer_mirrorVol(detName + "_outer_mirror", mirror_outer_cone_shape, mirrorMat);

      PlacedVolume mirror_outerPV = pfRICH_volume.placeVolume(outer_mirrorVol, Position(0, 0, 0));

      DetElement mirror_outerDE(sdet, "_outer_mirror_de", 0);
      mirror_outerDE.setPlacement(mirror_outerPV);

    } else {

      Cone mirror_inner_cone_shape(mlen / 2., mirror_r0[im], mirror_r0[im] + mirror_thickness,
                                   mirror_r1[im], mirror_r1[im] + mirror_thickness);

      SubtractionSolid mirror_inner_sub(mirror_inner_cone_shape, flange_final_shape);

      Volume inner_mirrorVol(detName + "_inner_mirror", mirror_inner_sub, mirrorMat);

      PlacedVolume mirror_innerPV = pfRICH_volume.placeVolume(inner_mirrorVol, Position(0, 0, 0));

      DetElement mirror_innerDE(sdet, "_inner_mirror_de", 0);
      mirror_innerDE.setPlacement(mirror_innerPV);
    }

  } //for im
#endif


  
  
#if _TODAY_
  /// tank solids

  //double boreDelta = vesselRmin1 - vesselRmin0;
  //Cone vesselTank(vesselLength / 2.0, vesselRmin1, vesselRmax1, vesselRmin0, vesselRmax0);
  //Cone gasvolTank(vesselLength / 2.0 - windowThickness, vesselRmin1 + wallThickness,
  //              vesselRmax1 - wallThickness, vesselRmin0 + wallThickness,
  //              vesselRmax0 - wallThickness);

  Cone vesselWall(vesselLength / 2.0, vesselRmax1 - 0.1, vesselRmax1, vesselRmax0 - 0.1,
                  vesselRmax0);

  //Box gasvolBox(1000, 1000, 1000);

  //Solid gasvolSolid;
  //gasvolSolid = gasvolTank;

  //Solid vesselSolid;
  //vesselSolid = vesselTank;

  //Solid mirrorSolid;
  //mirrorSolid = mirror_cone;

  Solid wallSolid;
  wallSolid = vesselWall;

  Volume vesselVol(detName + "_vesel_vol", wallSolid, vesselMat);
  vesselVol.setVisAttributes(vesselVis);

  PlacedVolume vesselPV = pfRICH_volume.placeVolume(vesselVol, Position(0, 0, 0));
  DetElement vesselDE(sdet, "vessel_de", 0);
  vesselDE.setPlacement(vesselPV);

  // BUILD RADIATOR //////////////////////////////////////

  // solid and volume: create aerogel and filter
  //Cone aerogelSolid(aerogelThickness / 2,
  //                radiatorRmin + boreDelta * aerogelThickness / vesselLength, /* at backplane */
  //                radiatorRmax, radiatorRmin,                                 /* at frontplane */
  //                radiatorRmax);
#if _UNUSED_
  Cone filterSolid(filterThickness / 2,
                   radiatorRmin + boreDelta *
                                      (aerogelThickness + airgapThickness + filterThickness) /
                                      vesselLength, /* at backplane */
                   radiatorRmax,
                   radiatorRmin + boreDelta * (aerogelThickness + airgapThickness) /
                                      vesselLength, /* at frontplane */
                   radiatorRmax);
  //Volume aerogelVol(detName + "_aerogel", aerogelSolid, aerogelMat);
  Volume filterVol(detName + "_filter", filterSolid, filterMat);
#endif
  
  // radiator material names
  description.add(Constant("PFRICH_aerogel_material", aerogelMat.ptr()->GetName(), "string"));
  description.add(Constant("PFRICH_filter_material", filterMat.ptr()->GetName(), "string"));
  description.add(Constant("PFRICH_gasvol_material", gasvolMat.ptr()->GetName(), "string"));

  Box sensorSolid(sensorSide / 2., sensorSide / 2., sensorThickness / 2.);
  Volume sensorVol(detName + "_sensor", sensorSolid, sensorMat);
  sensorVol.setVisAttributes(sensorVis);

  // -- Mirrors ---------------------------------------------------------------------------------
  // Some "standard" value applied to all mirrors;
  // At the downstream (sensor plane) location; upstream radii are calculated automatically;

  //Volume mirrorVol(detName, mirrorSolid, mirrorMat);
  //mirrorVol.setVisAttributes(mirrorVis);
#endif

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
    // auto _HRPPD_OPEN_AREA_SIZE_   = description.constant<double>("HRPPD_OPEN_AREA_SIZE");
    auto _HRPPD_ACTIVE_AREA_SIZE_ = description.constant<double>("HRPPD_ACTIVE_AREA_SIZE");
    //auto _HRPPD_CERAMIC_BODY_THICKNESS_ =
      description.constant<double>("HRPPD_CERAMIC_BODY_THICKNESS");
    //auto _HRPPD_BASEPLATE_THICKNESS_ = description.constant<double>("HRPPD_BASEPLATE_THICKNESS");
    //auto _HRPPD_PLATING_LAYER_THICKNESS_ =
      description.constant<double>("HRPPD_PLATING_LAYER_THICKNESS");
    //auto _EFFECTIVE_MCP_THICKNESS_ = description.constant<double>("EFFECTIVE_MCP_THICKNESS");
    
    //auto _READOUT_PCB_THICKNESS_ = description.constant<double>("READOUT_PCB_THICKNESS");
    //auto _READOUT_PCB_SIZE_      = description.constant<double>("READOUT_PCB_SIZE");
    //auto _ASIC_SIZE_XY_   = description.constant<double>("ASIC_SIZE_XY");
    //auto _ASIC_THICKNESS_ = description.constant<double>("ASIC_THICKNESS");
    
    double xysize = _HRPPD_TILE_SIZE_, wndthick = _HRPPD_WINDOW_THICKNESS_;

    // HRPPD assembly container volume;
    double hrppd_container_volume_thickness = _HRPPD_CONTAINER_VOLUME_HEIGHT_;

  
#if _LATER_
    // Ceramic body
    Box cerbox(xysize / 2, xysize / 2, certhick / 2);
    Box cut_box(xyopen / 2, xyopen / 2, certhick / 2);
    
    SubtractionSolid ceramic(cerbox, cut_box, Position(0, 0, -_HRPPD_BASEPLATE_THICKNESS_));
    
    Volume ceramicVol(detName + "_ceramic", ceramic, HRPPD_MPDMat);
    ceramicVol.setVisAttributes(gasvolVis);
    
    PlacedVolume ceramicPV =
      hrppdVol_air.placeVolume(ceramicVol, Position(0.0, 0.0, accu + certhick / 2));
    DetElement ceramicDE(sdet, "ceramic_de", 0);
    ceramicDE.setPlacement(ceramicPV);
    
    // Plating body
    
    Box plating_solid(xyopen / 2, xyopen / 2, _HRPPD_PLATING_LAYER_THICKNESS_ / 2);
    Volume platingVol(detName + "_plating", plating_solid, HRPPD_MPDMat);
    
    platingVol.setVisAttributes(gasvolVis);
    PlacedVolume platingPV =
      hrppdVol_air.placeVolume(platingVol, Position(0.0, 0.0, accu + certhick / 2));
    DetElement platingDE(sdet, "plating_de", 0);
    platingDE.setPlacement(platingPV);
    
    // MCP body
    
    Box mcp_solid(xyopen / 2, xyopen / 2, _EFFECTIVE_MCP_THICKNESS_ / 2);
    Volume mcpVol(detName + "_mcp", mcp_solid, HRPPD_MPDMat);
    
    mcpVol.setVisAttributes(gasvolVis);
    PlacedVolume mcpPV = hrppdVol_air.placeVolume(
						  mcpVol, Position(0.0, 0.0,
								   accu + certhick / 2 + _HRPPD_PLATING_LAYER_THICKNESS_ / 2 +
								   _EFFECTIVE_MCP_THICKNESS_ / 2));
    DetElement mcpDE(sdet, "mcp_de", 0);
    mcpDE.setPlacement(mcpPV);
#endif

    
#if _LATER_
    accu += certhick + 1 * mm;

    /// PCB Board
    Box pcb_solid(_READOUT_PCB_SIZE_ / 2, _READOUT_PCB_SIZE_ / 2, _READOUT_PCB_THICKNESS_ / 2);
    Volume pcbVol(detName + "_pcb", pcb_solid, HRPPD_PCBMat);
    
    pcbVol.setVisAttributes(gasvolVis);
    PlacedVolume pcbPV =
      hrppdVol_air.placeVolume(pcbVol, Position(0.0, 0.0, accu + _READOUT_PCB_THICKNESS_ / 2));
    
    DetElement pcbDE(sdet, "pcb_de", 0);
    pcbDE.setPlacement(pcbPV);
    
    accu += _READOUT_PCB_THICKNESS_ + 0.001;
    
    // ASIC Board
    Box asic_solid(_ASIC_SIZE_XY_ / 2, _ASIC_SIZE_XY_ / 2, _ASIC_THICKNESS_ / 2);
    Volume asicVol(detName + "_asic", asic_solid, HRPPD_ASICMat);
    asicVol.setVisAttributes(mirrorVis);
    
    double asic_pitch = _READOUT_PCB_SIZE_ / 2;
    
    imod = 0;

    for (unsigned ix = 0; ix < 2; ix++) {
      double xOffset = asic_pitch * (ix - (2 - 1) / 2.);
      
      for (unsigned iy = 0; iy < 2; iy++) {
	double yOffset = asic_pitch * (iy - (2 - 1) / 2.);
	
	auto asicPV = hrppdVol_air.placeVolume(
					       asicVol, Position(xOffset, yOffset, accu + _ASIC_THICKNESS_ / 2));
	
	DetElement asicDE(sdet, "asic_de_" + std::to_string(imod), 0);
	asicDE.setPlacement(asicPV);
	
	imod++;
	
      } //for iy
    } //for ix
    
    accu += _ASIC_THICKNESS_ + 0.01 * mm;
#endif
  // Loading the coordinates

    //#if _TODAY_
    std::vector<std::pair<TVector2, bool>> coord;
    {
      unsigned const hdim              = 9;
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
    
    /// Set sensors into the coordinates
    ///
#if 1//_NOW_
    for (auto xyptr : coord) {
      auto& xy = xyptr.first;

      double pdthick = 0.01*mm;//0.001;

      Box hrppd_Solid(xysize / 2, xysize / 2, hrppd_container_volume_thickness / 2);
      TString hrppdName; hrppdName.Form("%s-hrppd-%02d", detName.c_str(), imod);
      Volume hrppdVol_air(hrppdName.Data(), hrppd_Solid, air);
      
      //-hrppdVol_air.setSensitiveDetector(sens);
      hrppdVol_air.setVisAttributes(gasvolVis);
      //DetElement hrppdDE(sdet, "hrppd_de", 0);

    // Quartz Window
#if 1//_BACK_
      Box wnd_Solid(xysize / 2, xysize / 2, wndthick / 2);
      TString wndName; wndName.Form("%s-window-%02d", detName.c_str(), imod);
      Volume wndVol(wndName.Data(), wnd_Solid, HRPPD_WindowMat);
      wndVol.setVisAttributes(gasvolVis);
      
      double accu = -hrppd_container_volume_thickness / 2;
      
      PlacedVolume wndPV = hrppdVol_air.placeVolume(wndVol, Position(0, 0, accu + wndthick / 2));
      //DetElement wndDE(hrppdDE, "wnd_de", 0);
      //wndDE.setPlacement(wndPV);
      
      //#if _TODAY_
      //  double pitch = xysize + _HRPPD_INSTALLATION_GAP_;
      double xyactive = _HRPPD_ACTIVE_AREA_SIZE_;
      //+double xyopen   = _HRPPD_OPEN_AREA_SIZE_;
      //+double certhick = _HRPPD_CERAMIC_BODY_THICKNESS_; //, zcer = azOffset + wndthick + certhick/2;
      
      accu += wndthick;
#endif
    
      auto pcBox = Box(xyactive / 2, xyactive / 2, pdthick / 2);
      TString pcName; pcName.Form("%s-photocathode-%02d", detName.c_str(), imod);
      Volume pcVol(pcName.Data(), pcBox, HRPPD_pcMat);
      pcVol.setSensitiveDetector(sens);
      
      pcVol.setVisAttributes(gasvolVis);
      PlacedVolume pcPV = hrppdVol_air.placeVolume(pcVol, Position(0.0, 0.0, accu + pdthick / 2));
      {
	pcPV.addPhysVolID("hrppd", imod);
	    
	// sensor DetElement
	auto sensorID = encodeSensorID(pcPV.volIDs());
	//printf("@S@: %lu\n", sensorID);
	TString deName; deName.Form("%s-sensor-%02d", detName.c_str(), imod);//x, iy);
	DetElement pcDE(sdet, deName.Data(), sensorID);
	pcDE.setPlacement(pcPV);
      }
#if _BACK_
      //+DetElement pdboxDE(sdet, "pdbox_de", 0);
      //+pdboxDE.setPlacement(pdboxPV);
      accu += pdthick;
      
      Box qdbox_solid(xyactive / 2, xyactive / 2, pdthick / 2);
      TString absName; absName.Form("%s-absorber-%02d", detName.c_str(), imod);
      Volume qdboxVol(absName.Data(), qdbox_solid, HRPPD_MPDMat);
      
      // Yes, in case of ePIC implementation this layer is placed behind the actual photocathode
      // sensitive material in order to terminate optical photons escaping from photocathode
      // (which is somehow required to create hits);
      qdboxVol.setVisAttributes(gasvolVis);
      /*PlacedVolume qdboxPV =*/ hrppdVol_air.placeVolume(qdboxVol, Position(0.0, 0.0, accu + pdthick / 2));
      
      //+DetElement qdboxDE(sdet, "qdbox_de", 0);
      //+qdboxDE.setPlacement(qdboxPV);
#endif
      
#if _NOW_
      double sx = xy.X();
      double sy = xy.Y();
      
      // placement (note: transformations are in reverse order)
      auto sensorPlacement =
        Transform3D(Translation3D(sensorPlanePos.x(), sensorPlanePos.y(),
                                  sensorPlanePos.z() + 34) * // move to reference position
                    Translation3D(sx, sy, 0.)                // move to grid position
		    );
#endif
      double dz = _FIDUCIAL_VOLUME_LENGTH_/2 - _SENSOR_AREA_LENGTH_ + _HRPPD_CONTAINER_VOLUME_HEIGHT_/2;
      //-auto sensorPV = pfRICH_volume.placeVolume(hrppdVol_air, sensorPlacement);
      /*auto sensorPV =*/ pfRICH_volume.placeVolume(hrppdVol_air, Position(xy.X(), xy.Y(), dz));//sensorPlacement);
      
      // properties
#if _LATER_
      sensorPV.addPhysVolID(/*"module"*/"hrppd", imod); // NOTE: must be consistent with `sensorIDfields`
      auto imodEnc = encodeSensorID(sensorPV.volIDs());
      DetElement sensorDE(sdet, "sensor_de_" + std::to_string(imod), imodEnc);
      sensorDE.setPlacement(sensorPV);
#endif
#if _LATER_
      if (!debug_optics) {
	SkinSurface sensorSkin(description, sensorDE,
			       "sensor_optical_surface_" + std::to_string(imod), sensorSurf,
			       sensorVol);
	sensorSkin.isValid();
      };
#endif
      
      // increment sensor module number
      imod++;
    }//for coord
#endif 
  } 
  
  return sdet;
} // createDetector()

// -------------------------------------------------------------------------------------

// clang-format off
DECLARE_DETELEMENT(epic_PFRICH, createDetector)
