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

#include "TVector2.h"

using namespace dd4hep;

#ifdef WITH_IRT2_SUPPORT
#include "IRT2/CherenkovDetectorCollection.h"

#include "IRT2/ConicalSurface.h"

using namespace IRT2;
#endif

#define _SENSOR_PLANE_GEOMETRIC_EFFICIENCY_ (1.00)
#define _SAFETY_FACTOR_ (0.70)

// -------------------------------------------------------------------------------------

static UnionSolid FlangeCut(Detector& description, double length, double clearance) {
  // FIXME: not the most efficient way to recalculate them every time new,
  // but code is more readable, and overhead is negligible anyway;
  auto _FLANGE_EPIPE_DIAMETER_ = description.constant<double>("FLANGE_EPIPE_DIAMETER");
  auto _FLANGE_HPIPE_DIAMETER_ = description.constant<double>("FLANGE_HPIPE_DIAMETER");
  auto _FLANGE_HPIPE_OFFSET_   = description.constant<double>("FLANGE_HPIPE_OFFSET");

  // A wedge bridging two cylinders;
  Tube eflange(0.0, _FLANGE_EPIPE_DIAMETER_ / 2 + clearance, length / 2);
  Tube hflange(0.0, _FLANGE_HPIPE_DIAMETER_ / 2 + clearance, length / 2);

  double r0 = _FLANGE_EPIPE_DIAMETER_ / 2 + clearance;
  double r1 = _FLANGE_HPIPE_DIAMETER_ / 2 + clearance;
  double L  = _FLANGE_HPIPE_OFFSET_;
  double a  = r0 * L / (r0 - r1);
  double b  = r0 * r0 / a;
  double c  = r1 * (a - b) / r0;

  // GEANT variables to define G4Trap;
  double pDz = length / 2, pTheta = 0.0, pPhi = 0.0, pDy1 = (a - b - c) / 2, pDy2 = pDy1;
  double pDx1 = sqrt(r0 * r0 - b * b), pDx2 = pDx1 * r1 / r0, pDx3 = pDx1, pDx4 = pDx2, pAlp1 = 0.0,
         pAlp2 = 0.0;

  Trap wedge(pDz, pTheta, pPhi, pDy1, pDx1, pDx2, pAlp1, pDy2, pDx3, pDx4, pAlp2);

  UnionSolid flange_shape(eflange, hflange, Position(-_FLANGE_HPIPE_OFFSET_, 0.0, 0.0));
  Rotation3D rZ(RotationZYX(M_PI / 2.0, 0.0, 0.0));
  Transform3D transform_flange(rZ, Position(-b - pDy1, 0.0, 0.0));

  return UnionSolid(flange_shape, wedge, transform_flange);
} // FlangeCut()

// -------------------------------------------------------------------------------------

static Ref_t createDetector(Detector& description, xml_h e, SensitiveDetector sens) {
  xml::DetElement detElem = e;
  int det_id              = detElem.id();

  std::string det_name = detElem.nameStr();
  Material air         = description.air();

  DetElement sdet(det_name, det_id);

  sens.setType("tracker");
  description.invisible();

  std::string detName = detElem.nameStr();

  int id = detElem.hasAttr(_U(id)) ? detElem.id() : 0;

  // Start optical configuration if needed;
#ifdef WITH_IRT2_SUPPORT
  auto geometry = new CherenkovDetectorCollection();
  auto cdet     = geometry->AddNewDetector(detName.c_str());
#endif

  OpticalSurfaceManager surfMgr = description.surfaceManager();

  // - sensor module
  auto sensorElem      = detElem.child(_Unicode(sensors)).child(_Unicode(module));
  auto readoutName     = detElem.attr<std::string>(_Unicode(readout));
  auto HRPPD_WindowMat = description.material(sensorElem.attr<std::string>(_Unicode(windowmat)));
  auto HRPPD_pcMat     = description.material(sensorElem.attr<std::string>(_Unicode(photocathode)));
  auto HRPPD_MPDMat    = description.material(sensorElem.attr<std::string>(_Unicode(mpdmat)));
  auto HRPPD_PCBMat    = description.material(sensorElem.attr<std::string>(_Unicode(pcbmat)));

  auto gasvolVis = description.visAttributes("DRICH_gas_vis");
  auto sensorVis = description.visAttributes("PFRICH_sensor_vis");

  auto radiatorElem     = detElem.child(_Unicode(radiator));
  auto aerogelElem      = radiatorElem.child(_Unicode(aerogel));
  auto aerogelThickness = aerogelElem.attr<double>(_Unicode(thickness));
  auto aerogelMat       = description.material(aerogelElem.attr<std::string>(_Unicode(material)));

  // readout coder <-> unique sensor ID
  /* - `sensorIDfields` is a list of readout fields used to specify a unique sensor ID
     * - `cellMask` is defined such that a hit's `cellID & cellMask` is the corresponding sensor's unique ID
     * - this redundant generalization is for future flexibility, and consistency with dRICH
     */
  std::vector<std::string> sensorIDfields = {"hrppd"};
  const auto& readoutCoder                = *description.readout(readoutName).idSpec().decoder();
  // determine `cellMask` based on `sensorIDfields`
  uint64_t cellMask = 0;
  for (const auto& idField : sensorIDfields)
    cellMask |= readoutCoder[idField].mask();
  description.add(Constant("PFRICH_cell_mask", std::to_string(cellMask)));
#ifdef WITH_IRT2_SUPPORT
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

#ifdef WITH_IRT2_SUPPORT
  const bool flip   = true;
  const double sign = flip ? -1.0 : 1.0;
#endif

  auto mirrorElem = detElem.child(_Unicode(mirror));
  auto mirrorMat  = description.material(mirrorElem.attr<std::string>(_Unicode(material)));
  auto mirrorVis  = description.visAttributes(mirrorElem.attr<std::string>(_Unicode(vis)));
  auto vesselGas  = description.material(detElem.attr<std::string>(_Unicode(gas)));

  //
  // Fiducial volume; will be split into three subvolumes along the beam line: front wall, gas volume
  // and sensor compartment; FIXME: called 'air' but presently filled with nitrogen;
  //
  auto _FIDUCIAL_VOLUME_LENGTH_ = description.constant<double>("FIDUCIAL_VOLUME_LENGTH");
  auto _VESSEL_OUTER_RADIUS_    = description.constant<double>("VESSEL_OUTER_RADIUS");

  Tube pfRICH_air_volume(0.0, _VESSEL_OUTER_RADIUS_, _FIDUCIAL_VOLUME_LENGTH_ / 2);

  auto _FLANGE_CLEARANCE_ = description.constant<double>("FLANGE_CLEARANCE");
  SubtractionSolid pfRICH_volume_shape(
      pfRICH_air_volume,
      FlangeCut(description, _FIDUCIAL_VOLUME_LENGTH_ + 1 * mm, _FLANGE_CLEARANCE_));
  Volume pfRICH_volume(detName, pfRICH_volume_shape, vesselGas);

  Volume mother                 = description.pickMotherVolume(sdet);
  auto _FIDUCIAL_VOLUME_OFFSET_ = description.constant<double>("FIDUCIAL_VOLUME_OFFSET");
  Transform3D transform(RotationZYX(0, M_PI, 0), Position(0, 0, _FIDUCIAL_VOLUME_OFFSET_));
  PlacedVolume pv = mother.placeVolume(pfRICH_volume, transform);
  if (id != 0)
    pv.addPhysVolID("system", id);
  sdet.setPlacement(pv);

#ifdef WITH_IRT2_SUPPORT
  // FIXME: do it better later;
  double fvOffset = fabs(_FIDUCIAL_VOLUME_OFFSET_);
#endif

  //
  // Gas volume;
  //
  auto _VESSEL_FRONT_SIDE_THICKNESS_ = description.constant<double>("VESSEL_FRONT_SIDE_THICKNESS");
  auto _SENSOR_AREA_LENGTH_          = description.constant<double>("SENSOR_AREA_LENGTH");
  auto _VESSEL_OUTER_WALL_THICKNESS_ = description.constant<double>("VESSEL_OUTER_WALL_THICKNESS");

  double gas_volume_length =
      _FIDUCIAL_VOLUME_LENGTH_ - _VESSEL_FRONT_SIDE_THICKNESS_ - _SENSOR_AREA_LENGTH_;
  double gas_volume_radius = _VESSEL_OUTER_RADIUS_ - _VESSEL_OUTER_WALL_THICKNESS_;
  double gas_volume_offset = -(_SENSOR_AREA_LENGTH_ - _VESSEL_FRONT_SIDE_THICKNESS_) / 2;
#ifdef WITH_IRT2_SUPPORT
  double gvOffset = gas_volume_offset;
#endif

  Tube gasTube(0.0, gas_volume_radius, gas_volume_length / 2);
  SubtractionSolid gasSolid(gasTube,
                            FlangeCut(description, gas_volume_length + 1 * mm, _FLANGE_CLEARANCE_));
  Volume gasVolume(detName + "_GasVol", gasSolid, vesselGas);
  pfRICH_volume.placeVolume(gasVolume, Position(0, 0, gas_volume_offset));
#ifdef WITH_IRT2_SUPPORT
  {
    // FIXME: Z-location does not really matter here, right?;
    auto boundary =
        new FlatSurface(TVector3(0, 0, 0), sign * TVector3(1, 0, 0), TVector3(0, -1, 0));

    auto radiator =
        geometry->SetContainerVolume(cdet, "GasVolume", 0, (G4LogicalVolume*)(0x0), 0, boundary);
    radiator->SetAlternativeMaterialName("N2cherenkov");
  }
#endif

  //
  // Vessel walls -> cut'n'paste from pfRICH standalone code;;
  //

  auto _BUILDING_BLOCK_CLEARANCE_    = description.constant<double>("BUILDING_BLOCK_CLEARANCE");
  auto _VESSEL_INNER_WALL_THICKNESS_ = description.constant<double>("VESSEL_INNER_WALL_THICKNESS");

  // To be used in boolean operations in several places;
  auto flange =
      FlangeCut(description, gas_volume_length + 1 * mm,
                _FLANGE_CLEARANCE_ + _VESSEL_INNER_WALL_THICKNESS_ + _BUILDING_BLOCK_CLEARANCE_);

  double gzOffset = -gas_volume_length / 2 + _BUILDING_BLOCK_CLEARANCE_;

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
    auto _AEROGEL_SEPARATOR_WALL_THICKNESS_ =
        description.constant<double>("AEROGEL_SEPARATOR_WALL_THICKNESS");
    auto _AEROGEL_OUTER_WALL_THICKNESS_ =
        description.constant<double>("AEROGEL_OUTER_WALL_THICKNESS");

    const unsigned adim[_AEROGEL_BAND_COUNT_] = {9, 14, 20};
    double rheight =
        (m_r0max - m_r0min - (_AEROGEL_BAND_COUNT_ - 1) * _AEROGEL_SEPARATOR_WALL_THICKNESS_ -
         _AEROGEL_INNER_WALL_THICKNESS_ - _AEROGEL_OUTER_WALL_THICKNESS_) /
        _AEROGEL_BAND_COUNT_;

    double agthick = aerogelThickness;

    std::string aerogel_name = "a1040";
    int kkcounter            = 0;

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
      TString ag_name = "Tmp";

      if (ir)
        ag_name.Form("%s-%d-00", aerogel_name.c_str(), ir);

      Tube agtube(aerogel_r0, aerogel_r1, agthick / 2, 0 * degree, wd0);

      for (unsigned ia = 0; ia < adim[ir]; ia++) {

        Rotation3D r_aerogel_Z(RotationZYX(ia * apitch, 0.0, 0.0));
        Rotation3D r_aerogel_Zinv(RotationZYX(-1. * ia * apitch, 0.0, 0.0));

        if (ir) {
          ag_name.Form("%s-%d-%02d", "aerogel", ir, ia);

          Volume agtubeVol(ag_name.Data(), agtube, aerogelMat);
          auto aerogelTilePlacement =
              Transform3D(r_aerogel_Z, Position(0.0, 0.0, gzOffset + agthick / 2));
          auto aerogelTilePV = gasVolume.placeVolume(agtubeVol, aerogelTilePlacement);
          DetElement aerogelDE(sdet, "aerogel_de_" + std::to_string(kkcounter), 0);
          aerogelDE.setPlacement(aerogelTilePV);
        } else {
          ag_name.Form("%s-%d-%02d", "aerogel_inner", ir, ia);

          Tube agtube_inner(aerogel_r0, aerogel_r1, agthick / 2, 0 * degree + ia * apitch,
                            wd0 + ia * apitch);
          SubtractionSolid agsub(agtube_inner, flange);
          Volume agsubtubeVol(ag_name.Data(), agsub, aerogelMat);
          auto aerogelTilePlacement =
              Transform3D(RotationZYX(0.0, 0.0, 0.0), Position(0.0, 0.0, gzOffset + agthick / 2));
          auto agsubTilePV = gasVolume.placeVolume(agsubtubeVol, aerogelTilePlacement);

          DetElement aerogelDE(sdet, "agsubTile_de_" + std::to_string(counter), 0);
          aerogelDE.setPlacement(agsubTilePV);
        } //if

        counter++;
        kkcounter++;

      } //for ia
    } // for ir

#ifdef WITH_IRT2_SUPPORT
    {
      TVector3 nx(1 * sign, 0, 0), ny(0, -1, 0);

      auto surface = new FlatSurface(
          sign * (1 / mm) * TVector3(0, 0, fvOffset + gvOffset + gzOffset + agthick / 2), nx, ny);

      auto radiator = geometry->AddFlatRadiator(cdet, "Aerogel", CherenkovDetector::Upstream, 0,
                                                (G4LogicalVolume*)(0x1), 0, surface, agthick / mm);
      radiator->SetAlternativeMaterialName("Aerogel_PFRICH");
      // FIXME: what is it good for in ePIC IRT 2.0 implementation?;
      geometry->AddRadiatorLogicalVolume(radiator, (G4LogicalVolume*)(0x1));
    }
#endif

    // NB: there should be a small gap between aerogel and acrylic placed into the same
    // gas volume, otherwise IRT gets confused;
    gzOffset += agthick + _BUILDING_BLOCK_CLEARANCE_;
  }

  //
  // Acrylic filter;
  //
  double _ACRYLIC_THICKNESS_ = 3 * mm;
  {
    double acthick  = _ACRYLIC_THICKNESS_;
    auto acrylicMat = description.material("Acrylic_PFRICH");

    Tube ac_tube(m_r0min, m_r0max, acthick / 2, 0 * degree, 360 * degree);
    SubtractionSolid ac_shape(ac_tube, flange);
    Volume acVol(detName + "_ac", ac_shape, acrylicMat);

    gasVolume.placeVolume(acVol, Position(0, 0, gzOffset + acthick / 2));

#ifdef WITH_IRT2_SUPPORT
    {
      TVector3 nx(1 * sign, 0, 0), ny(0, -1, 0);

      auto surface = new FlatSurface(
          sign * (1 / mm) * TVector3(0, 0, fvOffset + gvOffset + gzOffset + acthick / 2), nx, ny);

      auto radiator = geometry->AddFlatRadiator(cdet, "Acrylic", CherenkovDetector::Upstream, 0,
                                                (G4LogicalVolume*)(0x2), 0, surface, acthick / mm);
      radiator->SetAlternativeMaterialName("Acrylic_PFRICH");
    }
#endif
  }

#ifdef WITH_IRT2_SUPPORT
  IRT2::OpticalBoundary* mboundaries[2] = {0, 0};
#endif

  //
  // Mirrors;
  //
  {
    auto mirrorSurf = surfMgr.opticalSurface("MirrorSurface_PFRICH");

    auto _CONICAL_MIRROR_INNER_RADIUS_ =
        description.constant<double>("CONICAL_MIRROR_INNER_RADIUS");
    auto _CONICAL_MIRROR_OUTER_RADIUS_ =
        description.constant<double>("CONICAL_MIRROR_OUTER_RADIUS");
    auto _INNER_MIRROR_THICKNESS_ = description.constant<double>("INNER_MIRROR_THICKNESS");
    auto _OUTER_MIRROR_THICKNESS_ = description.constant<double>("OUTER_MIRROR_THICKNESS");

    double mlen =
        gas_volume_length - 4 * _BUILDING_BLOCK_CLEARANCE_ - aerogelThickness - _ACRYLIC_THICKNESS_;
    double mzoffset = (aerogelThickness + _ACRYLIC_THICKNESS_ + 2 * _BUILDING_BLOCK_CLEARANCE_) / 2;

    double mirror_r0[2] = {m_r0min, m_r0max - _BUILDING_BLOCK_CLEARANCE_};
    double mirror_r1[2] = {_CONICAL_MIRROR_INNER_RADIUS_, _CONICAL_MIRROR_OUTER_RADIUS_};

    for (unsigned im = 0; im < 2; im++) {

      double mirror_thickness = im ? _OUTER_MIRROR_THICKNESS_ : _INNER_MIRROR_THICKNESS_;

      if (im) {

        Cone mirror_outer_cone_shape(mlen / 2.0, mirror_r0[im], mirror_r0[im] + mirror_thickness,
                                     mirror_r1[im], mirror_r1[im] + mirror_thickness);

        Volume outer_mirrorVol(detName + "_outer_mirror", mirror_outer_cone_shape, mirrorMat);
        outer_mirrorVol.setVisAttributes(mirrorVis);

        PlacedVolume mirror_outerPV =
            gasVolume.placeVolume(outer_mirrorVol, Position(0, 0, mzoffset));
        DetElement mirror_outerDE(sdet, "_outer_mirror_de", 0);
        mirror_outerDE.setPlacement(mirror_outerPV);
        SkinSurface mirrorSkin(description, mirror_outerDE, "outer_mirror_optical_surface_",
                               mirrorSurf, outer_mirrorVol);
        mirrorSkin.isValid();

#ifdef WITH_IRT2_SUPPORT
        auto msurface = new ConicalSurface(
            sign * (1 / mm) * TVector3(0, 0, fvOffset + gvOffset + mzoffset),
            sign * TVector3(0, 0, 1), mirror_r0[im] / mm, mirror_r1[im] / mm, mlen / mm);

        mboundaries[im] =
            new IRT2::OpticalBoundary(cdet->GetRadiator("GasVolume"), msurface, false);
        // Need to store it in a separate call (?), see a comment in CherenkovDetector.h;
        cdet->StoreOpticalBoundary(mboundaries[im]);
#endif
      } else {

        Cone mirror_inner_cone_shape(mlen / 2., mirror_r0[im], mirror_r0[im] + mirror_thickness,
                                     mirror_r1[im], mirror_r1[im] + mirror_thickness);

        SubtractionSolid mirror_inner_sub(mirror_inner_cone_shape, flange);

        Volume inner_mirrorVol(detName + "_inner_mirror", mirror_inner_sub, mirrorMat);
        inner_mirrorVol.setVisAttributes(mirrorVis);

        PlacedVolume mirror_innerPV =
            gasVolume.placeVolume(inner_mirrorVol, Position(0, 0, mzoffset));
        DetElement mirror_innerDE(sdet, "_inner_mirror_de", 0);
        mirror_innerDE.setPlacement(mirror_innerPV);
        SkinSurface mirrorSkin(description, mirror_innerDE, "inner_mirror_optical_surface_",
                               mirrorSurf, inner_mirrorVol);
        mirrorSkin.isValid();

#ifdef WITH_IRT2_SUPPORT
        auto msurface =
            new ConicalSurface(sign * (1 / mm) * TVector3(0, 0, fvOffset + gvOffset + mzoffset),
                               sign * TVector3(0, 0, 1), (mirror_r0[im] + mirror_thickness) / mm,
                               (mirror_r1[im] + mirror_thickness) / mm, mlen / mm);
        msurface->SetConvex();

        mboundaries[im] =
            new IRT2::OpticalBoundary(cdet->GetRadiator("GasVolume"), msurface, false);
        // Need to store it in a separate call (?), see a comment in CherenkovDetector.h;
        cdet->StoreOpticalBoundary(mboundaries[im]);
#endif
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

#ifdef WITH_IRT2_SUPPORT
    // [0,0]: have neither access to G4VSolid nor to G4Material; IRT code does not care; fine;
    auto pd = new IRT2::CherenkovPhotonDetector(0, 0);

    // FIXME: '0' stands for the unknown (and irrelevant) G4LogicalVolume;
    geometry->AddPhotonDetector(cdet, 0, pd);

    // Cannot access GEANT shapes in the reconstruction code -> store this value;
    pd->SetActiveAreaSize(_HRPPD_ACTIVE_AREA_SIZE_ / mm);

    pd->SetGeometricEfficiency(_SENSOR_PLANE_GEOMETRIC_EFFICIENCY_ * _SAFETY_FACTOR_);
#endif

    // HRPPD assembly container volume;
    double hrppd_container_volume_thickness = _HRPPD_CONTAINER_VOLUME_HEIGHT_;

    // For now assume it is a unique surface, same for all HRPPDs;
#ifdef WITH_IRT2_SUPPORT
    {
      TVector3 nx(1 * sign, 0, 0), ny(0, -1, 0);

      auto surface =
          new FlatSurface(sign * (1 / mm) *
                              TVector3(0, 0,
                                       fvOffset + _FIDUCIAL_VOLUME_LENGTH_ / 2 -
                                           _SENSOR_AREA_LENGTH_ + _HRPPD_WINDOW_THICKNESS_ / 2),
                          nx, ny);

      auto radiator = geometry->AddFlatRadiator(cdet, "QuartzWindow", CherenkovDetector::Downstream,
                                                0, (G4LogicalVolume*)(0x3), 0, surface,
                                                _HRPPD_WINDOW_THICKNESS_ / mm);
      radiator->SetAlternativeMaterialName("AirOptical");
    }
#endif

    // Create a vector of HRPPD XY-coordinates in a separate loop (it will be used place HRPPDs
    // but also create cutouts in the sensor plane later, etc);
    std::vector<std::pair<TVector2, bool>> coord;
    {
      unsigned const hdim              = 9;
      const unsigned flags[hdim][hdim] = {
          // clang-format off
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
      // clang-format on

      for (unsigned ix = 0; ix < hdim; ix++) {
        double xOffset = (_HRPPD_TILE_SIZE_ + _HRPPD_INSTALLATION_GAP_) * (ix - (hdim - 1) / 2.);

        for (unsigned iy = 0; iy < hdim; iy++) {
          double yOffset = (_HRPPD_TILE_SIZE_ + _HRPPD_INSTALLATION_GAP_) * (iy - (hdim - 1) / 2.);
          unsigned flag  = flags[hdim - iy - 1][ix];

          if (!flag)
            continue;

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

      double pdthick = 0.01 * mm;

      uint64_t sensorID = 0x0;

      // HRPPD container volume;
      //
      Box hrppd_Solid(xysize / 2, xysize / 2, hrppd_container_volume_thickness / 2);
      TString hrppdName;
      hrppdName.Form("%s-hrppd-%02d", detName.c_str(), imod);
      // FIXME: may want to use AirOptical here, but then return a dummy absorber
      // layer behind the actual phtocathode;
      Volume hrppdVol_air(hrppdName.Data(), hrppd_Solid, air);
      hrppdVol_air.setVisAttributes(gasvolVis);

      // A running variable to pack layers one after the other one;
      double accu = -hrppd_container_volume_thickness / 2;

      // Quartz Window;
      //
      Box wnd_Solid(xysize / 2, xysize / 2, wndthick / 2);
      TString wndName;
      wndName.Form("%s-window-%02d", detName.c_str(), imod);
      Volume wndVol(wndName.Data(), wnd_Solid, HRPPD_WindowMat);
      wndVol.setVisAttributes(gasvolVis);
      hrppdVol_air.placeVolume(wndVol, Position(0, 0, accu + wndthick / 2));

      double xyactive = _HRPPD_ACTIVE_AREA_SIZE_;

      accu += wndthick;

      // Photocathode layer (sensitive volume);
      //
      auto pcBox = Box(xyactive / 2, xyactive / 2, pdthick / 2);
      TString pcName;
      pcName.Form("%s-photocathode-%02d", detName.c_str(), imod);
      Volume pcVol(pcName.Data(), pcBox, HRPPD_pcMat);
      pcVol.setSensitiveDetector(sens);

      pcVol.setVisAttributes(sensorVis);
      PlacedVolume pcPV = hrppdVol_air.placeVolume(pcVol, Position(0.0, 0.0, accu + pdthick / 2));
      {
        pcPV.addPhysVolID("hrppd", imod);

        // sensor DetElement
        sensorID = encodeSensorID(pcPV.volIDs());
        TString deName;
        deName.Form("%s-sensor-%02d", detName.c_str(), imod);
        DetElement pcDE(sdet, deName.Data(), sensorID);
        pcDE.setPlacement(pcPV);
      }

      // A fake absorber layer behind the photocathode; FIXME: make sure that reflection
      // on the window and photocathode boundary still works correctly (no fake volume as
      // in a standalone code);
#if 1
      {
        TString absName;
        absName.Form("%s-absorber-%02d", detName.c_str(), imod);
        // Recycle the same pcBox shape;
        Volume absVol(absName.Data(), pcBox, HRPPD_PCBMat);
        hrppdVol_air.placeVolume(absVol, Position(0.0, 0.0, accu + pdthick + pdthick / 2));
      }
#endif

      double xyopen   = _HRPPD_OPEN_AREA_SIZE_;
      double certhick = _HRPPD_CERAMIC_BODY_THICKNESS_;

      // Ceramic body (sidewall and anode);
      //
      Box cerbox(xysize / 2, xysize / 2, certhick / 2);
      Box cut_box(xyopen / 2, xyopen / 2, certhick / 2);

      SubtractionSolid ceramic(cerbox, cut_box, Position(0, 0, -_HRPPD_BASEPLATE_THICKNESS_));

      TString cerName;
      cerName.Form("%s-ceramic-%02d", detName.c_str(), imod);
      Volume ceramicVol(cerName.Data(), ceramic, HRPPD_MPDMat);

      hrppdVol_air.placeVolume(ceramicVol, Position(0.0, 0.0, accu + certhick / 2));

      // Effective anode plating layer
      //
      Box plating_solid(xyopen / 2, xyopen / 2, _HRPPD_PLATING_LAYER_THICKNESS_ / 2);
      TString pltName;
      pltName.Form("%s-plating-%02d", detName.c_str(), imod);
      Volume platingVol(pltName.Data(), plating_solid, HRPPD_MPDMat);
      // Place somewhere in the middle of the ceramic body gap;
      hrppdVol_air.placeVolume(platingVol, Position(0.0, 0.0, accu + certhick / 2));

      // Effective MCP layer
      //
      Box mcp_solid(xyopen / 2, xyopen / 2, _EFFECTIVE_MCP_THICKNESS_ / 2);
      TString mcpName;
      mcpName.Form("%s-mcp-%02d", detName.c_str(), imod);
      Volume mcpVol(mcpName.Data(), mcp_solid, HRPPD_MPDMat);
      hrppdVol_air.placeVolume(mcpVol,
                               Position(0.0, 0.0,
                                        accu + certhick / 2 + _HRPPD_PLATING_LAYER_THICKNESS_ +
                                            _EFFECTIVE_MCP_THICKNESS_ / 2));

      accu += certhick;

      {
        auto _READOUT_PCB_THICKNESS_ = description.constant<double>("READOUT_PCB_THICKNESS");
        auto _READOUT_PCB_SIZE_      = description.constant<double>("READOUT_PCB_SIZE");

        // PCB Board
        //
        Box pcb_solid(_READOUT_PCB_SIZE_ / 2, _READOUT_PCB_SIZE_ / 2, _READOUT_PCB_THICKNESS_ / 2);
        TString pcbName;
        pcbName.Form("%s-pcb-%02d", detName.c_str(), imod);
        Volume pcbVol(pcbName.Data(), pcb_solid, HRPPD_PCBMat);
        hrppdVol_air.placeVolume(pcbVol, Position(0.0, 0.0, accu + _READOUT_PCB_THICKNESS_ / 2));
      }

      // Eventually place the whole HRPPD container volume;
      double dz =
          _FIDUCIAL_VOLUME_LENGTH_ / 2 - _SENSOR_AREA_LENGTH_ + _HRPPD_CONTAINER_VOLUME_HEIGHT_ / 2;
      pfRICH_volume.placeVolume(hrppdVol_air, Position(xy.X(), xy.Y(), dz));

#ifdef WITH_IRT2_SUPPORT
      {
        // Photocathode surface;
        double xOffset = xy.X(), yOffset = xy.Y();
        auto surface =
            new FlatSurface((1 / mm) * TVector3(sign * xOffset, yOffset,
                                                sign * (fvOffset + _FIDUCIAL_VOLUME_LENGTH_ / 2 -
                                                        _SENSOR_AREA_LENGTH_ +
                                                        _HRPPD_WINDOW_THICKNESS_ + pdthick / 2)),
                            TVector3(1 * sign, 0, 0), TVector3(0, -1, 0));

        {
          // '0': pfRICH has no division in sectors (unlike e.g. dRICH);
          unsigned sector = 0;

          for (unsigned iq = 0; iq < 4; iq++) {
            auto irt = pd->AllocateIRT(sector, sensorID);

            // Aerogel and acrylic;
            if (cdet->m_OpticalBoundaries[CherenkovDetector::Upstream].find(sector) !=
                cdet->m_OpticalBoundaries[CherenkovDetector::Upstream].end())
              for (auto boundary : cdet->m_OpticalBoundaries[CherenkovDetector::Upstream][sector])
                irt->AddOpticalBoundary(boundary);

            switch (iq) {
            case 0:
              // Direct hit;
              break;
            case 1:
            case 2:
              // Reflection on either inner or outer mirrors;
              irt->AddOpticalBoundary(mboundaries[iq - 1]);
              break;
            case 3:
              // Reflection on outer, then on inner mirror; happens at large angles; if the pyramids are
              // too high, these photons will undergo more reflections, and cannot be saved;
              irt->AddOpticalBoundary(mboundaries[1]);
              irt->AddOpticalBoundary(mboundaries[0]);
              break;
            } //switch

            // Fused silica windows;
            if (cdet->m_OpticalBoundaries[CherenkovDetector::Downstream].find(sector) !=
                cdet->m_OpticalBoundaries[CherenkovDetector::Downstream].end())
              for (auto boundary : cdet->m_OpticalBoundaries[CherenkovDetector::Downstream][sector])
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

  return sdet;
} // createDetector()

// -------------------------------------------------------------------------------------

// clang-format off
DECLARE_DETELEMENT(epic_PFRICH, createDetector)
