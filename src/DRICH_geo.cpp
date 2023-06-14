// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022, 2023 Christopher Dilks, Junhuai Xu

//==========================================================================
//  dRICH: Dual Ring Imaging Cherenkov Detector
//--------------------------------------------------------------------------
//
// Author: Christopher Dilks (Duke University)
//
// - Design Adapted from Standalone Fun4all and GEMC implementations
//   [ Evaristo Cisbani, Cristiano Fanelli, Alessio Del Dotto, et al. ]
//
//==========================================================================

#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/OpticalSurfaces.h"
#include "DD4hep/Printout.h"
#include "DDRec/DetectorData.h"
#include "DDRec/Surface.h"

#include <XML/Helper.h>

using namespace dd4hep;
using namespace dd4hep::rec;

// create the detector
static Ref_t createDetector(Detector& desc, xml::Handle_t handle, SensitiveDetector sens)
{

  xml::DetElement       detElem = handle;
  std::string           detName = detElem.nameStr();
  int                   detID   = detElem.id();
  xml::Component        dims    = detElem.dimensions();
  OpticalSurfaceManager surfMgr = desc.surfaceManager();
  DetElement            det(detName, detID);
  sens.setType("tracker");

  // attributes, from compact file =============================================
  // - vessel
  double vesselZmin      = dims.attr<double>(_Unicode(zmin));
  double vesselLength    = dims.attr<double>(_Unicode(length));
  double vesselRmin0     = dims.attr<double>(_Unicode(rmin0));
  double vesselRmin1     = dims.attr<double>(_Unicode(rmin1));
  double vesselRmax0     = dims.attr<double>(_Unicode(rmax0));
  double vesselRmax1     = dims.attr<double>(_Unicode(rmax1));
  double vesselRmax2     = dims.attr<double>(_Unicode(rmax2));
  double snoutLength     = dims.attr<double>(_Unicode(snout_length));
  int    nSectors        = dims.attr<int>(_Unicode(nsectors));
  double wallThickness   = dims.attr<double>(_Unicode(wall_thickness));
  double windowThickness = dims.attr<double>(_Unicode(window_thickness));
  auto   vesselMat       = desc.material(detElem.attr<std::string>(_Unicode(material)));
  auto   gasvolMat       = desc.material(detElem.attr<std::string>(_Unicode(gas)));
  auto   vesselVis       = desc.visAttributes(detElem.attr<std::string>(_Unicode(vis_vessel)));
  auto   gasvolVis       = desc.visAttributes(detElem.attr<std::string>(_Unicode(vis_gas)));
  // - radiator (applies to aerogel and filter)
  auto   radiatorElem       = detElem.child(_Unicode(radiator));
  double radiatorRmin       = radiatorElem.attr<double>(_Unicode(rmin));
  double radiatorRmax       = radiatorElem.attr<double>(_Unicode(rmax));
  double radiatorPitch      = radiatorElem.attr<double>(_Unicode(pitch));
  double radiatorFrontplane = radiatorElem.attr<double>(_Unicode(frontplane));
  // - aerogel
  auto   aerogelElem      = radiatorElem.child(_Unicode(aerogel));
  auto   aerogelMat       = desc.material(aerogelElem.attr<std::string>(_Unicode(material)));
  auto   aerogelVis       = desc.visAttributes(aerogelElem.attr<std::string>(_Unicode(vis)));
  double aerogelThickness = aerogelElem.attr<double>(_Unicode(thickness));
  // - filter
  auto   filterElem      = radiatorElem.child(_Unicode(filter));
  auto   filterMat       = desc.material(filterElem.attr<std::string>(_Unicode(material)));
  auto   filterVis       = desc.visAttributes(filterElem.attr<std::string>(_Unicode(vis)));
  double filterThickness = filterElem.attr<double>(_Unicode(thickness));
  // - airgap between filter and aerogel
  auto   airgapElem      = radiatorElem.child(_Unicode(airgap));
  auto   airgapMat       = desc.material(airgapElem.attr<std::string>(_Unicode(material)));
  auto   airgapVis       = desc.visAttributes(airgapElem.attr<std::string>(_Unicode(vis)));
  double airgapThickness = airgapElem.attr<double>(_Unicode(thickness));
  // - mirror
  auto   mirrorElem      = detElem.child(_Unicode(mirror));
  auto   mirrorMat       = desc.material(mirrorElem.attr<std::string>(_Unicode(material)));
  auto   mirrorVis       = desc.visAttributes(mirrorElem.attr<std::string>(_Unicode(vis)));
  auto   mirrorSurf      = surfMgr.opticalSurface(mirrorElem.attr<std::string>(_Unicode(surface)));
  double mirrorBackplane = mirrorElem.attr<double>(_Unicode(backplane));
  double mirrorThickness = mirrorElem.attr<double>(_Unicode(thickness));
  double mirrorRmin      = mirrorElem.attr<double>(_Unicode(rmin));
  double mirrorRmax      = mirrorElem.attr<double>(_Unicode(rmax));
  double mirrorPhiw      = mirrorElem.attr<double>(_Unicode(phiw));
  double focusTuneZ      = mirrorElem.attr<double>(_Unicode(focus_tune_z));
  double focusTuneX      = mirrorElem.attr<double>(_Unicode(focus_tune_x));
  // - sensor photosensitive surface (pss)
  auto   pssElem      = detElem.child(_Unicode(sensors)).child(_Unicode(pss));
  auto   pssMat       = desc.material(pssElem.attr<std::string>(_Unicode(material)));
  auto   pssVis       = desc.visAttributes(pssElem.attr<std::string>(_Unicode(vis)));
  auto   pssSurf      = surfMgr.opticalSurface(pssElem.attr<std::string>(_Unicode(surface)));
  double pssSide      = pssElem.attr<double>(_Unicode(side));
  double pssThickness = pssElem.attr<double>(_Unicode(thickness));
  // - sensor resin
  auto resinElem      = detElem.child(_Unicode(sensors)).child(_Unicode(resin));
  auto resinMat       = desc.material(resinElem.attr<std::string>(_Unicode(material)));
  auto resinVis       = desc.visAttributes(resinElem.attr<std::string>(_Unicode(vis)));
  auto resinSide      = resinElem.attr<double>(_Unicode(side));
  auto resinThickness = resinElem.attr<double>(_Unicode(thickness));
  // - photodetector unit (PDU)
  auto pduNumSensors = desc.constant<int>("DRICH_pdu_num_sensors");
  auto pduSensorGap  = desc.constant<double>("DRICH_pdu_sensor_gap");
  auto pduGap        = desc.constant<double>("DRICH_pdu_gap");
  // - sensor sphere
  auto   sensorSphElem    = detElem.child(_Unicode(sensors)).child(_Unicode(sphere));
  double sensorSphRadius  = sensorSphElem.attr<double>(_Unicode(radius));
  double sensorSphCenterX = sensorSphElem.attr<double>(_Unicode(centerx));
  double sensorSphCenterZ = sensorSphElem.attr<double>(_Unicode(centerz));
  // - sensor sphere patch cuts
  auto   sensorSphPatchElem = detElem.child(_Unicode(sensors)).child(_Unicode(sphericalpatch));
  double sensorSphPatchPhiw = sensorSphPatchElem.attr<double>(_Unicode(phiw));
  double sensorSphPatchRmin = sensorSphPatchElem.attr<double>(_Unicode(rmin));
  double sensorSphPatchRmax = sensorSphPatchElem.attr<double>(_Unicode(rmax));
  double sensorSphPatchZmin = sensorSphPatchElem.attr<double>(_Unicode(zmin));
  // - sensor readout
  auto readoutName = detElem.attr<std::string>(_Unicode(readout));
  // - settings and switches
  auto debugOpticsMode = desc.constant<int>("DRICH_debug_optics");
  bool debugSector     = desc.constant<int>("DRICH_debug_sector") == 1;
  bool debugMirror     = desc.constant<int>("DRICH_debug_mirror") == 1;
  bool debugSensors    = desc.constant<int>("DRICH_debug_sensors") == 1;

  // if debugging optics, override some settings
  bool debugOptics = debugOpticsMode > 0;
  if (debugOptics) {
    printout(WARNING, "DRICH_geo", "DEBUGGING DRICH OPTICS");
    switch (debugOpticsMode) {
    case 1:
      vesselMat = aerogelMat = filterMat = pssMat = gasvolMat = desc.material("VacuumOptical");
      break;
    case 2:
      vesselMat = aerogelMat = filterMat = pssMat = desc.material("VacuumOptical");
      break;
    case 3:
      vesselMat = aerogelMat = filterMat = gasvolMat = desc.material("VacuumOptical");
      break;
    default:
      printout(FATAL, "DRICH_geo", "UNKNOWN debugOpticsMode");
      return det;
    }
  }

  // if debugging anything, draw only one sector and adjust visibility
  if (debugOptics || debugMirror || debugSensors)
    debugSector = true;
  if (debugSector)
    gasvolVis = vesselVis = desc.invisible();

  // readout coder <-> unique sensor ID
  /* - `sensorIDfields` is a list of readout fields used to specify a unique sensor ID
   * - `cellMask` is defined such that a hit's `cellID & cellMask` is the corresponding sensor's unique ID
   */
  std::vector<std::string> sensorIDfields = {"module", "sector"};
  const auto&              readoutCoder   = *desc.readout(readoutName).idSpec().decoder();
  // determine `cellMask` based on `sensorIDfields`
  uint64_t cellMask = 0;
  for (const auto& idField : sensorIDfields)
    cellMask |= readoutCoder[idField].mask();
  desc.add(Constant("DRICH_cell_mask", std::to_string(cellMask)));
  // create a unique sensor ID from a sensor's PlacedVolume::volIDs
  auto encodeSensorID = [&readoutCoder](auto ids) {
    uint64_t enc = 0;
    for (const auto& [idField, idValue] : ids)
      enc |= uint64_t(idValue) << readoutCoder[idField].offset();
    return enc;
  };

  // BUILD VESSEL ====================================================================
  /* - `vessel`: aluminum enclosure, the mother volume of the dRICH
   * - `gasvol`: gas volume, which fills `vessel`; all other volumes defined below
   *   are children of `gasvol`
   * - the dRICH vessel geometry has two regions: the snout refers to the conic region
   *   in the front, housing the aerogel, while the tank refers to the cylindrical
   *   region, housing the rest of the detector components
   */

  // derived attributes
  double tankLength = vesselLength - snoutLength;
  double vesselZmax = vesselZmin + vesselLength;

  // snout solids
  double boreDelta  = vesselRmin1 - vesselRmin0;
  double snoutDelta = vesselRmax1 - vesselRmax0;
  Cone   vesselSnout(snoutLength / 2.0, vesselRmin0, vesselRmax0, vesselRmin0 + boreDelta * snoutLength / vesselLength,
                   vesselRmax1);
  Cone   gasvolSnout(
      /* note: `gasvolSnout` extends a bit into the tank, so it touches `gasvolTank`
       * - the extension distance is equal to the tank `windowThickness`, so the
       *   length of `gasvolSnout` == length of `vesselSnout`
       * - the extension backplane radius is calculated using similar triangles
       */
      snoutLength / 2.0, vesselRmin0 + wallThickness, vesselRmax0 - wallThickness,
      vesselRmin0 + boreDelta * (snoutLength - windowThickness) / vesselLength + wallThickness,
      vesselRmax1 - wallThickness + windowThickness * (vesselRmax1 - vesselRmax0) / snoutLength);

  // tank solids
  Cone vesselTank(tankLength / 2.0, vesselSnout.rMin2(), vesselRmax2, vesselRmin1, vesselRmax2);
  Cone gasvolTank(tankLength / 2.0 - windowThickness, gasvolSnout.rMin2(), vesselRmax2 - wallThickness,
                  vesselRmin1 + wallThickness, vesselRmax2 - wallThickness);

  // snout + tank solids
  UnionSolid vesselUnion(vesselTank, vesselSnout, Position(0., 0., -vesselLength / 2.));
  UnionSolid gasvolUnion(gasvolTank, gasvolSnout, Position(0., 0., -vesselLength / 2. + windowThickness));

  //  extra solids for `debugOptics` only
  Box vesselBox(1001, 1001, 1001);
  Box gasvolBox(1000, 1000, 1000);

  // choose vessel and gasvol solids (depending on `debugOpticsMode` (0=disabled))
  Solid vesselSolid, gasvolSolid;
  switch (debugOpticsMode) {
  case 0:
    vesselSolid = vesselUnion;
    gasvolSolid = gasvolUnion;
    break; // `!debugOptics`
  case 1:
  case 3:
    vesselSolid = vesselBox;
    gasvolSolid = gasvolBox;
    break;
  case 2:
    vesselSolid = vesselBox;
    gasvolSolid = gasvolUnion;
    break;
  }

  // volumes
  Volume vesselVol(detName, vesselSolid, vesselMat);
  Volume gasvolVol(detName + "_gas", gasvolSolid, gasvolMat);
  vesselVol.setVisAttributes(vesselVis);
  gasvolVol.setVisAttributes(gasvolVis);

  // reference positions
  // - the vessel is created such that the center of the cylindrical tank volume
  //   coincides with the origin; this is called the "origin position" of the vessel
  // - when the vessel (and its children volumes) is placed, it is translated in
  //   the z-direction to be in the proper EPIC-integration location
  // - these reference positions are for the frontplane and backplane of the vessel,
  //   with respect to the vessel origin position
  auto originFront = Position(0., 0., -tankLength / 2.0 - snoutLength);
  // auto originBack  = Position(0., 0., tankLength / 2.0);
  auto vesselPos = Position(0, 0, vesselZmin) - originFront;

  // place gas volume
  PlacedVolume gasvolPV = vesselVol.placeVolume(gasvolVol, Position(0, 0, 0));
  DetElement   gasvolDE(det, "gasvol_de", 0);
  gasvolDE.setPlacement(gasvolPV);

  // place mother volume (vessel)
  Volume       motherVol = desc.pickMotherVolume(det);
  PlacedVolume vesselPV  = motherVol.placeVolume(vesselVol, vesselPos);
  vesselPV.addPhysVolID("system", detID);
  det.setPlacement(vesselPV);

  // BUILD RADIATOR ====================================================================

  // solid and volume: create aerogel and filter
  Cone aerogelSolid(aerogelThickness / 2, radiatorRmin, radiatorRmax,
                    radiatorRmin + boreDelta * aerogelThickness / vesselLength,
                    radiatorRmax + snoutDelta * aerogelThickness / snoutLength);
  Cone airgapSolid(airgapThickness / 2, radiatorRmin + boreDelta * aerogelThickness / vesselLength,
                   radiatorRmax + snoutDelta * aerogelThickness / snoutLength,
                   radiatorRmin + boreDelta * (aerogelThickness + airgapThickness) / vesselLength,
                   radiatorRmax + snoutDelta * (aerogelThickness + airgapThickness) / snoutLength);
  Cone filterSolid(filterThickness / 2, radiatorRmin + boreDelta * (aerogelThickness + airgapThickness) / vesselLength,
                   radiatorRmax + snoutDelta * (aerogelThickness + airgapThickness) / snoutLength,
                   radiatorRmin + boreDelta * (aerogelThickness + airgapThickness + filterThickness) / vesselLength,
                   radiatorRmax + snoutDelta * (aerogelThickness + airgapThickness + filterThickness) / snoutLength);

  Volume aerogelVol(detName + "_aerogel", aerogelSolid, aerogelMat);
  Volume airgapVol(detName + "_airgap", airgapSolid, airgapMat);
  Volume filterVol(detName + "_filter", filterSolid, filterMat);
  aerogelVol.setVisAttributes(aerogelVis);
  airgapVol.setVisAttributes(airgapVis);
  filterVol.setVisAttributes(filterVis);

  // aerogel placement and surface properties
  // TODO [low-priority]: define skin properties for aerogel and filter
  // FIXME: radiatorPitch might not be working correctly (not yet used)
  auto radiatorPos      = Position(0., 0., radiatorFrontplane + 0.5 * aerogelThickness) + originFront;
  auto aerogelPlacement = Translation3D(radiatorPos) * // re-center to originFront
                          RotationY(radiatorPitch); // change polar angle to specified pitch
  auto       aerogelPV = gasvolVol.placeVolume(aerogelVol, aerogelPlacement);
  DetElement aerogelDE(det, "aerogel_de", 0);
  aerogelDE.setPlacement(aerogelPV);
  // SkinSurface aerogelSkin(desc, aerogelDE, "mirror_optical_surface", aerogelSurf, aerogelVol);
  // aerogelSkin.isValid();

  // airgap and filter placement and surface properties
  if (!debugOptics) {

    auto airgapPlacement =
        Translation3D(radiatorPos) *                                      // re-center to originFront
        RotationY(radiatorPitch) *                                        // change polar angle
        Translation3D(0., 0., (aerogelThickness + airgapThickness) / 2.); // move to aerogel backplane
    auto airgapPV = gasvolVol.placeVolume(airgapVol, airgapPlacement);
    DetElement airgapDE(det, "airgap_de", 0);
    airgapDE.setPlacement(airgapPV);

    auto filterPlacement =
        Translation3D(0., 0., airgapThickness) *                           // add an air gap
        Translation3D(radiatorPos) *                                       // re-center to originFront
        RotationY(radiatorPitch) *                                         // change polar angle
        Translation3D(0., 0., (aerogelThickness + filterThickness) / 2.);  // move to aerogel backplane
    auto filterPV = gasvolVol.placeVolume(filterVol, filterPlacement);
    DetElement filterDE(det, "filter_de", 0);
    filterDE.setPlacement(filterPV);
    // SkinSurface filterSkin(desc, filterDE, "mirror_optical_surface", filterSurf, filterVol);
    // filterSkin.isValid();

    // radiator z-positions (w.r.t. IP); only needed downstream if !debugOptics
    double aerogelZpos = vesselPos.z() + aerogelPV.position().z();
    double airgapZpos  = vesselPos.z() + airgapPV.position().z();
    double filterZpos  = vesselPos.z() + filterPV.position().z();
    desc.add(Constant("DRICH_aerogel_zpos", std::to_string(aerogelZpos)));
    desc.add(Constant("DRICH_airgap_zpos", std::to_string(airgapZpos)));
    desc.add(Constant("DRICH_filter_zpos", std::to_string(filterZpos)));
  }

  // radiator material names
  desc.add(Constant("DRICH_aerogel_material", aerogelMat.ptr()->GetName(), "string"));
  desc.add(Constant("DRICH_airgap_material", airgapMat.ptr()->GetName(), "string"));
  desc.add(Constant("DRICH_filter_material", filterMat.ptr()->GetName(), "string"));
  desc.add(Constant("DRICH_gasvol_material", gasvolMat.ptr()->GetName(), "string"));

  // SECTOR LOOP //////////////////////////////////////////////////////////////////////
  for (int isec = 0; isec < nSectors; isec++) {

    // debugging filters, limiting the number of sectors
    if (debugSector && isec != 0)
      continue;

    // sector rotation about z axis
    RotationZ   sectorRotation(isec * 2 * M_PI / nSectors);
    std::string secName = "sec" + std::to_string(isec);

    // BUILD MIRRORS ====================================================================

    // mirror positioning attributes
    // - sensor sphere center, w.r.t. IP
    double zS = sensorSphCenterZ + vesselZmin;
    double xS = sensorSphCenterX;
    // - distance between IP and mirror back plane
    double b = vesselZmax - mirrorBackplane;
    // - desired focal region: sensor sphere center, offset by focus-tune (z,x) parameters
    double zF = zS + focusTuneZ;
    double xF = xS + focusTuneX;

    // determine the mirror that focuses the IP to this desired region
    /* - uses point-to-point focusing to derive spherical mirror center
     *   `(mirrorCenterZ,mirrorCenterX)` and radius `mirrorRadius` for given
     *   image point coordinates `(zF,xF)` and `b`, defined as the z-distance
     *   between the object (IP) and the mirror surface
     * - all coordinates are specified w.r.t. the object point (IP)
     */
    double mirrorCenterZ = b * zF / (2 * b - zF);
    double mirrorCenterX = b * xF / (2 * b - zF);
    double mirrorRadius  = b - mirrorCenterZ;

    // translate mirror center to be w.r.t vessel front plane
    mirrorCenterZ -= vesselZmin;

    // spherical mirror patch cuts and rotation
    double mirrorThetaRot = std::asin(mirrorCenterX / mirrorRadius);
    double mirrorTheta1   = mirrorThetaRot - std::asin((mirrorCenterX - mirrorRmin) / mirrorRadius);
    double mirrorTheta2   = mirrorThetaRot + std::asin((mirrorRmax - mirrorCenterX) / mirrorRadius);

    // if debugging, draw full sphere
    if (debugMirror) {
      mirrorTheta1 = 0;
      mirrorTheta2 = M_PI; /*mirrorPhiw=2*M_PI;*/
    }

    // solid : create sphere at origin, with specified angular limits;
    // phi limits are increased to fill gaps (overlaps are cut away later)
    Sphere mirrorSolid1(mirrorRadius, mirrorRadius + mirrorThickness, mirrorTheta1, mirrorTheta2, -40 * degree,
                        40 * degree);

    // mirror placement transformation (note: transformations are in reverse order)
    auto mirrorPos = Position(mirrorCenterX, 0., mirrorCenterZ) + originFront;
    auto mirrorPlacement(Translation3D(mirrorPos) * // re-center to specified position
                         RotationY(-mirrorThetaRot) // rotate about vertical axis, to be within vessel radial walls
    );

    // cut overlaps with other sectors using "pie slice" wedges, to the extent specified
    // by `mirrorPhiw`
    Tube              pieSlice(0.01 * cm, vesselRmax2, tankLength / 2.0, -mirrorPhiw / 2.0, mirrorPhiw / 2.0);
    IntersectionSolid mirrorSolid2(pieSlice, mirrorSolid1, mirrorPlacement);

    // mirror volume, attributes, and placement
    Volume mirrorVol(detName + "_mirror_" + secName, mirrorSolid2, mirrorMat);
    mirrorVol.setVisAttributes(mirrorVis);
    auto mirrorSectorPlacement = Transform3D(sectorRotation); // rotate about beam axis to sector
    auto mirrorPV              = gasvolVol.placeVolume(mirrorVol, mirrorSectorPlacement);

    // properties
    DetElement mirrorDE(det, "mirror_de_" + secName, isec);
    mirrorDE.setPlacement(mirrorPV);
    SkinSurface mirrorSkin(desc, mirrorDE, "mirror_optical_surface_" + secName, mirrorSurf, mirrorVol);
    mirrorSkin.isValid();

    // reconstruction constants (w.r.t. IP)
    // - access sector center after `sectorRotation`
    auto mirrorFinalPlacement = mirrorSectorPlacement * mirrorPlacement;
    auto mirrorFinalCenter    = vesselPos + mirrorFinalPlacement.Translation().Vect();
    desc.add(Constant("DRICH_mirror_center_x_" + secName, std::to_string(mirrorFinalCenter.x())));
    desc.add(Constant("DRICH_mirror_center_y_" + secName, std::to_string(mirrorFinalCenter.y())));
    desc.add(Constant("DRICH_mirror_center_z_" + secName, std::to_string(mirrorFinalCenter.z())));
    if (isec == 0)
      desc.add(Constant("DRICH_mirror_radius", std::to_string(mirrorRadius)));

    // BUILD SENSORS ====================================================================

    // if debugging sphere properties, restrict number of sensors drawn
    if (debugSensors) {
      pssSide = 2 * M_PI * sensorSphRadius / 64;
    }

    // reconstruction constants
    auto sensorSphPos         = Position(sensorSphCenterX, 0., sensorSphCenterZ) + originFront;
    auto sensorSphFinalCenter = sectorRotation * Position(xS, 0.0, zS);
    desc.add(Constant("DRICH_sensor_sph_center_x_" + secName, std::to_string(sensorSphFinalCenter.x())));
    desc.add(Constant("DRICH_sensor_sph_center_y_" + secName, std::to_string(sensorSphFinalCenter.y())));
    desc.add(Constant("DRICH_sensor_sph_center_z_" + secName, std::to_string(sensorSphFinalCenter.z())));
    if (isec == 0)
      desc.add(Constant("DRICH_sensor_sph_radius", std::to_string(sensorSphRadius)));

    // SENSOR MODULE LOOP ------------------------
    /* ALGORITHM: generate sphere of positions
     * - NOTE: there are two coordinate systems here:
     *   - "global" the main EPIC coordinate system
     *   - "generator" (vars end in `Gen`) is a local coordinate system for
     *     generating points on a sphere; it is related to the global system by
     *     a rotation; we do this so the "patch" (subset of generated
     *     positions) of sensors we choose to build is near the equator, where
     *     point distribution is more uniform
     * - PROCEDURE: loop over `thetaGen`, with subloop over `phiGen`, each divided evenly
     *   - the number of points to generate depends how many PDUs
     *     can fit within each ring of constant `thetaGen` or `phiGen`
     *   - we divide the relevant circumference by the PDU size, and this
     *     number is allowed to be a fraction, because likely we don't care about
     *     generating a full sphere and don't mind a "seam" at the overlap point
     *   - if we pick a patch of the sphere near the equator, and not near
     *     the poles or seam, the sensor distribution will appear uniform
     */

    // initialize module number for this sector
    int imod = 0;

    // calculate PDU pitch: the distance between two adjacent PDUs
    double pduPitch = pduNumSensors * resinSide + (pduNumSensors + 1) * pduSensorGap + pduGap;

    // thetaGen loop: iterate less than "0.5 circumference / sensor size" times
    double nTheta = M_PI * sensorSphRadius / pduPitch;
    for (int t = 0; t < (int)(nTheta + 0.5); t++) {
      double thetaGen = t / ((double)nTheta) * M_PI;

      // phiGen loop: iterate less than "circumference at this latitude / sensor size" times
      double nPhi = 2 * M_PI * sensorSphRadius * std::sin(thetaGen) / pduPitch;
      for (int p = 0; p < (int)(nPhi + 0.5); p++) {
        double phiGen = p / ((double)nPhi) * 2 * M_PI - M_PI; // shift to [-pi,pi]

        // determine global phi and theta
        // - convert {radius,thetaGen,phiGen} -> {xGen,yGen,zGen}
        double xGen = sensorSphRadius * std::sin(thetaGen) * std::cos(phiGen);
        double yGen = sensorSphRadius * std::sin(thetaGen) * std::sin(phiGen);
        double zGen = sensorSphRadius * std::cos(thetaGen);
        // - convert {xGen,yGen,zGen} -> global {x,y,z} via rotation
        double x = zGen;
        double y = xGen;
        double z = yGen;
        // - convert global {x,y,z} -> global {phi,theta}
        // double phi   = std::atan2(y, x);
        // double theta = std::acos(z / sensorSphRadius);

        // shift global coordinates so we can apply spherical patch cuts
        double zCheck   = z + sensorSphCenterZ;
        double xCheck   = x + sensorSphCenterX;
        double yCheck   = y;
        double rCheck   = std::hypot(xCheck, yCheck);
        double phiCheck = std::atan2(yCheck, xCheck);

        // patch cut
        bool patchCut = std::fabs(phiCheck) < sensorSphPatchPhiw && zCheck > sensorSphPatchZmin &&
                        rCheck > sensorSphPatchRmin && rCheck < sensorSphPatchRmax;
        if (debugSensors)
          patchCut = std::fabs(phiCheck) < sensorSphPatchPhiw;
        if (patchCut) {

          // sensor assembly: collection of all objects for a single SiPM + services
          /* - coordinate system: the "origin" of the assembly will be the center of the
           *   outermost surface of the photosensitive surface (pss)
           *   - reconstruction can access the sensor surface position from the sensor
           *     assembly origin, which will ultimately have coordinates w.r.t. to the IP after
           *     placement in the dRICH vessel
           *   - the pss is segmented into SiPM pixels; gaps between the pixels
           *     are accounted for in reconstruction, and each pixel reads out as a unique `cellID`
           *   - `cellID` to postion conversion will give pixel centroids within the pss volume,
           *     (not exactly at the pss surface, but rather in the center of the pss volume,
           *     so keep in mind the very small offset)
           *
           * sensor assembly diagram:, where '0' denotes the origin:
           *
           *                                 axes:  z
           *    +-+--------0--------+-+             |
           *    | |       pss       | |             0--x
           *    | +-----------------+ |
           *    |        resin        |
           *    +---------------------+
           *
           */
          Assembly sensorAssembly(detName + "_sensor_" + secName);

          // photosensitive surface (pss) and resin solids
          Box pssSolid(pssSide / 2., pssSide / 2., pssThickness / 2.);
          Box resinSolid(resinSide / 2., resinSide / 2., resinThickness / 2.);

          // embed pss solid in resin solid, by subtracting `pssSolid` from `resinSolid`
          SubtractionSolid resinSolidEmbedded(resinSolid, pssSolid,
              Transform3D(Translation3D(0., 0., (resinThickness - pssThickness) / 2. )));

          /* NOTE:
           * Here we could add gaps (size=`DRICH_pixel_gap`) between the pixels
           * as additional resin volumes, but this would require several more
           * iterative boolean operations, which may cause significant
           * performance slow downs in the simulation. Alternatively, one can
           * create a pixel gap mask with several disjoint, thin `Box` volumes
           * just outside the pss surface (no booleans required), but this
           * would amount to a very large number of additional volumes. Instead,
           * we have decided to apply pixel gap masking to the digitization
           * algorithm, downstream in reconstruction.
           */

          // pss and resin volumes
          Volume pssVol(detName + "_pss_" + secName, pssSolid, pssMat);
          Volume resinVol(detName + "_resin_" + secName, resinSolidEmbedded, resinMat);
          pssVol.setVisAttributes(pssVis);
          resinVol.setVisAttributes(resinVis);

          // sensitivity
          if (!debugOptics || debugOpticsMode == 3)
            pssVol.setSensitiveDetector(sens);

          // placement of objects in `sensorAssembly`
          // clang-format off
          auto pssPV = sensorAssembly.placeVolume(
              pssVol,
              Transform3D(Translation3D(0., 0., -pssThickness / 2.0)) // set assembly origin to pss outermost surface centroid
              );
          sensorAssembly.placeVolume(
              resinVol,
              Transform3D(Translation3D(0., 0., -resinThickness / 2.0))
              );
          // clang-format on

          /* photodetector unit (PDU) assembly
           *
           *    Top view: 2x2 matrix of SiPMs (2 PDU units shown side-by-side)
           *    =============================
           *    
           *              ->:  :<- PDU gap size
           *                :  :
           *    +-----------+  +-----------+
           *    | +--+ +--+ |  | +--+ +--+ |
           *    | |  | |  | |  | |  | |  | |
           *    | +--+ +--+ |  | +--+ +--+ |
           *    | +--+ +--+ |  | +--+ +--+ |
           *    | |  | |  | |  | |  | |  | |
           *    | +--+ +--+ |  | +--+ +--+ |
           *    +-----------+  +-----------+
           *         : :  : :
           *         : :->: :<- sensor gap size
           *       ->: :<- sensor gap size (same)
           *
           *     Side view:
           *     ==========
           *
           *     +--------------------+
           *     |     SiPM Matrix    |
           *     +--------------------+
           *     |   Cooling, heat    |
           *     |   exchange, etc.   |
           *     +--------------------+
           *       ||  ||  ||  ||  ||
           *       ||  ||  ||  ||  || front-end and
           *       ||  ||  ||  ||  || readout boards
           *       ||  ||  ||  ||  ||
           *       ||  ||  ||  ||  ||
           *               :   :
           *             ->:   :<- offset of a board
           */

          Assembly pduAssembly(detName + "_pdu_" + secName);

          // generate matrix of sensors
          double pduSensorPitch     = resinSide + pduSensorGap;
          double pduSensorOffsetMax = pduSensorPitch * (pduNumSensors - 1) / 2.0;
          for(int sensorIx = 0; sensorIx < pduNumSensors; sensorIx++) {
            for(int sensorIy = 0; sensorIy < pduNumSensors; sensorIy++) {
              pduAssembly.placeVolume(
                  sensorAssembly,
                  Transform3D(Translation3D(
                      sensorIx * pduSensorPitch - pduSensorOffsetMax,
                      sensorIy * pduSensorPitch - pduSensorOffsetMax,
                      -resinThickness / 2.0
                      ))
                  );
            }
          }

          // place PDU assembly
          /* - transformations operate on global coordinates; the corresponding
           *   generator coordinates are provided in the comments
           * - transformations are applied in reverse order
           */
          // clang-format off
          auto pduAssemblyPlacement =
              sectorRotation *                               // rotate about beam axis to sector
              Translation3D(sensorSphPos) *                  // move sphere to reference position
              RotationX(phiGen) *                            // rotate about `zGen`
              RotationZ(thetaGen) *                          // rotate about `yGen`
              Translation3D(-resinThickness / 2.0, 0., 0.) * // pull back so sensor active surface is at spherical surface
              Translation3D(sensorSphRadius, 0., 0.) *       // push radially to spherical surface
              RotationY(M_PI / 2) *                          // rotate sensor to be compatible with generator coords
              RotationZ(-M_PI / 2);                          // correction for readout segmentation mapping
          auto pduPV = gasvolVol.placeVolume(pduAssembly, pduAssemblyPlacement);
          // clang-format on

          // generate LUT for module number -> sensor position, for readout mapping tests
          // if(isec==0) printf("%d %f %f\n",imod,pssPV.position().x(),pssPV.position().y());

          // sensor readout and DetElements
          pssPV.addPhysVolID("sector", isec).addPhysVolID("module", imod); // NOTE: follow `sensorIDfields`
          auto        imodsec    = encodeSensorID(pssPV.volIDs());
          std::string modsecName = secName + "_" + std::to_string(imod);
          DetElement pduDE(det, "sensor_de_" + modsecName, imodsec);
          pduDE.setPlacement(pduPV);

          // sensor surface properties
          if (!debugOptics || debugOpticsMode == 3) {
            DetElement pssDE(det, "pss_de_" + modsecName, imodsec);
            SkinSurface pssSkin(desc, pssDE, "sensor_optical_surface_" + modsecName, pssSurf, pssVol);
            pssSkin.isValid();
          }

          // increment sensor module number
          imod++;

        } // end patch cuts
      }   // end phiGen loop
    }     // end thetaGen loop

    // END SENSOR MODULE LOOP ------------------------

    // add constant for access to the number of modules per sector
    if (isec == 0)
      desc.add(Constant("DRICH_num_sensors", std::to_string(imod)));
    else if (imod != desc.constant<int>("DRICH_num_sensors"))
      printout(WARNING, "DRICH_geo", "number of sensors is not the same for each sector");

  } // END SECTOR LOOP //////////////////////////

  return det;
}

// clang-format off
DECLARE_DETELEMENT(epic_DRICH, createDetector)
