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

  // attributes -----------------------------------------------------------
  // - vessel
  double vesselZmin      = dims.attr<double>(_Unicode(zmin));
  double vesselLength    = dims.attr<double>(_Unicode(length));
  double vesselRmin0     = dims.attr<double>(_Unicode(rmin0));
  double vesselRmin1     = dims.attr<double>(_Unicode(rmin1));
  double vesselRmax0     = dims.attr<double>(_Unicode(rmax0));
  double vesselRmax1     = dims.attr<double>(_Unicode(rmax1));
  double proximityGap    = dims.attr<double>(_Unicode(proximity_gap));
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
  // - airgap between filter and aerogel // TODO: use these to place an airgap volume
  auto airgapElem = radiatorElem.child(_Unicode(airgap));
  // auto   airgapMat       = desc.material(airgapElem.attr<std::string>(_Unicode(material))); // TODO
  // auto   airgapVis       = desc.visAttributes(airgapElem.attr<std::string>(_Unicode(vis))); // TODO
  double airgapThickness = airgapElem.attr<double>(_Unicode(thickness));
  // - sensor module
  auto   sensorElem      = detElem.child(_Unicode(sensors)).child(_Unicode(module));
  auto   sensorMat       = desc.material(sensorElem.attr<std::string>(_Unicode(material)));
  auto   sensorVis       = desc.visAttributes(sensorElem.attr<std::string>(_Unicode(vis)));
  auto   sensorSurf      = surfMgr.opticalSurface(sensorElem.attr<std::string>(_Unicode(surface)));
  double sensorSide      = sensorElem.attr<double>(_Unicode(side));
  double sensorGap       = sensorElem.attr<double>(_Unicode(gap));
  double sensorThickness = sensorElem.attr<double>(_Unicode(thickness));
  auto   readoutName     = detElem.attr<std::string>(_Unicode(readout));
  // - sensor plane
  auto   sensorPlaneElem = detElem.child(_Unicode(sensors)).child(_Unicode(plane));
  double sensorPlaneRmin = sensorPlaneElem.attr<double>(_Unicode(rmin));
  double sensorPlaneRmax = sensorPlaneElem.attr<double>(_Unicode(rmax));
  // - debugging switches
  long debug_optics_mode = desc.constantAsLong("PFRICH_debug_optics");

  // if debugging optics, override some settings
  bool debug_optics = debug_optics_mode > 0;
  if (debug_optics) {
    printout(WARNING, "PFRICH_geo", "DEBUGGING PFRICH OPTICS");
    switch (debug_optics_mode) {
    case 1:
      vesselMat = aerogelMat = filterMat = sensorMat = gasvolMat = desc.material("VacuumOptical");
      break;
    case 2:
      vesselMat = aerogelMat = filterMat = sensorMat = desc.material("VacuumOptical");
      break;
    default:
      printout(FATAL, "PFRICH_geo", "UNKNOWN debug_optics_mode");
      return det;
    };
    aerogelVis = sensorVis;
    gasvolVis = vesselVis = desc.invisible();
  };

  // readout coder <-> unique sensor ID
  /* - `sensorIDfields` is a list of readout fields used to specify a unique sensor ID
   * - `cellMask` is defined such that a hit's `cellID & cellMask` is the corresponding sensor's unique ID
   * - this redundant generalization is for future flexibility, and consistency with dRICH
   */
  std::vector<std::string> sensorIDfields = {"module"};
  const auto&              readoutCoder   = *desc.readout(readoutName).idSpec().decoder();
  // determine `cellMask` based on `sensorIDfields`
  uint64_t cellMask = 0;
  for (const auto& idField : sensorIDfields)
    cellMask |= readoutCoder[idField].mask();
  // create a unique sensor ID from a sensor's PlacedVolume::volIDs
  auto encodeSensorID = [&readoutCoder](auto ids) {
    uint64_t enc = 0;
    for (const auto& [idField, idValue] : ids)
      enc |= uint64_t(idValue) << readoutCoder[idField].offset();
    return enc;
  };

  // define reconstruction geometry constants `PFRICH_RECON_*`
  /* - these are the numbers needed to rebuild the geometry in the
   *   reconstruction, in particular, the optical surfaces encountered by the
   *   Cherenkov photons
   * - positions are w.r.t. the IP
   * - check the values of all of the `PFRICH_RECON_*` constants after any change
   *   to the geometry
   * - some `PFRICH_RECON_*` constants are redundant, but are defined to make
   *   it clear that the reconstruction code depends on them
   */
  desc.add(Constant("PFRICH_RECON_zmin", std::to_string(vesselZmin)));
  desc.add(Constant("PFRICH_RECON_gasvolMaterial", gasvolMat.ptr()->GetName(), "string"));
  desc.add(Constant("PFRICH_RECON_cellMask", std::to_string(cellMask)));
  desc.add(Constant("PFRICH_RECON_sensorThickness", std::to_string(sensorThickness)));

  // BUILD VESSEL //////////////////////////////////////
  /* - `vessel`: aluminum enclosure, the mother volume of the pfRICH
   * - `gasvol`: gas volume, which fills `vessel`; all other volumes defined below
   *   are children of `gasvol`
   */

  // tank solids
  double boreDelta = vesselRmin1 - vesselRmin0;
  Cone   vesselTank(vesselLength / 2.0, vesselRmin1, vesselRmax1, vesselRmin0, vesselRmax0);
  Cone   gasvolTank(vesselLength / 2.0 - windowThickness, vesselRmin1 + wallThickness, vesselRmax1 - wallThickness,
                    vesselRmin0 + wallThickness, vesselRmax0 - wallThickness);

  //  extra solids for `debug_optics` only
  Box vesselBox(1001, 1001, 1001);
  Box gasvolBox(1000, 1000, 1000);

  // choose vessel and gasvol solids (depending on `debug_optics_mode` (0=disabled))
  Solid vesselSolid, gasvolSolid;
  switch (debug_optics_mode) {
  case 0:
    vesselSolid = vesselTank;
    gasvolSolid = gasvolTank;
    break; // `!debug_optics`
  case 1:
    vesselSolid = vesselBox;
    gasvolSolid = gasvolBox;
    break;
  case 2:
    vesselSolid = vesselBox;
    gasvolSolid = gasvolTank;
    break;
  };

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
  auto originFront = Position(0., 0., vesselLength / 2.0);
  // auto originBack = Position(0., 0., -vesselLength / 2.0);
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
  aerogelVol.setVisAttributes(aerogelVis);
  filterVol.setVisAttributes(filterVis);

  // aerogel placement and surface properties
  // TODO [low-priority]: define skin properties for aerogel and filter
  // FIXME: radiatorPitch might not be working correctly (not yet used)
  auto radiatorPos      = Position(0., 0., radiatorFrontplane - 0.5 * aerogelThickness) + originFront;
  auto aerogelPlacement = Translation3D(radiatorPos.x(), radiatorPos.y(), radiatorPos.z()) * // re-center to originFront
                          RotationY(radiatorPitch); // change polar angle to specified pitch
  auto       aerogelPV = gasvolVol.placeVolume(aerogelVol, aerogelPlacement);
  DetElement aerogelDE(det, "aerogel_de", 0);
  aerogelDE.setPlacement(aerogelPV);
  // SkinSurface aerogelSkin(desc, aerogelDE, "mirror_optical_surface", aerogelSurf, aerogelVol);
  // aerogelSkin.isValid();

  // filter placement and surface properties
  PlacedVolume filterPV;
  if (!debug_optics) {
    auto filterPlacement =
        Translation3D(0., 0., -airgapThickness) *                          // add an airgap
        Translation3D(radiatorPos.x(), radiatorPos.y(), radiatorPos.z()) * // re-center to originFront
        RotationY(radiatorPitch) *                                         // change polar angle
        Translation3D(0., 0., -(aerogelThickness + filterThickness) / 2.); // move to aerogel backplane
    filterPV = gasvolVol.placeVolume(filterVol, filterPlacement);
    DetElement filterDE(det, "filter_de", 0);
    filterDE.setPlacement(filterPV);
    // SkinSurface filterSkin(desc, filterDE, "mirror_optical_surface", filterSurf, filterVol);
    // filterSkin.isValid();
  };

  // reconstruction constants (w.r.t. IP)
  double aerogelZpos = vesselPos.z() + aerogelPV.position().z();
  double filterZpos  = vesselPos.z() + filterPV.position().z();
  desc.add(Constant("PFRICH_RECON_aerogelZpos", std::to_string(aerogelZpos)));
  desc.add(Constant("PFRICH_RECON_aerogelThickness", std::to_string(aerogelThickness)));
  desc.add(Constant("PFRICH_RECON_aerogelMaterial", aerogelMat.ptr()->GetName(), "string"));
  desc.add(Constant("PFRICH_RECON_filterZpos", std::to_string(filterZpos)));
  desc.add(Constant("PFRICH_RECON_filterThickness", std::to_string(filterThickness)));
  desc.add(Constant("PFRICH_RECON_filterMaterial", filterMat.ptr()->GetName(), "string"));

  // BUILD SENSORS ///////////////////////

  // solid and volume: single sensor module
  Box    sensorSolid(sensorSide / 2., sensorSide / 2., sensorThickness / 2.);
  Volume sensorVol(detName + "_sensor", sensorSolid, sensorMat);
  sensorVol.setVisAttributes(sensorVis);

  // sensitivity
  if (!debug_optics)
    sensorVol.setSensitiveDetector(sens);

  // sensor plane positioning: we want `proximityGap` to be the distance between the
  // aerogel backplane (i.e., aerogel/filter boundary) and the sensor active surface (e.g, photocathode)
  double sensorZpos     = radiatorFrontplane - aerogelThickness - proximityGap - 0.5 * sensorThickness;
  auto   sensorPlanePos = Position(0., 0., sensorZpos) + originFront; // reference position
  // miscellaneous
  int    imod    = 0;           // module number
  double tBoxMax = vesselRmax1; // sensors will be tiled in tBox, within annular limits

  // SENSOR MODULE LOOP ------------------------
  /* cartesian tiling loop
   * - start at (x=0,y=0), to center the grid
   * - loop over positive-x positions; for each, place the corresponding negative-x sensor too
   * - nested similar loop over y positions
   */
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
              Translation3D(sensorPlanePos.x(), sensorPlanePos.y(), sensorPlanePos.z()) * // move to reference position
              Translation3D(sx, sy, 0.)                                                   // move to grid position
          );
          auto sensorPV = gasvolVol.placeVolume(sensorVol, sensorPlacement);

          // generate LUT for module number -> sensor position, for readout mapping tests
          // printf("%d %f %f\n",imod,sensorPV.position().x(),sensorPV.position().y());

          // properties
          sensorPV.addPhysVolID("module", imod); // NOTE: must be consistent with `sensorIDfields`
          auto       imodEnc = encodeSensorID(sensorPV.volIDs());
          DetElement sensorDE(det, "sensor_de_" + std::to_string(imod), imodEnc);
          sensorDE.setPlacement(sensorPV);
          if (!debug_optics) {
            SkinSurface sensorSkin(desc, sensorDE, "sensor_optical_surface_" + std::to_string(imod), sensorSurf,
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
  //
  // Add service material if desired
  if (detElem.child("sensors").hasChild(_Unicode(services))) {
    xml_comp_t x_service = detElem.child("sensors").child(_Unicode(services));
    Assembly   service_vol("services");
    service_vol.setVisAttributes(desc, x_service.visStr());

    // Compute service total thickness from components
    double total_thickness = 0;
    for (xml_coll_t ci(x_service, _Unicode(component)); ci; ++ci) {
      total_thickness += xml_comp_t(ci).thickness();
    }

    int    ncomponents   = 0;
    double thickness_sum = -total_thickness / 2.0;
    for (xml_coll_t ci(x_service, _Unicode(component)); ci; ++ci, ncomponents++) {
      xml_comp_t x_comp    = ci;
      double     thickness = x_comp.thickness();
      Tube       c_tube{sensorPlaneRmin, sensorPlaneRmax, thickness / 2};
      Volume     c_vol{_toString(ncomponents, "component%d"), c_tube, desc.material(x_comp.materialStr())};
      c_vol.setVisAttributes(desc, x_comp.visStr());
      service_vol.placeVolume(c_vol, Position(0, 0, thickness_sum + thickness / 2.0));
      thickness_sum += thickness;
    }
    gasvolVol.placeVolume(service_vol,
                          Transform3D(Translation3D(sensorPlanePos.x(), sensorPlanePos.y(),
                                                    sensorPlanePos.z() - sensorThickness / 2 - total_thickness / 2)));
  }

  return det;
}

// clang-format off
DECLARE_DETELEMENT(epic_PFRICH, createDetector)
