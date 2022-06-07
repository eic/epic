//==========================================================================
//  dRICh: Dual Ring Imaging Cherenkov Detector
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
#include "GeometryHelpers.h"
#include "Math/Point2D.h"
#include "TMath.h"
#include "TString.h"
#include <XML/Helper.h>

using namespace dd4hep;
using namespace dd4hep::rec;

// create the detector
static Ref_t createDetector(Detector& desc, xml::Handle_t handle, SensitiveDetector sens)
{
  xml::DetElement detElem = handle;
  std::string     detName = detElem.nameStr();
  int             detID   = detElem.id();

  DetElement            det(detName, detID);
  xml::Component        dims    = detElem.dimensions();
  OpticalSurfaceManager surfMgr = desc.surfaceManager();

  // attributes -----------------------------------------------------------
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
  auto   radiatorElem = detElem.child(_Unicode(radiator));
  double radiatorRmin = radiatorElem.attr<double>(_Unicode(rmin));
  double radiatorRmax = radiatorElem.attr<double>(_Unicode(rmax));
  // FIXME unused
  // double radiatorPhiw       = radiatorElem.attr<double>(_Unicode(phiw));
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
  // - sensor module
  auto   sensorElem      = detElem.child(_Unicode(sensors)).child(_Unicode(module));
  auto   sensorMat       = desc.material(sensorElem.attr<std::string>(_Unicode(material)));
  auto   sensorVis       = desc.visAttributes(sensorElem.attr<std::string>(_Unicode(vis)));
  auto   sensorSurf      = surfMgr.opticalSurface(sensorElem.attr<std::string>(_Unicode(surface)));
  double sensorSide      = sensorElem.attr<double>(_Unicode(side));
  double sensorGap       = sensorElem.attr<double>(_Unicode(gap));
  double sensorThickness = sensorElem.attr<double>(_Unicode(thickness));
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
  // - debugging switches
  int  debug_optics_mode = detElem.attr<int>(_Unicode(debug_optics));
  bool debug_mirror      = mirrorElem.attr<bool>(_Unicode(debug));
  bool debug_sensors     = sensorSphElem.attr<bool>(_Unicode(debug));

  // if debugging optics, override some settings
  bool debug_optics = debug_optics_mode > 0;
  if (debug_optics) {
    printout(WARNING, "DRich_geo", "DEBUGGING DRICH OPTICS");
    switch (debug_optics_mode) {
    case 1:
      vesselMat = aerogelMat = filterMat = sensorMat = gasvolMat = desc.material("VacuumOptical");
      break;
    case 2:
      vesselMat = aerogelMat = filterMat = sensorMat = desc.material("VacuumOptical");
      break;
    default:
      printout(FATAL, "DRich_geo", "UNKNOWN debug_optics_mode");
      return det;
    };
    aerogelVis = sensorVis = mirrorVis;
    gasvolVis = vesselVis = desc.invisible();
  };

  // BUILD VESSEL ====================================================================
  /* - `vessel`: aluminum enclosure, the mother volume of the dRICh
   * - `gasvol`: gas volume, which fills `vessel`; all other volumes defined below
   *   are children of `gasvol`
   * - the dRICh vessel geometry has two regions: the snout refers to the conic region
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

  //  extra solids for `debug_optics` only
  Box vesselBox(1001, 1001, 1001);
  Box gasvolBox(1000, 1000, 1000);

  // choose vessel and gasvol solids (depending on `debug_optics_mode` (0=disabled))
  Solid vesselSolid, gasvolSolid;
  switch (debug_optics_mode) {
  case 0:
    vesselSolid = vesselUnion;
    gasvolSolid = gasvolUnion;
    break; // `!debug_optics`
  case 1:
    vesselSolid = vesselBox;
    gasvolSolid = gasvolBox;
    break;
  case 2:
    vesselSolid = vesselBox;
    gasvolSolid = gasvolUnion;
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
  //   the z-direction to be in the proper ATHENA-integration location
  // - these reference positions are for the frontplane and backplane of the vessel,
  //   with respect to the vessel origin position
  auto originFront = Position(0., 0., -tankLength / 2.0 - snoutLength);
  // FIXME unused
  // auto originBack  = Position(0., 0., tankLength / 2.0);

  // initialize sensor centroids (used for mirror parameterization below); this is
  // the average (x,y,z) of the placed sensors, w.r.t. originFront
  double sensorCentroidX = 0;
  double sensorCentroidZ = 0;
  int    sensorCount     = 0;

  // sensitive detector type
  sens.setType("tracker");

  // BUILD RADIATOR ====================================================================

  // attributes
  double airGap = 0.01 * mm; // air gap between aerogel and filter (FIXME? actually it's currently a gas gap)

  // solid and volume: create aerogel and filter
  Cone aerogelSolid(aerogelThickness / 2, radiatorRmin, radiatorRmax,
                    radiatorRmin + boreDelta * aerogelThickness / vesselLength,
                    radiatorRmax + snoutDelta * aerogelThickness / snoutLength);
  Cone filterSolid(filterThickness / 2, radiatorRmin + boreDelta * (aerogelThickness + airGap) / vesselLength,
                   radiatorRmax + snoutDelta * (aerogelThickness + airGap) / snoutLength,
                   radiatorRmin + boreDelta * (aerogelThickness + airGap + filterThickness) / vesselLength,
                   radiatorRmax + snoutDelta * (aerogelThickness + airGap + filterThickness) / snoutLength);

  Volume aerogelVol(detName + "_aerogel", aerogelSolid, aerogelMat);
  Volume filterVol(detName + "_filter", filterSolid, filterMat);
  aerogelVol.setVisAttributes(aerogelVis);
  filterVol.setVisAttributes(filterVis);

  // aerogel placement and surface properties
  // TODO [low-priority]: define skin properties for aerogel and filter
  auto radiatorPos = Position(0., 0., radiatorFrontplane) + originFront;
  auto aerogelPV   = gasvolVol.placeVolume(
      aerogelVol,
      Translation3D(radiatorPos.x(), radiatorPos.y(), radiatorPos.z()) // re-center to originFront
          *
          RotationY(
              radiatorPitch) // change polar angle to specified pitch // FIXME: probably broken (not yet in use anyway)
  );
  DetElement aerogelDE(det, "aerogel_de", 0);
  aerogelDE.setPlacement(aerogelPV);
  // SkinSurface aerogelSkin(desc, aerogelDE, "mirror_optical_surface", aerogelSurf, aerogelVol);
  // aerogelSkin.isValid();

  // filter placement and surface properties
  if (!debug_optics) {
    auto filterPV = gasvolVol.placeVolume(
        filterVol, Translation3D(0., 0., airGap)                                          // add an air gap
                       * Translation3D(radiatorPos.x(), radiatorPos.y(), radiatorPos.z()) // re-center to originFront
                       * RotationY(radiatorPitch)                                         // change polar angle
                       * Translation3D(0., 0., (aerogelThickness + filterThickness) / 2.) // move to aerogel backplane
    );
    DetElement filterDE(det, "filter_de", 0);
    filterDE.setPlacement(filterPV);
    // SkinSurface filterSkin(desc, filterDE, "mirror_optical_surface", filterSurf, filterVol);
    // filterSkin.isValid();
  };

  // SECTOR LOOP //////////////////////////////////
  for (int isec = 0; isec < nSectors; isec++) {

    // debugging filters, limiting the number of sectors
    if ((debug_mirror || debug_sensors || debug_optics) && isec != 0)
      continue;

    // sector rotation about z axis
    double      sectorRotation = isec * 360 / nSectors * degree;
    std::string secName        = "sec" + std::to_string(isec);

    // BUILD SENSORS ====================================================================

    // if debugging sphere properties, restrict number of sensors drawn
    if (debug_sensors) {
      sensorSide = 2 * M_PI * sensorSphRadius / 64;
    };

    // solid and volume: single sensor module
    Box    sensorSolid(sensorSide / 2., sensorSide / 2., sensorThickness / 2.);
    Volume sensorVol(detName + "_sensor_" + secName, sensorSolid, sensorMat);
    sensorVol.setVisAttributes(sensorVis);

    auto sensorSphPos = Position(sensorSphCenterX, 0., sensorSphCenterZ) + originFront;

    // sensitivity
    if (!debug_optics)
      sensorVol.setSensitiveDetector(sens);

    // SENSOR MODULE LOOP ------------------------
    /* ALGORITHM: generate sphere of positions
     * - NOTE: there are two coordinate systems here:
     *   - "global" the main ATHENA coordinate system
     *   - "generator" (vars end in `Gen`) is a local coordinate system for
     *     generating points on a sphere; it is related to the global system by
     *     a rotation; we do this so the "patch" (subset of generated
     *     positions) of sensors we choose to build is near the equator, where
     *     point distribution is more uniform
     * - PROCEDURE: loop over `thetaGen`, with subloop over `phiGen`, each divided evenly
     *   - the number of points to generate depends how many sensors (+`sensorGap`)
     *     can fit within each ring of constant `thetaGen` or `phiGen`
     *   - we divide the relevant circumference by the sensor
     *     size(+`sensorGap`), and this number is allowed to be a fraction,
     *     because likely we don't care about generating a full sphere and
     *     don't mind a "seam" at the overlap point
     *   - if we pick a patch of the sphere near the equator, and not near
     *     the poles or seam, the sensor distribution will appear uniform
     */

    // initialize module number for this sector
    int imod = 0;

    // thetaGen loop: iterate less than "0.5 circumference / sensor size" times
    double nTheta = M_PI * sensorSphRadius / (sensorSide + sensorGap);
    for (int t = 0; t < (int)(nTheta + 0.5); t++) {
      double thetaGen = t / ((double)nTheta) * M_PI;

      // phiGen loop: iterate less than "circumference at this latitude / sensor size" times
      double nPhi = 2 * M_PI * sensorSphRadius * std::sin(thetaGen) / (sensorSide + sensorGap);
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
        if (debug_sensors)
          patchCut = std::fabs(phiCheck) < sensorSphPatchPhiw;
        if (patchCut) {

          // append sensor position to centroid calculation
          if (isec == 0) {
            sensorCentroidX += xCheck;
            sensorCentroidZ += zCheck;
            sensorCount++;
          };

          // placement (note: transformations are in reverse order)
          // - transformations operate on global coordinates; the corresponding
          //   generator coordinates are provided in the comments
          auto sensorPV =
              gasvolVol.placeVolume(sensorVol,
                                    RotationZ(sectorRotation) // rotate about beam axis to sector
                                        * Translation3D(sensorSphPos.x(), sensorSphPos.y(),
                                                        sensorSphPos.z())        // move sphere to reference position
                                        * RotationX(phiGen)                      // rotate about `zGen`
                                        * RotationZ(thetaGen)                    // rotate about `yGen`
                                        * Translation3D(sensorSphRadius, 0., 0.) // push radially to spherical surface
                                        * RotationY(M_PI / 2)  // rotate sensor to be compatible with generator coords
                                        * RotationZ(-M_PI / 2) // correction for readout segmentation mapping
              );

          // generate LUT for module number -> sensor position, for readout mapping tests
          // if(isec==0) printf("%d %f %f\n",imod,sensorPV.position().x(),sensorPV.position().y());

          // properties
          sensorPV.addPhysVolID("sector", isec).addPhysVolID("module", imod);
          DetElement sensorDE(det, Form("sensor_de%d_%d", isec, imod),
                              (imod << 3) | isec); // id must match IRTAlgorithm usage
          sensorDE.setPlacement(sensorPV);
          if (!debug_optics) {
            SkinSurface sensorSkin(desc, sensorDE, Form("sensor_optical_surface%d", isec), sensorSurf, sensorVol);
            sensorSkin.isValid();
          };

          // increment sensor module number
          imod++;

        }; // end patch cuts
      };   // end phiGen loop
    };     // end thetaGen loop

    // calculate centroid sensor position
    if (isec == 0) {
      sensorCentroidX /= sensorCount;
      sensorCentroidZ /= sensorCount;
    };

    // END SENSOR MODULE LOOP ------------------------

    // BUILD MIRRORS ====================================================================

    // derive spherical mirror parameters `(zM,xM,rM)`, for given image point
    // coordinates `(zI,xI)` and `dO`, defined as the z-distance between the
    // object and the mirror surface
    // - all coordinates are specified w.r.t. the object point coordinates
    // - this is point-to-point focusing, but it can be used to effectively steer
    //   parallel-to-point focusing
    double zM, xM, rM;
    auto   FocusMirror = [&zM, &xM, &rM](double zI, double xI, double dO) {
      zM = dO * zI / (2 * dO - zI);
      xM = dO * xI / (2 * dO - zI);
      rM = dO - zM;
    };

    // attributes, re-defined w.r.t. IP, needed for mirror positioning
    double zS = sensorSphCenterZ + vesselZmin; // sensor sphere attributes
    double xS = sensorSphCenterX;
    // double rS = sensorSphRadius;
    double B = vesselZmax - mirrorBackplane; // distance between IP and mirror back plane

    // focus 1: set mirror to focus IP on center of sensor sphere `(zS,xS)`
    /*double zF = zS;
    double xF = xS;
    FocusMirror(zF,xF,B);*/

    // focus 2: move focal region along sensor sphere radius, according to `focusTuneLong`
    // - specifically, along the radial line which passes through the approximate centroid
    //   of the sensor region `(sensorCentroidZ,sensorCentroidX)`
    // - `focusTuneLong` is the distance to move, given as a fraction of `sensorSphRadius`
    // - `focusTuneLong==0` means `(zF,xF)==(zS,xS)`
    // - `focusTuneLong==1` means `(zF,xF)` will be on the sensor sphere, near the centroid
    /*
    double zC = sensorCentroidZ + vesselZmin;
    double xC = sensorCentroidX;
    double slopeF = (xC-xS) / (zC-zS);
    double thetaF = std::atan(std::fabs(slopeF));
    double zF = zS + focusTuneLong * sensorSphRadius * std::cos(thetaF);
    double xF = xS - focusTuneLong * sensorSphRadius * std::sin(thetaF);
    //FocusMirror(zF,xF,B);

    // focus 3: move along line perpendicular to focus 2's radial line,
    // according to `focusTunePerp`, with the same numerical scale as `focusTuneLong`
    zF += focusTunePerp * sensorSphRadius * std::cos(M_PI/2-thetaF);
    xF += focusTunePerp * sensorSphRadius * std::sin(M_PI/2-thetaF);
    FocusMirror(zF,xF,B);
    */

    // focus 4: use (z,x) coordinates for tune parameters
    double zF = zS + focusTuneZ;
    double xF = xS + focusTuneX;
    FocusMirror(zF, xF, B);

    // re-define mirror attributes to be w.r.t vessel front plane
    double mirrorCenterZ = zM - vesselZmin;
    double mirrorCenterX = xM;
    double mirrorRadius  = rM;

    // spherical mirror patch cuts and rotation
    double mirrorThetaRot = std::asin(mirrorCenterX / mirrorRadius);
    double mirrorTheta1   = mirrorThetaRot - std::asin((mirrorCenterX - mirrorRmin) / mirrorRadius);
    double mirrorTheta2   = mirrorThetaRot + std::asin((mirrorRmax - mirrorCenterX) / mirrorRadius);

    // if debugging, draw full sphere
    if (debug_mirror) {
      mirrorTheta1 = 0;
      mirrorTheta2 = M_PI; /*mirrorPhiw=2*M_PI;*/
    };

    // solid : create sphere at origin, with specified angular limits;
    // phi limits are increased to fill gaps (overlaps are cut away later)
    Sphere mirrorSolid1(mirrorRadius, mirrorRadius + mirrorThickness, mirrorTheta1, mirrorTheta2, -40 * degree,
                        40 * degree);

    // print mirror attributes for sector 0
    // if(isec==0) printf("dRICH mirror (zM, xM, rM) = (%f, %f, %f)\n",zM,xM,rM); // coords w.r.t. IP

    // mirror placement transformation (note: transformations are in reverse order)
    auto        mirrorPos = Position(mirrorCenterX, 0., mirrorCenterZ) + originFront;
    Transform3D mirrorPlacement(
        Translation3D(mirrorPos.x(), mirrorPos.y(), mirrorPos.z()) // re-center to specified position
        * RotationY(-mirrorThetaRot) // rotate about vertical axis, to be within vessel radial walls
    );

    // cut overlaps with other sectors using "pie slice" wedges, to the extent specified
    // by `mirrorPhiw`
    Tube              pieSlice(0.01 * cm, vesselRmax2, tankLength / 2.0, -mirrorPhiw / 2.0, mirrorPhiw / 2.0);
    IntersectionSolid mirrorSolid2(pieSlice, mirrorSolid1, mirrorPlacement);

    // mirror volume, attributes, and placement
    Volume mirrorVol(detName + "_mirror_" + secName, mirrorSolid2, mirrorMat);
    mirrorVol.setVisAttributes(mirrorVis);
    auto mirrorPV2 = gasvolVol.placeVolume(mirrorVol,
                                           RotationZ(sectorRotation) // rotate about beam axis to sector
                                               * Translation3D(0, 0, 0));

    // properties
    DetElement mirrorDE(det, Form("mirror_de%d", isec), isec);
    mirrorDE.setPlacement(mirrorPV2);
    SkinSurface mirrorSkin(desc, mirrorDE, Form("mirror_optical_surface%d", isec), mirrorSurf, mirrorVol);
    mirrorSkin.isValid();

  }; // END SECTOR LOOP //////////////////////////

  // place gas volume
  PlacedVolume gasvolPV = vesselVol.placeVolume(gasvolVol, Position(0, 0, 0));
  DetElement   gasvolDE(det, "gasvol_de", 0);
  gasvolDE.setPlacement(gasvolPV);

  // place mother volume (vessel)
  Volume       motherVol = desc.pickMotherVolume(det);
  PlacedVolume vesselPV  = motherVol.placeVolume(vesselVol, Position(0, 0, vesselZmin) - originFront);
  vesselPV.addPhysVolID("system", detID);
  det.setPlacement(vesselPV);

  return det;
};

// clang-format off
DECLARE_DETELEMENT(ecce_DRICH, createDetector)
