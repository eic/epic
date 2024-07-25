// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Whitney Armstrong

/** \addtogroup Trackers Trackers
 * \brief Type: **BarrelTrackerWithFrame**.
 * \author Nivedith Ramasubramanian
 *
 * \ingroup trackers
 *
 * @{
 */
#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/Printout.h"
#include "DD4hep/Shapes.h"
#include "DDRec/DetectorData.h"
#include "DDRec/Surface.h"
#include "XML/Layering.h"
#include "XML/Utilities.h"
#include <array>
#include "DD4hepDetectorHelper.h"

using namespace std;
using namespace dd4hep;
using namespace dd4hep::rec;
using namespace dd4hep::detail;

#include "Math/Vector2D.h"
using namespace ROOT::Math;

/** Micromegas Barrel Tracker with space frame
 *
 * - Designed to process "mpgd_barrel.xml" ("mpgd_barrel_ver1" as of 2024/02).
 *
 * - Derived from "BarrelTrackerWithFrame_geo.cpp".
 *
 * - "support" tag not addressed.
 *
 * - "frame" tag within the module element.
 *
 * - Several models of modules, with each a distinct radius of curvature
 *  but a single XML <module> and a single <layer>.
 *
 * \code
 * \endcode
 *
 *
 * @author Yann Bedfer
 */
static Ref_t create_MPGDCylinderBarrelTracker(Detector& description, xml_h e,
                                              SensitiveDetector sens) {
  xml_det_t x_det = e;
  Material air    = description.air();
  int det_id      = x_det.id();
  string det_name = x_det.nameStr();
  DetElement sdet(det_name, det_id);

  vector<Volume> volumes;
  vector<PlacedVolume> sensitives;
  vector<VolPlane> volplane_surfaces;

  PlacedVolume pv;

  //#define DEBUG_MPGDCylinderBarrelTracker
#ifdef DEBUG_MPGDCylinderBarrelTracker
  // TEMPORARILY INCREASE VERBOSITY level for debugging purposes
  PrintLevel priorPrintLevel = printLevel();
  setPrintLevel(DEBUG);
#endif

  // Set detector type flag
  dd4hep::xml::setDetectorTypeFlag(x_det, sdet);
  auto& params = DD4hepDetectorHelper::ensureExtension<dd4hep::rec::VariantParameters>(sdet);

  // Add the volume boundary material if configured
  for (xml_coll_t bmat(x_det, _Unicode(boundary_material)); bmat; ++bmat) {
    xml_comp_t x_boundary_material = bmat;
    DD4hepDetectorHelper::xmlToProtoSurfaceMaterial(x_boundary_material, params,
                                                    "boundary_material");
  }

  dd4hep::xml::Dimension dimensions(x_det.dimensions());
  Assembly assembly(det_name);

  sens.setType("tracker");

  // ********** MODULE
  // ***** ONE AND ONLY ONE MODULE
  xml_coll_t modules(x_det, _U(module));
  if (modules.size() != 1) {
    // Present detector constructor can only handle ONE <module> tag
    printout(ERROR, "MPGDCylinderBarrelTracker", "Number of modules = %u. Must be = 1",
             modules.size());
    throw runtime_error("Logics error in building modules.");
  }
  xml_comp_t x_mod   = modules;
  string m_nam       = x_mod.nameStr();
  double stave_width = dimensions.width(), det_length = dimensions.length();
  printout(DEBUG, "MPGDCylinderBarrelTracker", "Module \"%s\" width = %.2f cm length = %.2f cm",
           m_nam.c_str(), stave_width, det_length);

  // ********** LAYER
  // ***** ONE AND ONLY ONE LAYER
  xml_coll_t li(x_det, _U(layer));
  if (li.size() != 1) {
    printout(ERROR, "MPGDCylinderBarrelTracker", "Number of layers = %d. Must be = 1",
             (int)li.size());
    throw runtime_error("Logics error in building modules.");
  }
  // ***** RETRIEVE PARAMETERS
  xml_comp_t x_layer = li;
  int lay_id         = x_layer.id();
  if (x_layer.moduleStr() != m_nam) {
    printout(ERROR, "MPGDCylinderBarrelTracker", "Layer \"%s\" does not match module \"%s\"",
             x_layer.moduleStr().c_str(), m_nam.c_str());
    throw runtime_error("Logics error in building layer.");
  }
  xml_comp_t x_barrel  = x_layer.child(_U(barrel_envelope));
  double barrel_length = x_barrel.z_length();
  double barrel_z0     = getAttrOrDefault(x_barrel, _U(z0), 0.);
  // ***** LAYOUTS
  xml_comp_t x_layout = x_layer.child(_U(rphi_layout));
  double phi0         = x_layout.phi0();             // Starting phi of first module.
  xml_comp_t z_layout = x_layer.child(_U(z_layout)); // Get the <z_layout> element.
  double z_gap        = z_layout.gap();
  // ***** UNVALID LAYOUT PROPERTIES
  // # of staves (along phi) and sectors (along z) are now derived from stave
  // width and sector length. Used to be specified directly. In order to remind
  // the user this is no longer the case, let's forbid the use of the
  // corresponding (and a few more) tags.
  const int nUnvalids                     = 4;
  const xml::Tag_t unvalidTags[nUnvalids] = {_U(phi_tilt), _U(nphi), _U(rc), _U(dr)};
  for (int uv = 0; uv < nUnvalids; uv++) {
    if (x_barrel.hasChild(unvalidTags[uv])) {
      const string tag = _U(nphi);
      printout(ERROR, "MPGDCylinderBarrelTracker",
               "Layer \"%s\": Unvalid property \"%s\" in \"rphi_layout\"", m_nam.c_str(),
               tag.c_str());
      throw runtime_error("Logics error in building modules.");
    }
  }

  // ***** MODELS
  // Model = set of two radii of curvature used alternatingly for staves
  //       + an offset to be applied radially.
  // - The offset may be null if the two radii are different enough that
  //  staves adjacent along phi do not overlap.
  // - There needs to be 2, as of 2024/03: Inner and outer.
  // StaveModel = radius of curvature (=> "nphi", given "rmin")
  //            + offset
  // sector2Models[2][2]: Which model for [inner,outer][r1,r2]
  typedef struct {
    string name;
    double rmin;
    // "rsensor" is supposed to have been set equal to the radius specified in
    // the segmentation and is used here to double-check the consistency
    // between the stack of module components and that segmentation radius.
    double rsensor;
    int nphi;
    double offset;
    int io;
  } StaveModel;
  StaveModel staveModels[4];
  int sector2Models[2][2];
  xml_coll_t mi(x_mod, _U(model));
  if (mi.size() != 2) {
    printout(ERROR, "MPGDCylinderBarrelTracker", "Number of models = %d. Must be = 2",
             (int)mi.size());
    throw runtime_error("Logics error in building modules.");
  }
  int nStaveModels;
  for (nStaveModels = 0; mi; ++mi) {
    xml_comp_t x_model = mi;
    double rmin1 = x_model.rmin1(), rmin2 = x_model.rmin2();
    double rsensor = x_model.attr<double>(_Unicode(rsensor));
    // Determine "nphi" from stave width and radius
    int nphi  = int(2 * M_PI * rmin1 / stave_width + .5),
        nphi2 = int(2 * M_PI * rmin2 / stave_width + .5);
    if (nphi2 != nphi) {
      printout(ERROR, "MPGDCylinderBarrelTracker",
               "Model \"%s\": rmin1,2 = %.2f,%.2f cm are incompatible", x_model.nameStr().c_str(),
               rmin1, rmin2);
      throw runtime_error("Logics error in building modules.");
    }
    double offset          = x_model.offset();
    StaveModel& staveModel = staveModels[nStaveModels++];
    staveModel.rmin        = rmin1;
    staveModel.rsensor     = rsensor;
    staveModel.nphi        = nphi;
    staveModel.offset      = offset;
    int io;
    if (nStaveModels == 1) // 1st encountered "x_model" => ...
      io = 0;              // ...it's the inner one => io = 0.
    else
      io = 1;
    staveModel.io        = io;
    sector2Models[io][0] = nStaveModels - 1;
    if (fabs(rmin2 - rmin1) > 1.e6) {
      StaveModel& staveMode2 = staveModels[nStaveModels++];
      staveMode2.rmin        = rmin2;
      staveModel.rsensor     = rsensor + rmin2 - rmin1;
      staveModel.nphi        = nphi;
      staveMode2.offset      = offset;
      staveMode2.io          = io;
      staveModel.name        = x_model.nameStr() + string("_1");
      staveMode2.name        = x_model.nameStr() + string("_2");
    } else
      staveModel.name = x_model.nameStr();
    sector2Models[io][1] = nStaveModels - 1;
  }
  printout(DEBUG, "MPGDCylinderBarrelTracker", "%d Stave Models:", nStaveModels);
  for (int iSM = 0; iSM < nStaveModels; iSM++) {
    StaveModel& sM = staveModels[iSM];
    printout(DEBUG, "MPGDCylinderBarrelTracker",
             "Stave Model #%d \"%s\": rmin = %.2f cm offset = ±%.2f cm %s", iSM, sM.name.c_str(),
             sM.rmin, sM.offset / 2, sM.io ? "Outer" : "Inner");
  }

  // ***** FRAMES
  // There must be two:
  // - Outward frame (wider, because supporting connectors)
  // - Otherwise frame.
  // Widest is taken as outward frame.
  typedef struct {
    string name;
    double width;
    string material;
    string vis;
  } Frame;
  Frame frames[2];
  xml_coll_t fi(x_mod, _U(frame));
  if (fi.size() != 2) {
    printout(ERROR, "MPGDCylinderBarrelTracker", "Number of frames = %d. Must be = 2",
             (int)fi.size());
    throw runtime_error("Logics error in building modules.");
  }
  printout(DEBUG, "MPGDCylinderBarrelTracker", "2 Frames:");
  int iFr;
  for (iFr = 0; fi; ++fi, iFr++) {
    xml_comp_t x_frame = fi;
    string name        = x_frame.nameStr();
    double width       = x_frame.width();
    string material = x_frame.materialStr(), vis = x_frame.visStr();
    Frame& frame   = frames[iFr];
    frame.name     = name;
    frame.width    = width;
    frame.material = material;
    frame.vis      = vis;
  }
  if (frames[0].width < frames[1].width) // Outward = widest must be first
    swap(frames[0], frames[1]);
  for (iFr = 0; iFr < 2; iFr++) {
    Frame& frame = frames[iFr];
    printout(DEBUG, "MPGDCylinderBarrelTracker",
             "Frame #%d \"%s\": width = %.2f cm material = \"%s\" vis = \"%s\"", iFr,
             frame.name.c_str(), frame.width, frame.material.c_str(), frame.vis.c_str());
  }

  // ***** SERVICE
  xml_coll_t si(x_mod, _Unicode(service));
  if (si.size() != 1) {
    printout(ERROR, "MPGDCylinderBarrelTracker", "Number of services = %d. Must be = 1",
             (int)si.size());
    throw runtime_error("Logics error in building modules.");
  }
  xml_comp_t x_service     = si;
  double service_thickness = x_service.thickness();
  printout(DEBUG, "MPGDCylinderBarrelTracker",
           "1 Service \"%s\": thickness = %.4f cm, material = \"%s\"", x_service.nameStr().c_str(),
           service_thickness, x_service.materialStr().c_str());

  // ***** TOTAL THICKNESS from components (used to build frames)
  double total_thickness = 0;
  xml_coll_t ci(x_mod, _U(module_component));
  for (ci.reset(), total_thickness = 0.0; ci; ++ci) {
    const xml_comp_t x_comp = ci;
    printout(DEBUG, "MPGDCylinderBarrelTracker", "\"%s\": \t total_thickness %.4f cm",
             x_comp.nameStr().c_str(), total_thickness / cm);
    total_thickness += x_comp.thickness();
  }
  printout(DEBUG, "MPGDCylinderBarrelTracker", " => total_thickness %.4f cm", total_thickness / cm);

  // Now that we know total thickness, let's "printout" the characteristics
  // of the models (extrema, phi superposition).
  printout(DEBUG, "MPGDCylinderBarrelTracker", "2 Sector Models:");
  double outerPhiSuperpos = 0; // Later used to define phi range of service to inner
  for (int io = 0; io < 2; io++) {
    StaveModel& sM0 = staveModels[sector2Models[io][0]];
    StaveModel& sM1 = staveModels[sector2Models[io][1]];
    XYVector vOffset(sM0.offset / 2, 0);
    double phiEdge0 = stave_width / 2 / sM0.rmin;
    XYVector vEdge0(sM0.rmin * cos(phiEdge0), sM0.rmin * sin(phiEdge0));
    vEdge0 -= vOffset;
    // Module0: Stave model has smaller curvature than circle centred on
    // beam axis. => Minimum is in the middle.
    double RMn = sM0.rmin - sM0.offset / 2, Mag0 = sqrt(vEdge0.Mag2());
    double phiEdge1 = stave_width / 2 / sM1.rmin;
    XYVector vEdge1(sM1.rmin * cos(phiEdge1), sM1.rmin * sin(phiEdge1));
    vEdge1 += vOffset;
    // Module1: Stave model has larger curvature than circle centred on
    // beam axis. => Maximum is in the middle.
    double RMx = sM1.rmin + sM1.offset / 2, Mag1 = sqrt(vEdge1.Mag2());
    double dPhi = acos(vEdge0.Dot(vEdge1) / Mag0 / Mag1);
    RMx += total_thickness;
    if (io == 1) {
      RMn -= service_thickness; // Outer sector: acount for services to inner sector
      outerPhiSuperpos = dPhi;
    }
    printout(
        DEBUG, "MPGDCylinderBarrelTracker",
        "Sector Model #%d,\"%s\" = \"%s\"+\"%s\": phi overlap = %.1f deg, Extrema R = %.2f,%.2f",
        io, io ? "Outer" : "Inner", sM0.name.c_str(), sM1.name.c_str(), dPhi / M_PI * 180, RMn,
        RMx);
  }

  // ********** LOOP OVER STAVE MODELS
  double total_length = 0; // Total length including frames
  for (int iSM = 0; iSM < nStaveModels; iSM++) {
    // "sensor_number" = Bit field in cellID identifying the sensitive surface.
    // We intend to use it to set up two discrimination schemes:
    // i) Stave model w/ its distinctive radius.
    // ii) Readout coordinate: phi or Z.
    int sensor_number      = 2 * iSM;
    StaveModel& staveModel = staveModels[iSM];
    // phi range, when excluding frames
    double stave_rmin = staveModel.rmin;
    double dphi       = stave_width / stave_rmin / 2;
    double phi_start = -dphi, phi_end = dphi;

    // ***** ASSEMBLY VOLUME: ONE PER STAVE MODEL
    Assembly m_vol(staveModel.name);
    volumes.push_back(m_vol);
    m_vol.setVisAttributes(description.visAttributes(x_mod.visStr()));

    // Stave frames
    double zthickness, stave_rmax = stave_rmin + total_thickness;
    Frame &out_frame = frames[0], &frame = frames[1];
    total_length = det_length + out_frame.width + frame.width;
    // Outward frame
    zthickness = out_frame.width;
    Tube frame_tube_O(stave_rmin, stave_rmax, zthickness / 2, phi_start, phi_end);
    Volume frame_vol_O(out_frame.name, frame_tube_O, description.material(out_frame.material));
    m_vol.placeVolume(frame_vol_O, Position(0, 0, -(total_length - zthickness) / 2));
    // Inward frame
    zthickness = frame.width; // Update "zthickness"
    Tube frame_tube_I(stave_rmin, stave_rmax, zthickness / 2, phi_start, phi_end);
    Volume frame_vol_I(frame.name, frame_tube_I, description.material(frame.material));
    m_vol.placeVolume(frame_vol_I, Position(0, 0, (total_length - zthickness) / 2));

    // Start/stop frames
    double frame_dphi =
        zthickness / stave_rmin; //converting the thickness of the frame to angular radians.
    Tube frame_tube_3(stave_rmin, stave_rmax, total_length / 2, phi_start - frame_dphi, phi_start);
    const string start_frame_nam("StartFrame");
    Volume frame_vol_3(start_frame_nam, frame_tube_3, description.material(frame.material));
    m_vol.placeVolume(frame_vol_3);

    Tube frame_tube_4(stave_rmin, stave_rmax, total_length / 2, phi_end, phi_end + frame_dphi);
    const string stop_frame_nam("StopFrame");
    Volume frame_vol_4(stop_frame_nam, frame_tube_4, description.material(frame.material));
    m_vol.placeVolume(frame_vol_4);

    frame_vol_O.setVisAttributes(description, out_frame.vis);
    frame_vol_I.setVisAttributes(description, frame.vis);
    frame_vol_3.setVisAttributes(description, frame.vis);
    frame_vol_4.setVisAttributes(description, frame.vis);

    // ***** OUTER SECTORS: SERVICES TO INNER
    // Don't know what the purpose of the "(inner|outer)_thickness" arg.s
    // to the surface volume (see infra) and whether "inner_thickness" has to
    // include the extra thickness corresponding to the services.
    // Note: The phi range of the service volume is set so that it's smaller
    // than that of the rest of components. This, in order to avoid the
    // addition of a thickness of service to the already thick superposition
    // of 2 module thicknesses (or module+frame) on the edge.
    if (staveModel.io == 1) {
      // Superposition along Z
      double outerPos  = barrel_length / 2 - total_length / 2;
      double innerPos  = (total_length + z_gap) / 2;
      double zSuperpos = total_length - outerPos + innerPos;
      // Parameters
      double serv_thickness = service_thickness;
      double serv_rmin      = stave_rmin - service_thickness;
      double serv_length    = total_length - zSuperpos;
      // phi range: stay away from phi overlap and frame, add 1mm margin
      double dPhi       = outerPhiSuperpos / 2 + (frame.width + .1) / stave_rmin;
      double serv_start = phi_start + dPhi, serv_stop = phi_end - dPhi;
      Tube c_tube(serv_rmin, serv_rmin + serv_thickness, serv_length / 2, serv_start, serv_stop);
      Volume c_vol(x_service.nameStr(), c_tube, description.material(x_service.materialStr()));
      pv = m_vol.placeVolume(c_vol, Position(0, 0, -zSuperpos / 2));
      c_vol.setRegion(description, x_service.regionStr());
      c_vol.setLimitSet(description, x_service.limitsStr());
      c_vol.setVisAttributes(description, x_service.visStr());
    }

    // ********** LOOP OVER COMPONENTS
    double comp_rmin          = stave_rmin;
    double thickness_so_far   = 0;
    xml_comp_t* sensitiveComp = 0;
    for (xml_coll_t mci(x_mod, _U(module_component)); mci; ++mci) {
      xml_comp_t x_comp     = mci;
      const string c_nam    = x_comp.nameStr();
      double comp_thickness = x_comp.thickness();
      Tube c_tube(comp_rmin, comp_rmin + comp_thickness, det_length / 2, phi_start, phi_end);
      Volume c_vol(c_nam, c_tube, description.material(x_comp.materialStr()));
      pv = m_vol.placeVolume(c_vol, Position(0, 0, (out_frame.width - frame.width) / 2));
      c_vol.setRegion(description, x_comp.regionStr());
      c_vol.setLimitSet(description, x_comp.limitsStr());
      c_vol.setVisAttributes(description, x_comp.visStr());
      if (x_comp.isSensitive()) {
        // ***** SENSITIVE VOLUME
        if (sensitiveComp) {
          printout(ERROR, "MPGDCylinderBarrelTracker",
                   "Component \"%s\": 2nd sensitive component in module \"%s\" (1st is \"%s\"). "
                   "One only allowed.",
                   c_nam.c_str(), m_nam.c_str(), sensitiveComp->nameStr().c_str());
          throw runtime_error("Logics error in building modules.");
        }
        sensitiveComp = &x_comp; // TODO: Add second sensitive
        pv.addPhysVolID("sensor", sensor_number++);
        c_vol.setSensitiveDetector(sens);
        sensitives.push_back(pv);

        // -------- create a measurement plane for the tracking surface attached to the sensitive volume -----
        Vector3D u(-1., 0., 0.);
        Vector3D v(0., -1., 0.);
        Vector3D n(0., 0., 1.);

        // Compute the inner (i.e. thickness until mid-sensitive-volume) and
        //             outer (from mid-sensitive-volume to top)
        // thicknesses that need to be assigned to the tracking surface
        // depending on wether the support is above or below the sensor (!?)
        double inner_thickness = thickness_so_far + comp_thickness / 2;
        double outer_thickness = total_thickness - inner_thickness;
        // Consistency(+/-1um) check: segmentation = stack of module components
        double rXCheck = comp_rmin + comp_thickness / 2;
        if (fabs(staveModel.rsensor - rXCheck) > .0001 / cm) {
          printout(ERROR, "MPGDCylinderBarrelTracker",
                   "Sensitive Component \"%s\" of StaveModel #%d,\"%s\": rsensor(%.4f cm) != "
                   "radius @ sensitive surface(%.4f cm)",
                   iSM, c_nam.c_str(), staveModel.name.c_str(), staveModel.rsensor / cm,
                   rXCheck / cm);
          throw runtime_error("Logics error in building modules.");
        }
        printout(DEBUG, "MPGDCylinderBarrelTracker",
                 "Stave Model #%d,\"%s\": Sensitive surface @ R = %.4f (%.4f,%.4f) cm", iSM,
                 staveModel.name.c_str(), staveModel.rsensor / cm, inner_thickness / cm,
                 outer_thickness / cm);

        SurfaceType type(SurfaceType::Sensitive);
        VolPlane surf(c_vol, type, inner_thickness, outer_thickness, u, v, n); //,o ) ;
        volplane_surfaces.push_back(surf);
      }
      comp_rmin += comp_thickness;
      thickness_so_far += comp_thickness;
    } //end of module component loop
  } //end of stave model loop

  // ********** LAYER
  // ***** ENVELOPE
  string lay_nam = det_name + _toString(x_layer.id(), "_layer%d");
  Tube lay_tub(x_barrel.inner_r(), x_barrel.outer_r(), barrel_length / 2);
  Volume lay_vol(lay_nam, lay_tub, air); // Create the layer envelope volume.
  Position lay_pos(0, 0, barrel_z0);
  lay_vol.setVisAttributes(description.visAttributes(x_layer.visStr()));
  printout(DEBUG, "MPGDCylinderBarrelTracker",
           "Layer \"%s\": rmin,max = %.2f,%.2f cm 1/2length = %.2f cm", lay_nam.c_str(),
           x_barrel.inner_r(), x_barrel.outer_r(), barrel_length / 2);

  DetElement lay_elt(sdet, lay_nam, lay_id);

  // the local coordinate systems of modules in dd4hep and acts differ
  // see http://acts.web.cern.ch/ACTS/latest/doc/group__DD4hepPlugins.html
  auto& layerParams =
      DD4hepDetectorHelper::ensureExtension<dd4hep::rec::VariantParameters>(lay_elt);

  // ********** LOOP OVER THE SECTORS IN z
  // ***** SECTOR POSITIONS ALONG Z
  // These are the 4 central values in Z where the four sets of modules, called
  // sectors, will be placed.
  double modz_pos[4] = {-barrel_length / 2 + (total_length) / 2, -(total_length + z_gap) / 2,
                        +(total_length + z_gap) / 2, +barrel_length / 2 - (total_length) / 2};
  int nModules       = 0;
  for (int iz = 0; iz < 4; iz++) {
    int io                = (iz == 1 || iz == 2) ? 0 : 1;
    int iSMs[2]           = {sector2Models[io][0], sector2Models[io][1]};
    int iSM0              = iSMs[0];
    const StaveModel& sm0 = staveModels[iSM0];
    int nphi              = sm0.nphi;
    double offset         = sm0.offset;
    double phi_incr       = 2 * M_PI / nphi; // Phi increment for one module.
    // ***** LOOP OVER THE STAVES IN phi.
    int iphi;
    double phic;
    for (iphi = 0, phic = phi0; iphi < nphi; iphi++, phic += phi_incr, nModules++) {
      // Every other module...
      int jphi = iphi % 2, iV = iSMs[jphi];    // ...swap stave volume/sensitive
      double rc = (2 * jphi - 1) * offset / 2; // ...flip sign of offset
      if (iz >= 2)
        rc *= -1;    // Swap >0 and <0 offsets.
      double x1, y1; // Coordinates of the centre of curvature of
      x1                 = rc * std::cos(phic);
      y1                 = rc * std::sin(phic);
      string module_name = _toString(10 * iz + iphi, "module%02d");
      DetElement mod_elt(lay_elt, module_name, nModules);
      printout(DEBUG, "MPGDCylinderBarrelTracker",
               "System %d Layer \"%s\",id=%d Module \"%s\",id=%d: x,y,r: %7.4f,%7.4f,%7.4f cm",
               det_id, lay_nam.c_str(), lay_id, module_name.c_str(), nModules, x1 / cm, y1 / cm,
               sqrt(x1 * x1 + y1 * y1) / cm);
      RotationZYX rot;
      rot = RotationZYX(phic, 0, 0);
      if (iz >= 2) // Rotate so that outward-frame faces outwards.
        // Note: The product of rotations seems to have to be done in this
        // order (i.e. rotation by "phi" about Z times rotation by pi about Y.
        // This, for reason I don't fully understand. Anyway, a consequence of
        // that is that, in order to have iz=2,3 not out of phase w/ iz=0,1,
        // one has to swap >0 and <0 offsets (see instruction supra).
        rot *= RotationZYX(0, M_PI, 0);
      Transform3D tr(rot, Position(x1, y1, modz_pos[iz]));
      Volume& module_vol = volumes[iV];
      pv                 = lay_vol.placeVolume(module_vol, tr);
      pv.addPhysVolID("module", nModules);
      mod_elt.setPlacement(pv);
      // ***** SENSITIVE COMPONENT
      PlacedVolume& sens_pv = sensitives[iV];
      DetElement comp_de(mod_elt, std::string("de_") + sens_pv.volume().name(), nModules);
      comp_de.setPlacement(sens_pv);
      auto& comp_de_params =
          DD4hepDetectorHelper::ensureExtension<dd4hep::rec::VariantParameters>(comp_de);
      comp_de_params.set<string>("axis_definitions", "XYZ");
      volSurfaceList(comp_de)->push_back(volplane_surfaces[iV]);
    }
  }

  for (xml_coll_t lmat(x_layer, _Unicode(layer_material)); lmat; ++lmat) {
    xml_comp_t x_layer_material = lmat;
    DD4hepDetectorHelper::xmlToProtoSurfaceMaterial(x_layer_material, layerParams,
                                                    "layer_material");
  }

  // ***** CREATE THE PhysicalVolume FOR THE LAYER.
  pv = assembly.placeVolume(lay_vol, lay_pos); // Place layer in mother
  pv.addPhysVolID("layer", lay_id);            // Set the layer ID.
  lay_elt.setAttributes(description, lay_vol, x_layer.regionStr(), x_layer.limitsStr(),
                        x_layer.visStr());
  lay_elt.setPlacement(pv);

  sdet.setAttributes(description, assembly, x_det.regionStr(), x_det.limitsStr(), x_det.visStr());
  assembly.setVisAttributes(description.invisible());
  pv = description.pickMotherVolume(sdet).placeVolume(assembly);
  pv.addPhysVolID("system", det_id); // Set the subdetector system ID.
  sdet.setPlacement(pv);

#ifdef DEBUG_MPGDCylinderBarrelTracker
  // Reset initial print level before exiting
  setPrintLevel(priorPrintLevel);
#endif

  return sdet;
}

//@}
// clang-format off
DECLARE_DETELEMENT(MPGDCylinderBarrelTracker, create_MPGDCylinderBarrelTracker)
