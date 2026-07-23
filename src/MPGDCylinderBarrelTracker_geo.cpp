// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022-2024 Whitney Armstrong, Nivedith Ramasubramanian, Yann Bedfer

/** \addtogroup Trackers Trackers
 * \brief Type: **Barrel Cylinder MPGD with frames**.
 * \author Nivedith Ramasubramanian
 * Modified by M. Posik, Y. Bedfer
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
using ROOT::Math::XYVector;

/** Micromegas Barrel Tracker with fish scale layout
 *
 * - Designed to process "mpgd_barrel.xml".
 *
 * - Derived from "BarrelTrackerWithFrame_geo.cpp".
 *
 * - No provision for <support> nor <service> tags.
 *
 * - "frame" tag within the module element.
 *
 * - Single XML <module> and a single <layer>.
 *
 * - Two distinct versions of ".xml" are covered:
 *  I) Single Sensitive Volume (which name is then "DriftGap"): "eicrecon" is to
 *   be executed while setting bit 0x1 of option "MPGD:SiFactoryPattern".
 * II) Multiple Sensitive Volume: Bit 0x1 of above-mentioned option not set.
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

  // Sensitive volumes and associated surfaces
  // - There can be either one or five.
  int sensitiveVolumeSet = 5; // 1: single volume, 5: 5 volumes, -1: error
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

  Assembly assembly(det_name);

  sens.setType("tracker");

  // ********** MODULE
  // ***** ONE AND ONLY ONE <module>
  xml_coll_t modules(x_det, _U(module));
  if (modules.size() != 1) {
    // Present detector constructor can only handle ONE <module> tag
    printout(ERROR, "MPGDCylinderBarrelTracker", "Number of modules = %u. Must be = 1",
             modules.size());
    throw runtime_error("Logics error in building modules.");
  }
  // ***** RETRIEVE MODULE PARAMETERS
  // (Module is here the active piece, excluding frames and services.)
  xml_comp_t x_mod   = modules;
  string m_nam       = x_mod.nameStr();
  xml_comp_t x_dimensions = x_mod.dimensions();
  double mod_width = x_dimensions.width(), mod_length  = x_dimensions.length();
  double mod_r0    = x_dimensions.attr<double>(_Unicode(r0));
  double mod_rmin  = x_dimensions.rmin(),  mod_rsensor = x_dimensions.attr<double>(_Unicode(rsensor));
  // phi range, when excluding frames
  double dphi      = mod_width / mod_r0 / 2;
  double phi_start = -dphi, phi_end = dphi;
  printout(DEBUG, "MPGDCylinderBarrelTracker", "Module \"%s\": width,length %.2f,%.2f cm r0,rmin,rsensor %.2f,%.2f,%.2f cm phi +/-%.2f°",
           m_nam.c_str(), mod_width / cm, mod_length / cm,
	   mod_r0 / cm, mod_rmin / cm, mod_rsensor / cm, dphi / M_PI * 180);
  // ********** LAYER
  // ***** ONE AND ONLY ONE <layer>
  xml_coll_t li(x_det, _U(layer));
  if (li.size() != 1) {
    printout(ERROR, "MPGDCylinderBarrelTracker", "Number of layers = %d. Must be = 1",
             (int)li.size());
    throw runtime_error("Logics error in building modules.");
  }
  // ***** RETRIEVE LAYER PARAMETERS
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
  xml_comp_t rphi_layout = x_layer.child(_U(rphi_layout));
  double x_offset     = rphi_layout.x_offset();
  double y_offset     = rphi_layout.y_offset();
  double phi_tilt     = rphi_layout.phi_tilt();
  double delta        = rphi_layout.delta();
  int    nphi         = rphi_layout.nphi();
  double phi0         = rphi_layout.phi0(); // Overall shift in phi
  xml_comp_t z_layout = x_layer.child(_U(z_layout));
  double z0           = z_layout.z0();
  double z1           = z_layout.z1();
  double z2           = z_layout.z2();
  // ***** FRAMES
  // There must be two:
  // - Outward frame (wider, because supporting connectors)
  // - Otherwise frame (along phi or along z).
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
  // ***** RETRIEVE LAYER PARAMETERS
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
	     "Frame #%d \"%s\": width %.2f cm material \"%s\" vis \"%s\"", iFr,
	     frame.name.c_str(), frame.width, frame.material.c_str(), frame.vis.c_str());
  }

  // ***** TOTAL THICKNESS from components (later used to build frames)
  double total_thickness = 0;
  xml_coll_t ci(x_mod, _U(module_component));
  for (ci.reset(), total_thickness = 0.0; ci; ++ci) {
    const xml_comp_t x_comp = ci;
    total_thickness += x_comp.thickness();
    printout(DEBUG, "MPGDCylinderBarrelTracker", "\"%s\": \t thickness %.4f %.4f cm material \"%s\"",
             x_comp.nameStr().c_str(), x_comp.thickness() / cm, total_thickness / cm, x_comp.materialStr().c_str());
  }
  printout(DEBUG, "MPGDCylinderBarrelTracker", " => total_thickness %.4f cm", total_thickness / cm);

  // ***** ASSEMBLY VOLUME: A UNIQUE ONE
  Assembly m_vol(m_nam);
  m_vol.setVisAttributes(description.visAttributes(x_mod.visStr()));

  // ***** BUILD FRAMES
  double zthickness, mod_rmax = mod_rmin + total_thickness;
  Frame &out_frame = frames[0], &frame = frames[1];
  double total_length = // Total length including frames
    mod_length + out_frame.width + frame.width;
  // Outward frame
  zthickness = out_frame.width;
  Tube frame_tube_O(mod_rmin, mod_rmax, zthickness / 2, phi_start, phi_end);
  Volume frame_vol_O(out_frame.name, frame_tube_O, description.material(out_frame.material));
  m_vol.placeVolume(frame_vol_O, Position(0, 0, -(total_length - zthickness) / 2));
  // Inward frame
  zthickness = frame.width; // Update "zthickness"
  Tube frame_tube_I(mod_rmin, mod_rmax, zthickness / 2, phi_start, phi_end);
  Volume frame_vol_I(frame.name, frame_tube_I, description.material(frame.material));
  m_vol.placeVolume(frame_vol_I, Position(0, 0, (total_length - zthickness) / 2));
  // Start/stop frames
  double frame_dphi =
    zthickness / mod_r0; //converting the thickness of the frame to angular radians.
  Tube frame_tube_3(mod_rmin, mod_rmax, total_length / 2, phi_start - frame_dphi, phi_start);
  const string start_frame_nam("StartFrame");
  Volume frame_vol_3(start_frame_nam, frame_tube_3, description.material(frame.material));
  m_vol.placeVolume(frame_vol_3);
  Tube frame_tube_4(mod_rmin, mod_rmax, total_length / 2, phi_end, phi_end + frame_dphi);
  const string stop_frame_nam("StopFrame");
  Volume frame_vol_4(stop_frame_nam, frame_tube_4, description.material(frame.material));
  m_vol.placeVolume(frame_vol_4);

  frame_vol_O.setVisAttributes(description, out_frame.vis);
  frame_vol_I.setVisAttributes(description, frame.vis);
  frame_vol_3.setVisAttributes(description, frame.vis);
  frame_vol_4.setVisAttributes(description, frame.vis);

  // ********** LOOP OVER MODULE COMPONENTS
  // Initializations
  double comp_rmin        = mod_rmin;
  double thickness_so_far = 0;
  // Pattern of Multiple Sensitive Volumes
  // - The "DriftGap" may be subdivided into SUBVOLUMES (this, in order to
  //  have one sensitive component per strip coordinate and accessorily,
  //  some extras).
  int nSensitives = 0;
  for (xml_coll_t mci(x_mod, _U(module_component)); mci; ++mci) {
    xml_comp_t x_comp     = mci;
    const string c_nam    = x_comp.nameStr();
    double comp_thickness = x_comp.thickness();
    Tube c_tube(comp_rmin, comp_rmin + comp_thickness, mod_length / 2, phi_start, phi_end);
    Volume c_vol(c_nam, c_tube, description.material(x_comp.materialStr()));
    pv = m_vol.placeVolume(c_vol, Position(0, 0, (out_frame.width - frame.width) / 2));
    c_vol.setRegion(description, x_comp.regionStr());
    c_vol.setLimitSet(description, x_comp.limitsStr());
    c_vol.setVisAttributes(description, x_comp.visStr());
    if (x_comp.isSensitive()) {
      // ***** SENSITIVE VOLUME
      if (nSensitives >= 5) {
	sensitiveVolumeSet = -1;
	break;
      }
      pv.addPhysVolID("sensor", 0);
      // StripID. Single Sensitive Volume?
      if (c_nam == "DriftGap") {
	if (nSensitives != 0) {
	  sensitiveVolumeSet = -1;
	  break;
	}
	sensitiveVolumeSet = 1;
      }
      int strip_id = x_comp.key();
      pv.addPhysVolID("strip", strip_id);
      c_vol.setSensitiveDetector(sens);
      sensitives.push_back(pv);

      // -------- create a measurement plane for the tracking surface attached to the sensitive volume -----
      Vector3D u(-1., 0., 0.);
      Vector3D v(0., -1., 0.);
      Vector3D n(0., 0., 1.);

      if (strip_id == 0) {
	// Consistency(+/-1um) check:
	// "rsensor" is supposed to have been assigned to the radius specified
	// in the segmentation. Let's double-check the consistency between the
	// stack of module components and that radius.
	double rXCheck = comp_rmin + comp_thickness / 2;
	if (fabs(mod_rsensor - rXCheck) > .0001 / cm) {
	  printout(ERROR, "MPGDCylinderBarrelTracker",
		   "Sensitive Component \"%s\": rsensor(%.4f cm) != "
		   "radius @ sensitive surface(%.4f cm)",
		   c_nam.c_str(), mod_rsensor / cm, rXCheck / cm);
	  throw runtime_error("Logics error in building modules.");
	}
      }

      // Compute the inner (i.e. thickness until mid-sensitive-volume) and
      //             outer (from mid-sensitive-volume to top)
      // thicknesses that need to be assigned to the tracking surface
      // depending on whether the support is above or below the sensor.
      double inner_thickness, outer_thickness;
      if (sensitiveVolumeSet == 1) {
	inner_thickness = thickness_so_far + comp_thickness / 2;
	outer_thickness = total_thickness - inner_thickness;
      } else if (nSensitives == 0) {
	inner_thickness = thickness_so_far + comp_thickness / 2;
	outer_thickness = comp_thickness / 2;
      } else if (nSensitives == 4) {
	inner_thickness = comp_thickness / 2;
	outer_thickness = total_thickness - thickness_so_far - comp_thickness / 2;
      } else {
	inner_thickness = outer_thickness = comp_thickness / 2;
      }
      printout(DEBUG, "MPGDCylinderBarrelTracker",
	       "Sensitive surface @ R = %.4f (%.4f,%.4f) cm",
	       mod_rsensor, inner_thickness / cm, outer_thickness / cm);
      SurfaceType type(SurfaceType::Sensitive);
      VolPlane surf(c_vol, type, inner_thickness, outer_thickness, u, v, n); //,o ) ;
      volplane_surfaces.push_back(surf);
      nSensitives++;
    }
    comp_rmin += comp_thickness;
    thickness_so_far += comp_thickness;
  } //end of module component loop
  if (sensitiveVolumeSet < 0 || sensitiveVolumeSet != nSensitives) {
    printout(ERROR, "MPGDCylinderBarrelTracker",
	     "Invalid set of Sensitive Volumes: it's either one (named \"DriftGap\") or 5");
    throw runtime_error("Logics error in building modules.");
  }

  // ********** BUILD LAYER
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

  // ********** LOOP OVER THE FOUR SECTIONS ALONG Z
  for (int iz = 0; iz < 4; iz++) {
    // Z abscissa?
    bool isBackward = iz < 2, isExternal = iz == 0 || iz == 3;
    double zCorr = (out_frame.width - frame.width) / 2;
    double Zc = isExternal ? z2 : z1;
    Zc += zCorr; // Transform from sensitive area to overall module
    if (isBackward) Zc *= -1; // Apply symmetry w.r.t. overall offset
    Zc += z0;
    // ***** LOOP OVER THE STAVES IN phi.
    for (int iphi = 0; iphi < nphi; iphi++) {
      // Module index and name
      int moduleIdx = nphi * iz + iphi;
      string module_name = _toString(moduleIdx, "module%02d");
      // phi angle
      // Starting value? Sectors sit aside, as opposed to astride, the axes.
      // With still some "phi0" margin, to make for possible imperfections.
      double phi_incr = 2 * M_PI / nphi;
      double phi = phi0 +  phi_incr / 2 + iphi * phi_incr; 
      // Center of cylinder
      double Xc = x_offset, Yc = y_offset;
      if (isExternal) { // Add "delta" along tilt
	//double t = phi_tilt / 180 * M_PI;
	Xc += delta * cos(phi_tilt); Yc += delta * sin(phi_tilt);
      }
      // DetElement
      DetElement mod_elt(lay_elt, module_name, moduleIdx);
      // Plcament
      Translation3D origin2Centre(Xc,Yc,Zc);
      RotationZYX tilt(phi_tilt, 0, 0);
      RotationZYX rot(phi, 0, 0);
      Transform3D tr;
      if (isBackward)
	tr = rot * origin2Centre * tilt;
      else {
	// Reflection, so that outward-frame faces outwards.
	// Note: we need a reflection, as opposed to a rotation, so that fish
	// scales keep being oriented the same. This in turn, so that the set
	// of all four sections can be moved out in one go for maintenance.
	Translation3D restoreZ(0,0,2*Zc);
        tr = restoreZ * Rotation3D(1,0,0,0,1,0,0,0,-1) * rot * origin2Centre * tilt;
      }
      Volume& module_vol = m_vol;
      pv                 = lay_vol.placeVolume(module_vol, tr);
      pv.addPhysVolID("module", moduleIdx);
      mod_elt.setPlacement(pv);
      for (int iSensitive = 0; iSensitive < sensitiveVolumeSet; iSensitive++) {
        // ***** SENSITIVE COMPONENTS
        PlacedVolume& sens_pv = sensitives[iSensitive];
        int de_id             = nphi * 4 * iSensitive + moduleIdx;
        DetElement comp_de(mod_elt,
                           std::string("de_") + sens_pv.volume().name() + _toString(de_id, "%02d"),
                           de_id);
        comp_de.setPlacement(sens_pv);
        auto& comp_de_params =
            DD4hepDetectorHelper::ensureExtension<dd4hep::rec::VariantParameters>(comp_de);
        comp_de_params.set<string>("axis_definitions", "XYZ");
        volSurfaceList(comp_de)->push_back(volplane_surfaces[iSensitive]);
      }
      printout(DEBUG, "MPGDCylinderBarrelTracker",
               "System %d Layer %d Module \"%s\",id=%-2d: x,y,r,z: %7.4f,%7.4f,%6.4f,%8.4f cm",
               det_id, lay_id, module_name.c_str(), moduleIdx, pv.position().x() / cm, pv.position().y() / cm,
               sqrt(Xc * Xc + Yc * Yc) / cm, pv.position().z() / cm);

    }
  }

#ifdef DEBUG_MPGDCylinderBarrelTracker
  // Reset initial print level before exiting
  setPrintLevel(priorPrintLevel);
#endif

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

  return sdet;
}

//@}
// clang-format off
DECLARE_DETELEMENT(epic_CylinderMPGDBarrel, create_MPGDCylinderBarrelTracker)
