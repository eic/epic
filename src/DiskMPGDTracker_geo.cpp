// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2023 Wouter Deconinck, Matt Posik

/*
 * \brief Type: **Disk MPGD tracker segemented into two half disk modules
 * \author M. Posik
 * \ingroup Tracking
 *
 *
 * \code
 * \endcode
 *
 * @{
 */

#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/Printout.h"
#include "DD4hep/Shapes.h"
#include "DD4hepDetectorHelper.h"
#include "DDRec/DetectorData.h"
#include "DDRec/Surface.h"
#include "XML/Layering.h"
#include "XML/Utilities.h"
#include <array>

using namespace std;
using namespace dd4hep;
using namespace dd4hep::rec;

/** Disk MPGD detector with optional frame
 *
 * - Optional "frame" tag within the module element (frame eats part of the module width and length
 *   and surrounds the perimeter)
 * - Detector is setup as a "tracker" so we can use the hits
 *
 */
static Ref_t create_DiskMPGDTracker_geo(Detector& description, xml_h e, SensitiveDetector sens)
{
  typedef vector<PlacedVolume> Placements;
  xml_det_t                    x_det = e;
  int                                     det_id   = x_det.id();
  string                                  det_name = x_det.nameStr();
  DetElement                              sdet(det_name, det_id);
  map<string, Volume>                     volumes;
  map<string, Placements>                 sensitives;
  map<string, std::vector<rec::VolPlane>> volplane_surfaces;
  PlacedVolume                            pv;
  dd4hep::xml::Dimension                  dimensions(x_det.dimensions());
  Assembly                                assembly(det_name);

  // Set detector type flag
  dd4hep::xml::setDetectorTypeFlag(x_det, sdet);
  auto& params = DD4hepDetectorHelper::ensureExtension<dd4hep::rec::VariantParameters>(sdet);

  // Add the volume boundary material if configured
  for (xml_coll_t bmat(x_det, _Unicode(boundary_material)); bmat; ++bmat) {
    xml_comp_t x_boundary_material = bmat;
    DD4hepDetectorHelper::xmlToProtoSurfaceMaterial(x_boundary_material, params, "boundary_material");
  }

  map<string, std::array<double, 2>> module_thicknesses;
  sens.setType("tracker");

  // loop over the modules
  for (xml_coll_t mi(x_det, _U(module)); mi; ++mi) {
    xml_comp_t x_mod = mi;
    string     m_nam = x_mod.nameStr();

    if (volumes.find(m_nam) != volumes.end()) {
      printout(ERROR, "DiskMPGDTracker_geo",
               string((string("Module with named ") + m_nam + string(" already exists."))).c_str());
      throw runtime_error("Logics error in building modules.");
    }

    int    ncomponents     = 0;
    int    sensor_number   = 1;
    double total_thickness = 0;

    // Compute module total thickness from components
    xml_coll_t ci(x_mod, _U(module_component));
    for (ci.reset(), total_thickness = 0.0; ci; ++ci) {
      total_thickness += xml_comp_t(ci).thickness();
      std::cout << "Comp thickness[" << ci << "]: " << total_thickness << std::endl;
    
    }
    std::cout << "Detector Thickness: " << total_thickness << std::endl;

    // the module assembly volume
    Assembly m_vol(m_nam);
    volumes[m_nam] = m_vol;
    m_vol.setVisAttributes(description, x_mod.visStr());

    // Optional module frame.
    // frame consists of 4 parts:
    //  - around outer radial perimeter of disk
    //  - around the inner radial perimeter of the disk
    //  - 2 straight frame parts to complete the frame perimeter.
    //  The frame will eat the overlapping gas module area in the radial perimeter

    //set frame parameters to be used later to build the frame.
    string frame_mat = "";
    string frame_vis = "";
    double frame_width     = 0.0;
    double frame_thickness = 0.0;
    double frame_phi0      = 0.0;
    double frame_dphi      = 0.0;
    if (x_mod.hasChild(_U(frame))) {
      xml_comp_t m_frame = x_mod.child(_U(frame));
      frame_mat              = m_frame.materialStr();
      frame_vis              = m_frame.visStr();
      frame_width            = m_frame.width();
      frame_thickness        = m_frame.thickness();
      frame_phi0             = m_frame.phi0();
      frame_dphi             = m_frame.phi();
    }
    std::cout << "Detector Thickness + Frame Thickness: " << total_thickness + frame_thickness << std::endl;
    double thickness_so_far = 0.0;
    double thickness_sum    = -total_thickness / 2.0;
    double topframe_rmin    = 0.0; //top is upper radial frame
    double topframe_rmax    = 0.0;
    double botframe_rmin    = 0.0; //bot is lower radial frame
    double botframe_rmax    = 0.0;
    double gas_thickness    = 0.0;
    for (xml_coll_t mci(x_mod, _U(module_component)); mci; ++mci, ++ncomponents) {
      xml_comp_t x_comp    = mci;
      string     c_nam     = _toString(ncomponents, "component%d");
      string     comp_name = x_comp.nameStr();
      Tube    c_tube_tmp;
      Tube    beampipe_cutout;
      // Since MPGD frames are layed over the MPGD foils, the foil material is pressent under the frame.
      // The gas volumes are not present under the frames, so our frames must eat only the gas module areas
      // Example:
      //   ------------------- MPGD foil
      //   --               -- Frame
      //   --  gas volume   -- Frame
      //   --               -- Frame
      //   ------------------- MPGD foil

      // Look for gas modules to subtract radial frames from
      // FIXME: these module names are hard coded for now. Should find
      // a way to set a arribut via the moduel tag to flag what components
      // need to have frame thickness subtracted.
      printout(DEBUG, "DiskMPGDTracker_geo", "component: %s", comp_name.c_str());
      printout(DEBUG, "DiskMPGDTracker_geo", "x_comp.rmax(): %f", x_comp.rmax() );
      printout(DEBUG, "DiskMPGDTracker_geo", "x_comp.rmin(): %f", x_comp.rmin() );
      printout(DEBUG, "DiskMPGDTracker_geo", "x_comp.phi0(): %f", x_comp.phi0() );
      printout(DEBUG, "DiskMPGDTracker_geo", "x_comp.phi(): %f", x_comp.phi() );
      printout(DEBUG, "DiskMPGDTracker_geo", "x_comp.thickness(): %f", x_comp.thickness() );
      if ((comp_name == "DriftGap" || comp_name == "WindowGasGap")) {
        topframe_rmax  = x_comp.rmax();
        topframe_rmin  =  x_comp.rmax() - frame_width;
        botframe_rmax  = x_comp.rmin() + frame_width;
        botframe_rmin  =  x_comp.rmin();
        gas_thickness += x_comp.thickness();
        //create Tube shape with radial frame areas removed
        c_tube_tmp = { botframe_rmax, topframe_rmin, 0.5* x_comp.thickness(), -0.5*x_comp.phi() + x_comp.phi0(), 0.5*x_comp.phi() + x_comp.phi0()};
        beampipe_cutout = {0.0, x_comp.rmin()+5*cm, 0.5* x_comp.thickness(), -0.5*x_comp.phi() + x_comp.phi0(), 0.5*x_comp.phi() + x_comp.phi0()};
        std::cout << "Gas Volume, botframe_rmin: " << botframe_rmin << std::endl;
        std::cout << "Gas Volume, topframe_rmin: " << topframe_rmin << std::endl;
        std::cout << "Gas Volume, botframe_rmax: " << botframe_rmax << std::endl;
        std::cout << "Gas Volume, topframe_rmax: " << topframe_rmax << std::endl;
        std::cout << "Gas Volume, frame tube (rmin,rmax): " << botframe_rmin << " , " << topframe_rmin << std::endl;

        printout(DEBUG, "DiskMPGDTracker_geo", "topframe_rmax: %f", topframe_rmax);
        printout(DEBUG, "DiskMPGDTracker_geo", "topframe_rmin: %f", topframe_rmin);
        printout(DEBUG, "DiskMPGDTracker_geo", "botframe_rmax: %f", botframe_rmax);
        printout(DEBUG, "DiskMPGDTracker_geo", "botframe_rmin: %f", botframe_rmin);
      } else {
        //create Tube shape without radial frame areas removed
        c_tube_tmp = { x_comp.rmin(), x_comp.rmax(), 0.5* x_comp.thickness(), -0.5*x_comp.phi() + x_comp.phi0(), 0.5*x_comp.phi() + x_comp.phi0()};
        beampipe_cutout = {0.0, x_comp.rmin()+5*cm, 0.5* x_comp.thickness(), -0.5*x_comp.phi() + x_comp.phi0(), 0.5*x_comp.phi() + x_comp.phi0()};
        std::cout << "comp_rmin: " << x_comp.rmin() << std::endl;
        std::cout << "comp_rmax: " << x_comp.rmax() << std::endl;
      }
      //beampipe_cutout = {0.0, x_comp_rmin, 0.5* x_comp.thickness(), -0.5*x_comp.phi() + x_comp.phi0(), 0.5*x_comp.phi() + x_comp.phi0()};
      SubtractionSolid c_tube(c_tube_tmp, beampipe_cutout);
      Volume c_vol{c_nam, c_tube, description.material(x_comp.materialStr())};
      c_vol.setRegion(description, x_comp.regionStr());
      c_vol.setLimitSet(description, x_comp.limitsStr());
      c_vol.setVisAttributes(description, x_comp.visStr());

      pv = m_vol.placeVolume(c_vol, Position(0, 0, thickness_sum + 0.5*x_comp.thickness()));
      std::cout << "pv(Z) c_tube: " << comp_name.c_str() << " " << thickness_sum + 0.5*x_comp.thickness() << std::endl; 
      if (x_comp.isSensitive()) {
        pv.addPhysVolID("sensor", sensor_number++);
        c_vol.setSensitiveDetector(sens);
        sensitives[m_nam].push_back(pv);
        module_thicknesses[m_nam] = {thickness_so_far + 0.5*x_comp.thickness(),
                                     total_thickness - thickness_so_far - 0.5*x_comp.thickness()};
        // -------- create a measurement plane for the tracking surface attched to the sensitive volume -----
        Vector3D u(-1., 0., 0.);
        Vector3D v(0., -1., 0.);
        Vector3D n(0., 0., 1.);

        // compute the inner and outer thicknesses that need to be assigned to the tracking surface
        // depending on wether the support is above or below the sensor
        double inner_thickness = module_thicknesses[m_nam][0];
        double outer_thickness = module_thicknesses[m_nam][1];

        SurfaceType type(rec::SurfaceType::Sensitive);

        VolPlane surf(c_vol, type, inner_thickness, outer_thickness, u, v, n); //,o ) ;
        volplane_surfaces[m_nam].push_back(surf);
      }
      thickness_sum    += x_comp.thickness();
      thickness_so_far += x_comp.thickness();
    }
      std::cout << "total thickness: " << total_thickness  << std::endl;
      std::cout << "gas thickness: " << gas_thickness  << std::endl;
      std::cout << "frame thickness: " << frame_thickness  << std::endl;
      std::cout << "Detector Thickness + Frame Thickness - gas: " << total_thickness + frame_thickness - gas_thickness  << std::endl;
    //Build the frame: 2 radial segments + 2 straight segments
    if (x_mod.hasChild(_U(frame))) {
      double botframe_length = topframe_rmax - botframe_rmin;
      double frame_angle = 0.5*frame_dphi - frame_phi0;
      double x0 = botframe_rmin*TMath::Cos(frame_angle);
      double y0 = botframe_rmin*TMath::Sin(frame_angle);
      double xf = topframe_rmax*TMath::Cos(frame_angle); 
      double yf = topframe_rmax*TMath::Sin(frame_angle);
      double xp = xf - 0.5*botframe_length;
      double yp = xp*TMath::Tan(0.5*frame_dphi - frame_phi0) + 0.5*frame_width;
      //double yp = xp*TMath::Tan(0.5*frame_dphi - frame_phi0) + 0.5*frame_width;
      
      Tube topframe_tube{topframe_rmin, topframe_rmax, 0.5*frame_thickness, -0.5*frame_dphi + frame_phi0, 0.5*frame_dphi + frame_phi0};
      Tube botframe_tube{botframe_rmin, botframe_rmax, 0.5*frame_thickness, -0.5*frame_dphi + frame_phi0, 0.5*frame_dphi + frame_phi0};
      Box sideframe1_box{0.5*frame_width, 0.5*botframe_length, 0.5*frame_thickness};
      Box sideframe2_box{0.5*frame_width, 0.5*botframe_length, 0.5*frame_thickness};

      Volume topframe_vol{"topframe", topframe_tube, description.material(frame_mat)};
      Volume botframe_vol{"botframe", botframe_tube, description.material(frame_mat)};
      Volume sideframe1_vol{"sideframe1", sideframe1_box, description.material(frame_mat)};
      Volume sideframe2_vol{"sideframe2", sideframe2_box, description.material(frame_mat)};

      topframe_vol.setVisAttributes(description, frame_vis);
      botframe_vol.setVisAttributes(description, frame_vis);
      sideframe1_vol.setVisAttributes(description, frame_vis);
      sideframe2_vol.setVisAttributes(description, frame_vis);

      printout(DEBUG, "DiskMPGDTracker_geo", "frame_thickness: %f", frame_thickness);
      printout(DEBUG, "DiskMPGDTracker_geo", "total_thickness (no frames): %f", total_thickness);
      printout(DEBUG, "DiskMPGDTracker_geo", "frame_thickness + total_thickness: %f", frame_thickness + total_thickness);
      printout(DEBUG, "DiskMPGDTracker_geo", "botframe_length: %f", botframe_length);
      printout(DEBUG, "DiskMPGDTracker_geo", "frame_angle: %f", frame_angle);
      printout(DEBUG, "DiskMPGDTracker_geo", "(x0, y0): (%f, %f) ",x0,y0);
      printout(DEBUG, "DiskMPGDTracker_geo", "(xf, yf): (%f, %f)", xf,yf);
      printout(DEBUG, "DiskMPGDTracker_geo", "(xp, yp): (%f, %f)", xp,yp);

      std::cout << "frame_angle: " << frame_angle << std::endl;
      std::cout << "botframe_length: " << botframe_length << std::endl;
      std::cout << "x0,y0: " << x0 << ", " << y0 << std::endl;
      std::cout << "xf,yf: " << xf << ", " << yf << std::endl;
      std::cout << "xp,yp: " << xp << ", " << yp << std::endl;

      m_vol.placeVolume(topframe_vol, Position(0.0, 0.0, 0.5*(frame_thickness - total_thickness - gas_thickness)));
      m_vol.placeVolume(botframe_vol, Position(0.0, 0.0, 0.5*(frame_thickness - total_thickness - gas_thickness)));
      m_vol.placeVolume(sideframe1_vol, Transform3D(RotationZ(frame_angle + 270.0*deg),
                        Position(-xp, 
                                 -yp - 0.005*cm, 
                                  0.5*(frame_thickness - total_thickness - gas_thickness))));
      m_vol.placeVolume(sideframe2_vol, Transform3D(RotationZ(-frame_angle - 270.0*deg),
                        Position(xp, 
                                 -yp - 0.005*cm, 
                                  0.5*(frame_thickness - total_thickness - gas_thickness))));
    }
  }

  // build the layers to place the modules around
  for (xml_coll_t li(x_det, _U(layer)); li; ++li) {
    xml_comp_t x_layer            = li;
    xml_comp_t x_layout           = x_layer.child(_U(rphi_layout));
    xml_comp_t z_layout           = x_layer.child(_U(z_layout));
    int        lay_id             = x_layer.id();
    string     m_nam              = x_layer.moduleStr();
    string     lay_nam            = det_name + _toString(x_layer.id(), "_layer%d");
    xml_comp_t envelope_tolerance = x_layer.child(_Unicode(envelope_tolerance), false);
    double     envelope_r_min     = 0;
    double     envelope_r_max     = 0;
    double     envelope_z_min     = 0;
    double     envelope_z_max     = 0;
    if (envelope_tolerance) {
      envelope_r_min = getAttrOrDefault(envelope_tolerance, _Unicode(r_min), 0);
      envelope_r_max = getAttrOrDefault(envelope_tolerance, _Unicode(r_max), 0);
      envelope_z_min = getAttrOrDefault(envelope_tolerance, _Unicode(z_min), 0);
      envelope_z_max = getAttrOrDefault(envelope_tolerance, _Unicode(z_max), 0);
    }
    double phi0     = x_layout.phi0();     // starting phi of first module
    double phi_tilt = x_layout.phi_tilt(); // Phi tilit of module
    double rc       = x_layout.rc();       // Radius of the module
    int    nphi     = x_layout.nphi();     // Number of modules in phi
    double rphi_dr  = x_layout.dr();       // The delta radius of every other module
    double phi_incr = (2 * M_PI) / nphi;   // Phi increment for one module
    double phic     = phi0;                // Phi of the module
    double z_dr     = z_layout.dr();       // Radial offest of modules in z
    double z0       = z_layout.z0();       // Central Z position of module
    double z_offset = 0.5*z_layout.dz();   // Distance in z between nphi moduels
    Assembly    layer_assembly(lay_nam);
    Volume      module_env = volumes[m_nam];
    DetElement  lay_elt(sdet, lay_nam, lay_id);
    Placements& sensVols    = sensitives[m_nam];
    auto&       layerParams = DD4hepDetectorHelper::ensureExtension<dd4hep::rec::VariantParameters>(lay_elt);
    pv = assembly.placeVolume(layer_assembly);
    pv.addPhysVolID("layer", lay_id);
    lay_elt.setPlacement(pv);

    int module = 1;
    double module_sign = 1; //used to alternate half disk orientation: +1 = disk up, -1 = disk down
    // loop over the modules in phi
    for (int ii = 0; ii < nphi; ii++) {
      std::cout<< "module_sign: " << module_sign << std::endl; //MP
      double xc = rc * std::cos(phic);              // Basic x position of module
      double yc = rc * std::sin(phic);              // Basic y position of module
      double dx = z_dr * std::cos(phic + phi_tilt); // Deta x of module position
      double dy = z_dr * std::sin(phic + phi_tilt); // Deta y of module position
      string     module_name = _toString(module, "module%d");
      DetElement mod_elt(lay_elt, module_name, module);
      Transform3D tr(RotationZYX( ((M_PI/2) - phic -phi_tilt), 0.0, 0.0),
                     Position(xc, yc,  z0 + module_sign*z_offset)); // in x-y plane,
      std::cout << tr << std::endl;
      pv = layer_assembly.placeVolume(module_env, tr);
      pv.addPhysVolID("module", module);
      mod_elt.setPlacement(pv);
      for (size_t ic = 0; ic < sensVols.size(); ++ic) {
        PlacedVolume sens_pv = sensVols[ic];
        DetElement   comp_de(mod_elt, std::string("de_") + sens_pv.volume().name(), module);
        comp_de.setPlacement(sens_pv);
      }
      //increment counters and alternate relevant signs
      module++;
      xc += dx;
      yc += dy;
      module_sign *= -1.0;
      phic += phi_incr;
      rc += rphi_dr;
    }
    layer_assembly->GetShape()->ComputeBBox();
    layerParams.set<double>("envelope_r_min", envelope_r_min);
    layerParams.set<double>("envelope_r_max", envelope_r_max);
    layerParams.set<double>("envelope_z_min", envelope_z_min);
    layerParams.set<double>("envelope_z_max", envelope_z_max);

    for (xml_coll_t lmat(x_layer, _Unicode(layer_material)); lmat; ++lmat) {
      xml_comp_t x_layer_material = lmat;
      DD4hepDetectorHelper::xmlToProtoSurfaceMaterial(x_layer_material, layerParams, "layer_material");
    }
  }
  sdet.setAttributes(description, assembly, x_det.regionStr(), x_det.limitsStr(), x_det.visStr());
  assembly.setVisAttributes(description.invisible());
  pv = description.pickMotherVolume(sdet).placeVolume(assembly);
  pv.addPhysVolID("system", det_id); // Set the subdetector system ID
  sdet.setPlacement(pv);
  return sdet;
}

//@}
// clang-format off
DECLARE_DETELEMENT(epic_DiskMPGDTrackerN, create_DiskMPGDTracker_geo)
