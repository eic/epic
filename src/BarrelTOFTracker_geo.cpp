// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 - 2024, Whitney Armstrong, Chun Yuen Tsang

/** \addtogroup Trackers Trackers
 * \brief Type: **TOFBarrel**.
 * \author W. Armstrong
 * \modified by C.Y Tsang 3rd Aug, 2024
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
#include <iostream>
#include "DD4hepDetectorHelper.h"

using namespace std;
using namespace dd4hep;
using namespace dd4hep::rec;
using namespace dd4hep::detail;

/** Barrel Tracker with space frame.
 *
 * - Optional "support" tag within the detector element.
 *
 * The shapes are created using createShape which can be one of many basic geomtries.
 * See the examples Check_shape_*.xml in
 * [dd4hep's examples/ClientTests/compact](https://github.com/AIDASoft/DD4hep/tree/master/examples/ClientTests/compact)
 * directory.
 *
 *
 * - Optional "frame" tag within the module element.
 *
 * \ingroup trackers
 *
 * \code
 * \endcode
 *
 *
 * @author Whitney Armstrong
 */
static Ref_t create_TOFBarrel(Detector& description, xml_h e, SensitiveDetector sens) {
  typedef vector<PlacedVolume> Placements;
  xml_det_t x_det = e;
  Material air    = description.air();
  int det_id      = x_det.id();
  string det_name = x_det.nameStr();
  DetElement sdet(det_name, det_id);

  map<string, Volume> volumes;
  map<string, Placements> sensitives;
  map<string, std::vector<VolPlane>> volplane_surfaces;
  map<string, std::array<double, 2>> module_thicknesses;

  PlacedVolume pv;

  // Set detector type flag
  dd4hep::xml::setDetectorTypeFlag(x_det, sdet);
  auto& params = DD4hepDetectorHelper::ensureExtension<dd4hep::rec::VariantParameters>(sdet);

  // Add the volume boundary material if configured
  for (xml_coll_t bmat(x_det, _Unicode(boundary_material)); bmat; ++bmat) {
    xml_comp_t x_boundary_material = bmat;
    DD4hepDetectorHelper::xmlToProtoSurfaceMaterial(x_boundary_material, params,
                                                    "boundary_material");
  }

  // dd4hep::xml::Dimension dimensions(x_det.dimensions());
  // Tube topVolumeShape(dimensions.rmin(), dimensions.rmax(), dimensions.length() * 0.5);
  // Volume assembly(det_name,topVolumeShape,air);
  Assembly assembly(det_name);

  sens.setType("tracker");

  // Loop over the suports
  for (xml_coll_t su(x_det, _U(support)); su; ++su) {
    xml_comp_t x_support     = su;
    double support_thickness = getAttrOrDefault(x_support, _U(thickness), 2.0 * mm);
    double support_length    = getAttrOrDefault(x_support, _U(length), 2.0 * mm);
    double support_rmin      = getAttrOrDefault(x_support, _U(rmin), 2.0 * mm);
    double support_zstart    = getAttrOrDefault(x_support, _U(zstart), 2.0 * mm);
    std::string support_name =
        getAttrOrDefault<std::string>(x_support, _Unicode(name), "support_tube");
    std::string support_vis = getAttrOrDefault<std::string>(x_support, _Unicode(vis), "AnlRed");
    xml_dim_t pos(x_support.child(_U(position), false));
    xml_dim_t rot(x_support.child(_U(rotation), false));
    Solid support_solid;
    if (x_support.hasChild(_U(shape))) {
      xml_comp_t shape(x_support.child(_U(shape)));
      string shape_type = shape.typeStr();
      support_solid     = xml::createShape(description, shape_type, shape);
    } else {
      support_solid = Tube(support_rmin, support_rmin + support_thickness, support_length / 2);
    }
    Transform3D tr =
        Transform3D(Rotation3D(), Position(0, 0, (support_zstart + support_length / 2)));
    if (pos.ptr() && rot.ptr()) {
      Rotation3D rot3D(RotationZYX(rot.z(0), rot.y(0), rot.x(0)));
      Position pos3D(pos.x(0), pos.y(0), pos.z(0));
      tr = Transform3D(rot3D, pos3D);
    } else if (pos.ptr()) {
      tr = Transform3D(Rotation3D(), Position(pos.x(0), pos.y(0), pos.z(0)));
    } else if (rot.ptr()) {
      Rotation3D rot3D(RotationZYX(rot.z(0), rot.y(0), rot.x(0)));
      tr = Transform3D(rot3D, Position());
    }
    Material support_mat = description.material(x_support.materialStr());
    Volume support_vol(support_name, support_solid, support_mat);
    support_vol.setVisAttributes(description.visAttributes(support_vis));
    pv = assembly.placeVolume(support_vol, tr);
    // pv = assembly.placeVolume(support_vol, Position(0, 0, support_zstart + support_length / 2));
  }

  // loop over the modules
  for (xml_coll_t mi(x_det, _U(module)); mi; ++mi) {
    xml_comp_t x_mod = mi;
    string m_nam     = x_mod.nameStr();

    if (volumes.find(m_nam) != volumes.end()) {
      printout(ERROR, "TOFBarrel",
               string((string("Module with named ") + m_nam + string(" already exists."))).c_str());
      throw runtime_error("Logics error in building modules.");
    }

    int ncomponents        = 0;
    int sensor_number      = 1;
    double total_thickness = 0;

    // Compute module total thickness from components
    xml_coll_t ci(x_mod, _U(module_component));
    for (ci.reset(), total_thickness = 0.0; ci; ++ci) {
      total_thickness += xml_comp_t(ci).thickness();
    }
    // the module assembly volume
    Assembly m_vol(m_nam);
    volumes[m_nam] = m_vol;
    m_vol.setVisAttributes(description.visAttributes(x_mod.visStr()));

    // Optional module frame.
    if (x_mod.hasChild(_U(frame))) {
      xml_comp_t m_frame = x_mod.child(_U(frame));
      // xmleles[m_nam]  = x_mod;
      double frame_thickness = m_frame.thickness();
      double frame_width     = m_frame.width();
      double frame_height    = getAttrOrDefault<double>(m_frame, _U(height), 5.0 * mm);
      double tanth           = frame_height / (frame_width / 2.0);
      double costh           = 1. / sqrt(1 + tanth * tanth);
      double frame_height2   = frame_height - frame_thickness - frame_thickness / costh;
      double frame_width2    = 2.0 * frame_height2 / tanth;

      Trd1 moduleframe_part1(frame_width / 2, 0.001 * mm, m_frame.length() / 2, frame_height / 2);
      Trd1 moduleframe_part2(frame_width2 / 2, 0.001 * mm, m_frame.length() / 2 + 0.01 * mm,
                             frame_height2 / 2);

      SubtractionSolid moduleframe(moduleframe_part1, moduleframe_part2,
                                   Position(0.0, frame_thickness, 0.0));
      Volume v_moduleframe(m_nam + "_vol", moduleframe,
                           description.material(m_frame.materialStr()));
      v_moduleframe.setVisAttributes(description, m_frame.visStr());
      m_vol.placeVolume(v_moduleframe,
                        Position(0.0, 0.0, frame_height / 2 + total_thickness / 2.0));
    }

    double thickness_so_far = 0.0;
    double thickness_sum    = -total_thickness / 2.0;
    for (xml_coll_t mci(x_mod, _U(module_component)); mci; ++mci) {
      xml_comp_t x_comp = mci;
      xml_comp_t x_pos  = x_comp.position(false);
      xml_comp_t x_rot  = x_comp.rotation(false);
      auto make_box = [&](double width, double length, double thickness, 
		          double pos_x = 0, double pos_y = 0, double pos_z = 0, 
			  double rot_x = 0, double rot_y = 0, double rot_z = 0, bool z_stacking = true) {
        // Utility variable for the relative z-offset based off the previous components
        const double zoff = thickness_sum + thickness / 2.0;

        const string c_nam = _toString(ncomponents, "component%d");
        ++ncomponents;
        Box c_box(width / 2, length / 2, thickness / 2);
        Volume c_vol;
	
        xml_coll_t ci(x_comp, _Unicode(inner_tube));
	if(ci) {
	  double max_r = 0;
          for (; ci; ++ci) {
            // fill the hole with tube
	    xml_comp_t ct = ci;
	    max_r = std::max(max_r, ct.rmax());
       	    Tube c_tube(ct.rmin(), ct.rmax(), length / 2);
            Volume c_tubevol(c_nam + ct.nameStr(), c_tube, description.material(ct.materialStr()));
	    if(ct.visStr() != "") c_tubevol.setVisAttributes(description, ct.visStr());
            m_vol.placeVolume(c_tubevol, Transform3D(RotationZYX(0, 0, -M_PI / 2),
                                                        Position(pos_x, pos_y, pos_z + zoff)));
	  }

          Tube c_fbox(0, max_r, length / 2 + 1);
	  SubtractionSolid c_sbox(c_box, c_fbox, Transform3D(RotationZYX(0, 0, -M_PI / 2),
                                                  Position(0, 0, 0)));//pos_x, pos_y, pos_z + zoff)));

          c_vol = Volume(c_nam, c_sbox, description.material(x_comp.materialStr()));
	} else c_vol = Volume(c_nam, c_box, description.material(x_comp.materialStr()));

	Volume test;
	test = c_vol;


	// center if off by half the box length if box length is cut in half
        Position c_pos(pos_x, pos_y, pos_z + zoff);
        RotationZYX c_rot(rot_z, rot_y, rot_x);
        pv = m_vol.placeVolume(c_vol, Transform3D(c_rot, c_pos));


        c_vol.setRegion(description, x_comp.regionStr());
        c_vol.setLimitSet(description, x_comp.limitsStr());
        c_vol.setVisAttributes(description, x_comp.visStr());
        if (x_comp.isSensitive()) {
          pv.addPhysVolID("sensor", sensor_number++);
          c_vol.setSensitiveDetector(sens);
          sensitives[m_nam].push_back(pv);
          module_thicknesses[m_nam] = {thickness_so_far + thickness / 2.0,
                                       total_thickness - thickness_so_far -
                                           thickness / 2.0};

          // -------- create a measurement plane for the tracking surface attched to the sensitive volume -----
          Vector3D u(-1., 0., 0.);
          Vector3D v(0., -1., 0.);
          Vector3D n(0., 0., 1.);
          //    Vector3D o( 0. , 0. , 0. ) ;

          // compute the inner and outer thicknesses that need to be assigned to the tracking surface
          // depending on wether the support is above or below the sensor
          double inner_thickness = module_thicknesses[m_nam][0];
          double outer_thickness = module_thicknesses[m_nam][1];

          SurfaceType type(SurfaceType::Sensitive);

          // if( isStripDetector )
          //  type.setProperty( SurfaceType::Measurement1D , true ) ;

          VolPlane surf(c_vol, type, inner_thickness, outer_thickness, u, v, n); //,o ) ;
          volplane_surfaces[m_nam].push_back(surf);

          //--------------------------------------------
        }
        if (z_stacking) {
          thickness_sum += thickness;
          thickness_so_far += thickness;
          // apply relative offsets in z-position used to stack components side-by-side
          thickness_sum += pos_z;
          thickness_so_far += pos_z;
        }
      };

      double pos_x = 0, pos_y = 0, pos_z = 0;
      double rot_x = 0, rot_y = 0, rot_z = 0;
      if (x_rot) {
        rot_x = x_rot.x(0);
        rot_y = x_rot.y(0);
        rot_z = x_rot.z(0);
      }
      if (x_pos) {
        pos_x = x_pos.x(0);
        pos_y = x_pos.y(0);
        pos_z = x_pos.z(0);
      }
      double width = x_comp.width();
      double length = x_comp.length();
      double thickness = x_comp.thickness();

      if (x_comp.hasChild(_Unicode(GridSensors))) {
        auto x_comp_t = x_comp.child(_Unicode(GridSensors));
        // x-distance between centers of neighboring sensors
        double sensors_xdist = getAttrOrDefault<double>(x_comp_t, _Unicode(xdist), x_comp.width());
        // y-distance between centers of neighboring sensors
        double sensors_ydist = getAttrOrDefault<double>(x_comp_t, _Unicode(ydist), x_comp.length());
        // number of rows of sensors in a stave
        int nsensors_x = getAttrOrDefault<int>(x_comp_t, _Unicode(nx), 1);
        // number of column of sensors in a stave
        int nsensors_y = getAttrOrDefault<int>(x_comp_t, _Unicode(ny), 1);
        // x-location of the center of the leftmost sensor
        double start_x = getAttrOrDefault<double>(x_comp_t, _Unicode(start_x), 0);
        // y-location of the center of the uppermost sensor
        double start_y = getAttrOrDefault<double>(x_comp_t, _Unicode(start_y), 0);
        // z-locatino of the center of all sensors (All sensors appears at the same z-layer
        double start_z = getAttrOrDefault<double>(x_comp_t, _Unicode(start_z), 0);
        // central ring is located to the right of the ny_before_ring th sensor
        int ny_before_ring = getAttrOrDefault<int>(x_comp_t, _Unicode(ny_before_ring), 0);
        // Extra width caused by the ring
        // |<--sensors_ydist-->|<--sensors_ydist-->|<-----ring_extra_width------->|<--sensors_ydist-->|
        //   |<xcomp.width()>|   |<xcomp.width()>|                                   |<xcomp.width()>|
        //   ring_extra_width is the extra width between boundaries of the sensor boundaries (including dead space)
        double ring_extra_width = getAttrOrDefault<double>(x_comp_t, _Unicode(ring_extra_width), 0);
	auto half_length_str = getAttrOrDefault<std::string>(x_comp_t, _Unicode(half_length), "none");


        double current_x = start_x;
        for (int nx = 0; nx < nsensors_x; ++nx) {
          double current_y = start_y;
          for (int ny = 0; ny < nsensors_y; ++ny) {
            double sensor_length = length;
	    double tmp_sensors_ydist = sensors_ydist;
	    // when we draw half a sensor, the center has to be shifted by 0.25 times the length of a sensor
	    // distance between centers to the next sensor also has to be reduced by 0.25 times the length of a sensor
	    if(half_length_str == "left" && ny == 0) {
              sensor_length = 0.5*length;
	      current_y += 0.25*length;
	      tmp_sensors_ydist -= 0.25*length;
	    }
	    // same idea, but when you are drawing to the right, the right sensor center has to move in -y direction
	    if(half_length_str == "right" && ny == nsensors_y-1) {
              sensor_length = 0.5*length;
	      current_y -= 0.25*length;
	    }
            make_box(width, sensor_length, thickness, 
		     current_x, current_y, start_z, 
		     rot_x, rot_y, rot_z,
                     ((nx == nsensors_x - 1) && (ny == nsensors_y - 1))); // all sensors are located at the same z-layer
            // increment z-layers only at the end, after the last sensor is added
            current_y += tmp_sensors_ydist;
            if (ny + 1 == ny_before_ring)
              current_y += ring_extra_width;
          }
          current_x += sensors_xdist;
        }
      } else
        make_box(width, length, thickness, pos_x, pos_y, pos_z, rot_x, rot_y, rot_z);
    }
  }

  // now build the layers
  for (xml_coll_t li(x_det, _U(layer)); li; ++li) {
    xml_comp_t x_layer  = li;
    xml_comp_t x_barrel = x_layer.child(_U(barrel_envelope));
    xml_comp_t x_layout = x_layer.child(_U(rphi_layout));
    xml_comp_t z_layout = x_layer.child(_U(z_layout)); // Get the <z_layout> element.
    int lay_id          = x_layer.id();
    string m_nam        = x_layer.moduleStr();
    string lay_nam      = det_name + _toString(x_layer.id(), "_layer%d");
    Tube lay_tub(x_barrel.inner_r(), x_barrel.outer_r(), x_barrel.z_length() / 2.0);
    Volume lay_vol(lay_nam, lay_tub, air); // Create the layer envelope volume.
    Position lay_pos(0, 0, getAttrOrDefault(x_barrel, _U(z0), 0.));
    lay_vol.setVisAttributes(description.visAttributes(x_layer.visStr()));

    double phi0     = x_layout.phi0();     // Starting phi of first module.
    double phi_tilt = x_layout.phi_tilt(); // Phi tilt of a module.
    double rc       = x_layout.rc();       // Radius of the module center.
    int nphi        = x_layout.nphi();     // Number of modules in phi.
    double rphi_dr  = x_layout.dr();       // The delta radius of every other module.
    double phi_incr = (M_PI * 2) / nphi;   // Phi increment for one module.
    double phic     = phi0;                // Phi of the module center.
    double z0       = z_layout.z0();       // Z position of first module in phi.
    double nz       = z_layout.nz();       // Number of modules to place in z.
    double z_dr     = z_layout.dr();       // Radial displacement parameter, of every other module.

    Volume module_env = volumes[m_nam];
    DetElement lay_elt(sdet, lay_nam, lay_id);
    Placements& sensVols = sensitives[m_nam];

    // the local coordinate systems of modules in dd4hep and acts differ
    // see http://acts.web.cern.ch/ACTS/latest/doc/group__DD4hepPlugins.html
    auto& layerParams =
        DD4hepDetectorHelper::ensureExtension<dd4hep::rec::VariantParameters>(lay_elt);

    for (xml_coll_t lmat(x_layer, _Unicode(layer_material)); lmat; ++lmat) {
      xml_comp_t x_layer_material = lmat;
      DD4hepDetectorHelper::xmlToProtoSurfaceMaterial(x_layer_material, layerParams,
                                                      "layer_material");
    }

    // Z increment for module placement along Z axis.
    // Adjust for z0 at center of module rather than
    // the end of cylindrical envelope.
    double z_incr = nz > 1 ? (2.0 * z0) / (nz - 1) : 0.0;
    // Starting z for module placement along Z axis.
    double module_z = -z0;
    int module      = 1;

    // Loop over the number of modules in phi.
    for (int ii = 0; ii < nphi; ii++) {
      double dx = z_dr * std::cos(phic + phi_tilt); // Delta x of module position.
      double dy = z_dr * std::sin(phic + phi_tilt); // Delta y of module position.
      double x  = rc * std::cos(phic);              // Basic x module position.
      double y  = rc * std::sin(phic);              // Basic y module position.

      // Loop over the number of modules in z.
      for (int j = 0; j < nz; j++) {
        string module_name = _toString(module, "module%d");
        DetElement mod_elt(lay_elt, module_name, module);

        Transform3D tr(RotationZYX(0, ((M_PI / 2) - phic - phi_tilt), -M_PI / 2),
                       Position(x, y, module_z));

        pv = lay_vol.placeVolume(module_env, tr);
        pv.addPhysVolID("module", module);
        mod_elt.setPlacement(pv);
        for (size_t ic = 0; ic < sensVols.size(); ++ic) {
          PlacedVolume sens_pv = sensVols[ic];
          DetElement comp_de(mod_elt, std::string("de_") + sens_pv.volume().name(), module);
          comp_de.setPlacement(sens_pv);

          auto& comp_de_params =
              DD4hepDetectorHelper::ensureExtension<dd4hep::rec::VariantParameters>(comp_de);
          comp_de_params.set<string>("axis_definitions", "XYZ");
          // comp_de.setAttributes(description, sens_pv.volume(), x_layer.regionStr(), x_layer.limitsStr(),
          //                       xml_det_t(xmleles[m_nam]).visStr());
          //

          volSurfaceList(comp_de)->push_back(volplane_surfaces[m_nam][ic]);
        }

        /// Increase counters etc.
        module++;
        // Adjust the x and y coordinates of the module.
        x += dx;
        y += dy;
        // Flip sign of x and y adjustments.
        dx *= -1;
        dy *= -1;
        // Add z increment to get next z placement pos.
        module_z += z_incr;
      }
      phic += phi_incr; // Increment the phi placement of module.
      rc += rphi_dr;    // Increment the center radius according to dr parameter.
      rphi_dr *= -1;    // Flip sign of dr parameter.
      module_z = -z0;   // Reset the Z placement parameter for module.
    }
    // Create the PhysicalVolume for the layer.
    pv = assembly.placeVolume(lay_vol, lay_pos); // Place layer in mother
    pv.addPhysVolID("layer", lay_id);            // Set the layer ID.
    lay_elt.setAttributes(description, lay_vol, x_layer.regionStr(), x_layer.limitsStr(),
                          x_layer.visStr());
    lay_elt.setPlacement(pv);
  }
  sdet.setAttributes(description, assembly, x_det.regionStr(), x_det.limitsStr(), x_det.visStr());
  assembly.setVisAttributes(description.invisible());
  pv = description.pickMotherVolume(sdet).placeVolume(assembly);
  pv.addPhysVolID("system", det_id); // Set the subdetector system ID.
  sdet.setPlacement(pv);
  return sdet;
}

//@}
// clang-format off
DECLARE_DETELEMENT(epic_TOFBarrel,       create_TOFBarrel)
