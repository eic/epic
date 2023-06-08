// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Nicolas Schmidt

/** \addtogroup Trackers Trackers
 * \brief Type: **Endcap Tracker with TOF**.
 * \author N. Schmidt
 *
 * \ingroup trackers
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
#include <map>

using namespace std;
using namespace dd4hep;
using namespace dd4hep::rec;
using namespace dd4hep::detail;

static Ref_t create_detector(Detector& description, xml_h e, SensitiveDetector sens)
{
  typedef vector<PlacedVolume> Placements;

  xml_det_t    x_det    = e;
  int          det_id   = x_det.id();
  std::string  det_name = x_det.nameStr();
  DetElement   sdet(det_name, det_id);
  Material     air = description.material("Air");
  PlacedVolume pv;

  map<string, Volume>                volumes;
  map<string, Placements>            sensitives;
  map<string, std::vector<VolPlane>> volplane_surfaces;
  map<string, std::array<double, 2>> module_thicknesses;

  // Set detector type flag
  dd4hep::xml::setDetectorTypeFlag(x_det, sdet);
  auto& params = DD4hepDetectorHelper::ensureExtension<dd4hep::rec::VariantParameters>(sdet);

  // Add the volume boundary material if configured
  for (xml_coll_t bmat(x_det, _Unicode(boundary_material)); bmat; ++bmat) {
    xml_comp_t x_boundary_material = bmat;
    DD4hepDetectorHelper::xmlToProtoSurfaceMaterial(x_boundary_material, params, "boundary_material");
  }

  // xml_comp_t  x_mod         = x_det.child(_Unicode(module));
  // std::string m_nam         = x_mod.nameStr();
  // xml_comp_t  diskdimension = x_mod.dimensions();

  // xml_comp_t  envelope      = x_det.child(_Unicode(envelope), false);
  // double envelope_r_min     = 0;
  // double envelope_r_max     = 0;
  // double envelope_z_min     = 0;
  // double envelope_z_max     = 0;
  // if(envelope){
  //   envelope_r_min = getAttrOrDefault(envelope, _Unicode(r_min), 0);
  //   envelope_r_max = getAttrOrDefault(envelope, _Unicode(r_max), 0);
  //   envelope_z_min = getAttrOrDefault(envelope, _Unicode(z_min), 0);
  //   envelope_z_max = getAttrOrDefault(envelope, _Unicode(z_max), 0);
  // }

  // // load all information from the diskdimension definitions
  // double disk_zPos              = getAttrOrDefault(diskdimension, _Unicode(zPos), 0.);
  // double disk_rMin              = getAttrOrDefault(diskdimension, _Unicode(rMin), 0.);
  // double disk_rMax              = getAttrOrDefault(diskdimension, _Unicode(rMax), 100.);
  // double disk_xOffset           = getAttrOrDefault(diskdimension, _Unicode(xOffset), 0.);
  // double disk_det_height        = getAttrOrDefault(diskdimension, _Unicode(det_height), 0.);
  // double cooling_plate_height   = getAttrOrDefault(diskdimension, _Unicode(cooling_plate_height), 0.);
  // double cooling_tube_thickness = getAttrOrDefault(diskdimension, _Unicode(wallthickness_coolingtube), 0.);
  // double cooling_tube_diameter  = getAttrOrDefault(diskdimension, _Unicode(diameter_coolingtube), 0.);

  Assembly assembly(det_name);
  assembly.setVisAttributes(description.invisible());
  sens.setType("tracker");

  // loop over the modules
  for (xml_coll_t mi(x_det, _U(module)); mi; ++mi) {
    xml_comp_t x_mod = mi;
    string     m_nam = x_mod.nameStr();

    if (volumes.find(m_nam) != volumes.end()) {
      printout(ERROR, "EndcapTOF", string((string("Module with named ") + m_nam + string(" already exists."))).c_str());
      throw runtime_error("Logics error in building modules.");
    }

    int    ncomponents     = 0;
    int    sensor_number   = 1;
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

    double thickness_so_far = 0.0;
    double thickness_sum    = -total_thickness / 2.0;
    for (xml_coll_t mci(x_mod, _U(module_component)); mci; ++mci, ++ncomponents) {
      xml_comp_t   x_comp = mci;
      xml_comp_t   x_pos  = x_comp.position(false);
      xml_comp_t   x_rot  = x_comp.rotation(false);
      const string c_nam  = _toString(ncomponents, "component%d");
      Box          c_box(x_comp.width() / 2, x_comp.length() / 2, x_comp.thickness() / 2);
      Volume       c_vol(c_nam, c_box, description.material(x_comp.materialStr()));

      // Utility variable for the relative z-offset based off the previous components
      const double zoff = thickness_sum + x_comp.thickness() / 2.0;
      if (x_pos && x_rot) {
        Position    c_pos(x_pos.x(0), x_pos.y(0), x_pos.z(0) + zoff);
        RotationZYX c_rot(x_rot.z(0), x_rot.y(0), x_rot.x(0));
        pv = m_vol.placeVolume(c_vol, Transform3D(c_rot, c_pos));
      } else if (x_rot) {
        Position c_pos(0, 0, zoff);
        pv = m_vol.placeVolume(c_vol, Transform3D(RotationZYX(x_rot.z(0), x_rot.y(0), x_rot.x(0)), c_pos));
      } else if (x_pos) {
        pv = m_vol.placeVolume(c_vol, Position(x_pos.x(0), x_pos.y(0), x_pos.z(0) + zoff));
      } else {
        pv = m_vol.placeVolume(c_vol, Position(0, 0, zoff));
      }
      c_vol.setRegion(description, x_comp.regionStr());
      c_vol.setLimitSet(description, x_comp.limitsStr());
      c_vol.setVisAttributes(description, x_comp.visStr());
      if (x_comp.isSensitive()) {
        pv.addPhysVolID("sensor", sensor_number++);
        c_vol.setSensitiveDetector(sens);
        sensitives[m_nam].push_back(pv);
        module_thicknesses[m_nam] = {thickness_so_far + x_comp.thickness() / 2.0,
                                     total_thickness - thickness_so_far - x_comp.thickness() / 2.0};

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
      thickness_sum += x_comp.thickness();
      thickness_so_far += x_comp.thickness();
      // apply relative offsets in z-position used to stack components side-by-side
      if (x_pos) {
        thickness_sum += x_pos.z(0);
        thickness_so_far += x_pos.z(0);
      }
    }
  }
  float zPos = 0;
  // now build the layers
  for (xml_coll_t li(x_det, _U(layer)); li; ++li) {
    xml_comp_t x_layer = li;
    // xml_comp_t  envelope      = x_det.child(_Unicode(envelope), false);

    xml_comp_t x_endcap = x_layer.child(_Unicode(endcap_envelope), false);
    // xml_comp_t x_layout = x_layer.child(_Unicode(x_layout), false); // was rphi
    // xml_comp_t y_layout = x_layer.child(_Unicode(y_layout), false); // was z
    xml_comp_t z_layout = x_layer.child(_Unicode(z_layout), false); // was z
    int        lay_id   = x_layer.id();
    string     m_nam    = x_layer.moduleStr();
    string     lay_nam  = det_name + _toString(x_layer.id(), "_layer%d");
    Tube       lay_tub(x_endcap.inner_r(), x_endcap.outer_r(), x_endcap.z_length() / 2.0);
    Volume     lay_vol(lay_nam, lay_tub, air); // Create the layer envelope volume.
    // Position   lay_pos(0, 0, z_layout.z0());
    zPos = z_layout.z0();
    Position   lay_pos(0, 0, 0);
    lay_vol.setVisAttributes(description.visAttributes(x_layer.visStr()));

    Volume      module_env = volumes[m_nam];
    DetElement  lay_elt(sdet, lay_nam, lay_id);
    Placements& sensVols = sensitives[m_nam];

    // the local coordinate systems of modules in dd4hep and acts differ
    // see http://acts.web.cern.ch/ACTS/latest/doc/group__DD4hepPlugins.html
    auto& layerParams = DD4hepDetectorHelper::ensureExtension<dd4hep::rec::VariantParameters>(lay_elt);

    // layerParams.set<double>("envelope_r_min", x_endcap.inner_r()/dd4hep::cm);
    // layerParams.set<double>("envelope_r_max", x_endcap.outer_r()/dd4hep::cm);
    // layerParams.set<double>("envelope_z_min", (z_layout.z0()-x_endcap.z_length() / 2.0)/dd4hep::cm);
    // layerParams.set<double>("envelope_z_max", (z_layout.z0()+x_endcap.z_length() / 2.0)/dd4hep::cm);
    for (xml_coll_t lmat(x_layer, _Unicode(layer_material)); lmat; ++lmat) {
      xml_comp_t x_layer_material = lmat;
      DD4hepDetectorHelper::xmlToProtoSurfaceMaterial(x_layer_material, layerParams, "layer_material");
    }

    int    module   = 1;

    int ny = 16;
    double y = -35*cm;
    double x = 0;
    double module_z = 0;
    double dy = 5.5*cm;
    for (int ii = 0; ii < ny; ii++) {
      if(abs(y)<15*cm){
        y+=dy;
        continue;
      }
      string     module_name = _toString(module, "module%d");
      DetElement mod_elt(lay_elt, module_name, module);

      Transform3D tr(RotationZYX(M_PI / 2, 0, 0), Position(x, y, module_z));

      pv = lay_vol.placeVolume(module_env, tr);
      pv.addPhysVolID("module", module);
      mod_elt.setPlacement(pv);
      for (size_t ic = 0; ic < sensVols.size(); ++ic) {
        PlacedVolume sens_pv = sensVols[ic];
        DetElement   comp_de(mod_elt, std::string("de_") + sens_pv.volume().name(), module);
        comp_de.setPlacement(sens_pv);

        auto& comp_de_params = DD4hepDetectorHelper::ensureExtension<dd4hep::rec::VariantParameters>(comp_de);
        comp_de_params.set<string>("axis_definitions", "XYZ");
        // comp_de.setAttributes(description, sens_pv.volume(), x_layer.regionStr(), x_layer.limitsStr(),
        //                       xml_det_t(xmleles[m_nam]).visStr());
        //

        volSurfaceList(comp_de)->push_back(volplane_surfaces[m_nam][ic]);
      }

      /// Increase counters etc.
      module++;
      // Adjust the x and y coordinates of the module.
      y += dy;
      // Flip sign of x and y adjustments.
      // dx *= -1;
      // dy *= -1;
      // Add z increment to get next z placement pos.
    }
    // Create the PhysicalVolume for the layer.
    pv = assembly.placeVolume(lay_vol, lay_pos); // Place layer in mother
    pv.addPhysVolID("layer", lay_id);            // Set the layer ID.
    lay_elt.setAttributes(description, lay_vol, x_layer.regionStr(), x_layer.limitsStr(), x_layer.visStr());
    lay_elt.setPlacement(pv);
  }
  // assembly->GetShape()->ComputeBBox();

  // pv = description.pickMotherVolume(sdet).placeVolume(assembly, Position(0, 0, disk_zPos));
  pv = description.pickMotherVolume(sdet).placeVolume(assembly, Position(0, 0, zPos));
  pv.addPhysVolID("system", det_id);
  sdet.setPlacement(pv);

  return sdet;
}

//@}
// clang-format off
DECLARE_DETELEMENT(epic_TOFEndcap, create_detector)
