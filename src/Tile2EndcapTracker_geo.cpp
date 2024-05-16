// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Whitney Armstrong

/** \addtogroup Trackers Trackers
 * \brief Type: **BarrelTrackerWithFrame**.
 * \author W. Armstrong
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
#include "DD4hepDetectorHelper.h"
#include <array>
#include <map>

using namespace std;
using namespace dd4hep;
using namespace dd4hep::rec;
using namespace dd4hep::detail;

/** Endcap Trapezoidal Tracker.
 *
 * @author Whitney Armstrong
 *
 */
static Ref_t create_detector(Detector& description, xml_h e, SensitiveDetector sens) {
  typedef vector<PlacedVolume> Placements;
  xml_det_t x_det = e;
  Material vacuum = description.vacuum();
  int det_id      = x_det.id();
  string det_name = x_det.nameStr();
  bool reflect    = x_det.reflect(false);
  DetElement sdet(det_name, det_id);
  Assembly assembly(det_name);

  Material air     = description.material("Air");
  Volume motherVol = description.pickMotherVolume(sdet);
  int m_id = 0, c_id = 0, n_sensor = 0;
  map<string, Volume> modules;
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

  assembly.setVisAttributes(description.invisible());
  sens.setType("tracker");

  for (xml_coll_t mi(x_det, _U(module)); mi; ++mi, ++m_id) {
    xml_comp_t x_mod = mi;
    string m_nam     = x_mod.nameStr();
    xml_comp_t trd   = x_mod.trd();

    double posY;
    double y              = trd.y();
    double x              = trd.x();
    double m_dz              = trd.dz();
    // double z               = trd.z();
    double total_thickness = 0.;
    xml_coll_t ci(x_mod, _U(module_component));
    for (ci.reset(), total_thickness = 0.0; ci; ++ci)
      total_thickness += xml_comp_t(ci).thickness();

    double thickness_so_far = 0.0;
    double y1               = total_thickness / 2e0;

    std::cout<<"module dimension xyz: "<< x/2.0<<"  "<<y/2.0<<"  "<<m_dz/2.0<<endl;
    Box m_solid(x/2e0, y/2e0, m_dz/2e0);
    Volume m_volume(m_nam, m_solid, vacuum);

    // Assembly m_volume(m_nam);
    m_volume.setVisAttributes(description.visAttributes(x_mod.visStr()));


    for (ci.reset(), n_sensor = 1, c_id = 0, posY = -y1; ci; ++ci, ++c_id) {
      xml_comp_t c     = ci;
      double c_thick   = c.thickness();
      auto comp_x     = getAttrOrDefault(c, _Unicode(x), x);
      auto comp_y     = getAttrOrDefault(c, _Unicode(y), y);
      // auto comp_height = getAttrOrDefault(c, _Unicode(height), z);

      Material c_mat = description.material(c.materialStr());
      string c_name  = _toString(c_id, "component%d");

      Box comp_s1(comp_x/2.0, comp_y/2.0, c_thick / 2e0);
      Solid comp_shape = comp_s1;

      Volume c_vol(c_name, comp_shape, c_mat);

      c_vol.setVisAttributes(description.visAttributes(c.visStr()));
      pv = m_volume.placeVolume(c_vol, Position(0, 0, posY + c_thick / 2));
      if (c.isSensitive()) {
        module_thicknesses[m_nam] = {thickness_so_far + c_thick / 2.0,
                                     total_thickness - thickness_so_far - c_thick / 2.0};
        // std::cout << " adding sensitive volume" << c_name << "\n";
        sdet.check(n_sensor > 2,
                   "SiTrackerEndcap2::fromCompact: " + c_name + " Max of 2 modules allowed!");
        pv.addPhysVolID("sensor", n_sensor);
        c_vol.setSensitiveDetector(sens);
        sensitives[m_nam].push_back(pv);
        ++n_sensor;
        // -------- create a measurement plane for the tracking surface attched to the sensitive volume -----
        Vector3D u(1., 0., 0.);
        Vector3D v(0. , 1., 0.);
        Vector3D n(0., 0., 1.);
        Vector3D o( 0. , 0. , 0. ) ;

        // compute the inner and outer thicknesses that need to be assigned to the tracking surface
        // depending on wether the support is above or below the sensor
        double inner_thickness = module_thicknesses[m_nam][0];
        double outer_thickness = module_thicknesses[m_nam][1];

        SurfaceType type(SurfaceType::Sensitive);

        // if( isStripDetector )
        //  type.setProperty( SurfaceType::Measurement1D , true ) ;

        VolPlane surf(c_vol, type, inner_thickness, outer_thickness, u, v, n ,o ) ;
        volplane_surfaces[m_nam].push_back(surf);

        //--------------------------------------------
      }
      posY += c_thick;
      thickness_so_far += c_thick;
    }
    modules[m_nam] = m_volume;
  }

  for (xml_coll_t li(x_det, _U(layer)); li; ++li) {
    xml_comp_t x_layer(li);
    int l_id    = x_layer.id();
    int mod_num = 1;

    xml_comp_t l_env  = x_layer.child(_U(envelope));
    string layer_name = det_name + std::string("_layer") + std::to_string(l_id);

    std::string layer_vis = l_env.attr<std::string>(_Unicode(vis));
    double layer_rmin     = l_env.attr<double>(_Unicode(rmin));
    double layer_rmax     = l_env.attr<double>(_Unicode(rmax));
    double layer_length   = l_env.attr<double>(_Unicode(length));
    double layer_zstart   = l_env.attr<double>(_Unicode(zstart));
    double layer_center_z = layer_zstart + layer_length / 2.0;
    std::cout<<"Tile2 layer rmin/max:"<<layer_rmin<<"  "<<layer_rmax<<std::endl;

    // printout(INFO,"ROOTGDMLParse","+++ Read geometry from GDML file file:%s",input.c_str());
    // std::cout << "SiTracker Endcap layer " << l_id << " zstart = " << layer_zstart/dd4hep::mm << "mm ( " <<
    // layer_length/dd4hep::mm << " mm thick )\n";

    // Assembly    layer_assembly(layer_name);
    // assembly.placeVolume(layer_assembly);
    Tube layer_tub(layer_rmin, layer_rmax, layer_length / 2);
    Volume layer_vol(layer_name, layer_tub, air); // Create the layer envelope volume.
    layer_vol.setVisAttributes(description.visAttributes(layer_vis));

    PlacedVolume layer_pv;
    if (reflect) {
      layer_pv = assembly.placeVolume(
          layer_vol, Transform3D(RotationZYX(0.0, -M_PI, 0.0), Position(0, 0, -layer_center_z)));
      layer_pv.addPhysVolID("layer", l_id);
      layer_name += "_N";
    } else {
      layer_pv = assembly.placeVolume(layer_vol, Position(0, 0, layer_center_z));
      layer_pv.addPhysVolID("layer", l_id);
      layer_name += "_P";
    }
    DetElement layer_element(sdet, layer_name, l_id);
    layer_element.setPlacement(layer_pv);

    auto& layerParams =
        DD4hepDetectorHelper::ensureExtension<dd4hep::rec::VariantParameters>(layer_element);
    
    for (xml_coll_t lmat(x_layer, _Unicode(layer_material)); lmat; ++lmat) {
      xml_comp_t x_layer_material = lmat;
      DD4hepDetectorHelper::xmlToProtoSurfaceMaterial(x_layer_material, layerParams,
                                                      "layer_material");
    }

    for (xml_coll_t ri(x_layer, _U(ring)); ri; ++ri) {
      xml_comp_t x_ring    = ri;
      // double r             = x_ring.r();
      double x             = x_ring.x();
      double y             = x_ring.y();
      // double phi0          = x_ring.phi0(0);
      double zstart        = x_ring.zstart();
      double dz            = x_ring.dz(0);
      int nmodules         = x_ring.nmodules();
      string m_nam         = x_ring.moduleStr();
      Volume m_vol         = modules[m_nam];
      // double iphi          = 2 * M_PI / nmodules;
      // double phi           = phi0;
      Placements& sensVols = sensitives[m_nam];

      for (int k = 0; k < nmodules; ++k) {
      // for (int k = 0; k < 2; ++k) {

        string m_base = _toString(l_id, "layer%d") + _toString(mod_num, "_module%d");
        string pos;
        pos = (reflect) ? "_neg" : "_pos";
        double flip = (reflect) ? -1.0:1.0;
        DetElement r_module(layer_element, m_base + pos, det_id);
        std::cout<<"layer_vol xyz center:"<<x<<"  "<<y<<"  "<<flip*(zstart + dz)<<std::endl;
        pv = layer_vol.placeVolume(m_vol, Transform3D(RotationZYX(0, 0, 0),
                                                Position(x, y, flip*(zstart + dz))));
        pv.addPhysVolID("module", mod_num);
        r_module.setPlacement(pv);
        for (size_t ic = 0; ic < sensVols.size(); ++ic) {
          PlacedVolume sens_pv = sensVols[ic];
          DetElement comp_elt(r_module, sens_pv.volume().name(), mod_num);
          auto& comp_elt_params =
          DD4hepDetectorHelper::ensureExtension<dd4hep::rec::VariantParameters>(comp_elt);
          comp_elt_params.set<string>("axis_definitions", "XYZ");
          comp_elt.setPlacement(sens_pv);
          volSurfaceList(comp_elt)->push_back(volplane_surfaces[m_nam][ic]);
        }
      mod_num++;
      y = -(y >= 0 ? 1 : -1) * (std::abs(y) + 0.5);

      }

        // if (!reflect) {
        //   DetElement r_module(layer_element, m_base + "_pos", det_id);
        //   // pv = layer_vol.placeVolume(m_vol, Transform3D(RotationZYX(0, -M_PI / 2 - phi, -M_PI / 2),
        //   pv = layer_vol.placeVolume(m_vol, Transform3D(RotationZYX(0, 0, 0),
        //                                                 Position(x, y, zstart + dz)));
        //   pv.addPhysVolID("module", mod_num);
        //   r_module.setPlacement(pv);
        //   for (size_t ic = 0; ic < sensVols.size(); ++ic) {
        //     PlacedVolume sens_pv = sensVols[ic];
        //     DetElement comp_elt(r_module, sens_pv.volume().name(), mod_num);
        //     auto& comp_elt_params =
        //         DD4hepDetectorHelper::ensureExtension<dd4hep::rec::VariantParameters>(comp_elt);
        //     comp_elt_params.set<string>("axis_definitions", "XYZ");
        //     comp_elt.setPlacement(sens_pv);
        //     volSurfaceList(comp_elt)->push_back(volplane_surfaces[m_nam][ic]);
        //   }
        // } else {
        //   // pv = layer_vol.placeVolume(m_vol, Transform3D(RotationZYX(0, -M_PI / 2 - phi, -M_PI / 2),
        //   pv = layer_vol.placeVolume(m_vol, Transform3D(RotationZYX(0, 0, 0),
        //                                                 Position(x, y, -zstart - dz)));
        //   pv.addPhysVolID("module", mod_num);
        //   DetElement r_module(layer_element, m_base + "_neg", det_id);
        //   r_module.setPlacement(pv);
        //   for (size_t ic = 0; ic < sensVols.size(); ++ic) {
        //     PlacedVolume sens_pv = sensVols[ic];
        //     DetElement comp_elt(r_module, sens_pv.volume().name(), mod_num);
        //     auto& comp_elt_params =
        //         DD4hepDetectorHelper::ensureExtension<dd4hep::rec::VariantParameters>(comp_elt);
        //     comp_elt_params.set<string>("axis_definitions", "XYZ");
        //     comp_elt.setPlacement(sens_pv);
        //     volSurfaceList(comp_elt)->push_back(volplane_surfaces[m_nam][ic]);
        //   }
        // }
        // dz = -dz;
        // phi += iphi;
        // ++mod_num;
      // }
    }
  }
  pv = motherVol.placeVolume(assembly, Position(0, 0, (reflect ? -1.0e-9 : 1.0e-9)));
  pv.addPhysVolID("system", det_id);
  sdet.setPlacement(pv);
  return sdet;
}

DECLARE_DETELEMENT(epic_Tile2EndcapTracker, create_detector)
