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
  map<string, Volume> modules_overlap; //M.S.
  map<string, Placements> sensitives;
  map<string, Placements> sensitives_overlap; //M.S.
  map<string, std::vector<VolPlane>> volplane_surfaces;
  map<string, std::array<double, 2>> module_thicknesses;
  PlacedVolume pv; 
  PlacedVolume pv_overlap; //M.S. pv_overlap to automatize some parts

  double xb = 0, zbox = 0, Zbox = 0; // M.S. variable to store width of overlap

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
        Transform3D(Rotation3D(),
                    Position(0, 0, (reflect ? -1.0 : 1.0) * (support_zstart + support_length / 2)));
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

  for (xml_coll_t mi(x_det, _U(module)); mi; ++mi, ++m_id) {
    xml_comp_t x_mod = mi;
    string m_nam     = x_mod.nameStr();
    string m_nam_overlap = m_nam + "_overlap"; //M.S.
    xml_comp_t trd   = x_mod.trd();
    xml_comp_t box = x_mod.child(_U(box)); //M.S.

    double posY;
    double x1              = trd.x1();
    double x2              = trd.x2();
    double z               = trd.z();
    double total_thickness = 0.;
    xml_coll_t ci(x_mod, _U(module_component));
    for (ci.reset(), total_thickness = 0.0; ci; ++ci)
      total_thickness += xml_comp_t(ci).thickness();

    double thickness_so_far = 0.0;
    double y1               = total_thickness / 2;
    double y2               = total_thickness / 2; 
    Trapezoid m_solid(x1, x2, y1, y2, z);
    Volume m_volume(m_nam, m_solid, vacuum);
    m_volume.setVisAttributes(description.visAttributes(x_mod.visStr()));

    //M.S. define the two boxes for overlap
    xb = box.x();
    zbox = box.z(); 
    Zbox = getAttrOrDefault(box, _Unicode(Z), zbox);
    Box m_box_overlap(xb, y1, Zbox);

    Volume m_volume_overlap(m_nam_overlap, m_box_overlap, vacuum); //M.S.
    m_volume_overlap.setVisAttributes(description.visAttributes(x_mod.visStr())); //M.S.
    //end of M.S. 

    Solid frame_s;
    if (x_mod.hasChild(_U(frame))) {
      // build frame from trd (assumed to be smaller)
      xml_comp_t m_frame     = x_mod.child(_U(frame));
      xml_comp_t f_pos       = m_frame.child(_U(position));
      xml_comp_t frame_trd   = m_frame.trd();
      double frame_thickness = getAttrOrDefault(m_frame, _U(thickness), total_thickness);
      double frame_x1        = frame_trd.x1();
      double frame_x2        = frame_trd.x2();
      double frame_z         = frame_trd.z();
      // make the frame match the total thickness if thickness attribute is not given
      Trapezoid f_solid1(x1, x2, frame_thickness / 2.0, frame_thickness / 2.0, z);
      Trapezoid f_solid(frame_x1, frame_x2, frame_thickness / 2.0, frame_thickness / 2.0, frame_z);
      SubtractionSolid frame_shape(f_solid1, f_solid);
      frame_s = frame_shape;

      Material f_mat = description.material(m_frame.materialStr());
      Volume f_vol(m_nam + "_frame", frame_shape, f_mat);
      f_vol.setVisAttributes(description.visAttributes(m_frame.visStr()));

      // figure out how to best place
      pv = m_volume.placeVolume(f_vol, Position(f_pos.x(), f_pos.y(), f_pos.z()));
    }

    for (ci.reset(), n_sensor = 1, c_id = 0, posY = -y1; ci; ++ci, ++c_id) {
      xml_comp_t c     = ci;
      double c_thick   = c.thickness();
      auto comp_x1     = getAttrOrDefault(c, _Unicode(x1), x1);
      auto comp_x2     = getAttrOrDefault(c, _Unicode(x2), x2);
      auto comp_xb     = getAttrOrDefault(c, _Unicode(xb), xb); // M.S.
      auto comp_height = getAttrOrDefault(c, _Unicode(height), z);
      //auto comp_zb     = getAttrOrDefault(c, _Unicode(z), Zbox); // M.S.
      auto comp_zb = Zbox; 
      if(comp_height < z) comp_zb = zbox; 

      std::cout << Zbox << " " << zbox << " " << comp_zb << std::endl; // M.S. printout to check the correct assignment of comp_zb

      Material c_mat = description.material(c.materialStr());
      string c_name  = _toString(c_id, "component%d");
      string c_name_overlap = c_name + "_overlap"; //M.S.

      Trapezoid comp_s1(comp_x1, comp_x2, c_thick / 2e0, c_thick / 2e0, comp_height);
      Box comp_s_overlap(comp_xb, c_thick / 2e0, comp_zb); //M.S.
      Solid comp_shape = comp_s1;
      if (frame_s.isValid()) {
        comp_shape = SubtractionSolid(comp_s1, frame_s);
      }
      Volume c_vol(c_name, comp_shape, c_mat);
      Volume c_vol_overlap(c_name_overlap, comp_s_overlap, c_mat); //M.S.

      c_vol.setVisAttributes(description.visAttributes(c.visStr()));
      c_vol_overlap.setVisAttributes(description.visAttributes(c.visStr())); //M.S.
      pv = m_volume.placeVolume(c_vol, Position(0, posY + c_thick / 2, comp_height-z));
      pv_overlap = m_volume_overlap.placeVolume(c_vol_overlap, Position(0, posY + c_thick / 2, comp_zb-Zbox)); //M.S.
      cout << zbox << " " << Zbox << " " << comp_zb << " " << comp_zb-zbox << endl; 
      if (c.isSensitive()) {
        module_thicknesses[m_nam] = {thickness_so_far + c_thick / 2.0,
                                     total_thickness - thickness_so_far - c_thick / 2.0};
        // std::cout << " adding sensitive volume" << c_name << "\n";
        sdet.check(n_sensor > 2,
                   "SiTrackerEndcap2::fromCompact: " + c_name + " Max of 2 modules allowed!");
        pv.addPhysVolID("sensor", n_sensor);
        pv_overlap.addPhysVolID("sensor", n_sensor); //M.S.
        c_vol.setSensitiveDetector(sens);
        c_vol_overlap.setSensitiveDetector(sens); //M.S.
        sensitives[m_nam].push_back(pv);
        sensitives_overlap[m_nam].push_back(pv_overlap); //M.S.
        ++n_sensor;
        // -------- create a measurement plane for the tracking surface attched to the sensitive volume -----
        Vector3D u(0., 0., -1.);
        Vector3D v(-1., 0., 0.);
        Vector3D n(0., 1., 0.);
        // Vector3D o( 0. , 0. , 0. ) ;

        // compute the inner and outer thicknesses that need to be assigned to the tracking surface
        // depending on wether the support is above or below the sensor
        double inner_thickness = module_thicknesses[m_nam][0];
        double outer_thickness = module_thicknesses[m_nam][1];

        SurfaceType type(SurfaceType::Sensitive);

        // if( isStripDetector )
        //  type.setProperty( SurfaceType::Measurement1D , true ) ;

        VolPlane surf(c_vol, type, inner_thickness, outer_thickness, u, v, n); //,o ) ;
        VolPlane surf_overlap(c_vol_overlap, type, inner_thickness, outer_thickness, u, v, n); //M.S.
        volplane_surfaces[m_nam].push_back(surf);
        volplane_surfaces[m_nam].push_back(surf_overlap); //M.S.
        //--------------------------------------------
      }
      posY += c_thick;
      thickness_so_far += c_thick;
    }
    modules[m_nam] = m_volume;
    modules_overlap[m_nam] = m_volume_overlap; //M.S.
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
    //std::cout << layer_length << std::endl;
    double layer_zstart   = l_env.attr<double>(_Unicode(zstart));
    double layer_center_z = layer_zstart + layer_length / 2.0;
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
      double r             = x_ring.r();
      double R            = x_ring.R(); //M.S. 
      double phi0          = x_ring.phi0(0);
      double zstart        = x_ring.zstart();
      double dz            = x_ring.dz(0);
      double dz_offset     = getAttrOrDefault(x_ring, _Unicode(dz_offset), 0.0); //M.S.
      int nmodules         = x_ring.nmodules();
      int nmodules_per_quadrant = (int)(nmodules / 4);
      string m_nam         = x_ring.moduleStr();
      Volume m_vol         = modules[m_nam];
      Volume m_vol_overlap = modules_overlap[m_nam]; //M.S.
      double iphi          = 2 * M_PI / nmodules;
      double phi           = phi0;
      Placements& sensVols = sensitives[m_nam];
      Placements& sensVols_overlap = sensitives_overlap[m_nam];

      for (int k = 0; k < nmodules; ++k) {
        string m_base = _toString(l_id, "layer%d") + _toString(mod_num, "_module%d");
        double x      = -r * std::cos(phi);
        double y      = -r * std::sin(phi);

        // M.S. stagger quadrants 
        double dz_final = dz;
        if(k < nmodules_per_quadrant  || (k >= 2*nmodules_per_quadrant && k < 3*nmodules_per_quadrant)) dz_final = dz + dz_offset;

        if (!reflect) {
          DetElement module(layer_element, m_base + "_pos", det_id);
          pv = layer_vol.placeVolume(m_vol, Transform3D(RotationZYX(0, -M_PI / 2 - phi, -M_PI / 2),
                                                        Position(x, y, zstart + dz_final)));
          pv.addPhysVolID("module", mod_num);
          module.setPlacement(pv);
          for (size_t ic = 0; ic < sensVols.size(); ++ic) {
            PlacedVolume sens_pv = sensVols[ic];
            DetElement comp_elt(module, sens_pv.volume().name(), mod_num);
            auto& comp_elt_params =
                DD4hepDetectorHelper::ensureExtension<dd4hep::rec::VariantParameters>(comp_elt);
            comp_elt_params.set<string>("axis_definitions", "XZY");
            comp_elt.setPlacement(sens_pv);
            volSurfaceList(comp_elt)->push_back(volplane_surfaces[m_nam][ic]);
          }
        } else {
          pv = layer_vol.placeVolume(m_vol, Transform3D(RotationZYX(0, -M_PI / 2 - phi, -M_PI / 2),
                                                        Position(x, y, -zstart - dz_final)));
          pv.addPhysVolID("module", mod_num);
          DetElement r_module(layer_element, m_base + "_neg", det_id);
          r_module.setPlacement(pv);
          for (size_t ic = 0; ic < sensVols.size(); ++ic) {
            PlacedVolume sens_pv = sensVols[ic];
            DetElement comp_elt(r_module, sens_pv.volume().name(), mod_num);
            auto& comp_elt_params =
                DD4hepDetectorHelper::ensureExtension<dd4hep::rec::VariantParameters>(comp_elt);
            comp_elt_params.set<string>("axis_definitions", "XZY");
            comp_elt.setPlacement(sens_pv);
            volSurfaceList(comp_elt)->push_back(volplane_surfaces[m_nam][ic]);
          }
        }
        dz = -dz;
        phi += iphi;
        ++mod_num;

        //M.S. add two overlaps volumes for each quadrant 
        if(k%nmodules_per_quadrant == 0 || k%(nmodules_per_quadrant) == nmodules_per_quadrant-1){
        //if(-1 == 1){
        double xo = 0, yo = 0, roto = 0; 

        if(k == 0)                              xo = -R,  yo = xb,  roto = -M_PI/2;
        else if (k == nmodules_per_quadrant-1)  xo = xb,  yo = -R,  roto = M_PI;
        else if(k == nmodules_per_quadrant)     xo = -xb,  yo = -R,  roto = M_PI;
        else if(k == 2*nmodules_per_quadrant-1) xo = R,   yo = xb, roto = M_PI/2;
        else if(k == 2*nmodules_per_quadrant)   xo = R,   yo = -xb,  roto = M_PI/2;
        else if(k == 3*nmodules_per_quadrant-1) xo = -xb,  yo = R,   roto = 0;
        else if(k == 3*nmodules_per_quadrant)   xo = xb, yo = R,   roto = 0;
        else if(k == 4*nmodules_per_quadrant-1) xo = -R,  yo = -xb, roto = -M_PI/2;

        if(!reflect){
          DetElement module(layer_element, m_base + "_overlap_pos", det_id);
          pv = layer_vol.placeVolume(m_vol_overlap, Transform3D(RotationZYX(0, roto, -M_PI / 2),
                                                        Position(xo, yo, zstart + dz_final)));
          pv.addPhysVolID("module", mod_num);
          module.setPlacement(pv);
          for (size_t ic = 0; ic < sensVols_overlap.size(); ++ic) {
            PlacedVolume sens_pv = sensVols_overlap[ic];
            DetElement comp_elt(module, sens_pv.volume().name(), mod_num);
            auto& comp_elt_params =
                DD4hepDetectorHelper::ensureExtension<dd4hep::rec::VariantParameters>(comp_elt);
            comp_elt_params.set<string>("axis_definitions", "XZY");
            comp_elt.setPlacement(sens_pv);
            volSurfaceList(comp_elt)->push_back(volplane_surfaces[m_nam][ic]);
          } 
        }else{
          DetElement r_module(layer_element, m_base + "_overlap_neg", det_id);
          pv = layer_vol.placeVolume(m_vol_overlap, Transform3D(RotationZYX(0, roto, -M_PI / 2),
                                                        Position(xo, yo, -zstart - dz_final)));
          pv.addPhysVolID("module", mod_num);
          r_module.setPlacement(pv);
          for (size_t ic = 0; ic < sensVols_overlap.size(); ++ic) {
            PlacedVolume sens_pv = sensVols_overlap[ic];
            DetElement comp_elt(r_module, sens_pv.volume().name(), mod_num);
            auto& comp_elt_params =
                DD4hepDetectorHelper::ensureExtension<dd4hep::rec::VariantParameters>(comp_elt);
            comp_elt_params.set<string>("axis_definitions", "XZY");
            comp_elt.setPlacement(sens_pv);
            volSurfaceList(comp_elt)->push_back(volplane_surfaces[m_nam][ic]);
        }


        ++mod_num;
      }
    }
      }
    }
  }
  pv = motherVol.placeVolume(assembly, Position(0, 0, (reflect ? -1.0e-9 : 1.0e-9)));
  pv.addPhysVolID("system", det_id);
  sdet.setPlacement(pv);
  return sdet;
}

DECLARE_DETELEMENT(epic_MPGDTrapEndcapTracker, create_detector)
