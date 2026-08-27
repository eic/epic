// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Marco Spreafico

/** \addtogroup Trackers Trackers
 * \brief Type: **MPGDEndcap**.
 * \author M. Spreafico
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
  map<string, Volume> modules_N, modules_S;
  map<string, Volume> modules_overlap_N, modules_overlap_S; //M.S.
  map<string, Placements> sensitives_N, sensitives_S;
  map<string, Placements> sensitives_overlap_N, sensitives_overlap_S; //M.S.
  map<string, std::vector<VolPlane>> volplane_surfaces_N, volplane_surfaces_S;
  map<string, std::array<double, 2>> module_thicknesses;
  PlacedVolume pv_N, pv_S, pv; 
  PlacedVolume pv_overlap_N, pv_overlap_S; //M.S. pv_overlap to automatize some parts

  //Parameters to store 
  double xb = 0;; // overlap size 
  std::array<double, 3> IRN = {0.0, 0.0, 0.0};
  std::array<double, 3> IRS = {0.0, 0.0, 0.0};
  std::array<double, 4> OR = {0.0, 0.0, 0.0, 0.0};
  double theta_trd = 0.0; 

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

  //Removed because no support 
  /*for (xml_coll_t su(x_det, _U(support)); su; ++su) {
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
  }*/

  for (xml_coll_t mi(x_det, _U(module)); mi; ++mi, ++m_id) {
    xml_comp_t x_mod = mi;
    string m_nam     = x_mod.nameStr();

    double posY;

    // get and store parameters 
    xml_comp_t trd   = x_mod.trd(); 
    xml_comp_t box = x_mod.child(_U(box)); 

    IRN = {
        getAttrOrDefault<double>(trd, _Unicode(IR1N), 0.0),
        getAttrOrDefault<double>(trd, _Unicode(IR2N), 0.0),
        getAttrOrDefault<double>(trd, _Unicode(IR3N), 0.0)};
    IRS = {
        getAttrOrDefault<double>(trd, _Unicode(IR1S), 0.0),
        getAttrOrDefault<double>(trd, _Unicode(IR2S), 0.0),
        getAttrOrDefault<double>(trd, _Unicode(IR3S), 0.0)};
    OR = {
        getAttrOrDefault<double>(trd, _Unicode(OR1), 0.0),
        getAttrOrDefault<double>(trd, _Unicode(OR2), 0.0), 
        getAttrOrDefault<double>(trd, _Unicode(OR3), 0.0),
        getAttrOrDefault<double>(trd, _Unicode(OR4), 0.0)};

    theta_trd = getAttrOrDefault<double>(trd, _Unicode(theta), 0.0);

    xb = box.x();
    // end of getting and storing parameters 
    
    double total_thickness = 0.;
    xml_coll_t ci(x_mod, _U(module_component));
    for (ci.reset(), total_thickness = 0.0; ci; ++ci)
      total_thickness += xml_comp_t(ci).thickness();

    double thickness_so_far = 0.0;

    double y = total_thickness / 2;

    double xN = IRN[0] * sin(theta_trd);
    double XN = OR[0] * sin(theta_trd);
    double zN = (OR[0] - IRN[0]) * cos(theta_trd) / 2.0;

    double xS = IRS[0] * sin(theta_trd);
    double XS = OR[0] * sin(theta_trd);
    double zS = (OR[0] - IRS[0]) * cos(theta_trd) / 2.0;

    // Define volumes for north and south trapezoids
    Trapezoid m_solid_N(xN, XN, y, y, zN);
    Trapezoid m_solid_S(xS, XS, y, y, zS);

    Volume m_volume_N(m_nam + "_N", m_solid_N, vacuum);
    Volume m_volume_S(m_nam + "_S", m_solid_S, vacuum);

    m_volume_N.setVisAttributes(description.visAttributes(x_mod.visStr()));
    m_volume_S.setVisAttributes(description.visAttributes(x_mod.visStr()));

    double zbox_N = (OR[0] - IRN[0])/2;
    double zbox_S = (OR[0] - IRS[0])/2;
    // Define volume for overlap box
    Box m_box_overlap_N(xb, y, zbox_N);
    Box m_box_overlap_S(xb, y, zbox_S);

    Volume m_volume_overlap_N(m_nam + "_overlap_N", m_box_overlap_N, vacuum);
    Volume m_volume_overlap_S(m_nam + "_overlap_S", m_box_overlap_S, vacuum);

    m_volume_overlap_N.setVisAttributes(description.visAttributes(x_mod.visStr()));
    m_volume_overlap_S.setVisAttributes(description.visAttributes(x_mod.visStr()));
    
    //Commented frame part because frames are not present
    /*
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
    */

    for (ci.reset(), n_sensor = 1, c_id = 0, posY = -y; ci; ++ci, ++c_id) {
      xml_comp_t c     = ci;
      double c_thick   = c.thickness();
      auto comp_IRN = IRN[getAttrOrDefault(c, _Unicode(IR), 1)-1];
      auto comp_IRS = IRS[getAttrOrDefault(c, _Unicode(IR), 1)-1];
      auto comp_OR = OR[getAttrOrDefault(c, _Unicode(OR), 1)-1];

      double comp_x1_N = comp_IRN * sin(theta_trd);
      double comp_x2_N = comp_OR * sin(theta_trd);
      double comp_z_N = (comp_OR - comp_IRN) * cos(theta_trd)/2;

      double comp_x1_S = comp_IRS * sin(theta_trd);
      double comp_x2_S = comp_OR * sin(theta_trd);
      double comp_z_S = (comp_OR - comp_IRS) * cos(theta_trd)/2;

      double comp_zb_N = (comp_OR - comp_IRN)/2.0; 
      double comp_zb_S = (comp_OR - comp_IRS)/2.0;

      Material c_mat = description.material(c.materialStr());

      string c_name  = _toString(c_id, "component%d");
      string c_name_overlap = c_name + "_overlap"; //M.S.

      //Define volumes for north and south sectors 
      Trapezoid comp_s_N(comp_x1_N, comp_x2_N, c_thick / 2e0, c_thick / 2e0, comp_z_N);
      Trapezoid comp_s_S(comp_x1_S, comp_x2_S, c_thick / 2e0, c_thick / 2e0, comp_z_S);

      Box comp_o_N(xb, c_thick / 2e0, comp_zb_N); 
      Box comp_o_S(xb, c_thick / 2e0, comp_zb_S); 

      Solid comp_shape_N = comp_s_N;
      Solid comp_shape_S = comp_s_S;
      Solid comp_overlap_N = comp_o_N;
      Solid comp_overlap_S = comp_o_S;
      
      /*f (frame_s.isValid()) {
        comp_shape_N = SubtractionSolid(comp_s_N, frame_s);
        comp_shape_S = SubtractionSolid(comp_s_S, frame_s);
      }*/

      Volume c_vol_N(_toString(c_id, "%d_N"), comp_shape_N, c_mat);
      Volume c_vol_S(_toString(c_id, "%d_S"), comp_shape_S, c_mat);

      Volume c_vol_overlap_N(_toString(c_id, "%d_overlap_N"), comp_overlap_N, c_mat); //M.S.
      Volume c_vol_overlap_S(_toString(c_id, "%d_overlap_S"), comp_overlap_S, c_mat); //M.S.

      c_vol_N.setVisAttributes(description.visAttributes(c.visStr()));
      c_vol_S.setVisAttributes(description.visAttributes(c.visStr()));

      c_vol_overlap_N.setVisAttributes(description.visAttributes(c.visStr())); //M.S.
      c_vol_overlap_S.setVisAttributes(description.visAttributes(c.visStr())); //M.S.

      double comp_pos_Z_N = (comp_IRN - IRN[0]) * cos(theta_trd) + comp_z_N; 
      double comp_pos_Z_S = (comp_IRS - IRS[0]) * cos(theta_trd) + comp_z_S; 

      pv_N = m_volume_N.placeVolume(c_vol_N, Position(0, posY + c_thick / 2, comp_pos_Z_N-zN));
      pv_S = m_volume_S.placeVolume(c_vol_S, Position(0, posY + c_thick / 2, comp_pos_Z_S-zS));

      double comp_pos_Zb_N = (comp_IRN - IRN[0]) + comp_zb_N;
      double comp_pos_Zb_S = (comp_IRS - IRS[0]) + comp_zb_S;
      
      pv_overlap_N = m_volume_overlap_N.placeVolume(c_vol_overlap_N, Position(0, posY + c_thick / 2, comp_pos_Zb_N-zbox_N)); //M.S.
      pv_overlap_S = m_volume_overlap_S.placeVolume(c_vol_overlap_S, Position(0, posY + c_thick / 2, comp_pos_Zb_S-zbox_S)); //M.S.

      if (c.isSensitive()) {
        module_thicknesses[m_nam] = {thickness_so_far + c_thick / 2.0, total_thickness - thickness_so_far - c_thick / 2.0};
                                     
        sdet.check(n_sensor > 2, "SiTrackerEndcap2::fromCompact: " + c_name + " Max of 2 modules allowed!");

        pv_N.addPhysVolID("sensor", n_sensor);
        pv_S.addPhysVolID("sensor", n_sensor);
        pv_overlap_N.addPhysVolID("sensor", n_sensor); //M.S.
        pv_overlap_S.addPhysVolID("sensor", n_sensor); //M.S.

        c_vol_N.setSensitiveDetector(sens);
        c_vol_S.setSensitiveDetector(sens);
        c_vol_overlap_N.setSensitiveDetector(sens); //M.S.
        c_vol_overlap_S.setSensitiveDetector(sens); //M.S.

        sensitives_N[m_nam].push_back(pv_N);
        sensitives_S[m_nam].push_back(pv_S);
        sensitives_overlap_N[m_nam].push_back(pv_overlap_N); //M.S.
        sensitives_overlap_S[m_nam].push_back(pv_overlap_S); //M.S.

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

        VolPlane surf_N(c_vol_N, type, inner_thickness, outer_thickness, u, v, n); //,o ) ;
        VolPlane surf_S(c_vol_S, type, inner_thickness, outer_thickness, u, v, n); //M.S.
        VolPlane surf_overlap_N(c_vol_overlap_N, type, inner_thickness, outer_thickness, u, v, n); //M.S.
        VolPlane surf_overlap_S(c_vol_overlap_S, type, inner_thickness, outer_thickness, u, v, n); //M.S.
        volplane_surfaces_N[m_nam].push_back(surf_N);
        volplane_surfaces_S[m_nam].push_back(surf_S);
        volplane_surfaces_N[m_nam].push_back(surf_overlap_N); //M.S.
        volplane_surfaces_S[m_nam].push_back(surf_overlap_S); //M.S.
        //--------------------------------------------
      }
      posY += c_thick;
      thickness_so_far += c_thick;
    }
    modules_N[m_nam] = m_volume_N; //M.S.
    modules_S[m_nam] = m_volume_S; //M.S.
    modules_overlap_N[m_nam] = m_volume_overlap_N; //M.S.
    modules_overlap_S[m_nam] = m_volume_overlap_S; //M.S.
  }

  for (xml_coll_t li(x_det, _U(layer)); li; ++li) {
    xml_comp_t x_layer(li);
    int l_id    = x_layer.id();
    int mod_num = 1;

    xml_comp_t l_env  = x_layer.child(_U(envelope));
    string layer_name = det_name + std::string("_layer") + std::to_string(l_id);

    std::string layer_vis = l_env.attr<std::string>(_Unicode(vis));
    double layer_rmin_N = IRN[0] * cos(theta_trd); 
    double layer_rmin_S = IRS[0] * cos(theta_trd);
    double layer_rmax = sqrt(4*xb*xb+OR[0]*OR[0]);

    double layer_length   = l_env.attr<double>(_Unicode(length));
    double layer_zstart   = l_env.attr<double>(_Unicode(zstart));
    double layer_center_z = layer_zstart + layer_length / 2.0;

    // printout(INFO,"ROOTGDMLParse","+++ Read geometry from GDML file file:%s",input.c_str());
    // std::cout << "SiTracker Endcap layer " << l_id << " zstart = " << layer_zstart/dd4hep::mm << "mm ( " <<
    // layer_length/dd4hep::mm << " mm thick )\n";

    // Assembly    layer_assembly(layer_name);
    // assembly.placeVolume(layer_assembly);
    Tube layer_tub(std::min(layer_rmin_N, layer_rmin_S), layer_rmax, layer_length / 2);
    Volume layer_vol(layer_name, layer_tub, air); // Create the layer envelope volume.
    layer_vol.setVisAttributes(description.visAttributes(layer_vis));

    

    PlacedVolume layer_pv;
    if (reflect) {
      layer_pv = assembly.placeVolume(
          layer_vol, Transform3D(RotationZYX(0.0, -M_PI, 0.0), Position(0, 0, -layer_center_z)));
      layer_pv.addPhysVolID("layer", l_id);
      layer_name += "_B";
    } else {
      layer_pv = assembly.placeVolume(layer_vol, Position(0, 0, layer_center_z));
      layer_pv.addPhysVolID("layer", l_id);
      layer_name += "_F";
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
      double r             = (IRS[0]+OR[0])*cos(theta_trd)/2;
      double R             = (IRS[0]+OR[0])/2; 
      double phi0          = x_ring.phi0(0);
      double zstart        = x_ring.zstart();
      double dz            = x_ring.dz(0);
      double dz_offset     = getAttrOrDefault(x_ring, _Unicode(dz_offset), 0.0); //M.S.
      int nmodules         = x_ring.nmodules();
      int nmodules_per_quadrant = (int)(nmodules / 4);
      string m_nam         = x_ring.moduleStr();
      Volume m_vol         = modules_S[m_nam];
      Volume m_vol_overlap = modules_overlap_S[m_nam]; //M.S.
      double iphi          = 2 * M_PI / nmodules;
      double phi           = phi0;
      Placements& sensVols = sensitives_S[m_nam];
      Placements& sensVols_overlap = sensitives_overlap_S[m_nam];

      for (int k = 0; k < nmodules; ++k) {
        mod_num = k; 
        string m_base = _toString(l_id, "layer%d") + _toString(mod_num, "_module%d");
 
        bool isN = false; // default = S
        if(k >= nmodules_per_quadrant && k < 3*nmodules_per_quadrant) {isN = true;}

        // M.S. stagger quadrants 
        double dz_final = dz;
        if(k < nmodules_per_quadrant  || (k >= 2*nmodules_per_quadrant && k < 3*nmodules_per_quadrant)) dz_final = dz + dz_offset;

        if(isN) {
          m_vol = modules_N[m_nam];
          m_vol_overlap = modules_overlap_N[m_nam];
          sensVols = sensitives_N[m_nam];
          sensVols_overlap = sensitives_overlap_N[m_nam];
          r             = (IRN[0]+ OR[0]) * cos(theta_trd) /2;
          R             = (IRN[0]+OR[0] * cos(theta_trd))/2; 
        }else{
          m_vol = modules_S[m_nam];
          m_vol_overlap = modules_overlap_S[m_nam];
          sensVols = sensitives_S[m_nam];
          sensVols_overlap = sensitives_overlap_S[m_nam];
          r             = (IRS[0]+OR[0]) * cos(theta_trd) /2;
          R             = (IRS[0]+OR[0] * cos(theta_trd))/2; 
        }

        double x      = -r * std::cos(phi);
        double y      = -r * std::sin(phi);

        if (!reflect) {
          DetElement module(layer_element, m_base + "_pos", det_id);
          pv = layer_vol.placeVolume(m_vol, mod_num, Transform3D(RotationZYX(0, -M_PI / 2 - phi, -M_PI / 2),
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
            if(isN) {
            volSurfaceList(comp_elt)->push_back(volplane_surfaces_N[m_nam][ic]);
            } else {
            volSurfaceList(comp_elt)->push_back(volplane_surfaces_S[m_nam][ ic]);
            }
          }
        } else {
          pv = layer_vol.placeVolume(m_vol, mod_num, Transform3D(RotationZYX(0, -M_PI / 2 - phi, -M_PI / 2),
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
            if(isN) {
            volSurfaceList(comp_elt)->push_back(volplane_surfaces_N[m_nam][ic]);
            } else {
            volSurfaceList(comp_elt)->push_back(volplane_surfaces_S[m_nam][ ic]);
            }
          }
        }
        dz = -dz;
        phi += iphi;

        if (isN) {
          R = (IRN[0] + OR[0]) / 2.0;
        } else {
          R = (IRS[0] + OR[0]) / 2.0;
        }

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
          pv = layer_vol.placeVolume(m_vol_overlap, mod_num, Transform3D(RotationZYX(0, roto, -M_PI / 2),
                                                        Position(xo, yo, zstart + dz_final)));
          pv.addPhysVolID("module_overlap", mod_num);
          module.setPlacement(pv);
          for (size_t ic = 0; ic < sensVols_overlap.size(); ++ic) {
            PlacedVolume sens_pv = sensVols_overlap[ic];
            DetElement comp_elt(module, sens_pv.volume().name(), mod_num);
            auto& comp_elt_params =
                DD4hepDetectorHelper::ensureExtension<dd4hep::rec::VariantParameters>(comp_elt);
            comp_elt_params.set<string>("axis_definitions", "XZY");
            comp_elt.setPlacement(sens_pv);
            if(k >= nmodules_per_quadrant && k < 3*nmodules_per_quadrant) {
            volSurfaceList(comp_elt)->push_back(volplane_surfaces_N[m_nam][ic]);
            } else {
            volSurfaceList(comp_elt)->push_back(volplane_surfaces_S[m_nam][ ic]);
            }
          } 
        }else{
          DetElement r_module(layer_element, m_base + "_overlap_neg", det_id);
          pv = layer_vol.placeVolume(m_vol_overlap, mod_num, Transform3D(RotationZYX(0, roto, -M_PI / 2),
                                                        Position(xo, yo, -zstart - dz_final)));
          pv.addPhysVolID("module_overlap", mod_num);
          r_module.setPlacement(pv);
          for (size_t ic = 0; ic < sensVols_overlap.size(); ++ic) {
            PlacedVolume sens_pv = sensVols_overlap[ic];
            DetElement comp_elt(r_module, sens_pv.volume().name(), mod_num);
            auto& comp_elt_params =
                DD4hepDetectorHelper::ensureExtension<dd4hep::rec::VariantParameters>(comp_elt);
            comp_elt_params.set<string>("axis_definitions", "XZY");
            comp_elt.setPlacement(sens_pv);
            if(k >= nmodules_per_quadrant && k < 3*nmodules_per_quadrant) {
            volSurfaceList(comp_elt)->push_back(volplane_surfaces_N[m_nam][ic]);
            } else {
            volSurfaceList(comp_elt)->push_back(volplane_surfaces_S[m_nam][ ic]);
            }
        }
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
