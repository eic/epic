#include <map>
#include "DD4hep/DetFactoryHelper.h"
#include "XML/Utilities.h"
#include "XML/Layering.h"
#include "Acts/Plugins/DD4hep/ActsExtension.hpp"
#include "Acts/Definitions/Units.hpp"

using namespace std;
using namespace dd4hep;
using namespace dd4hep::detail;

/**  Forward off-momentum tracker.
 *
 *
 * \ingroup tracking
 *
 * @author Whitney Armstrong
 */
static Ref_t create_OffMomentumTracker(Detector& description, xml_h e, SensitiveDetector sens)
{
  typedef vector<PlacedVolume> Placements;
  xml_det_t                    x_det    = e;
  Material                     vacuum   = description.vacuum();
  int                          det_id   = x_det.id();
  string                       det_name = x_det.nameStr();
  DetElement                   sdet(det_name, det_id);
  Assembly                     assembly(det_name);
  xml::Component pos   = x_det.position();
  xml::Component rot   = x_det.rotation();


  Material  air  = description.material("Air");
  // Volume      assembly    (det_name,Box(10000,10000,10000),vacuum);
  Volume                  motherVol = description.pickMotherVolume(sdet);
  int                     m_id = 0, c_id = 0, n_sensor = 0;
  map<string, Volume>     modules;
  map<string, Placements> sensitives;
  PlacedVolume            pv;


  Acts::ActsExtension* detWorldExt = new Acts::ActsExtension();
  detWorldExt->addType("endcap", "detector");
  sdet.addExtension<Acts::ActsExtension>(detWorldExt);

  assembly.setVisAttributes(description.invisible());
  sens.setType("tracker");

  for (xml_coll_t su(x_det, _U(support)); su; ++su) {
    xml_comp_t x_support = su;
    double      support_thickness = getAttrOrDefault(x_support, _U(thickness), 2.0 * mm);
    double      support_length    = getAttrOrDefault(x_support, _U(length), 2.0 * mm);
    double      support_rmin      = getAttrOrDefault(x_support, _U(rmin), 2.0 * mm);
    double      support_zstart    = getAttrOrDefault(x_support, _U(zstart), 2.0 * mm);
    std::string support_name      = getAttrOrDefault<std::string>(x_support, _Unicode(name), "support_tube");
    std::string support_vis       = getAttrOrDefault<std::string>(x_support, _Unicode(vis), "AnlRed");
    xml_dim_t  pos        (x_support.child(_U(position), false));
    xml_dim_t  rot        (x_support.child(_U(rotation), false));
    Solid support_solid;
    if(x_support.hasChild("shape")){
      xml_comp_t shape(x_support.child(_U(shape)));
      string     shape_type = shape.typeStr();
      support_solid  = xml::createShape(description, shape_type, shape);
    } else {
      support_solid = Tube(support_rmin, support_rmin + support_thickness, support_length / 2);
    }
    Transform3D tr = Transform3D(Rotation3D(),Position(0,0,(support_zstart + support_length / 2)));
    if ( pos.ptr() && rot.ptr() )  {
      Rotation3D  rot3D(RotationZYX(rot.z(0),rot.y(0),rot.x(0)));
      Position    pos3D(pos.x(0),pos.y(0),pos.z(0));
      tr = Transform3D(rot3D, pos3D);
    }
    else if ( pos.ptr() )  {
      tr = Transform3D(Rotation3D(),Position(pos.x(0),pos.y(0),pos.z(0)));
    }
    else if ( rot.ptr() )  {
      Rotation3D rot3D(RotationZYX(rot.z(0),rot.y(0),rot.x(0)));
      tr = Transform3D(rot3D,Position());
    }
    Material    support_mat       = description.material(x_support.materialStr());
    Volume      support_vol(support_name, support_solid, support_mat);
    support_vol.setVisAttributes(description.visAttributes(support_vis));
    pv = assembly.placeVolume(support_vol, tr);
    // pv = assembly.placeVolume(support_vol, Position(0, 0, support_zstart + support_length / 2));
  }

  for (xml_coll_t mi(x_det, _U(module)); mi; ++mi, ++m_id) {
    xml_comp_t x_mod = mi;
    string     m_nam = x_mod.nameStr();
    xml_comp_t x_box = x_mod.shape();

    double     x1              = x_box.x();
    double     y1              = x_box.y();
    //double     z1              = x_box.z();
    double     total_thickness = 0.;

    xml_coll_t ci(x_mod, _U(module_component));
    for (ci.reset(), total_thickness = 0.0; ci; ++ci) {
      total_thickness += xml_comp_t(ci).thickness();
    }

    Box m_solid(x1/2.0, y1/2.0, total_thickness/2.0);
    Volume m_volume(m_nam, m_solid, vacuum);
    m_volume.setVisAttributes(description.visAttributes(x_mod.visStr()));

    Solid  frame_s;
    if(x_mod.hasChild("frame")){
      // todo
      // build frame from trd (assumed to be smaller)
      //xml_comp_t m_frame         = x_mod.child(_U(frame));
      //xml_comp_t f_pos           = m_frame.child(_U(position));
      //xml_comp_t frame_trd       = m_frame.trd();
      //double     frame_thickness = getAttrOrDefault(m_frame, _U(thickness), total_thickness);
      //double     frame_x1        = frame_trd.x1();
      //double     frame_x2        = frame_trd.x2();
      //double     frame_z         = frame_trd.z();
      //// make the frame match the total thickness if thickness attribute is not given
      //Trapezoid        f_solid1(x1, x2,frame_thickness / 2.0, frame_thickness / 2.0, z);
      //Trapezoid        f_solid(frame_x1, frame_x2, frame_thickness / 2.0, frame_thickness / 2.0, frame_z) ;
      //SubtractionSolid frame_shape(f_solid1, f_solid);
      //frame_s = frame_shape;

      //Material f_mat  = description.material(m_frame.materialStr());
      //Volume f_vol(m_nam + "_frame", frame_shape, f_mat);
      //f_vol.setVisAttributes(description.visAttributes(m_frame.visStr()));

      //// figure out how to best place
      //pv = m_volume.placeVolume(f_vol, Position(f_pos.x(), f_pos.y(),  f_pos.z()));
    }

    double     posZ            = -total_thickness/2.0;
    for (ci.reset(), n_sensor = 1, c_id = 0, posZ = -total_thickness/2.0; ci; ++ci, ++c_id) {
      xml_comp_t c           = ci;
      double     c_thick     = c.thickness();
      auto       comp_x     = getAttrOrDefault(c, _Unicode(x), x1);
      auto       comp_y     = getAttrOrDefault(c, _Unicode(y), y1);

      Material c_mat  = description.material(c.materialStr());
      string   c_name = _toString(c_id, "component%d");

      Box comp_s1(comp_x/2.0, comp_y/2.0, c_thick / 2e0);
      Solid  comp_shape = comp_s1;
      //if(frame_s.isValid()) {
      //  comp_shape = SubtractionSolid( comp_s1, frame_s); 
      //}
      Volume   c_vol(c_name, comp_shape, c_mat);

      c_vol.setVisAttributes(description.visAttributes(c.visStr()));
      pv = m_volume.placeVolume(c_vol, Position(0, 0,posZ + c_thick / 2.0));
      if (c.isSensitive()) {
        //std::cout << " adding sensitive volume" << c_name << "\n";
        sdet.check(n_sensor > 2, "SiTrackerEndcap2::fromCompact: " + c_name + " Max of 2 modules allowed!");
        pv.addPhysVolID("sensor", n_sensor);
        sens.setType("tracker");
        c_vol.setSensitiveDetector(sens);
        sensitives[m_nam].push_back(pv);
        ++n_sensor;
      }
      posZ += c_thick;
    }
    modules[m_nam] = m_volume;
  }

  for (xml_coll_t li(x_det, _U(layer)); li; ++li) {
    xml_comp_t x_layer(li);
    int        l_id    = x_layer.id();
    int        mod_num = 1;

    xml_comp_t l_env      = x_layer.child(_U(envelope));
    string     layer_name = det_name + std::string("_layer") + std::to_string(l_id);

    std::string layer_vis    = l_env.attr<std::string>(_Unicode(vis));
    double      layer_x  = l_env.attr<double>(_Unicode(x));
    double      layer_y   = l_env.attr<double>(_Unicode(y));
    double      layer_length = l_env.attr<double>(_Unicode(length));
    double      layer_zstart = l_env.attr<double>(_Unicode(zstart));
    double      layer_center_z =  layer_zstart + layer_length/2.0;
    //printout(INFO,"ROOTGDMLParse","+++ Read geometry from GDML file file:%s",input.c_str());
    //std::cout << "SiTracker Endcap layer " << l_id << " zstart = " << layer_zstart/dd4hep::mm << "mm ( " << layer_length/dd4hep::mm << " mm thick )\n";

    Assembly    layer_vol(layer_name);
    //assembly.placeVolume(layer_assembly);
    //Tube       layer_tub(layer_rmin, layer_rmax, layer_length / 2);
    //Volume     layer_vol(layer_name, layer_tub, air); // Create the layer envelope volume.
    layer_vol.setVisAttributes(description.visAttributes(layer_vis));

    PlacedVolume layer_pv;
    //if (reflect) {
    //  layer_pv =
    //      assembly.placeVolume(layer_vol, Transform3D(RotationZYX(0.0, -M_PI, 0.0), Position(0, 0, -layer_center_z)));
    //  layer_pv.addPhysVolID("barrel", 3).addPhysVolID("layer", l_id);
    //  layer_name += "_N";
    //} else {
      layer_pv = assembly.placeVolume(layer_vol, Position(0, 0, layer_center_z));
      layer_pv.addPhysVolID("layer", l_id);
    //}
    DetElement layer_element(sdet, layer_name, l_id);
    layer_element.setPlacement(layer_pv);
    Acts::ActsExtension* layerExtension = new Acts::ActsExtension();
    layerExtension->addType("layer", "layer");
    //layerExtension->addType("axes", "definitions", "XZY");
    //layerExtension->addType("sensitive disk", "layer");
    //layerExtension->addType("axes", "definitions", "XZY");
    layer_element.addExtension<Acts::ActsExtension>(layerExtension);

    string      m_nam    = x_layer.moduleStr();
    Volume      m_vol    = modules[m_nam];
    Placements& sensVols = sensitives[m_nam];

    DetElement module(layer_element, "module_", l_id);
    pv = layer_vol.placeVolume(m_vol, Position(0, 0, 0));
    pv.addPhysVolID("layer", l_id).addPhysVolID("module", mod_num);
    module.setPlacement(pv);
    for (size_t ic = 0; ic < sensVols.size(); ++ic) {
      PlacedVolume sens_pv = sensVols[ic];
      DetElement   comp_elt(module, sens_pv.volume().name(), mod_num);
      comp_elt.setPlacement(sens_pv);
      Acts::ActsExtension* moduleExtension = new Acts::ActsExtension();
      comp_elt.addExtension<Acts::ActsExtension>(moduleExtension);
    }

    // for (xml_coll_t ri(x_layer, _U(ring)); ri; ++ri) {
    //  xml_comp_t  x_ring   = ri;
    //  double      r        = x_ring.r();
    //  double      phi0     = x_ring.phi0(0);
    //  double      zstart   = x_ring.zstart();
    //  double      dz       = x_ring.dz(0);
    //  int         nmodules = x_ring.nmodules();
    //  string      m_nam    = x_ring.moduleStr();
    //  Volume      m_vol    = modules[m_nam];
    //  double      iphi     = 2 * M_PI / nmodules;
    //  double      dphi     = dd4hep::getAttrOrDefault(x_ring,_Unicode(dphi),iphi);
    //  double      phi      = phi0;
    //  Placements& sensVols = sensitives[m_nam];

    //  for (int k = 0; k < nmodules; ++k) {
    //    string     m_base = _toString(l_id, "layer%d") + _toString(mod_num, "_module%d");
    //    double     x      = -r * std::cos(phi);
    //    double     y      = -r * std::sin(phi);

    //    //if (!reflect) {
    //      DetElement module(layer_element, m_base + "_pos", det_id);
    //      pv = layer_vol.placeVolume(
    //          m_vol, Transform3D(RotationZYX(0, -M_PI / 2 - phi, -M_PI / 2), Position(x, y, zstart + dz)));
    //      pv.addPhysVolID("layer", l_id).addPhysVolID("module", mod_num);
    //      module.setPlacement(pv);
    //      for (size_t ic = 0; ic < sensVols.size(); ++ic) {
    //        PlacedVolume sens_pv = sensVols[ic];
    //        DetElement   comp_elt(module, sens_pv.volume().name(), mod_num);
    //        comp_elt.setPlacement(sens_pv);
    //    //std::cout << " adding ACTS extension" << "\n";
    //        Acts::ActsExtension* moduleExtension = new Acts::ActsExtension();
    //        comp_elt.addExtension<Acts::ActsExtension>(moduleExtension);
    //      }
    //    //} else {
    //    //  pv = layer_vol.placeVolume(
    //    //      m_vol, Transform3D(RotationZYX(0, -M_PI / 2 - phi, -M_PI / 2), Position(x, y, -zstart - dz)));
    //    //  pv.addPhysVolID("barrel", 2).addPhysVolID("layer", l_id).addPhysVolID("module", mod_num);
    //    //  DetElement r_module(layer_element, m_base + "_neg", det_id);
    //    //  r_module.setPlacement(pv);
    //    //  for (size_t ic = 0; ic < sensVols.size(); ++ic) {
    //    //    PlacedVolume sens_pv = sensVols[ic];
    //    //    DetElement   comp_elt(r_module, sens_pv.volume().name(), mod_num);
    //    //    comp_elt.setPlacement(sens_pv);
    //    ////std::cout << " adding ACTS extension" << "\n";
    //    //    Acts::ActsExtension* moduleExtension = new Acts::ActsExtension("XZY");
    //    //    comp_elt.addExtension<Acts::ActsExtension>(moduleExtension);
    //    //  }
    //    //}
    //    dz = -dz;
    //    phi += dphi;
    //    ++mod_num;
    //  }
    //}
    ++mod_num;
  }
  Transform3D posAndRot(RotationZYX(rot.z(), rot.y(), rot.x()), Position(pos.x(), pos.y(), pos.z()));
  pv = motherVol.placeVolume(assembly,posAndRot );
  pv.addPhysVolID("system", det_id);
  sdet.setPlacement(pv);
  return sdet;
}

// clang-format off
DECLARE_DETELEMENT(ip6_OffMomentumTracker, create_OffMomentumTracker)
