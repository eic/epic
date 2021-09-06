#include <map>
#include "DD4hep/DetFactoryHelper.h"
#include "Acts/Plugins/DD4hep/ActsExtension.hpp"
#include "Acts/Definitions/Units.hpp"

using namespace std;
using namespace dd4hep;
using namespace dd4hep::detail;

/*! B0 Tracker.
 *
 * @author Whitney Armstrong
 *
 */
static Ref_t create_B0Preshower(Detector& description, xml_h e, SensitiveDetector sens)
{
  typedef vector<PlacedVolume> Placements;
  xml_det_t                    x_det    = e;
  Material                     vacuum   = description.vacuum();
  int                          det_id   = x_det.id();
  string                       det_name = x_det.nameStr();
  bool                         reflect  = x_det.reflect(false);
  DetElement                   sdet(det_name, det_id);
  Assembly                     assembly(det_name);
  xml::Component pos   = x_det.position();
  xml::Component rot   = x_det.rotation();


  //Material  air  = description.material("Air");
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

  for (xml_coll_t mi(x_det, _U(module)); mi; ++mi, ++m_id) {
    xml_comp_t x_mod = mi;
    string     m_nam = x_mod.nameStr();
    xml_comp_t trd   = x_mod.trd();

    double     posY;
    double     x1 = trd.x1();
    double     x2 = trd.x2();
    double     z  = trd.z();
    double     y1, y2, total_thickness = 0.;
    xml_coll_t ci(x_mod, _U(module_component));
    for (ci.reset(), total_thickness = 0.0; ci; ++ci)
      total_thickness += xml_comp_t(ci).thickness();

    y1 = y2 = total_thickness / 2;
    Trapezoid m_solid(x1, x2, y1, y2, z);
    Volume m_volume(m_nam, m_solid, vacuum);
    m_volume.setVisAttributes(description.visAttributes(x_mod.visStr()));

    Solid  frame_s;
    if(x_mod.hasChild("frame")){
      // build frame from trd (assumed to be smaller)
      xml_comp_t m_frame         = x_mod.child(_U(frame));
      xml_comp_t f_pos           = m_frame.child(_U(position));
      xml_comp_t frame_trd       = m_frame.trd();
      double     frame_thickness = getAttrOrDefault(m_frame, _U(thickness), total_thickness);
      double     frame_x1        = frame_trd.x1();
      double     frame_x2        = frame_trd.x2();
      double     frame_z         = frame_trd.z();
      // make the frame match the total thickness if thickness attribute is not given
      Trapezoid        f_solid1(x1, x2,frame_thickness / 2.0, frame_thickness / 2.0, z);
      Trapezoid        f_solid(frame_x1, frame_x2, frame_thickness / 2.0, frame_thickness / 2.0, frame_z) ;
      SubtractionSolid frame_shape(f_solid1, f_solid);
      frame_s = frame_shape;

      Material f_mat  = description.material(m_frame.materialStr());
      Volume f_vol(m_nam + "_frame", frame_shape, f_mat);
      f_vol.setVisAttributes(description.visAttributes(m_frame.visStr()));

      // figure out how to best place
      pv = m_volume.placeVolume(f_vol, Position(f_pos.x(), f_pos.y(),  f_pos.z()));
    }

    for (ci.reset(), n_sensor = 1, c_id = 0, posY = -y1; ci; ++ci, ++c_id) {
      xml_comp_t c           = ci;
      double     c_thick     = c.thickness();
      auto       comp_x1     = getAttrOrDefault(c, _Unicode(x1), x1);
      auto       comp_x2     = getAttrOrDefault(c, _Unicode(x2), x2);
      auto       comp_height = getAttrOrDefault(c, _Unicode(height), z);

      Material c_mat  = description.material(c.materialStr());
      string   c_name = _toString(c_id, "component%d");

      Trapezoid comp_s1(comp_x1, comp_x2, c_thick / 2e0, c_thick / 2e0, comp_height);
      Solid  comp_shape = comp_s1;
      if(frame_s.isValid()) {
        comp_shape = SubtractionSolid( comp_s1, frame_s); 
      }
      Volume   c_vol(c_name, comp_shape, c_mat);

      c_vol.setVisAttributes(description.visAttributes(c.visStr()));
      pv = m_volume.placeVolume(c_vol, Position(0, posY + c_thick / 2, 0));
      if (c.isSensitive()) {
        //std::cout << " adding sensitive volume" << c_name << "\n";
        sdet.check(n_sensor > 2, "SiTrackerEndcap2::fromCompact: " + c_name + " Max of 2 modules allowed!");
        pv.addPhysVolID("sensor", n_sensor);
        sens.setType("tracker");
        c_vol.setSensitiveDetector(sens);
        sensitives[m_nam].push_back(pv);
        ++n_sensor;
      }
      posY += c_thick;
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
    //double      layer_rmin   = l_env.attr<double>(_Unicode(rmin));
    //double      layer_rmax   = l_env.attr<double>(_Unicode(rmax));
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
    if (reflect) {
      layer_pv =
          assembly.placeVolume(layer_vol, Transform3D(RotationZYX(0.0, -M_PI, 0.0), Position(0, 0, -layer_center_z)));
      layer_pv.addPhysVolID("barrel", 3).addPhysVolID("layer", l_id);
      layer_name += "_N";
    } else {
      layer_pv = assembly.placeVolume(layer_vol, Position(0, 0, layer_center_z));
      layer_pv.addPhysVolID("barrel", 2).addPhysVolID("layer", l_id);
      layer_name += "_P";
    }
    DetElement layer_element(sdet, layer_name, l_id);
    layer_element.setPlacement(layer_pv);
    Acts::ActsExtension* layerExtension = new Acts::ActsExtension();
    layerExtension->addType("layer", "layer");
    //layerExtension->addType("axes", "definitions", "XZY");
    //layerExtension->addType("sensitive disk", "layer");
    //layerExtension->addType("axes", "definitions", "XZY");
    layer_element.addExtension<Acts::ActsExtension>(layerExtension);

    for (xml_coll_t ri(x_layer, _U(ring)); ri; ++ri) {
      xml_comp_t  x_ring   = ri;
      double      r        = x_ring.r();
      double      phi0     = x_ring.phi0(0);
      double      zstart   = x_ring.zstart();
      double      dz       = x_ring.dz(0);
      int         nmodules = x_ring.nmodules();
      string      m_nam    = x_ring.moduleStr();
      Volume      m_vol    = modules[m_nam];
      double      iphi     = 2 * M_PI / nmodules;
      double      dphi     = dd4hep::getAttrOrDefault(x_ring,_Unicode(dphi),iphi);
      double      phi      = phi0;
      Placements& sensVols = sensitives[m_nam];

      for (int k = 0; k < nmodules; ++k) {
        string     m_base = _toString(l_id, "layer%d") + _toString(mod_num, "_module%d");
        double     x      = -r * std::cos(phi);
        double     y      = -r * std::sin(phi);

        if (!reflect) {
          DetElement module(layer_element, m_base + "_pos", det_id);
          pv = layer_vol.placeVolume(
              m_vol, Transform3D(RotationZYX(0, -M_PI / 2 - phi, -M_PI / 2), Position(x, y, zstart + dz)));
          pv.addPhysVolID("barrel", 1).addPhysVolID("layer", l_id).addPhysVolID("module", mod_num);
          module.setPlacement(pv);
          for (size_t ic = 0; ic < sensVols.size(); ++ic) {
            PlacedVolume sens_pv = sensVols[ic];
            DetElement   comp_elt(module, sens_pv.volume().name(), mod_num);
            comp_elt.setPlacement(sens_pv);
        //std::cout << " adding ACTS extension" << "\n";
            Acts::ActsExtension* moduleExtension = new Acts::ActsExtension("XZY");
            comp_elt.addExtension<Acts::ActsExtension>(moduleExtension);
          }
        } else {
          pv = layer_vol.placeVolume(
              m_vol, Transform3D(RotationZYX(0, -M_PI / 2 - phi, -M_PI / 2), Position(x, y, -zstart - dz)));
          pv.addPhysVolID("barrel", 2).addPhysVolID("layer", l_id).addPhysVolID("module", mod_num);
          DetElement r_module(layer_element, m_base + "_neg", det_id);
          r_module.setPlacement(pv);
          for (size_t ic = 0; ic < sensVols.size(); ++ic) {
            PlacedVolume sens_pv = sensVols[ic];
            DetElement   comp_elt(r_module, sens_pv.volume().name(), mod_num);
            comp_elt.setPlacement(sens_pv);
        //std::cout << " adding ACTS extension" << "\n";
            Acts::ActsExtension* moduleExtension = new Acts::ActsExtension("XZY");
            comp_elt.addExtension<Acts::ActsExtension>(moduleExtension);
          }
        }
        dz = -dz;
        phi += dphi;
        ++mod_num;
      }
    }
  }
  Transform3D posAndRot(RotationZYX(rot.z(), rot.y(), rot.x()), Position(pos.x(), pos.y(), pos.z()));
  pv = motherVol.placeVolume(assembly,posAndRot );
  pv.addPhysVolID("system", det_id);
  sdet.setPlacement(pv);
  return sdet;
}

// clang-format off
DECLARE_DETELEMENT(ip6_B0Preshower, create_B0Preshower)
