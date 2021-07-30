#include "DD4hep/DetFactoryHelper.h"
#include <map>

using namespace std;
using namespace dd4hep;
using namespace dd4hep::detail;

using Placements = vector<PlacedVolume>;

static Ref_t create_detector(Detector& description, xml_h e, SensitiveDetector sens)  {
  xml_det_t      x_det    = e;
  Material       air      = description.air();
  Material       vacuum   = description.vacuum();
  string         det_name = x_det.nameStr();
  xml::Component pos      = x_det.position();
  xml::Component rot      = x_det.rotation();
  DetElement     sdet(det_name, x_det.id());
  Assembly       assembly(det_name);
  sens.setType("tracker");

  PlacedVolume pv;

  map<string, Volume>     modules;
  map<string, Placements> sensitives;
  map<string, Volume>     module_assemblies;

  int m_id = 0;
  // mi ~ module iterator
  for (xml_coll_t mi(x_det, _U(module)); mi; ++mi, ++m_id) {
    xml_comp_t x_mod = mi;
    string     m_nam = x_mod.nameStr();
    double     mod_width  = getAttrOrDefault(x_mod, _U(width), 3.2*cm);
    double     mod_height = getAttrOrDefault(x_mod, _U(height), 3.2*cm);
    double     mod_total_thickness = 0.;

    xml_coll_t ci(x_mod, _U(module_component));
    for (ci.reset(), mod_total_thickness = 0.0; ci; ++ci)
      mod_total_thickness += xml_comp_t(ci).thickness();

    Box m_solid(mod_width/2.0, mod_height/2.0, mod_total_thickness/2.0);
    Volume m_volume(m_nam, m_solid, vacuum);
    m_volume.setVisAttributes(description.visAttributes(x_mod.visStr()));

    double comp_z_pos = -mod_total_thickness/2.0;
    int n_sensor = 1;
    int c_id;
    for (ci.reset(), n_sensor = 1, c_id = 0; ci; ++ci, ++c_id) {
      xml_comp_t c       = ci;
      double     c_thick = c.thickness();
      double     comp_x  = getAttrOrDefault(c, _Unicode(width), mod_width);
      double     comp_y  = getAttrOrDefault(c, _Unicode(height), mod_height);

      Material c_mat  = description.material(c.materialStr());
      string   c_name = _toString(c_id, "component%d");

      Box comp_s1(comp_x/2.0, comp_y/2.0, c_thick / 2.0);
      Solid  comp_shape = comp_s1;
      Volume c_vol(c_name, comp_shape, c_mat);
      c_vol.setVisAttributes(description.visAttributes(c.visStr()));

      Position c_position(0, 0.0,  comp_z_pos +c_thick/2.0);
      if(c.hasChild("position")){
        xml_comp_t c_pos = c.child(_U(position));
        c_position       = Position(c_pos.x(), c_pos.y(), c_pos.z());
      }
      pv = m_volume.placeVolume(c_vol, c_position);
      if (c.isSensitive()) {
        //sdet.check(n_sensor > 2, "SiTrackerEndcap2::fromCompact: " + c_name + " Max of 2 modules allowed!");
        pv.addPhysVolID("sensor", n_sensor);
        c_vol.setSensitiveDetector(sens);
        sensitives[m_nam].push_back(pv);
        ++n_sensor;
      }
      comp_z_pos += c_thick;
    }
    modules[m_nam] = m_volume;
  }

  std::map<std::string,DetElement> module_assembly_delements;
  // module assemblies
  for (xml_coll_t ma(x_det, _Unicode(module_assembly)); ma; ++ma) {
    xml_comp_t x_ma    = ma;
    string     ma_name = x_ma.nameStr();
    Assembly ma_vol(ma_name);
    DetElement ma_de(ma_name,x_det.id());
    module_assemblies[ma_name] = ma_vol;
    module_assembly_delements[ma_name] = ma_de;

    int i_mod =0;
    // array of modules
    for (xml_coll_t ai(x_ma, _Unicode(array)); ai; ++ai) {
      xml_comp_t x_array = ai;
      double nx      = getAttrOrDefault(x_array, _Unicode(nx), 1);
      double ny      = getAttrOrDefault(x_array, _Unicode(ny), 1);
      double dz      = getAttrOrDefault(x_array, _Unicode(dz), 0*mm);
      double arr_width       = getAttrOrDefault(x_array, _Unicode(width ), 3.2*cm);
      double arr_height      = getAttrOrDefault(x_array, _Unicode(height), 3.2*cm);
      std::string arr_module      = getAttrOrDefault(x_array, _Unicode(module), "");
      // TODO: add check here
      auto arr_vol = modules[arr_module];
      Placements& sensVols = sensitives[arr_module];

      double arr_x_delta = arr_width/double(nx);
      double arr_y_delta = arr_height/double(ny);

      xml_comp_t x_pos  = x_array.position(false);

      for (int ix = 0; ix < nx; ix++) {
        for (int iy = 0; iy < ny; iy++) {
          double z_pos = dz; //(ix%2)? dz : -dz;
          i_mod++;
          Position arr_pos(-arr_width / 2.0 + arr_x_delta / 2.0 + ix * arr_x_delta,
                           -arr_height / 2.0 + arr_y_delta / 2.0 + iy * arr_y_delta, z_pos);
          if (x_pos) {
            arr_pos += Position(x_pos.x(0), x_pos.y(0), x_pos.z(0));
          }
          DetElement mod_de(ma_de, ma_name + std::string("_mod") + std::to_string(i_mod), i_mod);
          pv = ma_vol.placeVolume(arr_vol, arr_pos);
          pv.addPhysVolID("module", i_mod);
          mod_de.setPlacement(pv);
          for (size_t ic = 0; ic < sensVols.size(); ++ic) {
            PlacedVolume sens_pv = sensVols[ic];
            DetElement   comp_de(mod_de, std::string("de_") + sens_pv.volume().name(), ic+1);
            comp_de.setPlacement(sens_pv);
            //Acts::ActsExtension* sensorExtension = new Acts::ActsExtension();
            //// sensorExtension->addType("sensor", "detector");
            //comp_de.addExtension<Acts::ActsExtension>(sensorExtension);
            //// comp_de.setAttributes(description, sens_pv.volume(), x_layer.regionStr(),
            //// x_layer.limitsStr(),
            ////                      xml_det_t(xmleles[m_nam]).visStr());
          }
        }
      }
    }
  }

  int l_num = 0;
  for (xml_coll_t i(x_det, _U(layer)); i; ++i, ++l_num) {
    xml_comp_t x_layer    = i;
    string     l_nam      = det_name + _toString(l_num, "_layer%d");
    xml_comp_t l_pos      = x_layer.position(false);
    Assembly l_vol(l_nam); //(l_nam, l_box, air);

    Position layer_pos(0, 0,0);
    if(l_pos) {
      layer_pos = Position(l_pos.x(),l_pos.y(),l_pos.z());
    }
    DetElement layer(sdet, l_nam + "_pos", l_num);

    int i_assembly = 1;
    xml_coll_t ci(x_layer, _U(component));
    for (ci.reset(); ci; ++ci) {
      xml_comp_t x_comp  = ci;
      xml_comp_t c_pos   = x_comp.position(false);

      // string     ma_name = x_comp.nameStr();
      string comp_assembly = getAttrOrDefault(x_comp, _Unicode(assembly), "");


      auto comp_vol = module_assemblies[comp_assembly];
      //auto de       = ;
      auto comp_de  = module_assembly_delements[comp_assembly].clone(comp_assembly + std::to_string(l_num));
      if (c_pos) {
        pv = l_vol.placeVolume(comp_vol, Position(c_pos.x(), c_pos.y(), c_pos.z()));
      } else {
        pv = l_vol.placeVolume(comp_vol);
      }
      pv.addPhysVolID("assembly", i_assembly);
      comp_de.setPlacement(pv);
      layer.add(comp_de);
      i_assembly++;
      // DetElement det = module > 1 ? stave.clone(_toString(module,"stave%d")) : stave;
      // Transform3D trafo(RotationZYX(0, rotY, rotX), Translation3D(-posX, -posY, 0));
      // PlacedVolume pv = envelopeVolume.placeVolume(sectVolume,trafo);
      //// Not a valid volID: pv.addPhysVolID("stave", 0);
      // pv.addPhysVolID("module", module);
      // det.setPlacement(pv);
      // parent.add(det);
    }
    pv = assembly.placeVolume(l_vol, l_pos);
    pv.addPhysVolID("layer", l_num);
  }

  //pv = description.pickMotherVolume(sdet).placeVolume(assembly, Position(pos.x(), pos.y(), pos.z()));
  Transform3D posAndRot(RotationZYX(rot.z(), rot.y(), rot.x()), Position(pos.x(), pos.y(), pos.z()));
  //pv = description.pickMotherVolume(sdet).placeVolume(assembly, Position(pos.x(), pos.y(), pos.z()));
  pv = description.pickMotherVolume(sdet).placeVolume(assembly, posAndRot);
  pv.addPhysVolID("system", x_det.id()); // Set the subdetector system ID.
  sdet.setPlacement(pv);
  return sdet;
}

DECLARE_DETELEMENT(ip6_ForwardRomanPot, create_detector)
