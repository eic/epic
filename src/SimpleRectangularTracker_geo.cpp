#include "DD4hep/DetFactoryHelper.h"

using namespace std;
using namespace dd4hep;
using namespace dd4hep::detail;

static Ref_t create_detector(Detector& description, xml_h e, SensitiveDetector sens)  {
  xml_det_t  x_det     = e;
  Material   air       = description.air();
  string     det_name  = x_det.nameStr();
  bool       reflect   = x_det.reflect();
  DetElement sdet(det_name,x_det.id());
  Assembly   assembly(det_name);
  sens.setType("tracker");

  PlacedVolume pv;
  int l_num = 0;
  xml::Component  pos  = x_det.position();

  for(xml_coll_t i(x_det,_U(layer)); i; ++i, ++l_num)  {
    xml_comp_t x_layer    = i;
    string     l_nam      = det_name + _toString(l_num, "_layer%d");
    double     x_lay      = x_layer.x();
    double     y_lay      = x_layer.y();
    double     z          = 0;
    double     zmin       = 0;
    double     layerWidth = 0.;
    int        m_num      = 0;
    for(xml_coll_t j(x_layer,_U(slice)); j; ++j)  {
      double thickness = xml_comp_t(j).thickness();
      layerWidth += thickness;
    }
    Box    l_box(x_lay/2.0, y_lay/2.0, layerWidth/2.0 );
    Volume l_vol(l_nam, l_box, air);
    l_vol.setVisAttributes(description, x_layer.visStr());
    
    for (xml_coll_t j(x_layer, _U(module)); j; ++j, ++m_num) {
      xml_comp_t x_module = j;
      string     m_nam   = l_nam + _toString(m_num, "_module%d");
      double     x_mod      = x_module.x();
      double     y_mod      = x_module.y();
      Volume     m_vol(m_nam, Box(x_mod/2.0, y_mod/2.0, layerWidth/2.0), air);
      xml::Component  m_pos  = x_module.position();
      int         s_num      = 0;

      for (xml_coll_t k(x_module, _U(slice)); k; ++k, ++s_num) {
        xml_comp_t x_slice = k;
        double     thick   = x_slice.thickness();
        Material   mat     = description.material(x_slice.materialStr());
        string     s_nam   = m_nam + _toString(s_num, "_slice%d");
        Volume     s_vol(s_nam, Box(x_lay/2.0, y_lay/2.0, thick/2.0), mat);
        if (x_slice.isSensitive()) {
         sens.setType("tracker");
         s_vol.setSensitiveDetector(sens);
        }
        s_vol.setAttributes(description, x_slice.regionStr(), x_slice.limitsStr(), x_slice.visStr());
        pv = s_vol.placeVolume(s_vol, Position(0, 0, z - zmin - layerWidth / 2 + thick / 2));
        pv.addPhysVolID("slice", s_num);
      }

      m_vol.setAttributes(description, x_module.regionStr(), x_module.limitsStr(), x_module.visStr());
      pv = m_vol.placeVolume(m_vol, Position(m_pos.x(), m_pos.y(), m_pos.y()));
      pv.addPhysVolID("module", m_num);
    }
    
    DetElement layer(sdet,l_nam+"_pos",l_num);
    pv = assembly.placeVolume(l_vol,Position(0,0,zmin+layerWidth/2.));
    pv.addPhysVolID("layer",l_num);
    layer.setPlacement(pv);
    if ( reflect )  {
      pv = assembly.placeVolume(l_vol,Transform3D(RotationY(M_PI),Position(0,0,-zmin-layerWidth/2)));
      pv.addPhysVolID("layer",l_num);
      DetElement layerR = layer.clone(l_nam+"_neg");
      sdet.add(layerR.setPlacement(pv));
    }
  }
  if ( x_det.hasAttr(_U(combineHits)) ) {
    sdet.setCombineHits(x_det.attr<bool>(_U(combineHits)),sens);
  }
  pv = description.pickMotherVolume(sdet).placeVolume(assembly,Position(pos.x(),pos.y(),pos.z()));
  pv.addPhysVolID("system", x_det.id());      // Set the subdetector system ID.
  sdet.setPlacement(pv);
  return sdet;
}

DECLARE_DETELEMENT(ip6_RectangularTracker,create_detector)
