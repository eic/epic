//==========================================================================
//  AIDA Detector description implementation
//--------------------------------------------------------------------------
// Copyright (C) Organisation europeenne pour la Recherche nucleaire (CERN)
// All rights reserved.
//
// For the licensing terms see $DD4hepINSTALL/LICENSE.
// For the list of contributors see $DD4hepINSTALL/doc/CREDITS.
//
// Author     : M.Frank
//
//==========================================================================
//
// Modified for ATHENA detector
//
//==========================================================================
#include "DD4hep/DetFactoryHelper.h"
#include "XML/Layering.h"

using namespace std;
using namespace dd4hep;
using namespace dd4hep::detail;

static Ref_t create_detector(Detector& description, xml_h e, SensitiveDetector sens)
{
  xml_det_t      x_det    = e;
  xml_dim_t      dim      = x_det.dimensions();
  int            det_id   = x_det.id();
  bool           reflect  = x_det.reflect(true);
  string         det_name = x_det.nameStr();
  Material       air      = description.air();
  int            numsides = dim.numsides();
  xml::Component pos      = x_det.position();
  double         rmin     = dim.rmin();
  double         rmax     = dim.rmax();
  double         zmin     = dim.zmin();
  Layering       layering(x_det);
  double         totalThickness = layering.totalThickness();
  Volume         endcapVol("endcap", PolyhedraRegular(numsides, rmin, rmax, totalThickness), air);
  DetElement     endcap("endcap", det_id);

  //std::cout << "totalThickness = " << totalThickness << "\n";
  //std::cout << "zmin = " << zmin << "\n";
  //std::cout << "rmin = " << rmin << "\n";
  //std::cout << "rmax = " << rmax << "\n";
  //std::cout << "nlayers = " << std::size(layering.layers()) << "\n";
  int    l_num     = 1;
  int    layerType = 0;
  double layerZ    = -totalThickness / 2;

  endcapVol.setAttributes(description, x_det.regionStr(), x_det.limitsStr(), x_det.visStr());

  for (xml_coll_t xc(x_det, _U(layer)); xc; ++xc) {
    //std::cout << "l_num = " << l_num << "\n";
    //std::cout << "xc = " << xc << "\n";
    xml_comp_t x_layer = xc;
    double     l_thick = layering.layer(l_num - 1)->thickness();
    //std::cout << "xc = " << xc << "\n";
    string               l_name   = _toString(layerType, "layer%d");
    int                  l_repeat = x_layer.repeat();
    Volume               l_vol(l_name, PolyhedraRegular(numsides, rmin, rmax, l_thick), air);
    vector<PlacedVolume> sensitives;

    int    s_num  = 1;
    double sliceZ = -l_thick / 2;
    for (xml_coll_t xs(x_layer, _U(slice)); xs; ++xs) {
      xml_comp_t x_slice = xs;
      string     s_name  = _toString(s_num, "slice%d");
      double     s_thick = x_slice.thickness();
      Material   s_mat   = description.material(x_slice.materialStr());
      Volume     s_vol(s_name, PolyhedraRegular(numsides, rmin, rmax, s_thick), s_mat);

      s_vol.setVisAttributes(description.visAttributes(x_slice.visStr()));
      sliceZ += s_thick / 2;
      PlacedVolume s_phv = l_vol.placeVolume(s_vol, Position(0, 0, sliceZ));
      s_phv.addPhysVolID("slice", s_num);
      if (x_slice.isSensitive()) {
        sens.setType("calorimeter");
        s_vol.setSensitiveDetector(sens);
        sensitives.push_back(s_phv);
      }
      sliceZ += s_thick / 2;
      s_num++;
    }
    l_vol.setVisAttributes(description.visAttributes(x_layer.visStr()));
    if (l_repeat <= 0)
      throw std::runtime_error(x_det.nameStr() + "> Invalid repeat value");
    for (int j = 0; j < l_repeat; ++j) {
      string phys_lay = _toString(l_num, "layer%d");
      layerZ += l_thick / 2;
      DetElement   layer_elt(endcap, phys_lay, l_num);
      PlacedVolume pv = endcapVol.placeVolume(l_vol, Position(0, 0, layerZ));
      pv.addPhysVolID("layer", l_num);
      layer_elt.setPlacement(pv);
      for (size_t ic = 0; ic < sensitives.size(); ++ic) {
        PlacedVolume sens_pv = sensitives[ic];
        DetElement   comp_elt(layer_elt, sens_pv.volume().name(), l_num);
        comp_elt.setPlacement(sens_pv);
      }
      layerZ += l_thick / 2;
      ++l_num;
    }
    ++layerType;
  }

  double       z_pos = zmin + totalThickness / 2;
  PlacedVolume pv;
  // Reflect it.
  Assembly   assembly(det_name);
  DetElement endcapAssyDE(det_name, det_id);
  Volume     motherVol = description.pickMotherVolume(endcapAssyDE);
  if (reflect) {
    pv = assembly.placeVolume(endcapVol, Transform3D(RotationZYX(M_PI / numsides, M_PI, 0), Position(0, 0, -z_pos)));
    pv.addPhysVolID("barrel", 2);
    Ref_t(endcap)->SetName((det_name + "_backward").c_str());
    endcap.setPlacement(pv);
  } else {
    pv = assembly.placeVolume(endcapVol, Transform3D(RotationZYX(M_PI / numsides, 0, 0), Position(0, 0, z_pos)));
    pv.addPhysVolID("barrel", 1);
    Ref_t(endcap)->SetName((det_name + "_forward").c_str());
    endcap.setPlacement(pv);
  }
  endcapAssyDE.add(endcap);
  pv = motherVol.placeVolume(assembly,Position(pos.x(),pos.y(),pos.z()));
  pv.addPhysVolID("system", det_id);
  endcapAssyDE.setPlacement(pv);
  return endcapAssyDE;
}

// clang-format off
DECLARE_DETELEMENT(athena_PolyhedraEndcapCalorimeter2, create_detector)
DECLARE_DETELEMENT(athena_PolyhedraEndcapCalorimeter, create_detector)
