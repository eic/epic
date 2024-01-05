#include "DD4hep/DetFactoryHelper.h"
#include "TVector3.h"
#include "XML/Layering.h"
using namespace std;
using namespace dd4hep;
using namespace dd4hep::detail;


static Ref_t create_detector(Detector& description, xml_h e, [[maybe_unused]] SensitiveDetector sens)
{
  xml_det_t x_det    = e;
  //xml_dim_t dim      = x_det.dimensions();
  int       det_id   = x_det.id();
  string    det_name = x_det.nameStr();
  Material  air      = description.air();
  xml::Component x_pos  = x_det.position();


  Assembly     assembly(det_name);
  DetElement   sdet(det_name, det_id);
  PlacedVolume pv;

  double disksGap = 0.0;

  Layering layering(x_det);
  //double   totalThickness = layering.totalThickness();

  // Looping through all the different layer sections
  for (xml_coll_t xc(x_det, _U(layer)); xc; ++xc) {
    xml_comp_t x_layer         = xc;

    int     layer_id = x_layer.id();
    double     layer_rmax = x_layer.rmax();
    double     layer_rmin = x_layer.rmin();
    double     layer_thickness = x_layer.thickness();
    double     layer_zpos = x_layer.zpos();

    Material   l_mat   = description.material(x_layer.materialStr());


    DetElement disk_ele("disk_ele", layer_id);
    Volume disk(x_layer.nameStr(), Tube(layer_rmin, layer_rmax, layer_thickness / 2, 0.0, 2.0 * M_PI), air);
    disk.setVisAttributes(description.visAttributes(x_layer.visStr()));

    Volume halfdisk("halfdisk", Tube(layer_rmin, layer_rmax, layer_thickness / 2, M_PI / 2, M_PI * 3 / 2), l_mat);
    halfdisk.setVisAttributes(description.visAttributes(x_layer.visStr()));
    PlacedVolume s_phv1 = disk.placeVolume(halfdisk, Position(-disksGap / 2, 0, 0));
    s_phv1.addPhysVolID("halfdisk", 0);
    PlacedVolume s_phv2 = disk.placeVolume(halfdisk, Transform3D(RotationZYX(M_PI, 0, 0), Position(+disksGap / 2, 0, 0)));
    s_phv2.addPhysVolID("halfdisk", 1);


    pv = assembly.placeVolume(disk, Position(0, 0, -layer_thickness/2-layer_zpos));
    pv.addPhysVolID(x_layer.nameStr(), layer_id);
    disk_ele.setPlacement(pv);
  }

  // Get position and place volume
  Position       pos(x_pos.x(), x_pos.y(), x_pos.z());
  pv = description.pickMotherVolume(sdet).placeVolume(assembly, pos);
  pv.addPhysVolID("system", det_id).addPhysVolID("barrel", 0);
  sdet.setPlacement(pv);
  return sdet;
}

// clang-format off
DECLARE_DETELEMENT(epic_EndcapFluxReturnN, create_detector)
