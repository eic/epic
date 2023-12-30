#include "DD4hep/DetFactoryHelper.h"
#include "TVector3.h"
#include "XML/Layering.h"
using namespace std;
using namespace dd4hep;
using namespace dd4hep::detail;


static Ref_t create_detector(Detector& description, xml_h e, [[maybe_unused]] SensitiveDetector sens)
{
  xml_det_t x_det    = e;
  xml_dim_t dim      = x_det.dimensions();
  int       det_id   = x_det.id();
  string    det_name = x_det.nameStr();
  Material  air      = description.air();
  // int            numsides = dim.numsides();
  xml::Component pos  = x_det.position();
  double         rmin = dim.rmin();
  double         rmax = dim.rmax();
  double         zmin = dim.zmin();

  double disksGap = 0.0;

  Layering layering(x_det);
  double   totalThickness = layering.totalThickness();

  Volume endcapVol("endcapFlux", Tube(rmin, rmax + 2 * cm, totalThickness), air); // rmax + 2 * cm -to account for disk gap
  DetElement endcap("endcapFlux", det_id);

  endcapVol.setAttributes(description, x_det.regionStr(), x_det.limitsStr(), x_det.visStr());

  //int    iLayer  = 0;
  //double globalZ = zmin;
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

    //disk_ele.setPlacement(s_phv1);
    //disk_ele.setPlacement(s_phv2);


    PlacedVolume s_phv = endcapVol.placeVolume(disk, Position(0, 0, -layer_thickness/2-layer_zpos));
    s_phv.addPhysVolID(x_layer.nameStr(), layer_id);
    disk_ele.setPlacement(s_phv);

  }

  double       z_pos = zmin + totalThickness / 2;
  PlacedVolume pv;
  // Reflect it.
  Assembly   assembly(det_name);
  DetElement endcapAssyDE(det_name, det_id);
  Volume     motherVol = description.pickMotherVolume(endcapAssyDE);

	pv = assembly.placeVolume(endcapVol, Transform3D(RotationZYX(0, 0, 0), Position(0, 0, -z_pos)));
	pv.addPhysVolID("barrel", 1);
	Ref_t(endcap)->SetName((det_name + "_backward").c_str());
	endcap.setPlacement(pv);

  endcapAssyDE.add(endcap);
  pv = motherVol.placeVolume(assembly, Position(pos.x(), pos.y(), pos.z()));
  pv.addPhysVolID("system", det_id);
  endcapAssyDE.setPlacement(pv);
  return endcapAssyDE;
}

// clang-format off
DECLARE_DETELEMENT(epic_EndcapFluxReturnN, create_detector)
