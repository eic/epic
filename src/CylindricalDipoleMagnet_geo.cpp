#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/Printout.h"
#include "DD4hep/Shapes.h"
#include "DDRec/DetectorData.h"
#include "DDRec/Surface.h"
#include "TMath.h"
#include "XML/Layering.h"

using namespace std;
using namespace dd4hep;
using namespace dd4hep::rec;
using namespace ROOT::Math;

static Ref_t build_magnet(Detector& dtor, xml_h e, SensitiveDetector /* sens */)
{
  xml_det_t x_det     = e;
  int       det_id    = x_det.id();
  string    det_name  = x_det.nameStr();
  xml_dim_t pos       = x_det.child(_U(placement));
  double    pos_x     = pos.x();
  double    pos_y     = pos.y();
  double    pos_z     = pos.z();
  double    pos_theta = pos.attr<double>(_U(theta));
  xml_dim_t dims      = x_det.dimensions();
  double    dim_r     = dims.r();
  double    dim_z     = dims.z();
  xml_dim_t apperture = x_det.child(_Unicode(apperture));
  double    app_r     = apperture.r();
  xml_dim_t coil      = x_det.child(_Unicode(coil));
  double    coil_x    = coil.dx();
  double    coil_y    = coil.dy();
  Material  iron      = dtor.material("Iron");
  Material  niobium   = dtor.material("Niobium");

  // std::cout << det_name << " positioned at z=" << pos.z() << ", x=" << pos.x() << "\n";

  DetElement sdet(det_name, det_id);
  Assembly   assembly(det_name + "_assembly");

  const string module_name = "Quad_magnet";

  const string yoke_vis = dd4hep::getAttrOrDefault<std::string>(x_det, _Unicode(vis), "FFMagnetVis");
  const string coil_vis = dd4hep::getAttrOrDefault<std::string>(coil,  _Unicode(vis), "FFMagnetCoilVis");

  sdet.setAttributes(dtor, assembly, x_det.regionStr(), x_det.limitsStr(), yoke_vis);

  // -- yoke
  Tube   yoke_tube(app_r + coil_y, dim_r, 0.5 * dim_z);
  Volume yoke_vol("yoke_vol", yoke_tube, iron);
  auto   yoke_pv = assembly.placeVolume(yoke_vol);
  yoke_pv.addPhysVolID("element", 1);
  DetElement yoke_de(sdet, "yoke_de", 1);
  yoke_de.setPlacement(yoke_pv);
  yoke_de.setAttributes(dtor, yoke_vol, x_det.regionStr(), x_det.limitsStr(), yoke_vis);

  // -- coils
  double offset       = 1.5 * coil_x;
  double appc_r       = app_r + 0.5 * coil_y;
  double offset_angle = atan(offset / appc_r);

  Tube longrod(app_r, app_r + coil_y, 0.5 * dim_z, atan(coil_x / app_r));
  Tube connector(app_r, app_r + coil_y, coil_y, 0.5 * M_PI - offset_angle);

  UnionSolid coil1(longrod, longrod, Transform3D(RotationZ(-offset_angle)));
  UnionSolid coil2(coil1, longrod, Transform3D(RotationZ(0.5 * M_PI)));
  UnionSolid coil3(coil2, longrod, Transform3D(RotationZ(M_PI)));
  UnionSolid coil4(coil3, longrod, Transform3D(RotationZ(1.5 * M_PI)));

  UnionSolid coil5(coil4, longrod, Transform3D(RotationZ(0.5 * M_PI - offset_angle)));
  UnionSolid coil6(coil5, longrod, Transform3D(RotationZ(M_PI - offset_angle)));
  UnionSolid coil7(coil6, longrod, Transform3D(RotationZ(1.5 * M_PI - offset_angle)));

  UnionSolid coil8(coil7, connector, Transform3D(Translation3D(0.0, 0.0, -0.5 * dim_z + coil_y)));
  UnionSolid coil9(coil8, connector,
                   Transform3D(Translation3D(0.0, 0.0, -0.5 * dim_z + coil_y) * RotationZ(0.5 * M_PI)));
  UnionSolid coil10(coil9, connector, Transform3D(Translation3D(0.0, 0.0, -0.5 * dim_z + coil_y) * RotationZ(M_PI)));
  UnionSolid coil11(coil10, connector,
                    Transform3D(Translation3D(0.0, 0.0, -0.5 * dim_z + coil_y) * RotationZ(1.5 * M_PI)));
  UnionSolid coil12(coil11, connector, Transform3D(Translation3D(0.0, 0.0, 0.5 * dim_z - coil_y)));
  UnionSolid coil13(coil12, connector,
                    Transform3D(Translation3D(0.0, 0.0, 0.5 * dim_z - coil_y) * RotationZ(0.5 * M_PI)));
  UnionSolid coil14(coil13, connector, Transform3D(Translation3D(0.0, 0.0, 0.5 * dim_z - coil_y) * RotationZ(M_PI)));
  UnionSolid coil15(coil14, connector,
                    Transform3D(Translation3D(0.0, 0.0, 0.5 * dim_z - coil_y) * RotationZ(1.5 * M_PI)));

  Volume coil_vol("coil_vol", coil14, niobium);

  auto       coil_pos = Transform3D(RotationZ(0.5 * M_PI - offset_angle));
  auto       coil_pv  = assembly.placeVolume(coil_vol, coil_pos);
  DetElement coil_de(sdet, "coil_de", 2);
  coil_de.setAttributes(dtor, coil_vol, x_det.regionStr(), x_det.limitsStr(), coil_vis);
  coil_de.setPlacement(coil_pv);

  // -- finishing steps
  auto final_pos = Transform3D(Translation3D(pos_x, pos_y, pos_z) * RotationY(pos_theta));
  auto pv        = dtor.pickMotherVolume(sdet).placeVolume(assembly, final_pos);
  pv.addPhysVolID("system", det_id);
  sdet.setPlacement(pv);

  assembly->GetShape()->ComputeBBox();
  return sdet;
}

DECLARE_DETELEMENT(ip6_CylindricalDipoleMagnet, build_magnet)
