
#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/Printout.h"
#include "TMath.h"
#include <XML/Helper.h>

using namespace std;
using namespace dd4hep;

static Ref_t createDetector(Detector& lccdd, xml_h e, SensitiveDetector /*sens*/)
{

  xml_det_t x_det = e;

  // EightPointSolid -> TGeoArb8 -> G4GenericTrap
  // generic trapezoid native coordinates: 4 xy points plane at -dz, 4 xy points plane at +dz, both clockwise
  // rotation by +pi/2 about x from generic trapezoid coordinates to detector frame: y -> z,  z -> y

  // double ysiz = 200*mm; // full size in y

  // full size in y
  double ysiz = dd4hep::getAttrOrDefault<double>(x_det, _Unicode(ysiz), 0);

  xml::Component points = x_det.child(_Unicode(points));

  double zQT = points.attr<double>(_Unicode(zQT));
  double xQT = points.attr<double>(_Unicode(xQT));
  double zQB = points.attr<double>(_Unicode(zQB));
  double xQB = points.attr<double>(_Unicode(xQB));

  printout(DEBUG, "BackwardsVacDrift", "zQT: %f", zQT);
  printout(DEBUG, "BackwardsVacDrift", "xQT: %f", xQT);
  printout(DEBUG, "BackwardsVacDrift", "zQB: %f", zQB);
  printout(DEBUG, "BackwardsVacDrift", "xQB: %f", xQB);

  double zW = points.attr<double>(_Unicode(zW));
  double xW = points.attr<double>(_Unicode(xW));
  double xA = points.attr<double>(_Unicode(xA));

  printout(DEBUG, "BackwardsVacDrift", "zW: %f", zW);
  printout(DEBUG, "BackwardsVacDrift", "xW: %f", xW);
  printout(DEBUG, "BackwardsVacDrift", "xA: %f", xA);

  // vertices for the trapezoid
  double ver[16];

  // plane at upper y
  ver[0] = xQB;
  ver[1] = zQB;

  ver[2] = xQT;
  ver[3] = zQT;

  ver[4] = xW;
  ver[5] = zW;

  ver[6] = xA;
  ver[7] = zW;

  // plane at lower y
  for (int i = 8; i < 16; i++) {
    ver[i] = ver[i - 8];
  }

  EightPointSolid shape(ysiz / 2, ver);

  Volume vol(x_det.nameStr() + "_vol", shape, lccdd.material("Vacuum"));
  vol.setVisAttributes(x_det.visStr());

  DetElement det(x_det.nameStr(), x_det.id());

  Transform3D  pos(RotationZYX(0, 0, M_PI / 2), Position(0, 0, 0));
  PlacedVolume pv = lccdd.pickMotherVolume(det).placeVolume(vol, pos);
  det.setPlacement(pv);

  return det;
}

DECLARE_DETELEMENT(BackwardsVacDrift, createDetector)
