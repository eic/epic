
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

  double zB = points.attr<double>(_Unicode(zB));
  double xB = points.attr<double>(_Unicode(xB));

  double zT  = points.attr<double>(_Unicode(zT));
  double xT  = points.attr<double>(_Unicode(xT));
  double zTB = points.attr<double>(_Unicode(zTB));
  double xTB = points.attr<double>(_Unicode(xTB));

  printout(DEBUG, "BackwardsTagWin", "zB: %f", zB);
  printout(DEBUG, "BackwardsTagWin", "xB: %f", xB);
  printout(DEBUG, "BackwardsTagWin", "zT: %f", zT);
  printout(DEBUG, "BackwardsTagWin", "xT: %f", xT);
  printout(DEBUG, "BackwardsTagWin", "zTB: %f", zTB);
  printout(DEBUG, "BackwardsTagWin", "xTB: %f", xTB);

  // vertices for the trapezoid
  double ver[16];

  // plane at upper y
  ver[0] = xTB;
  ver[1] = zTB;

  ver[2] = xB;
  ver[3] = zB;

  ver[4] = xB;
  ver[5] = zB;

  ver[6] = xT;
  ver[7] = zT;

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

DECLARE_DETELEMENT(BackwardsTagWin, createDetector)
