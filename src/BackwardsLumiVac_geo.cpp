
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

  // start and end in z
  double z0 = dd4hep::getAttrOrDefault<double>(x_det, _Unicode(z0), 0);
  double z1 = dd4hep::getAttrOrDefault<double>(x_det, _Unicode(z1), 0);

  // full length in z and center position along z
  double zsiz = z0 - z1;
  double zcen = (z0 + z1) / 2;

  xml::Component points = x_det.child(_Unicode(points));

  double dX0 = points.attr<double>(_Unicode(dX0));
  double dY0 = points.attr<double>(_Unicode(dY0));
  double dX1 = points.attr<double>(_Unicode(dX1));
  double dY1 = points.attr<double>(_Unicode(dY1));

  printout(DEBUG, "BackwardsLumiVac", "dX0: %f", dX0);
  printout(DEBUG, "BackwardsLumiVac", "dY0: %f", dY0);
  printout(DEBUG, "BackwardsLumiVac", "dX1: %f", dX1);
  printout(DEBUG, "BackwardsLumiVac", "dY1: %f", dY1);

  // vertices for the trapezoid
  double ver[16];

  // plane at lower z
  ver[0] = -dX1;
  ver[1] = -dY1;

  ver[2] = -dX1;
  ver[3] = dY1;

  ver[4] = dX1;
  ver[5] = dY1;

  ver[6] = dX1;
  ver[7] = -dY1;

  // plane at higher z
  ver[8] = -dX0;
  ver[9] = -dY0;

  ver[10] = -dX0;
  ver[11] = dY0;

  ver[12] = dX0;
  ver[13] = dY0;

  ver[14] = dX0;
  ver[15] = -dY0;

  EightPointSolid shape(zsiz / 2, ver);

  Volume vol(x_det.nameStr() + "_vol", shape, lccdd.material("Vacuum"));
  vol.setVisAttributes(x_det.visStr());

  DetElement det(x_det.nameStr(), x_det.id());

  Transform3D  pos(RotationZYX(0, 0, 0), Position(0, 0, zcen));
  PlacedVolume pv = lccdd.pickMotherVolume(det).placeVolume(vol, pos);
  det.setPlacement(pv);

  return det;
}

DECLARE_DETELEMENT(BackwardsLumiVac, createDetector)
