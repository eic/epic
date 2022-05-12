
#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/Printout.h"
#include <XML/Helper.h>
#include "TMath.h"

using namespace std;
using namespace dd4hep;

static Ref_t createDetector(Detector& lccdd, xml_h e, SensitiveDetector /*sens*/) {

  xml_det_t x_det = e;

  //EightPointSolid -> TGeoArb8 -> G4GenericTrap
  //generic trapezoid native coordinates: 4 xy points plane at -dz, 4 xy points plane at +dz, both clockwise
  //rotation by +pi/2 about x from generic trapezoid coordinates to detector frame: y -> z,  z -> y

  //double exit_win_x_beam = -20*mm; // x of beam pipe at the front of exit window
  //double exit_win_z = -100*mm;

  xml::Component points = x_det.child(_Unicode(points));

  //front of exit window volume
  //double win_xmin = -250*mm; // given by beam pipe to tagger

  // given by beam pipe to tagger
  double win_xmin = points.attr<double>( _Unicode(win_xmin) );
  printout(DEBUG, "BackwardsVacBox", "win_xmin: %f", win_xmin);

  //double zQT = points.attr<double>( _Unicode(zQT) );
  //double xQT = points.attr<double>( _Unicode(xQT) );
  //double zQB = points.attr<double>( _Unicode(zQB) );
  //double xQB = points.attr<double>( _Unicode(xQB) );

  //printout(DEBUG, "BackwardsVacBox", "zQT: %f", zQT);
  //printout(DEBUG, "BackwardsVacBox", "xQT: %f", xQT);
  //printout(DEBUG, "BackwardsVacBox", "zQB: %f", zQB);
  //printout(DEBUG, "BackwardsVacBox", "xQB: %f", xQB);

  //double zB = points.attr<double>( _Unicode(zB) );
  //double xB = points.attr<double>( _Unicode(xB) );

  //printout(DEBUG, "BackwardsVacBox", "zB: %f", zB);
  //printout(DEBUG, "BackwardsVacBox", "xB: %f", xB);

  //auxiliary point from direct calculation
  //double xA = xQB + (zB-win_z)*(xB-xQB)/(zB-zQB);
  //double xA = xQB + (win_z-zQB)*(xB-xQB)/(zB-zQB);
  //win_xmin = xA;
  //printout(DEBUG, "BackwardsVacBox", "xA: %f", xA);

  double win_z = -18500*mm; // z of front of exit window volume
  double win_xmax = 100*mm; // given by angular apperture to exit window
  //double win_ysiz = 200*mm; // full size in y at the front of exit window volume

  // full size in y at the front of exit window volume
  double win_ysiz = dd4hep::getAttrOrDefault<double>(x_det, _Unicode(ysiz), 0);

  //end of B2BeR magnet
  double b2b_end_xy = 196*mm; // full size in x and y at the end of B2eR, equal to inner diameter
  double b2b_end_z = -14865*mm; // z position of B2BeR end

  //derived quantities
  double zsiz = TMath::Abs(win_z - b2b_end_z); //vacuum element full size in z
  double zpos = (win_z + b2b_end_z)/2; // z position of the vacuum element

  //vertices for the trapezoid
  double ver[16];

  //plane at the front of exit window
  ver[0] = win_xmin;
  ver[1] = -win_ysiz/2;

  ver[2] = win_xmin;
  ver[3] = win_ysiz/2;

  ver[4] = win_xmax;
  ver[5] = win_ysiz/2;

  ver[6] = win_xmax;
  ver[7] = -win_ysiz/2;

  //plane at the end of B2BeR
  ver[8] = -b2b_end_xy/2;
  //ver[9] = -b2b_end_xy/2;
  ver[9] = -win_ysiz/2;

  ver[10] = -b2b_end_xy/2;
  //ver[11] = b2b_end_xy/2;
  ver[11] = win_ysiz/2;

  ver[12] = b2b_end_xy/2;
  //ver[13] = b2b_end_xy/2;
  ver[13] = win_ysiz/2;

  ver[14] = b2b_end_xy/2;
  //ver[15] = -b2b_end_xy/2;
  ver[15] = -win_ysiz/2;

  //double dz = 25*mm;
  //double vertices[] = {-30*mm, -30*mm, -30*mm, 30*mm, 30*mm, 30*mm, 30*mm, -30*mm, -5*mm, -20*mm, -20*mm, 20*mm, 20*mm, 20*mm, 20*mm, -20*mm};

  //for(int i=1; i<16; i+=2) {
  //  vertices[i] -= 40*mm;
  //}

  //EightPointSolid shape(dz, vertices);
  EightPointSolid shape(zsiz/2, ver);

  Volume vol(x_det.nameStr() + "_vol", shape, lccdd.material("Vacuum"));
  vol.setVisAttributes(x_det.visStr());

  DetElement det(x_det.nameStr(), x_det.id());

  //Transform3D pos(RotationZYX(0, 0, M_PI/2), Position(0, 0, 0));
  Transform3D pos(RotationZYX(0, 0, 0), Position(0, 0, zpos));
  PlacedVolume pv = lccdd.pickMotherVolume(det).placeVolume(vol, pos);
  det.setPlacement(pv);

  return det;

}

DECLARE_DETELEMENT(BackwardsVacBox, createDetector)


