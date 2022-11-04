#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/OpticalSurfaces.h"
#include "DD4hep/Printout.h"
#include "DDRec/DetectorData.h"
#include "DDRec/Surface.h"
#include <XML/Helper.h>
#include <algorithm>
#include <iostream>
#include <math.h>
#include <tuple>
//////////////////////////////////////////////////
// Far Forward Ion Zero Degree Calorimeter - Ecal
// Reference from ATHENA ScFiCalorimeter_geo.cpp
//////////////////////////////////////////////////

using namespace std;
using namespace dd4hep;

// main
static Ref_t create_detector(Detector& desc, xml_h handle, SensitiveDetector sens)
{
  xml::DetElement detElem = handle;
  std::string     detName = detElem.nameStr();
  int             detID   = detElem.id();
  DetElement      det(detName, detID);
  sens.setType("calorimeter");
  auto      dim    = detElem.dimensions();
  auto      xwidth = dim.x();
  auto      ywidth = dim.y();
  auto      length = dim.z();
  xml_dim_t pos    = detElem.position();
  xml_dim_t rot    = detElem.rotation();

  // envelope
  Box    envShape(xwidth * 0.5, ywidth * 0.5, length * 0.5);
  Volume env(detName + "_envelope", envShape, desc.material("Air"));
  env.setVisAttributes(desc.visAttributes(detElem.visStr()));

  int    layerid = 0;
  double zpos_0  = -length / 2.;

  for (xml_coll_t li(detElem, _Unicode(layer)); li; ++li) {

    xml_comp_t x_lyr = li;
    auto       nlyr  = x_lyr.attr<int>(_Unicode(nlayer));
    auto       gap_z = x_lyr.attr<double>(_Unicode(gapspace));

    map<int, string>    v_sl_name;
    map<string, Volume> slices;
    map<string, double> sl_thickness;

    int        nsl = 0;
    xml_coll_t ci(x_lyr, _Unicode(slice));
    for (ci.reset(); ci; ++ci) {
      xml_comp_t x_sl    = ci;
      Material   sl_mat  = desc.material(x_sl.materialStr());
      string     sl_name = x_sl.nameStr();
      double     sl_z    = x_sl.thickness();

      Box    sl_Shape(xwidth / 2., ywidth / 2., sl_z / 2.);
      Volume sl_Vol("slice_vol", sl_Shape, sl_mat);
      sl_Vol.setVisAttributes(desc.visAttributes(x_sl.visStr()));
      if (x_sl.isSensitive())
        sl_Vol.setSensitiveDetector(sens);

      nsl++;
      v_sl_name[nsl]        = sl_name;
      slices[sl_name]       = sl_Vol;
      sl_thickness[sl_name] = sl_z;
    }

    for (int ilyr = 0; ilyr < nlyr; ilyr++) {
      layerid++;
      for (int isl = 0; isl < nsl; isl++) {
        string sl_name = v_sl_name[isl + 1];

        double       zpos = zpos_0 + sl_thickness[sl_name] / 2.;
        Position     sl_pos(0, 0, zpos);
        PlacedVolume pv = env.placeVolume(slices[sl_name], sl_pos);
        if (slices[sl_name].isSensitive())
          pv.addPhysVolID(sl_name, layerid);

        zpos_0 += sl_thickness[sl_name];
      }

      zpos_0 += gap_z;
    }
  }

  // detector position and rotation
  Volume       motherVol = desc.pickMotherVolume(det);
  Transform3D  tr(RotationZYX(rot.z(), -rot.y(), rot.x()), Position(pos.x(), pos.y(), pos.z()));
  PlacedVolume envPV = motherVol.placeVolume(env, tr);
  envPV.addPhysVolID("system", detID);
  det.setPlacement(envPV);
  return det;
}

DECLARE_DETELEMENT(ZDC_ImagingCal, create_detector)
