//==========================================================================
//  Implementation for shashlik calorimeter modules
//  it supports disk placements with (rmin, rmax), and (phimin, phimax)
//--------------------------------------------------------------------------
//  Author: Chao Peng (ANL)
//  Date: 06/22/2021
//==========================================================================

#include "DD4hep/DetFactoryHelper.h"
#include "GeometryHelpers.h"
#include <XML/Helper.h>
#include <XML/Layering.h>
#include <algorithm>
#include <iostream>
#include <math.h>
#include <tuple>

using namespace dd4hep;

static void add_disk_shashlik(Detector& desc, Assembly& env, xml::Collection_t& plm, SensitiveDetector& sens, int id);

// helper function to get x, y, z if defined in a xml component
template <class XmlComp>
Position get_xml_xyz(XmlComp& comp, dd4hep::xml::Strng_t name)
{
  Position pos(0., 0., 0.);
  if (comp.hasChild(name)) {
    auto child = comp.child(name);
    pos.SetX(dd4hep::getAttrOrDefault<double>(child, _Unicode(x), 0.));
    pos.SetY(dd4hep::getAttrOrDefault<double>(child, _Unicode(y), 0.));
    pos.SetZ(dd4hep::getAttrOrDefault<double>(child, _Unicode(z), 0.));
  }
  return pos;
}

static Ref_t create_detector(Detector& desc, xml::Handle_t handle, SensitiveDetector sens)
{
  static const std::string func    = "ShashlikCalorimeter";
  xml::DetElement          detElem = handle;
  std::string              detName = detElem.nameStr();
  int                      detID   = detElem.id();
  DetElement               det(detName, detID);
  sens.setType("calorimeter");
  // envelope
  Assembly assembly(detName);

  // module placement
  xml::Component plm    = detElem.child(_Unicode(placements));
  int            sector = 1;
  for (xml::Collection_t mod(plm, _Unicode(disk)); mod; ++mod) {
    add_disk_shashlik(desc, assembly, mod, sens, sector++);
  }

  // detector position and rotation
  auto         pos       = get_xml_xyz(detElem, _Unicode(position));
  auto         rot       = get_xml_xyz(detElem, _Unicode(rotation));
  Volume       motherVol = desc.pickMotherVolume(det);
  Transform3D  tr        = Translation3D(pos.x(), pos.y(), pos.z()) * RotationZYX(rot.z(), rot.y(), rot.x());
  PlacedVolume envPV     = motherVol.placeVolume(assembly, tr);
  envPV.addPhysVolID("system", detID);
  det.setPlacement(envPV);
  return det;
}

// helper function to build module with or w/o wrapper
std::tuple<Volume, int, double, double> build_shashlik(Detector& desc, xml::Collection_t& plm, SensitiveDetector& sens)
{
  auto mod = plm.child(_Unicode(module));
  // a modular volume
  std::string shape = dd4hep::getAttrOrDefault<std::string>(mod, _Unicode(shape), "square");
  std::transform(shape.begin(), shape.end(), shape.begin(), [](char c) { return std::tolower(c); });
  int nsides = 4;
  if (shape == "hexagon") {
    nsides = 6;
  } else if (shape != "square") {
    std::cerr << "ShashlikCalorimeter Error: Unsupported shape of module " << shape
              << ". Please choose from (square, hexagon). Proceed with square shape." << std::endl;
  }
  double   slen = mod.attr<double>(_Unicode(side_length));
  double   rmax = slen / 2. / std::sin(M_PI / nsides);
  Layering layering(mod);
  auto     len = layering.totalThickness();

  // wrapper info
  PolyhedraRegular mpoly(nsides, 0., rmax, len);
  Volume           mvol("shashlik_module_vol", mpoly, desc.air());
  mvol.setVisAttributes(desc.visAttributes(dd4hep::getAttrOrDefault<std::string>(mod, _Unicode(vis), "GreenVis")));

  double gap = 0.;
  Volume wvol("shashlik_wrapper_vol");
  if (plm.hasChild(_Unicode(wrapper))) {
    auto wrap = plm.child(_Unicode(wrapper));
    gap       = wrap.attr<double>(_Unicode(thickness));
    if (gap > 1e-6 * mm) {
      wvol.setSolid(PolyhedraRegular(nsides, 0., rmax + gap, len));
      wvol.setMaterial(desc.material(wrap.attr<std::string>(_Unicode(material))));
      wvol.setVisAttributes(desc.visAttributes(dd4hep::getAttrOrDefault<std::string>(wrap, _Unicode(vis), "WhiteVis")));
      wvol.placeVolume(mvol, Position{0., 0., 0.});
    }
  }

  // layer start point
  double lz   = -len / 2.;
  int    lnum = 1;
  // Loop over the sets of layer elements in the detector.
  for (xml_coll_t li(mod, _U(layer)); li; ++li) {
    int repeat = li.attr<int>(_Unicode(repeat));
    // Loop over number of repeats for this layer.
    for (int j = 0; j < repeat; j++) {
      std::string      lname  = Form("layer%d", lnum);
      double           lthick = layering.layer(lnum - 1)->thickness(); // Layer's thickness.
      PolyhedraRegular lpoly(nsides, 0., rmax, lthick);
      Volume           lvol(lname, lpoly, desc.air());

      // Loop over the sublayers or slices for this layer.
      int    snum = 1;
      double sz   = -lthick / 2.;
      for (xml_coll_t si(li, _U(slice)); si; ++si) {
        std::string      sname  = Form("slice%d", snum);
        double           sthick = si.attr<double>(_Unicode(thickness));
        PolyhedraRegular spoly(nsides, 0., rmax, sthick);
        Volume           svol(sname, spoly, desc.material(si.attr<std::string>(_Unicode(material))));

        std::string issens = dd4hep::getAttrOrDefault<std::string>(si, _Unicode(sensitive), "no");
        std::transform(issens.begin(), issens.end(), issens.begin(), [](char c) { return std::tolower(c); });
        if ((issens == "yes") || (issens == "y") || (issens == "true")) {
          svol.setSensitiveDetector(sens);
        }
        svol.setAttributes(desc, dd4hep::getAttrOrDefault<std::string>(si, _Unicode(region), ""),
                           dd4hep::getAttrOrDefault<std::string>(si, _Unicode(limits), ""),
                           dd4hep::getAttrOrDefault<std::string>(si, _Unicode(vis), "InvisibleNoDaughters"));

        // Slice placement.
        auto slicePV = lvol.placeVolume(svol, Position(0, 0, sz + sthick / 2.));
        slicePV.addPhysVolID("slice", snum++);
        // Increment Z position of slice.
        sz += sthick;
      }

      // Set region, limitset, and vis of layer.
      lvol.setAttributes(desc, dd4hep::getAttrOrDefault<std::string>(li, _Unicode(region), ""),
                         dd4hep::getAttrOrDefault<std::string>(li, _Unicode(limits), ""),
                         dd4hep::getAttrOrDefault<std::string>(li, _Unicode(vis), "InvisibleNoDaughters"));

      auto layerPV = mvol.placeVolume(lvol, Position(0, 0, lz + lthick / 2));
      layerPV.addPhysVolID("layer", lnum++);
      // Increment to next layer Z position.
      lz += lthick;
    }
  }

  if (gap > 1e-6 * mm) {
    return std::make_tuple(wvol, nsides, 2. * std::sin(M_PI / nsides) * (rmax + gap), len);
  } else {
    return std::make_tuple(mvol, nsides, slen, len);
  }
}

// place disk of modules
static void add_disk_shashlik(Detector& desc, Assembly& env, xml::Collection_t& plm, SensitiveDetector& sens, int sid)
{
  auto [mvol, nsides, sidelen, len] = build_shashlik(desc, plm, sens);
  int    sector_id                  = dd4hep::getAttrOrDefault<int>(plm, _Unicode(sector), sid);
  int    id_begin                   = dd4hep::getAttrOrDefault<int>(plm, _Unicode(id_begin), 1);
  double rmin                       = plm.attr<double>(_Unicode(rmin));
  double rmax                       = plm.attr<double>(_Unicode(rmax));
  double phimin                     = dd4hep::getAttrOrDefault<double>(plm, _Unicode(phimin), 0.);
  double phimax                     = dd4hep::getAttrOrDefault<double>(plm, _Unicode(phimax), 2. * M_PI);

  auto points = (nsides == 6) ? ecce::geo::fillHexagons({0., 0.}, sidelen, rmin, rmax, phimin, phimax)
                              : ecce::geo::fillSquares({0., 0.}, sidelen * 1.414, rmin, rmax, phimin, phimax);
  // placement to mother
  auto pos = get_xml_xyz(plm, _Unicode(position));
  auto rot = get_xml_xyz(plm, _Unicode(rotation));
  int  mid = 0;
  for (auto& p : points) {
    Transform3D tr = RotationZYX(rot.z(), rot.y(), rot.x()) *
                     Translation3D(pos.x() + p.x(), pos.y() + p.y(), pos.z() + len / 2.) *
                     RotationZ((nsides == 4) ? 45 * degree : 0);
    auto modPV = env.placeVolume(mvol, tr);
    modPV.addPhysVolID("sector", sector_id).addPhysVolID("module", id_begin + mid++);
  }
}

DECLARE_DETELEMENT(ecce_ShashlikCalorimeter, create_detector)
