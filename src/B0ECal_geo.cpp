#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/OpticalSurfaces.h"
#include "DD4hep/Printout.h"
#include "DDRec/DetectorData.h"
#include "DDRec/Surface.h"
#include "GeometryHelper.h"
#include "Math/Point2D.h"
#include <XML/Helper.h>
#include <vector>

//////////////////////////////////////////////////
// Far Forward B0 Electromagnetic Calorimeter
//////////////////////////////////////////////////

using std::make_tuple;
using std::map;
using std::string;
using std::tuple;
using std::vector;
using namespace dd4hep;

static tuple<Volume, Position> build_module(Detector& desc, xml_coll_t& plm, SensitiveDetector& sens);

static Ref_t createDetector(Detector& desc, xml_h e, SensitiveDetector sens)
{
  xml_det_t  x_det   = e;
  string     detName = x_det.nameStr();
  int        detID   = x_det.id();
  DetElement det(detName, detID);
  sens.setType("calorimeter");

  // assembly
  Assembly detVol(detName);

  xml_dim_t pos = x_det.position();
  xml_dim_t rot = x_det.rotation();

  // module placement
  xml_comp_t    plm = x_det.child(_Unicode(placements));
  map<int, int> sectorModuleNumbers;
  auto          addModuleNumbers = [&sectorModuleNumbers](int sector, int nmod) {
    auto it = sectorModuleNumbers.find(sector);
    if (it != sectorModuleNumbers.end()) {
      it->second += nmod;
    } else {
      sectorModuleNumbers[sector] = nmod;
    }
  };

  int sector_id = 1;

  for (xml_coll_t mod(plm, _Unicode(individuals)); mod; ++mod) {
    auto [sector, nmod] = ip6::geo::add_individuals(build_module, desc, detVol, mod, sens, sector_id++);
    addModuleNumbers(sector, nmod);
  }

  for (xml_coll_t disk(plm, _Unicode(disk)); disk; ++disk) {
    auto [sector, nmod] = ip6::geo::add_disk(build_module, desc, detVol, disk, sens, sector_id++);
    addModuleNumbers(sector, nmod);
  }

  // position and rotation of parent volume
  Volume       motherVol = desc.pickMotherVolume(det);
  Transform3D  tr(RotationZYX(rot.z(), rot.y(), rot.x()), Position(pos.x(), pos.y(), pos.z()));
  PlacedVolume detPV = motherVol.placeVolume(detVol, tr);
  detPV.addPhysVolID("system", detID);
  det.setPlacement(detPV);
  return det;
}

// helper function to build module with or w/o wrapper
static tuple<Volume, Position> build_module(Detector& desc, xml::Collection_t& plm, SensitiveDetector& sens)
{
  auto   mod = plm.child(_Unicode(module));
  auto   sx  = mod.attr<double>(_Unicode(sizex));
  auto   sy  = mod.attr<double>(_Unicode(sizey));
  auto   sz  = mod.attr<double>(_Unicode(sizez));
  Box    modShape(sx / 2., sy / 2., sz / 2.);
  auto   modMat = desc.material(mod.attr<string>(_Unicode(material)));
  Volume modVol("module_vol", modShape, modMat);
  modVol.setSensitiveDetector(sens);
  modVol.setVisAttributes(desc.visAttributes(mod.attr<string>(_Unicode(vis))));

  // no wrapper
  if (!plm.hasChild(_Unicode(wrapper))) {
    return make_tuple(modVol, Position{sx, sy, sz});
    // build wrapper
  } else {
    auto wrp       = plm.child(_Unicode(wrapper));
    auto thickness = wrp.attr<double>(_Unicode(thickness));
    if (thickness < 1e-12 * mm) {
      return make_tuple(modVol, Position{sx, sy, sz});
    }
    auto   wrpMat = desc.material(wrp.attr<string>(_Unicode(material)));
    Box    wrpShape((sx + thickness) / 2., (sy + thickness) / 2., sz / 2.);
    Volume wrpVol("wrapper_vol", wrpShape, wrpMat);
    wrpVol.placeVolume(modVol, Position(0., 0., 0.));
    wrpVol.setVisAttributes(desc.visAttributes(wrp.attr<string>(_Unicode(vis))));
    return make_tuple(wrpVol, Position{sx + thickness, sy + thickness, sz});
  }
}

DECLARE_DETELEMENT(B0_ECAL, createDetector)
