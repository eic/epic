//==========================================================================
//  Scintillating fiber calorimeter with tower shape blocks
//  reference: https://github.com/adamjaro/lmon/blob/master/calo/src/WScFiZXv3.cxx
//  Support disk placement
//--------------------------------------------------------------------------
//  Author: Chao Peng (ANL)
//  Date: 07/19/2021
//==========================================================================

#include "GeometryHelpers.h"
#include "DD4hep/DetFactoryHelper.h"
#include <XML/Helper.h>
#include <iostream>
#include <algorithm>
#include <tuple>
#include <math.h>

using namespace dd4hep;
using Point =  ROOT::Math::XYPoint;

std::tuple<Volume, Position> build_module(const Detector &desc, const xml::Component &mod_x, SensitiveDetector &sens);

// helper function to get x, y, z if defined in a xml component
template<class XmlComp>
Position get_xml_xyz(const XmlComp &comp, dd4hep::xml::Strng_t name)
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

// main
static Ref_t create_detector(Detector& desc, xml::Handle_t handle, SensitiveDetector sens)
{
    xml::DetElement detElem = handle;
    std::string detName = detElem.nameStr();
    int detID = detElem.id();
    DetElement det(detName, detID);
    sens.setType("calorimeter");
    auto dim = detElem.dimensions();
    auto rmin = dim.rmin();
    auto rmax = dim.rmax();
    auto length = dim.length();
    auto phimin = dd4hep::getAttrOrDefault<double>(dim, _Unicode(phimin), 0.);
    auto phimax = dd4hep::getAttrOrDefault<double>(dim, _Unicode(phimax), 2.*M_PI);
    // envelope
    Tube envShape(rmin, rmax, length/2., phimin, phimax);
    Volume env(detName + "_envelope", envShape, desc.material("Air"));
    env.setVisAttributes(desc.visAttributes(detElem.visStr()));

    // build module
    auto [modVol, modSize] = build_module(desc, detElem.child(_Unicode(module)), sens);
    double modSizeR = std::sqrt(modSize.x() * modSize.x() + modSize.y() * modSize.y());
    double assembly_rwidth = modSizeR*2.;
    int nas = int((rmax - rmin) / assembly_rwidth) + 1;
    std::vector<Assembly> assemblies;
    // calorimeter block z-offsets (as blocks are shorter than the volume length)
    const double block_offset = -0.5*(length - modSize.z());
    for (int i = 0; i < nas; ++i) {
        Assembly assembly(detName + Form("_ring%d", i + 1));
        auto assemblyPV = env.placeVolume(assembly, Position{0., 0., block_offset});
        assemblyPV.addPhysVolID("ring", i + 1);
        assemblies.emplace_back(std::move(assembly));
    }
    // std::cout << assemblies.size() << std::endl;

    int modid = 1;
    for (int ix = 0; ix < int(2.*rmax / modSize.x()) + 1; ++ix) {
        double mx = modSize.x() * ix - rmax;
        for (int iy = 0; iy < int(2.*rmax / modSize.y()) + 1; ++iy) {
            double my = modSize.y() * iy - rmax;
            double mr = std::sqrt(mx*mx + my*my);
            if (mr - modSizeR >= rmin && mr + modSizeR <= rmax) {
                int ias = int((mr - rmin) / assembly_rwidth);
                auto &assembly = assemblies[ias];
                auto modPV = assembly.placeVolume(modVol, Position(mx, my, 0.));
                modPV.addPhysVolID("module", modid++);
            }
        }
    }

    desc.add(Constant(detName + "_NModules", std::to_string(modid - 1)));

    for (auto &assembly : assemblies) {
        assembly.ptr()->Voxelize("");
    }

    // detector position and rotation
    auto pos = get_xml_xyz(detElem, _Unicode(position));
    auto rot = get_xml_xyz(detElem, _Unicode(rotation));
    Volume motherVol = desc.pickMotherVolume(det);
    Transform3D tr = Translation3D(pos.x(), pos.y(), pos.z()) * RotationZYX(rot.z(), rot.y(), rot.x());
    PlacedVolume envPV = motherVol.placeVolume(env, tr);
    envPV.addPhysVolID("system", detID);
    det.setPlacement(envPV);
    return det;
}

// helper function to build module with scintillating fibers
std::tuple<Volume, Position> build_module(const Detector &desc, const xml::Component &mod_x, SensitiveDetector &sens)
{
    auto sx = mod_x.attr<double>(_Unicode(sizex));
    auto sy = mod_x.attr<double>(_Unicode(sizey));
    auto sz = mod_x.attr<double>(_Unicode(sizez));

    Box modShape(sx/2., sy/2., sz/2.);
    auto modMat = desc.material(mod_x.attr<std::string>(_Unicode(material)));
    Volume modVol("module_vol", modShape, modMat);
    if (mod_x.hasAttr(_Unicode(vis))) {
        modVol.setVisAttributes(desc.visAttributes(mod_x.attr<std::string>(_Unicode(vis))));
    }

    if (mod_x.hasChild("fiber")) {
      auto fiber_x  = mod_x.child(_Unicode(fiber));
      auto fr       = fiber_x.attr<double>(_Unicode(radius));
      auto fsx      = fiber_x.attr<double>(_Unicode(spacex));
      auto fsy      = fiber_x.attr<double>(_Unicode(spacey));
      auto foff     = dd4hep::getAttrOrDefault<double>(fiber_x, _Unicode(offset), 0.5*mm);
      auto fiberMat = desc.material(fiber_x.attr<std::string>(_Unicode(material)));
      Tube fiberShape(0., fr, sz/2.);
      Volume fiberVol("fiber_vol", fiberShape, fiberMat);
      fiberVol.setSensitiveDetector(sens);

      // Fibers are placed in a honeycomb with the radius = sqrt(3)/2. * hexagon side length
      // So each fiber is fully contained in a regular hexagon, which are placed as
      //           ______________________________________
      //           |          ____        ____          |
      // even:     |         /    \      /    \         |
      //           |    ____/      \____/      \____    |
      //           |   /    \      /    \      /    \   |
      // odd:      |  /      \____/      \____/      \  |
      //           |  \      /    \      /    \      /  |
      //           |   \____/      \____/      \____/   |
      // even:     |        \      /    \      /        |
      //           |         \____/      \____/      ___|___
      //           |____________________________________|___offset
      //                                              | |
      //                                              |offset
      // the parameters space x and space y are used to add additional spaces between the hexagons
      double fside  = 2. / std::sqrt(3.) * fr;
      double fdistx = 2. * fside + fsx;
      double fdisty = 2. * fr + fsy;

      // maximum numbers of the fibers, help narrow the loop range
      int nx = int(sx / (2.*fr)) + 1;
      int ny = int(sy / (2.*fr)) + 1;

      // std::cout << sx << ", " << sy << ", " << fr << ", " << nx << ", " << ny << std::endl;

      // place the fibers
      double y0 = (foff + fside);
      int nfibers = 0;
      for (int iy = 0; iy < ny; ++iy) {
          double y = y0 + fdisty * iy;
          // about to touch the boundary
          if ((sy - y) < y0) { break; }
          double x0 = (iy % 2) ? (foff + fside) : (foff + fside + fdistx / 2.);
          for (int ix = 0; ix < nx; ++ix) {
              double x = x0 + fdistx * ix;
              // about to touch the boundary
              if ((sx - x) < x0) { break; }
              auto fiberPV = modVol.placeVolume(fiberVol, nfibers++, Position{x - sx/2., y - sy/2., 0});
              //std::cout << "(" << ix << ", " << iy << ", " << x - sx/2. << ", " << y - sy/2. << ", " << fr << "),\n";
              fiberPV.addPhysVolID("fiber_x", ix + 1).addPhysVolID("fiber_y", iy + 1);
          }
      }
    // if no fibers we make the module itself sensitive
    } else {
      modVol.setSensitiveDetector(sens);
    }


    return std::make_tuple(modVol, Position{sx, sy, sz});
}

DECLARE_DETELEMENT(ScFiCalorimeter, create_detector)

