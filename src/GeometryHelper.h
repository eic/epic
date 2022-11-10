#pragma once
#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/OpticalSurfaces.h"
#include "DD4hep/Printout.h"
#include "DDRec/DetectorData.h"
#include "DDRec/Surface.h"
#include "Math/Point2D.h"
#include <XML/Helper.h>
#include <functional>
#include <vector>

// some utility functions that can be shared
namespace ip6::geo {

  using std::function;
  using std::make_tuple;
  using std::map;
  using std::string;
  using std::tuple;
  using std::vector;
  using Volume            = dd4hep::Volume;
  using Position          = dd4hep::Position;
  using Assembly          = dd4hep::Assembly;
  using Detector          = dd4hep::Detector;
  using SensitiveDetector = dd4hep::SensitiveDetector;
  using Point             = ROOT::Math::XYPoint;
  using Transform3D       = dd4hep::Transform3D;
  using Translation3D     = dd4hep::Translation3D;
  using RotationZYX       = dd4hep::RotationZYX;
  using UnionSolid        = dd4hep::UnionSolid;
  using Material          = dd4hep::Material;
  using Tube              = dd4hep::Tube;

  /* Fill rectangles in a pacman disk
   *
   * @param ref  2D reference point.
   * @param sx   x side length
   * @param sy   y side length
   * @param rmin inner radius of disk
   * @param rintermediate outer radius of inner full azimuth disk
   * @param rmax  outer radius of outer disk section ranging from phi min to phi max
   * @param phmin  phi min
   * @param phmax phi max
   * N.B. Current implementation only guarantees that inner radial bound will be strictly observed.
   * Outer radial bunds can be penetrated by internal modules on or near the boundary.
   */

  vector<Point> fillRectangles(Point ref, double sx, double sy, double rmin, double rintermediate, double rmax,
                               double phmin = -M_PI, double phmax = M_PI);

  tuple<int, int>
  add_individuals(function<tuple<Volume, Position>(Detector&, xml_coll_t&, SensitiveDetector&)> build_module,
                  Detector& desc, Assembly& env, xml_coll_t& plm, SensitiveDetector& sens, int id);

  tuple<int, int> add_disk(function<tuple<Volume, Position>(Detector&, xml_coll_t&, SensitiveDetector&)> build_module,
                           Detector& desc, Assembly& env, xml_coll_t& plm, SensitiveDetector& sens, int id);

  Position get_xml_xyz(xml_comp_t& comp, dd4hep::xml::Strng_t name);

} // namespace ip6::geo
