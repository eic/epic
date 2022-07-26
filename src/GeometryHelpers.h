#pragma once
#include "Math/Point2D.h"
#include <vector>

// some utility functions that can be shared
namespace epic::geo {

  using Point = ROOT::Math::XYPoint;

  /** Fill rectangles in a ring (disk).
   *
   * @param ref  2D reference point.
   * @param sx   x side length
   * @param sy   y side length
   * @param rmin inner radius of disk
   * @param rmax  outer radius of disk to fill
   * @param phmin  phi min
   * @param phmax phi max
   */
  std::vector<Point> fillRectangles(Point ref, double sx, double sy, double rmin, double rmax, double phmin = -M_PI,
                                    double phmax = M_PI);
  // fill squares in a ring
  inline std::vector<Point> fillSquares(Point ref, double size, double rmin, double rmax, double phmin = -M_PI,
                                        double phmax = M_PI)
  {
    return fillRectangles(ref, size, size, rmin, rmax, phmin, phmax);
  }

  std::vector<Point> fillHexagons(Point ref, double lside, double rmin, double rmax, double phmin = -M_PI,
                                  double phmax = M_PI);

} // namespace epic::geo
