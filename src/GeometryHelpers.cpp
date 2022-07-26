#include "GeometryHelpers.h"

// some utility functions that can be shared
namespace epic::geo {

  typedef ROOT::Math::XYPoint Point;

  // check if a 2d point is already in the container
  bool already_placed(const Point& p, const std::vector<Point>& vec, double xs = 1.0, double ys = 1.0,
                      double tol = 1e-6)
  {
    for (auto& pt : vec) {
      if ((std::abs(pt.x() - p.x()) / xs < tol) && std::abs(pt.y() - p.y()) / ys < tol) {
        return true;
      }
    }
    return false;
  }

  // check if a square in a ring
  inline bool rec_in_ring(const Point& pt, double sx, double sy, double rmin, double rmax, double phmin, double phmax)
  {
    if (pt.r() > rmax || pt.r() < rmin) {
      return false;
    }

    // check four corners
    std::vector<Point> pts{
        Point(pt.x() - sx / 2., pt.y() - sy / 2.),
        Point(pt.x() - sx / 2., pt.y() + sy / 2.),
        Point(pt.x() + sx / 2., pt.y() - sy / 2.),
        Point(pt.x() + sx / 2., pt.y() + sy / 2.),
    };
    for (auto& p : pts) {
      if (p.r() > rmax || p.r() < rmin || p.phi() > phmax || p.phi() < phmin) {
        return false;
      }
    }
    return true;
  }

  // a helper function to recursively fill square in a ring
  void add_rectangle(Point p, std::vector<Point>& res, double sx, double sy, double rmin, double rmax, double phmin,
                     double phmax, int max_depth = 20, int depth = 0)
  {
    // std::cout << depth << "/" << max_depth << std::endl;
    // exceeds the maximum depth in searching or already placed
    if ((depth > max_depth) || (already_placed(p, res, sx, sy))) {
      return;
    }

    bool in_ring = rec_in_ring(p, sx, sy, rmin, rmax, phmin, phmax);
    if (in_ring) {
      res.emplace_back(p);
    }

    // continue search for a good placement or if no placement found yet
    if (in_ring || res.empty()) {
      // check adjacent squares
      add_rectangle(Point(p.x() + sx, p.y()), res, sx, sy, rmin, rmax, phmin, phmax, max_depth, depth + 1);
      add_rectangle(Point(p.x() - sx, p.y()), res, sx, sy, rmin, rmax, phmin, phmax, max_depth, depth + 1);
      add_rectangle(Point(p.x(), p.y() + sy), res, sx, sy, rmin, rmax, phmin, phmax, max_depth, depth + 1);
      add_rectangle(Point(p.x(), p.y() - sy), res, sx, sy, rmin, rmax, phmin, phmax, max_depth, depth + 1);
    }
  }

  // fill squares
  std::vector<Point> fillRectangles(Point ref, double sx, double sy, double rmin, double rmax, double phmin,
                                    double phmax)
  {
    // convert (0, 2pi) to (-pi, pi)
    if (phmax > M_PI) {
      phmin -= M_PI;
      phmax -= M_PI;
    }
    // start with a seed square and find one in the ring
    // move to center
    ref = ref - Point(int(ref.x() / sx) * sx, int(ref.y() / sy) * sy);

    std::vector<Point> res;
    add_rectangle(ref, res, sx, sy, rmin, rmax, phmin, phmax, (int(rmax / sx) + 1) * (int(rmax / sy) + 1) * 2);
    return res;
  }

  // check if a regular polygon is inside a ring
  bool poly_in_ring(const Point& p, int nsides, double lside, double rmin, double rmax, double phmin, double phmax)
  {
    // outer radius is contained
    if ((p.r() + lside <= rmax) && (p.r() - lside >= rmin)) {
      return true;
    }

    // inner radius is not contained
    double rin = std::cos(M_PI / nsides) * lside;
    if ((p.r() + rin > rmax) || (p.r() - rin < rmin)) {
      return false;
    }

    // in between, check every corner
    for (int i = 0; i < nsides; ++i) {
      double phi = (i + 0.5) * 2. * M_PI / static_cast<double>(nsides);
      Point  p2(p.x() + 2. * lside * std::sin(phi), p.y() + 2. * lside * std::cos(phi));
      if ((p2.r() > rmax) || (p2.r() < rmin) || p.phi() > phmax || p.phi() < phmin) {
        return false;
      }
    }
    return true;
  }

  // recursively fill square (nside=4) or hexagon (nside=6) in a ring, other polygons won't work
  void add_poly(Point p, std::vector<Point>& res, int nsides, double lside, double rmin, double rmax, double phmin,
                double phmax, int max_depth = 20, int depth = 0)
  {
    // std::cout << depth << "/" << max_depth << std::endl;
    // exceeds the maximum depth in searching or already placed
    if ((depth > max_depth) || (already_placed(p, res, lside, lside))) {
      return;
    }

    bool in_ring = poly_in_ring(p, nsides, lside, rmin, rmax, phmin, phmax);
    if (in_ring) {
      res.emplace_back(p);
    }

    // recursively add neigbors, continue if it was a good placement or no placement found yet
    if (in_ring || res.empty()) {
      for (int i = 0; i < nsides; ++i) {
        double phi = i * 2. * M_PI / static_cast<double>(nsides);
        add_poly(Point(p.x() + 2. * lside * std::sin(phi), p.y() + 2. * lside * std::cos(phi)), res, nsides, lside,
                 rmin, rmax, phmin, phmax, max_depth, depth + 1);
      }
    }
  }

  std::vector<Point> fillHexagons(Point ref, double lside, double rmin, double rmax, double phmin, double phmax)
  {
    // convert (0, 2pi) to (-pi, pi)
    if (phmax > M_PI) {
      phmin -= M_PI;
      phmax -= M_PI;
    }
    // start with a seed and find one in the ring
    // move to center
    ref = ref - Point(int(ref.x() / lside) * lside, int(ref.y() / lside) * lside);

    std::vector<Point> res;
    add_poly(ref, res, 6, lside, rmin, rmax, phmin, phmax, std::pow(int(rmax / lside) + 1, 2) * 2);
    return res;
  }

} // namespace epic::geo
