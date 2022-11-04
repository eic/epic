#include "GeometryHelper.h"

namespace ip6::geo {

  Position get_xml_xyz(xml_coll_t& comp, dd4hep::xml::Strng_t name)
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

  // place modules, id must be provided
  tuple<int, int>
  add_individuals(function<tuple<Volume, Position>(Detector&, xml_coll_t&, SensitiveDetector&)> build_module,
                  Detector& desc, Assembly& env, xml_coll_t& plm, SensitiveDetector& sens, int sid)
  {
    auto [modVol, modSize] = build_module(desc, plm, sens);
    int sector_id          = dd4hep::getAttrOrDefault<int>(plm, _Unicode(sector), sid);
    int nmodules           = 0;
    for (xml_coll_t pl(plm, _Unicode(placement)); pl; ++pl) {
      Position    pos(dd4hep::getAttrOrDefault<double>(pl, _Unicode(x), 0.),
                      dd4hep::getAttrOrDefault<double>(pl, _Unicode(y), 0.),
                      dd4hep::getAttrOrDefault<double>(pl, _Unicode(z), 0.));
      Position    rot(dd4hep::getAttrOrDefault<double>(pl, _Unicode(rotx), 0.),
                      dd4hep::getAttrOrDefault<double>(pl, _Unicode(roty), 0.),
                      dd4hep::getAttrOrDefault<double>(pl, _Unicode(rotz), 0.));
      auto        mid   = pl.attr<int>(_Unicode(id));
      Transform3D tr    = Translation3D(pos.x(), pos.y(), pos.z()) * RotationZYX(rot.z(), rot.y(), rot.x());
      auto        modPV = env.placeVolume(modVol, tr);
      modPV.addPhysVolID("sector", sector_id).addPhysVolID("module", mid);
      nmodules++;
    }

    return {sector_id, nmodules};
  }

  // place disk of modules
  tuple<int, int> add_disk(function<tuple<Volume, Position>(Detector&, xml_coll_t&, SensitiveDetector&)> build_module,
                           Detector& desc, Assembly& env, xml_coll_t& plm, SensitiveDetector& sens, int sid)
  {
    auto [modVol, modSize]       = build_module(desc, plm, sens);
    int    sector_id             = dd4hep::getAttrOrDefault<int>(plm, _Unicode(sector), sid);
    int    id_begin              = dd4hep::getAttrOrDefault<int>(plm, _Unicode(id_begin), 1);
    double rmin                  = plm.attr<double>(_Unicode(rmin));
    double rintermediate         = plm.attr<double>(_Unicode(rintermediate));
    double r_envelopeclearance   = dd4hep::getAttrOrDefault<int>(plm, _Unicode(r_envelopeclearance), 0);
    double phi_envelopeclearance = dd4hep::getAttrOrDefault<int>(plm, _Unicode(phi_envelopeclearance), 0);
    double rmax                  = plm.attr<double>(_Unicode(rmax));
    double phimin                = dd4hep::getAttrOrDefault<double>(plm, _Unicode(phimin), 0.);
    double phimax                = dd4hep::getAttrOrDefault<double>(plm, _Unicode(phimax), 2. * M_PI);

    // placement inside mother
    auto pos = get_xml_xyz(plm, _Unicode(position));
    auto rot = get_xml_xyz(plm, _Unicode(rotation));

    // optional envelope volume
    bool       has_envelope = dd4hep::getAttrOrDefault<bool>(plm, _Unicode(envelope), false);
    Material   material     = desc.material(dd4hep::getAttrOrDefault<string>(plm, _U(material), "Air"));
    Tube       inner_solid(rmin, rintermediate + r_envelopeclearance, modSize.z() / 2.0, 0, 2. * M_PI);
    Tube       outer_solid(rintermediate, rmax + r_envelopeclearance, modSize.z() / 2.0, phimin - phi_envelopeclearance,
                           phimax + phi_envelopeclearance);
    UnionSolid solid(inner_solid, outer_solid);
    Volume     env_vol(string(env.name()) + "_envelope", solid, material);
    Transform3D tr_global = RotationZYX(rot.z(), rot.y(), rot.x()) * Translation3D(pos.x(), pos.y(), pos.z());
    if (has_envelope) {
      env.placeVolume(env_vol, tr_global);
    }

    // local placement of modules
    int  mid    = 0;
    auto points = fillRectangles({0., 0.}, modSize.x(), modSize.y(), rmin, rintermediate, rmax, phimin, phimax);
    for (auto& p : points) {
      Transform3D tr_local = RotationZYX(0.0, 0.0, 0.0) * Translation3D(p.x(), p.y(), 0.0);
      auto        modPV =
          (has_envelope ? env_vol.placeVolume(modVol, tr_local) : env.placeVolume(modVol, tr_global * tr_local));
      modPV.addPhysVolID("sector", sector_id).addPhysVolID("module", id_begin + mid++);
    }
    return {sector_id, mid};
  }

  // check if a 2d point is already in the container
  bool already_placed(const Point& p, const vector<Point>& vec, double xs = 1.0, double ys = 1.0, double tol = 1e-6)
  {
    for (auto& pt : vec) {
      if ((std::abs(pt.x() - p.x()) / xs < tol) && std::abs(pt.y() - p.y()) / ys < tol) {
        return true;
      }
    }
    return false;
  }

  // check if a point is in a ring
  inline bool rec_in_ring(const Point& pt, double sx, double sy, double rmin, double rintermediate, double rmax,
                          double phmin, double phmax)
  {

    // check four corners
    vector<Point> pts{
        Point(pt.x() - sx / 2., pt.y() - sy / 2.),
        Point(pt.x() - sx / 2., pt.y() + sy / 2.),
        Point(pt.x() + sx / 2., pt.y() - sy / 2.),
        Point(pt.x() + sx / 2., pt.y() + sy / 2.),
    };

    bool inside   = false;
    int  minindex = 0;
    int  i        = 0;

    for (auto& p : pts) {
      minindex = (p.r() < pts[minindex].r()) ? i : minindex;
      i++;
    }

    double rmax_pacman = (pts[minindex].phi() < phmin || pts[minindex].phi() > phmax) ? rintermediate : rmax;
    inside             = pts[minindex].r() <= rmax_pacman && pts[minindex].r() >= rmin;

    return inside;
  }

  // a helper function to recursively fill square in a ring
  void add_rectangle(Point p, vector<Point>& res, double sx, double sy, double rmin, double rintermediate, double rmax,
                     double phmin, double phmax, int max_depth = 20, int depth = 0)
  {
    // exceeds the maximum depth in searching or already placed
    if ((depth > max_depth) || (already_placed(p, res, sx, sy))) {
      return;
    }

    bool in_ring = rec_in_ring(p, sx, sy, rmin, rintermediate, rmax, phmin, phmax);
    if (in_ring) {
      res.emplace_back(p);
    }
    // continue search for a good placement or if no placement found yet
    if (in_ring || res.empty()) {
      // check adjacent squares
      add_rectangle(Point(p.x() + sx, p.y()), res, sx, sy, rmin, rintermediate, rmax, phmin, phmax, max_depth,
                    depth + 1);
      add_rectangle(Point(p.x() - sx, p.y()), res, sx, sy, rmin, rintermediate, rmax, phmin, phmax, max_depth,
                    depth + 1);
      add_rectangle(Point(p.x(), p.y() + sy), res, sx, sy, rmin, rintermediate, rmax, phmin, phmax, max_depth,
                    depth + 1);
      add_rectangle(Point(p.x(), p.y() - sy), res, sx, sy, rmin, rintermediate, rmax, phmin, phmax, max_depth,
                    depth + 1);
    }
  }

  // fill squares
  vector<Point> fillRectangles(Point ref, double sx, double sy, double rmin, double rintermediate, double rmax,
                               double phmin, double phmax)
  {
    // convert (0, 2pi) to (-pi, pi)
    if (phmax > M_PI) {
      phmin -= M_PI;
      phmax -= M_PI;
    }
    // start with a seed square and find one in the ring
    // move to center
    ref = ref - Point(int(ref.x() / sx) * sx, int(ref.y() / sy) * sy);
    vector<Point> res;
    add_rectangle(ref, res, sx, sy, rmin, rintermediate, rmax, phmin, phmax,
                  (int(rmax / sx) + 1) * (int(rmax / sy) + 1) * 2);
    return res;
  }
} // namespace ip6::geo
