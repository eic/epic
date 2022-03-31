//==========================================================================
//  Gaseous Ring Imaging Cherenkov Detector
//--------------------------------------------------------------------------
//
// Author: C. Peng (ANL)
// Date: 09/30/2020
//
//==========================================================================

#include <XML/Helper.h>
#include "TMath.h"
#include "TString.h"
#include "GeometryHelpers.h"
#include "Math/Point2D.h"
#include "DDRec/Surface.h"
#include "DDRec/DetectorData.h"
#include "DD4hep/OpticalSurfaces.h"
#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/Printout.h"

using namespace std;
using namespace dd4hep;
using namespace dd4hep::rec;

// headers
void build_radiator(Detector &desc, Volume &env, xml::Component plm, const Position &offset);
void build_mirrors(Detector &desc, DetElement &sdet, Volume &env, xml::Component plm, const Position &offset);
void build_sensors(Detector &desc, Volume &env, xml::Component plm, const Position &offset, SensitiveDetector &sens);

// helper function to get x, y, z if defined in a xml component
template<class XmlComp>
Position get_xml_xyz(XmlComp &comp, dd4hep::xml::Strng_t name)
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

// create the detector
static Ref_t createDetector(Detector& desc, xml::Handle_t handle, SensitiveDetector sens)
{
    xml::DetElement detElem = handle;
    std::string detName = detElem.nameStr();
    int detID = detElem.id();

    DetElement det(detName, detID);
    xml::Component dims = detElem.dimensions();

    // build a big envelope for all the components, filled with optical material to ensure transport of optical photons
    double z0 = dims.z0();
    double length = dims.length();
    double rmin = dims.rmin();
    double rmax0 = dims.attr<double>(_Unicode(rmax0));
    double rmax1 = dims.attr<double>(_Unicode(rmax1));
    double rmax2 = dims.attr<double>(_Unicode(rmax2));
    double snout_length = dims.attr<double>(_Unicode(snout_length));
    // fill envelope with radiator materials (need optical property for optical photons)
    auto gasMat = desc.material(dd4hep::getAttrOrDefault(detElem, _Unicode(gas), "AirOptical"));

    Cone snout(snout_length/2.0, rmin, rmax0, rmin, rmax1);
    Tube tank(rmin, rmax2, (length - snout_length)/2., 0., 2*M_PI);
    // shift the snout to the left side of tank
    UnionSolid envShape(tank, snout, Position(0., 0., -length/2.));
    // some reference points
    auto snout_front = Position(0., 0., -(length + snout_length)/2.);
    // auto snout_end = Position(0., 0., -(length - snout_length)/2.);  // tank_front
    // auto tank_end = Position(0., 0., (length - snout_length)/2.);

    Volume envVol(detName + "_envelope", envShape, gasMat);
    envVol.setVisAttributes(desc.visAttributes(detElem.visStr()));

    // sensitive detector type
    sens.setType("tracker");

    // @TODO: place a radiator
    // build_radiator(desc, envVol, detElement.child(_Unicode(radiator)), snout_front);

    // place mirrors
    build_mirrors(desc, det, envVol, detElem.child(_Unicode(mirrors)), snout_front);

    // place photo-sensors
    build_sensors(desc, envVol, detElem.child(_Unicode(sensors)), snout_front, sens);

    Volume motherVol = desc.pickMotherVolume(det);
    PlacedVolume envPV = motherVol.placeVolume(envVol, Position(0, 0, z0) - snout_front);
    envPV.addPhysVolID("system", detID);
    det.setPlacement(envPV);

    return det;
}

// @TODO: implement a radiator, now the envelope serves as the radiator
void build_radiator(Detector &desc, Volume &env, xml::Component plm, const Position &offset)
{
    // place holder
}

// place mirrors
void build_mirrors(Detector &desc, DetElement &sdet, Volume &env, xml::Component plm, const Position &offset)
{
    double thickness = dd4hep::getAttrOrDefault<double>(plm, _Unicode(dz), 1.*dd4hep::mm);
    auto mat = desc.material(plm.attr<std::string>(_Unicode(material)));
    auto vis = desc.visAttributes(plm.attr<std::string>(_Unicode(vis)));

    // optical surface
    OpticalSurfaceManager surfMgr = desc.surfaceManager();
    auto surf = surfMgr.opticalSurface(dd4hep::getAttrOrDefault(plm, _Unicode(surface), "MirrorOpticalSurface"));

    // placements
    auto gpos = get_xml_xyz(plm, _Unicode(position)) + offset;
    auto grot = get_xml_xyz(plm, _Unicode(position));
    int imod = 1;
    for (xml::Collection_t mir(plm, _Unicode(mirror)); mir; ++mir, ++imod) {
        double rmin = mir.attr<double>(_Unicode(rmin));
        double rmax = mir.attr<double>(_Unicode(rmax));
        double phiw = dd4hep::getAttrOrDefault<double>(mir, _Unicode(phiw), 2.*M_PI);
        double rotz = dd4hep::getAttrOrDefault<double>(mir, _Unicode(rotz), 0.);
        double roty = dd4hep::getAttrOrDefault<double>(mir, _Unicode(roty), 0.);
        double rotx = dd4hep::getAttrOrDefault<double>(mir, _Unicode(rotx), 0.);

        Volume vol(Form("mirror_v_dummy%d", imod));
        vol.setMaterial(mat);
        vol.setVisAttributes(vis);
        // mirror curvature
        double curve = dd4hep::getAttrOrDefault<double>(mir, _Unicode(curve), 0.);
        // spherical mirror
        if (curve > 0.) {
            double th1 = std::asin(rmin/curve);
            double th2 = std::asin(rmax/curve);
            vol.setSolid(Sphere(curve, curve + thickness, th1, th2, 0., phiw));
        // plane mirror
        } else {
            vol.setSolid(Tube(rmin, rmax, thickness/2., 0., phiw));
        }
        // transforms are in a reverse order
        Transform3D tr = Translation3D(gpos.x(), gpos.y(), gpos.z())
                       * RotationZYX(grot.z(), grot.y(), grot.x())
                       * RotationZ(rotz)                // rotation of the piece
                       * RotationY(roty)                // rotation of the piece
                       * RotationX(rotx)                // rotation of the piece
                       * Translation3D(0., 0., -curve)  // move spherical shell to origin (no move for planes)
                       * RotationZ(-phiw/2.);           // center phi angle to 0. (-phiw/2., phiw/2.)
        auto pv = env.placeVolume(vol, tr);
        DetElement de(sdet, Form("mirror_de%d", imod), imod);
        de.setPlacement(pv);
        SkinSurface skin(desc, de, Form("mirror_optical_surface%d", imod), surf, vol);
        skin.isValid();
    }
}

// place photo-sensors
void build_sensors(Detector &desc, Volume &env, xml::Component plm, const Position &offset, SensitiveDetector &sens)
{
    // build sensor unit geometry
    auto mod = plm.child(_Unicode(module));
    double sx = mod.attr<double>(_Unicode(sx));
    double sy = mod.attr<double>(_Unicode(sy));
    double sz = mod.attr<double>(_Unicode(sz));
    double gap = mod.attr<double>(_Unicode(gap));
    auto mat = desc.material(mod.attr<std::string>(_Unicode(material)));
    auto vis = desc.visAttributes(mod.attr<std::string>(_Unicode(vis)));

    Box sensor(sx/2., sy/2., sz/2.);
    Volume svol("sensor_v", sensor, mat);
    svol.setVisAttributes(vis);

    // a thin layer of cherenkov gas for accepting optical photons
    auto opt = plm.child(_Unicode(optical));
    double opthick = opt.attr<double>(_Unicode(thickness));
    auto opmat = desc.material(opt.attr<std::string>(_Unicode(material)));

    Box opshape(sx/2., sy/2., sz/2. + opthick/2.);
    Volume opvol("sensor_v_optical", opshape, opmat);
    opvol.placeVolume(svol, Position(0., 0., 0.));
    opvol.setSensitiveDetector(sens);

    // photo-detector plane envelope
    auto gpos = get_xml_xyz(plm, _Unicode(position)) + offset;
    auto grot = get_xml_xyz(plm, _Unicode(position));
    int isec = 1;
    for (xml::Collection_t sec(plm, _Unicode(sector)); sec; ++sec, ++isec) {
        double rmin = sec.attr<double>(_Unicode(rmin));
        double rmax = sec.attr<double>(_Unicode(rmax));
        double phiw = dd4hep::getAttrOrDefault<double>(sec, _Unicode(phiw), 2.*M_PI);
        double rotz = dd4hep::getAttrOrDefault<double>(sec, _Unicode(rotz), 0.);
        double roty = dd4hep::getAttrOrDefault<double>(sec, _Unicode(roty), 0.);
        double rotx = dd4hep::getAttrOrDefault<double>(sec, _Unicode(rotx), 0.);

        // fill sensors to the piece
        auto points = athena::geo::fillRectangles({0., 0.}, sx + gap, sy + gap, rmin - gap, rmax + gap, -phiw/2., phiw/2.);
        int imod = 1;
        for (auto &p : points) {
            // transofrms are in a reversed order
            Transform3D tr = Translation3D(gpos.x(), gpos.y(), gpos.z())
                           * RotationZYX(grot.z(), grot.y(), grot.x())
                           * RotationZ(rotz)                    // rotation of the sector
                           * RotationY(roty)                    // rotation of the sector
                           * RotationX(rotx)                    // rotation of the sector
                           * Translation3D(p.x(), p.y(), 0.);   // place modules in each sector
            auto pv = env.placeVolume(opvol, tr);
            pv.addPhysVolID("sector", isec).addPhysVolID("module", imod++);
        }
    }
}

//@}

// clang-format off
DECLARE_DETELEMENT(athena_GaseousRICH, createDetector)

