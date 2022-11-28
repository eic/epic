#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/OpticalSurfaces.h"
#include "DD4hep/Printout.h"
#include "DDRec/DetectorData.h"
#include "DDRec/Surface.h"
#include <XML/Helper.h>

//////////////////////////////////////////////////
// Low Q2 tagger calorimeters based on HomogeniousCalorimeter_geo.cpp
//////////////////////////////////////////////////

using namespace std;
using namespace dd4hep;

// Helper function to make the calorimeter detectors
//static void Make_Sampling_Calorimeter(Detector& desc, xml_coll_t& mod, Assembly& env, SensitiveDetector& sens);
static void Make_PbW_Calorimeter(Detector& desc, xml_coll_t& mod, Assembly& env, SensitiveDetector& sens);
static Volume build_crystal(Detector& desc, xml_coll_t& mod, SensitiveDetector& sens);

static Ref_t create_detector(Detector& desc, xml_h e, SensitiveDetector sens)
{

  xml_det_t x_det   = e;
  string    detName = x_det.nameStr();
  int       detID   = x_det.id();

  string vis_name = dd4hep::getAttrOrDefault<std::string>(x_det, _Unicode(vis), "BackwardsCalorimeter");
  sens.setType("calorimeter");

  // Central focal point of the geometry
  xml::Component pos = x_det.child(_Unicode(focus));

  // Beamline rotation
  xml_dim_t rot = x_det.rotation();

  // Beampipe thickness
  double wall = dd4hep::getAttrOrDefault<double>(x_det, _Unicode(wall), 1 * mm);

  Assembly DetAssemblyAir(detName + "_assembly");

  DetElement det(detName, detID);

  //-----------------------------------------------------------------
  // Add Tagger box containers and vacuum box extension for modules
  //-----------------------------------------------------------------
  for (xml_coll_t mod(x_det, _Unicode(module)); mod; ++mod) {

    int    moduleID   = dd4hep::getAttrOrDefault<int>(mod, _Unicode(id), 0);
    string moduleName = dd4hep::getAttrOrDefault<std::string>(mod, _Unicode(modname), "Tagger0");

    // Offset from the electron beam
    double tagoff = dd4hep::getAttrOrDefault<double>(mod, _Unicode(offset_min), 50.0 * mm);

    // Overlap left beyond theta setting
    double overlap = dd4hep::getAttrOrDefault<double>(mod, _Unicode(overlap), 0.0 * mm);

    // Theta coverage expected
    double thetamin = dd4hep::getAttrOrDefault<double>(mod, _Unicode(theta_min), 0.030 * rad) - rot.theta();
    double thetamax = dd4hep::getAttrOrDefault<double>(mod, _Unicode(theta_max), 0.030 * rad) - rot.theta();

    // Align box to max or minimum theta expected at the tagger from focal point
    bool max_align = dd4hep::getAttrOrDefault<bool>(mod, _Unicode(max_align), false);


    // Size f the actual tagger box, replicated in BackwardsTagger
    xml_dim_t moddim  = mod.child(_Unicode(dimensions));
    double    w       = moddim.x();
    double    tagboxL = moddim.z();

    // Width and height of box volume
    auto vac_w = w;
    auto box_w = w + wall;

    auto theta      = thetamin;
    auto offsetx    = -(box_w - wall) * (cos(theta));
    auto vacoffsetx = -vac_w * (cos(theta));
    auto vacoffsetz = vac_w * (sin(theta));
    auto l          = (tagoff) / (sin(theta));
    //    auto l          = (tagoff)/(sin(theta));

    auto tagoffsetx = vacoffsetx - (l + tagboxL) * sin(theta);
    auto tagoffsetz = vacoffsetz - (l + tagboxL) * cos(theta);
    //     auto tagoffsetx = vacoffsetx-(l+tagboxL/2)*sin(theta);
    //     auto tagoffsetz = vacoffsetz-(l+tagboxL/2)*cos(theta);

    if (max_align) {
      theta      = thetamax;
      offsetx    = (overlap+box_w - wall) * (cos(theta));
      vacoffsetx = (overlap+vac_w) * (cos(theta));
      vacoffsetz = -(overlap+vac_w) * (sin(theta));
      l          = (2 * offsetx + tagoff) / sin(theta);
      tagoffsetx = vacoffsetx - (l + tagboxL) * sin(theta);
      tagoffsetz = vacoffsetz - (l + tagboxL) * cos(theta);
      //       tagoffsetx = -wall+vacoffsetx-(l+tagboxL/2)*sin(theta);
      //       tagoffsetz = vacoffsetz-(l+tagboxL/2)*cos(theta);
    }


    RotationY rotate(theta);

    Assembly TaggerAssembly("Calorimeter_module_assembly");

    Make_PbW_Calorimeter(desc,mod,TaggerAssembly,sens);

    PlacedVolume pv_mod2 = DetAssemblyAir.placeVolume(
        TaggerAssembly,
        Transform3D(rotate, Position(tagoffsetx, 0,
                                     tagoffsetz))); // Very strange y is not centered and offset needs correcting for...

    DetElement moddet(moduleName, moduleID);
    pv_mod2.addPhysVolID("module", moduleID);
    moddet.setPlacement(pv_mod2);
    det.add(moddet);

  }

  // placement in mother volume
  Transform3D  tr(RotationY(rot.theta()), Position(pos.x(), pos.y(), pos.z()));
  PlacedVolume detPV = desc.pickMotherVolume(det).placeVolume(DetAssemblyAir, tr);
  detPV.addPhysVolID("system", detID);

  det.setPlacement(detPV);

  return det;
}


// place array of modules
static void Make_PbW_Calorimeter(Detector& desc, xml_coll_t& mod, Assembly& env, SensitiveDetector& sens)
{

  // placement inside mother
  xml_dim_t pos = mod.child(_Unicode(offset));

  

  xml_dim_t modSize =  mod.child(_Unicode(dimensions));
  int nrow               = mod.attr<int>(_Unicode(nrow));
  int ncol               = mod.attr<int>(_Unicode(ncol));

  double componentX = (2*modSize.x()+pos.x()) / ncol;
  double componentY = 2*modSize.y() / nrow;


  // compute array position
  double begx = -modSize.x()-pos.x() + componentX;
  double begy = -modSize.y() + componentY;


  // optional envelope volume
  Volume      crystal      = build_crystal(desc, mod, sens); //Currently without wrapper
  Assembly    env_vol(std::string(env.name()) + "_envelope");
  Transform3D tr_global = RotationZYX(0.0, 0.0, 0.0) * Translation3D(pos.x(), pos.y(), pos.z());

  env.placeVolume(env_vol, tr_global);

  // local placement of modules
  int ncrystals = 0;
  for (int i = 0; i < nrow; ++i) {
    for (int j = 0; j < ncol; ++j) {
      double      px       = begx + componentX * j;
      double      py       = begy + componentY * i;
      Transform3D tr_local = RotationZYX(0.0, 0.0, 0.0) * Translation3D(px, py, 0.0);
      auto        modPV = env_vol.placeVolume(crystal, tr_local);
      modPV.addPhysVolID("crystal", i * ncol + j);
      ncrystals++;
    }
  }
}

// helper function to build module with or w/o wrapper
static Volume build_crystal(Detector& desc, xml_coll_t& mod, SensitiveDetector& sens)
{
  xml_dim_t modSize =  mod.child(_Unicode(dimensions));
  int nrow               = mod.attr<int>(_Unicode(nrow));
  int ncol               = mod.attr<int>(_Unicode(ncol));
  double X = modSize.x() / ncol;
  double Y = modSize.y() / nrow;
  double Z = modSize.z();

  Box crystalshape(X, Y, Z);

  Material   crystalMat    = desc.material(getAttrOrDefault<std::string>(mod, _U(material), "leadtungsten_optical"));

  Volume crystalVol("crystal_vol", crystalshape, crystalMat);
  crystalVol.setVisAttributes(desc.visAttributes(mod.attr<std::string>(_Unicode(cryvis))));
  crystalVol.setSensitiveDetector(sens);

  return crystalVol;

//   if (!plm.hasChild(_Unicode(wrapper)))    // no wrapper
//     {
//       printout(DEBUG, "HomogeneousCalorimeter", "without wrapper");

//       return std::make_tuple(modVol, Position{sx, sy, sz});

//     }
//   else // build wrapper
//     {
//       auto wrp       = plm.child(_Unicode(wrapper)); //Read all the content in the wrapper block
//       auto carbon_thickness = wrp.attr<double>(_Unicode(carbon_thickness));
//       auto wrap_thickness = wrp.attr<double>(_Unicode(wrap_thickness));
//       auto length_wrapper = wrp.attr<double>(_Unicode(length));
//       auto carbonMat = desc.material(wrp.attr<std::string>(_Unicode(material_carbon)));
//       auto wrpMat = desc.material(wrp.attr<std::string>(_Unicode(material_wrap)));
//       auto gapMat = desc.material(wrp.attr<std::string>(_Unicode(material_gap)));

//       if (carbon_thickness < 1e-12 * mm)
//         return std::make_tuple(modVol, Position{sx, sy, sz});

//       Box carbonShape(sx / 2., sy / 2., length_wrapper / 2.);
//       Box carbonShape_sub((sx - 2.*carbon_thickness) / 2., (sy - 2.*carbon_thickness) / 2., length_wrapper / 2.);
//       SubtractionSolid carbon_subtract(carbonShape, carbonShape_sub, Position(0., 0., 0.));

//       Box wrpShape((sx - 2.*carbon_thickness) / 2., (sy - 2.*carbon_thickness) / 2., (sz + wrap_thickness) / 2.);
//       Box wrpShape_sub((sx - 2.*carbon_thickness - 2.*wrap_thickness) / 2., (sy - 2.*carbon_thickness - 2.*wrap_thickness) / 2., sz / 2.);
//       SubtractionSolid wrp_subtract(wrpShape, wrpShape_sub, Position(0., 0., -wrap_thickness));

//       Box gapShape(sx / 2., sy / 2., (sz - 2. * length_wrapper) / 2.);
//       Box gapShape_sub((sx - 2.*carbon_thickness) / 2., (sy - 2.*carbon_thickness) / 2., (sz - 2. * length_wrapper) / 2.);
//       SubtractionSolid gap_subtract(gapShape, gapShape_sub, Position(0., 0., 0.));

//       Volume carbonVol("carbon_vol", carbon_subtract, carbonMat);
//       Volume wrpVol("wrapper_vol", wrp_subtract, wrpMat);
//       Volume gapVol("gap_vol", gap_subtract, gapMat);

//       modVol.placeVolume(carbonVol, Position(0., 0., (sz - length_wrapper - wrap_thickness) / 2.)); //put the wrap in the both ends of crystal
//       modVol.placeVolume(carbonVol, Position(0., 0., (length_wrapper - sz - wrap_thickness) / 2.));
//       modVol.placeVolume(wrpVol, Position(0., 0., wrap_thickness / 2.));
//       modVol.placeVolume(gapVol, Position(0., 0., -wrap_thickness / 2.)); //put the gap in the middle of crystal

//       carbonVol.setVisAttributes(desc.visAttributes(wrp.attr<std::string>(_Unicode(vis_carbon))));
//       wrpVol.setVisAttributes(desc.visAttributes(wrp.attr<std::string>(_Unicode(vis_wrap))));
//       gapVol.setVisAttributes(desc.visAttributes(wrp.attr<std::string>(_Unicode(vis_gap))));

//       printout(DEBUG, "HomogeneousCalorimeter", "with wrapper");

//       return std::make_tuple(modVol, Position{sx, sy, sz});
//     }
}

DECLARE_DETELEMENT(BackwardsTaggerCalorimeter, create_detector)
