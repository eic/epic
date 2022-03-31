//
// Author     : Whit Armstrong (warmstrong@anl.gov)
//
#include <XML/Helper.h>
#include "TMath.h"
#include "TString.h"
#include "DDRec/Surface.h"
#include "DDRec/DetectorData.h"
#include "DD4hep/OpticalSurfaces.h"
#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/Printout.h"
#include "GeometryHelpers.h"
#include "Math/Vector3D.h"
#include "Math/AxisAngle.h"
#include "Math/VectorUtil.h"

using namespace std;
using namespace dd4hep;
using namespace dd4hep::rec;

using Placements = vector<PlacedVolume>;



static Ref_t createDetector(Detector& description, xml::Handle_t e, SensitiveDetector sens){
  xml_det_t      x_det    = e;
  Material       air      = description.material("AirOptical");
  Material       vacuum   = description.vacuum();
  string         det_name = x_det.nameStr();
  DetElement     sdet(det_name, x_det.id());
  Assembly       assembly(det_name);
  sens.setType("tracker");
  OpticalSurfaceManager surfMgr = description.surfaceManager();

  bool projective = getAttrOrDefault(x_det, _Unicode(projective), false);
  bool reflect    = x_det.reflect(true);

  PlacedVolume pv;

  map<string, Volume>     modules;
  map<string, Placements> sensitives;
  map<string, Volume>     module_assemblies;
  std::map<std::string,DetElement> module_assembly_delements;

  int                     n_sensor = 1;

  // dimensions
  xml::Component dims   = x_det.dimensions();
  auto           rmin   = dims.rmin();
  auto           rmax   = dims.rmax();
  auto           length = dims.length();
  auto           zmin   = dims.zmin();
  auto           zpos   = zmin + length / 2;

  // envelope
  Tube   envShape(rmin, rmax, length / 2., 0., 2 * M_PI);
  Volume envVol("MRICH_Envelope", envShape, air);
  envVol.setVisAttributes(description.visAttributes(x_det.visStr()));
  if (x_det.hasChild(_Unicode(envelope))) {
    xml_comp_t x_envelope = x_det.child(_Unicode(envelope));
    double thickness  = x_envelope.thickness();
    Material material = description.material(x_envelope.materialStr());
    Tube   envInsideShape(rmin + thickness, rmax - thickness, length / 2. - thickness);
    SubtractionSolid envShellShape(envShape, envInsideShape);
    Volume envShell("MRICH_Envelope_Inside", envShellShape, material);
    envVol.placeVolume(envShell);
  }

  // expect only one module (for now)
  xml_comp_t x_mod = x_det.child(_U(module));
  string     mod_name            = x_mod.nameStr();
  double     mod_width           = getAttrOrDefault(x_mod, _U(width), 130.0 * mm);
  double     mod_height          = getAttrOrDefault(x_mod, _U(height), 130.0 * mm);
  double     mod_length          = getAttrOrDefault(x_mod, _U(length), 130.0 * mm);

  // module
  Box    m_solid(mod_width / 2.0, mod_height / 2.0, mod_length / 2.0);
  Volume m_volume(mod_name, m_solid, air);
  m_volume.setVisAttributes(description.visAttributes(x_mod.visStr()));
  DetElement mod_de( mod_name + std::string("_mod_") + std::to_string(1), 1);
  double z_placement = - mod_length / 2.0;

  // todo module frame
  if (x_mod.hasChild(_Unicode(frame))) {
    xml_comp_t x_frame    = x_mod.child(_Unicode(frame));
    double     frame_thickness = getAttrOrDefault(x_frame, _U(thickness), 2.0 * mm);
    Box        frame_inside(mod_width / 2.0 - frame_thickness, mod_height / 2.0 - frame_thickness, mod_length / 2.0 - frame_thickness);
    SubtractionSolid frame_solid(m_solid, frame_inside);
    Material   frame_mat       = description.material(x_frame.materialStr());
    Volume     frame_vol(mod_name+"_frame", frame_solid, frame_mat);
    auto       frame_vis       = getAttrOrDefault<std::string>(x_frame, _U(vis), std::string("GrayVis"));
    frame_vol.setVisAttributes(description.visAttributes(frame_vis));
    // update position
    z_placement += frame_thickness / 2.0;
    // place volume
    m_volume.placeVolume(frame_vol);
    // update position
    z_placement += frame_thickness / 2.0;
  }

  // aerogel box
  if (x_mod.hasChild(_Unicode(aerogel))) {
    xml_comp_t x_aerogel  = x_mod.child(_Unicode(aerogel));
    double     aerogel_width   = getAttrOrDefault(x_aerogel, _U(width), 130.0 * mm);
    double     aerogel_length  = getAttrOrDefault(x_aerogel, _U(length), 130.0 * mm);
    Material   aerogel_mat     = description.material(x_aerogel.materialStr());
    auto       aerogel_vis     = getAttrOrDefault<std::string>(x_aerogel, _U(vis), std::string("InvisibleWithDaughters"));

    xml_comp_t x_aerogel_frame = x_aerogel.child(_Unicode(frame));
    double     foam_thickness  = getAttrOrDefault(x_aerogel_frame, _U(thickness), 2.0 * mm);
    Material   foam_mat        = description.material(x_aerogel_frame.materialStr());
    auto       foam_vis        = getAttrOrDefault<std::string>(x_aerogel_frame, _U(vis), std::string("RedVis"));

    // foam frame
    Box foam_box(aerogel_width / 2.0 + foam_thickness, aerogel_width / 2.0 + foam_thickness, (aerogel_length + foam_thickness) / 2.0);
    Box foam_sub_box(aerogel_width / 2.0, aerogel_width / 2.0, (aerogel_length + foam_thickness) / 2.0);
    SubtractionSolid foam_frame_solid(foam_box, foam_sub_box, Position(0, 0, foam_thickness));
    Volume           foam_vol(mod_name+"_aerogel_frame", foam_frame_solid, foam_mat);
    foam_vol.setVisAttributes(description.visAttributes(foam_vis));

    // aerogel
    Box              aerogel_box(aerogel_width / 2.0, aerogel_width / 2.0, (aerogel_length) / 2.0);
    Volume           aerogel_vol(mod_name+"_aerogel", aerogel_box, aerogel_mat);
    aerogel_vol.setVisAttributes(description.visAttributes(aerogel_vis));

    // update position
    z_placement += (aerogel_length + foam_thickness) / 2.0;
    // place foam frame
    pv = m_volume.placeVolume(foam_vol,Position(0,0,z_placement));
    // place aerogel
    z_placement += foam_thickness / 2.0;
    pv = m_volume.placeVolume(aerogel_vol,Position(0,0,z_placement));
    DetElement aerogel_de(mod_de, mod_name + std::string("_aerogel_de") + std::to_string(1), 1);
    aerogel_de.setPlacement(pv);
    // update position
    z_placement += aerogel_length / 2.0;

    // optical surfaces
    auto aerogel_surf = surfMgr.opticalSurface(dd4hep::getAttrOrDefault(x_aerogel, _Unicode(surface), "MRICH_AerogelOpticalSurface"));
    SkinSurface skin_surf(description, aerogel_de, Form("MRICH_aerogel_skin_surface_%d", 1), aerogel_surf, aerogel_vol);
    skin_surf.isValid();
  }

  // Fresnel Lens
  if (x_mod.hasChild(_Unicode(lens))) {
    xml_comp_t x_lens     = x_mod.child(_Unicode(lens));

    //  - The lens has a constant groove pitch (delta r) as opposed to fixing the groove height.
    //  - The lens area outside of the effective diamtere is flat.
    //  - The grooves are not curved, rather they are polycone shaped, ie a flat approximating the curvature.
    auto   lens_vis       = getAttrOrDefault<std::string>(x_lens, _U(vis), std::string("AnlBlue"));
    double groove_pitch   = getAttrOrDefault(x_lens, _Unicode(pitch), 0.2 * mm);// 0.5 * mm);
    double lens_f         = getAttrOrDefault(x_lens, _Unicode(focal_length), 6.0*2.54*cm);
    double eff_diameter   = getAttrOrDefault(x_lens, _Unicode(effective_diameter), 152.4 * mm);
    double lens_width     = getAttrOrDefault(x_lens, _Unicode(width), 6.7*2.54*cm);
    double center_thickness = getAttrOrDefault(x_lens, _U(thickness), 0.068 * 2.54 * cm);//2.0 * mm);

    double n_acrylic        = 1.49;
    double lens_curvature   = 1.0 / (lens_f*(n_acrylic - 1.0)); //confirmed
    double full_ring_rmax   = std::min(eff_diameter / 2.0, lens_width/2.0);

    double N_grooves        = std::ceil((full_ring_rmax) / groove_pitch);
    double groove_last_rmin = (N_grooves - 1) * groove_pitch;
    double groove_last_rmax = N_grooves * groove_pitch;

    auto   groove_sagitta = [&](double r) { return lens_curvature * std::pow(r, 2) / (1.0 + 1.0); };
    double lens_thickness = groove_sagitta(groove_last_rmax) - groove_sagitta(groove_last_rmin) + center_thickness;

    Material         lens_mat = description.material(x_lens.materialStr());
    Box              lens_box(lens_width / 2.0, lens_width / 2.0, (center_thickness) / 2.0);
    SubtractionSolid flat_lens(lens_box, Tube(0.0, full_ring_rmax, 2 * center_thickness));

    Assembly lens_vol(mod_name + "_lens");
    Volume   flatpart_lens_vol( "flatpart_lens", flat_lens, lens_mat);
    lens_vol.placeVolume(flatpart_lens_vol);

    int    i_groove           = 0;
    double groove_rmax        = groove_pitch;
    double groove_rmin        = 0;

    while ( groove_rmax <= full_ring_rmax ) {
      double   dZ = groove_sagitta(groove_rmax) - groove_sagitta(groove_rmin);
      Polycone groove_solid(0, 2.0 * M_PI,
                            {groove_rmin, groove_rmin, groove_rmin},
                            {groove_rmax, groove_rmax, groove_rmin},
                            {-lens_thickness/2.0, lens_thickness/2.0-dZ, lens_thickness/2.0});
      Volume   lens_groove_vol("lens_groove_" + std::to_string(i_groove), groove_solid, lens_mat);
      lens_vol.placeVolume(lens_groove_vol);

      i_groove++;
      groove_rmin = (i_groove  )*groove_pitch;
      groove_rmax = (i_groove+1)*groove_pitch;
    }

    lens_vol.setVisAttributes(description.visAttributes(lens_vis));

    // update position
    z_placement += lens_thickness/2.0;
    // place volume
    pv = m_volume.placeVolume(lens_vol,Position(0,0,z_placement));
    DetElement lens_de(mod_de, mod_name + std::string("_lens_de") + std::to_string(1), 1);
    lens_de.setPlacement(pv);
    // update position
    z_placement += lens_thickness/2.0;

    // optical surfaces
    auto lens_surf = surfMgr.opticalSurface(dd4hep::getAttrOrDefault(x_lens, _Unicode(surface), "MRICH_LensOpticalSurface"));
    SkinSurface skin_surf(description, lens_de, Form("MRichFresnelLens_skin_surface_%d", 1), lens_surf, lens_vol);
    skin_surf.isValid();
  }

  // mirror
  if (x_mod.hasChild(_Unicode(space))) {
    xml_comp_t x_space = x_mod.child(_Unicode(space));
    z_placement += getAttrOrDefault(x_space, _U(thickness), 0.0 * mm);
  }

  // mirror
  if (x_mod.hasChild(_Unicode(mirror))) {
    xml_comp_t x_mirror   = x_mod.child(_Unicode(mirror));
    auto   mirror_vis     = getAttrOrDefault<std::string>(x_mirror, _U(vis), std::string("AnlGray"));
    double mirror_x1      = getAttrOrDefault(x_mirror, _U(x1), 100.0 * mm);
    double mirror_x2      = getAttrOrDefault(x_mirror, _U(x2), 80.0 * mm);
    double mirror_length  = getAttrOrDefault(x_mirror, _U(length), 130.0 * mm);
    double mirror_thickness  = getAttrOrDefault(x_mirror, _U(thickness), 2.0 * mm);
    double outer_x1 = (mirror_x1+mirror_thickness)/2.0;
    double outer_x2 = (mirror_x2+mirror_thickness)/2.0;
    Trd2   outer_mirror_trd(outer_x1, outer_x2, outer_x1,  outer_x2, mirror_length/2.0);
    Trd2   inner_mirror_trd(mirror_x1 / 2.0,  mirror_x2 / 2.0, mirror_x1 / 2.0,mirror_x2 / 2.0, mirror_length/2.0+0.1*mm);
    SubtractionSolid mirror_solid(outer_mirror_trd, inner_mirror_trd);
    Material mirror_mat        = description.material(x_mirror.materialStr());
    Volume   mirror_vol(mod_name+"_mirror", mirror_solid, mirror_mat);

    // update position
    z_placement += mirror_length/2.0;
    // place volume
    pv = m_volume.placeVolume(mirror_vol,Position(0,0,z_placement));
    DetElement mirror_de(mod_de, mod_name + std::string("_mirror_de") + std::to_string(1), 1);
    mirror_de.setPlacement(pv);
    // update position
    z_placement += mirror_length/2.0;

    // optical surfaces
    auto mirror_surf = surfMgr.opticalSurface(dd4hep::getAttrOrDefault(x_mirror, _Unicode(surface), "MRICH_MirrorOpticalSurface"));
    SkinSurface skin_surf(description, mirror_de, Form("MRICH_mirror_skin_surface_%d", 1), mirror_surf, mirror_vol);
    skin_surf.isValid();
  }

  // photon detector
  if (x_mod.hasChild(_Unicode(photodet))) {
    xml_comp_t x_photodet = x_mod.child(_Unicode(photodet));
    auto       photodet_vis       = getAttrOrDefault<std::string>(x_photodet, _U(vis), std::string("AnlRed"));
    double     photodet_width     = getAttrOrDefault(x_photodet, _U(width), 130.0 * mm);
    double     photodet_thickness = getAttrOrDefault(x_photodet, _U(thickness), 2.0 * mm);
    Material   photodet_mat       = description.material(x_photodet.materialStr());
    Box        window_box(photodet_width/2.0,photodet_width/2.0,photodet_thickness/2.0);
    Volume     window_vol(mod_name+"_window", window_box, photodet_mat);

    // update position
    z_placement += photodet_thickness/2.0;
    // place volume
    pv = m_volume.placeVolume(window_vol,Position(0,0,z_placement));
    DetElement   comp_de(mod_de, mod_name + std::string("_sensor_de_") + std::to_string(1) ,  1);
    comp_de.setPlacement(pv);
    // update position
    z_placement += photodet_thickness/2.0;

    // sensitive
    pv.addPhysVolID("sensor", n_sensor);
    window_vol.setSensitiveDetector(sens);
    sensitives[mod_name].push_back(pv);
    ++n_sensor;

    // sensor
    xml_comp_t x_sensor  = x_photodet.child(_Unicode(sensor));
    double     sensor_thickness   = getAttrOrDefault(x_sensor, _U(thickness), 2.0 * mm);
    Material   sensor_mat         = description.material(x_sensor.materialStr());
    int        sensor_nx          = getAttrOrDefault(x_sensor, _Unicode(nx), 2);
    int        sensor_ny          = getAttrOrDefault(x_sensor, _Unicode(ny), 2);

    // layers
    int i_layer = 1;
    for (xml_coll_t li(x_photodet, _Unicode(layer)); li; ++li) {
      xml_comp_t x_layer = li;
      Material layer_mat = description.material(x_layer.materialStr());
      double   layer_thickness = x_layer.thickness();
      Box      layer_box(photodet_width/2.0,photodet_width/2.0,layer_thickness/2.0);
      Volume   layer_vol(mod_name + "_layer_" + std::to_string(i_layer), layer_box, layer_mat);

      // update position
      z_placement += layer_thickness / 2.0;
      // place volume
      pv = m_volume.placeVolume(layer_vol,Position(0,0,z_placement));
      DetElement layer_de(mod_de, mod_name + std::string("_layer_de_") + std::to_string(i_layer),  1);
      layer_de.setPlacement(pv);
      // update position
      z_placement += layer_thickness / 2.0;

      i_layer++;
    }
  }

  //for (size_t ic = 0; ic < sensVols.size(); ++ic) {
  //  PlacedVolume sens_pv = sensVols[ic];
  //  DetElement   comp_de(mod_de, std::string("de_") + sens_pv.volume().name(), ic + 1);
  //  comp_de.setPlacement(sens_pv);
  //  // Acts::ActsExtension* sensorExtension = new Acts::ActsExtension();
  //  //// sensorExtension->addType("sensor", "detector");
  //  // comp_de.addExtension<Acts::ActsExtension>(sensorExtension);
  //  //// comp_de.setAttributes(description, sens_pv.volume(), x_layer.regionStr(),
  //  //// x_layer.limitsStr(),
  //  ////                      xml_det_t(xmleles[m_nam]).visStr());
  //}
  //DetElement window_de(sdet, mod_name + std::string("_window_de") + std::to_string(1), 1);
  //window_de.setPlacement(pv);

  modules[mod_name]                   = m_volume;
  module_assembly_delements[mod_name] = mod_de;
  // end module

  // place modules in the sectors (disk)
  auto points = athena::geo::fillSquares({0., 0.}, mod_width, rmin, rmax);

  // mod_name = ...
  Placements& sensVols = sensitives[mod_name];
  auto        mod_v    = modules[mod_name];
  // determine module direction, always facing z = 0
  double roty = dims.zmin() < 0. ? -M_PI : 0 ;

  // read module positions
  std::vector<std::tuple<double,double,double>> positions;
  for (xml_coll_t x_positions_i(x_det, _Unicode(positions)); x_positions_i; ++x_positions_i) {
    xml_comp_t x_positions = x_positions_i;
    for (xml_coll_t x_position_i(x_positions, _U(position)); x_position_i; ++x_position_i) {
      xml_comp_t x_position = x_position_i;
      positions.push_back(
        std::make_tuple(x_positions.scale() * x_position.x() * mm,
                        x_positions.scale() * x_position.y() * mm,
                        -x_positions.z0()));
    }
  }
  // if no positions, then autoplacement
  if (positions.empty()) {
    for (double x = mod_width / 2.0; x < rmax - mod_width / 2.0; x += mod_width) {
      for (double y = mod_width / 2.0; y < rmax - mod_width / 2.0; y += mod_width) {
        if (pow(x + mod_width / 2.0,2) + pow(y + mod_width / 2.0,2) > rmax*rmax) continue;
        if (pow(x - mod_width / 2.0,2) + pow(y - mod_width / 2.0,2) < rmin*rmin) continue;
        positions.push_back(std::make_tuple(x, y, 0));
      }
    }
  }

  // place modules
  int i_mod = 1; // starts at 1
  for (auto& p: positions) {

    // get positions in one quadrant
    double x = std::get<0>(p);
    double y = std::get<1>(p);
    double z0 = std::get<2>(p);

    // and place in all quadrants (intentional shadowing)
    for (auto& p: decltype(positions){{x,y,z0}, {y,-x,z0}, {-x,-y,z0}, {-y,x,z0}}) {

      // get positions (intentional shadowing)
      double x = std::get<0>(p);
      double y = std::get<1>(p);
      double z0 = std::get<2>(p);

      Transform3D tr;
      if(projective) {
        double rotAngX = atan(y/z0);
        double rotAngY = -1.*atan(x/z0);
        tr = Translation3D(x, y, 0) * RotationX(rotAngX) * RotationY(rotAngY);
      } else {
        tr = Translation3D(x, y, 0) * RotationX(0);
      }

      // mod placement
      pv = envVol.placeVolume(mod_v, tr);
      pv.addPhysVolID("module", i_mod);

      auto mod_det_element  = module_assembly_delements[mod_name].clone(mod_name + "__" + std::to_string(i_mod));
      mod_det_element.setPlacement(pv);
      sdet.add(mod_det_element);

      i_mod++;
    }
  }

  // additional layers
  if (x_det.hasChild(_Unicode(layer))) {
    xml_comp_t x_layer = x_det.child(_Unicode(layer));
    double   layer_thickness = x_layer.thickness();
    Material layer_mat = description.material(x_layer.materialStr());
    Tube   frameShape(rmin, rmax, layer_thickness / 2., 0., 2 * M_PI);
    Volume frameVol("MRICH_Frame", frameShape, layer_mat);
    pv = envVol.placeVolume(frameVol, Position(0, 0, (length - layer_thickness) / 2.0));
  }

  // place envelope
  Volume motherVol = description.pickMotherVolume(sdet);
  if (reflect) {
    pv = motherVol.placeVolume(envVol, Transform3D(RotationZYX(0, M_PI, 0), Position(0, 0, -zpos)));
  } else {
    pv = motherVol.placeVolume(envVol, Transform3D(RotationZYX(0,    0, 0), Position(0, 0, +zpos)));
  }
  pv.addPhysVolID("system", x_det.id());
  sdet.setPlacement(pv);
  return sdet;
}


//void addModules(Volume &mother, xml::DetElement &detElem, Detector &description, SensitiveDetector &sens)
//{
//    xml::Component dims = detElem.dimensions();
//    xml::Component mods = detElem.child(_Unicode(modules));
//
//    auto rmin = dims.rmin();
//    auto rmax = dims.rmax();
//
//    auto mThick = mods.attr<double>(_Unicode(thickness));
//    auto mWidth = mods.attr<double>(_Unicode(width));
//    auto mGap = mods.attr<double>(_Unicode(gap));
//
//    auto modMat = description.material(mods.materialStr());
//    auto gasMat = description.material("AirOptical");
//
//    // single module
//    Box mShape(mWidth/2., mWidth/2., mThick/2. - 0.1*mm);
//    Volume mVol("ce_MRICH_mod_Solid", mShape, modMat);
//
//    // a thin gas layer to detect optical photons
//    Box modShape(mWidth/2., mWidth/2., mThick/2.);
//    Volume modVol("ce_MRICH_mod_Solid_v", modShape, gasMat);
//    // thin gas layer is on top (+z) of the material
//    modVol.placeVolume(mVol, Position(0., 0., -0.1*mm));
//
//    modVol.setVisAttributes(description.visAttributes(mods.visStr()));
//    sens.setType("tracker");
//    modVol.setSensitiveDetector(sens);
//
//    // place modules in the sectors (disk)
//    auto points = ref::utils::fillSquares({0., 0.}, mWidth + mGap, rmin - mGap, rmax + mGap);
//
//    // determine module direction, always facing z = 0
//    double roty = dims.z() > 0. ? M_PI/2. : -M_PI/2.;
//    int imod = 1;
//    for (auto &p : points) {
//        // operations are inversely ordered
//        Transform3D tr = Translation3D(p.x(), p.y(), 0.)        // move to position
//                       * RotationY(roty);                       // facing z = 0.
//        auto modPV = mother.placeVolume(modVol, tr);
//        modPV.addPhysVolID("sector", 1).addPhysVolID("module", imod ++);
//    }
//}

// clang-format off
DECLARE_DETELEMENT(athena_MRICH, createDetector)

