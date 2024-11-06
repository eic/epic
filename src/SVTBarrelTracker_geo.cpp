// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022-2024 Whitney Armstrong, Nivedith Ramasubramanian, Yann Bedfer, Shujie Li, Jonathan Witte

/** \addtogroup Trackers Trackers
 * \brief Type: **Curved Silicon Vertex Tracker Barrel with inactive area**.
 * \author Shujie Li, Jonathan Witte
 *
 *
 * \ingroup trackers
 *
 * @{
 */
#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/Printout.h"
#include "DD4hep/Shapes.h"
#include "DDRec/DetectorData.h"
#include "DDRec/Surface.h"
#include "XML/Layering.h"
#include "XML/Utilities.h"
#include <array>
#include "DD4hepDetectorHelper.h"

using namespace std;
using namespace dd4hep;
using namespace dd4hep::rec;
using namespace dd4hep::detail;

#include "Math/Vector2D.h"
using ROOT::Math::XYVector;

/** Curved Silicon Vertex Tracker Barrel with inactive area
 *
 * - Designed to process "vertex_barrel_curved.xml" 
 *
 * - Derived from "MPGDCylinderBarrelTracker_geo.cpp".
 * - Build-in RSU structure with four tiles
 *
 * \code
 * \endcode
 *
 *
 * @author Yann Bedfer
 */
static Ref_t create_SVTBarrelTracker(Detector& description, xml_h e, SensitiveDetector sens) {
  xml_det_t x_det = e;
  Material air    = description.air();
  int det_id      = x_det.id();
  string det_name = x_det.nameStr();
  DetElement sdet(det_name, det_id);

  // vector<Volume> volumes;
  map<string, Volume> volumes;
  typedef vector<PlacedVolume> Placements;
  map<string, Placements> sensitives;
  map<string, std::vector<VolPlane>> volplane_surfaces;

  PlacedVolume pv;
  // Set detector type flag
  dd4hep::xml::setDetectorTypeFlag(x_det, sdet);
  auto& params = DD4hepDetectorHelper::ensureExtension<dd4hep::rec::VariantParameters>(sdet);

  // Add the volume boundary material if configured
  for (xml_coll_t bmat(x_det, _Unicode(boundary_material)); bmat; ++bmat) {
    xml_comp_t x_boundary_material = bmat;
    DD4hepDetectorHelper::xmlToProtoSurfaceMaterial(x_boundary_material, params,
                                                    "boundary_material");
  }

  dd4hep::xml::Dimension dimensions(x_det.dimensions());
  Assembly assembly(det_name);
  sens.setType("tracker");
  double det_length = dimensions.length();
  double det_z0     = getAttrOrDefault(dimensions, _U(z0), 0.);
  double det_start  = det_z0 - (det_length / 2);
  double det_end    = det_z0 + (det_length / 2);

  // ********** MODULE
  struct StaveModule {
    std::string name;
    Double_t rmin;
    Double_t width = 0.0; 
    Double_t length = 0.0;
    bool   has_rsu = false; 
    Double_t uthickness = 0.0; // for the RSU
    string uname;
    string uvis;
    string umaterial;
  };
  map<string, StaveModule> staveModules;

  for (xml_coll_t mi(x_det, _U(module)); mi; ++mi) {
    xml_comp_t x_mod = mi;
    StaveModule sM;
    sM.rmin  = x_mod.rmin();
    sM.name  = x_mod.nameStr();
    sM.width = x_mod.width();
    sM.length = getAttrOrDefault(x_mod, _U(length), det_length); 

    // read RSU dimension from frame, and save in sM
    if (x_mod.hasChild(_U(frame))) {
      // xml_comp_t m_frame = x_mod.child(_U(frame));
      xml_coll_t fi(x_mod, _U(frame));
      if (fi.size() > 1) {
        printout(ERROR, "SVTBarrelTracker", "Number of frames = %d. Must be 0 or 1",
                 (int)fi.size());
        throw runtime_error("Logics error in building modules.");
      }
      xml_comp_t x_frame = fi;
      sM.uthickness      = x_frame.thickness();
      sM.has_rsu         = 1;
      sM.uname           = x_frame.nameStr();
      sM.uvis            = x_frame.visStr();
      sM.umaterial       = x_frame.materialStr();
    }
    if (staveModules.find(sM.name) != staveModules.end()) {
      printout(ERROR, "SVTBarrelTracker: ",
               string((string("staveModule with named ") + sM.name + string(" already exists.")))
                   .c_str());
      throw runtime_error("Logics error in building modules.");
    }
    staveModules[sM.name] = sM;

    // ***** TOTAL THICKNESS from components
    double total_thickness = 0;
    xml_coll_t ci(x_mod, _U(module_component));
    for (ci.reset(), total_thickness = 0.0; ci; ++ci) {
      const xml_comp_t x_comp = ci;
      double cthickness       = x_comp.thickness();
      printout(DEBUG, "SVTBarrelTracker", "\"%s\": \t comp_thickness %.4f cm",
               x_comp.nameStr().c_str(), cthickness / cm);
      total_thickness += cthickness;
    }
    if(sM.has_rsu){
      total_thickness += sM.uthickness;
      printout(DEBUG, "SVTBarrelTracker", "\"%s\": \t comp_thickness %.4f cm", sM.name.c_str(),
             sM.uthickness / cm);
    }
    printout(DEBUG, "SVTBarrelTracker", " => total_thickness %.4f cm", total_thickness / cm);

    // ***** ASSEMBLY VOLUME
    Assembly m_vol(sM.name);
    if (volumes.find(sM.name) != volumes.end()) {
      printout(ERROR, "SVTBarrelTracker: ",
               string((string("Module volumes with named ") + sM.name + string(" already exists.")))
                   .c_str());
      throw runtime_error("Logics error in building modules.");
    }
    volumes[sM.name] = m_vol;
    m_vol.setVisAttributes(description.visAttributes(x_mod.visStr()));


    // creat the total RSU volume, then subtract the tiles to get frame.
    // this method doesn't work well, has strange residuals
    // Solid rsu_frame;
    // IntersectionSolid rsu_frame;
    // if (sM.has_rsu) {
      // double dphi = sM.width / sM.rmin / 2;
      // rsu_frame=Tube(sM.rmin, sM.rmin + sM.uthickness, sM.length / 2, -dphi, dphi);
    // }
 

    // ********** LOOP OVER COMPONENTS 
    // always put RSU at innermost surface 
    double comp_rmin      = sM.rmin; 
    double thickness_so_far   = 0;
    if (sM.has_rsu){ 
      thickness_so_far+=sM.uthickness;
      comp_rmin+=sM.uthickness; 
    }
    // xml_comp_t* sensitiveComp = 0;
    int sensor_number       = 1;
    for (xml_coll_t mci(x_mod, _U(module_component)); mci; ++mci) {
      xml_comp_t x_comp     = mci;
      const string c_nam    = x_comp.nameStr();
      double comp_thickness = x_comp.thickness();
      double comp_length = sM.length; // default value for a regular component (not a frame)
      double comp_width  = sM.width; // default value for a regular component (not a frame)

      Volume c_vol;
      PlacedVolume pv_sens;

      bool is_tile=0;
      // place four tiles. Other components are set at center
      if ((c_nam == "UpperRightTile")||(c_nam == "UpperLeftTile")||(c_nam == "LowerRightTile")||(c_nam == "LowerLeftTile")){
        if (sM.has_rsu) { // ignore tiles if no frame provided
          is_tile=1;
          double rphi_lo=0, rphi_hi=0, z_left=0, z_right=0, rphi_hihi=0;
          double drphi = x_comp.width();   //tile_width along rphi
          double dz    = x_comp.length();  //tile_length along z  
          // **** hard-coded RSU design with 4 tiles, plus backbones, readout pads, biasing
          // "|" = backbone:
          //
          // | ------readout-------- | -------readout--------
          // | tile                  | tile
          // | ------biasing-------- | -------biasing--------
          // | ------biasing-------- | -------biasing--------
          // | tile                  | tile
          // | ------readout-------- | -------readout--------
          // please keep, RSU dimensions are calculated from the number below in mm:
          const double tile_width          = drphi;//9.197;
          const double tile_length         = dz;//10.773; // along the stave (z)
          const double BiasingWidth        = 0.06*mm;   // need to x2 for two sets
          // const double BiasingLength       = tile_length; // along the stave (z)
          // const double ReadoutPadsLength   = tile_length;
          const double ReadoutPadsWidth    = sM.width/2-BiasingWidth-tile_width;
          const double BackboneLength      = sM.length/2-tile_length;//0.06;
          const double BackboneWidth       = sM.width/2;

          if ((c_nam == "UpperRightTile")) {
            rphi_lo = BiasingWidth;
            rphi_hi = rphi_lo + drphi;
            rphi_hihi = sM.width/2;
            z_left  = BackboneLength;
            z_right = z_left + dz;
          } else if (c_nam == "UpperLeftTile") {
            rphi_lo = BiasingWidth;
            rphi_hihi = sM.width/2;
            rphi_hi = rphi_lo + drphi;
            z_right  = 0;
            z_left   = z_right - dz;
          } else if (c_nam == "LowerRightTile") {
            rphi_lo = -BiasingWidth;
            rphi_hi = rphi_lo - drphi; //hi: farther away from phi=0. Signed
            rphi_hihi = -sM.width/2;
            z_left  = BackboneLength;
            z_right = z_left + dz;
          } else if (c_nam == "LowerLeftTile") {
            rphi_lo = -BiasingWidth;
            rphi_hi = rphi_lo - drphi;
            rphi_hihi = -sM.width/2;
            z_right  = 0;
            z_left   = z_right - dz;
          } 
          double dphi  = drphi/sM.rmin;
          double pz    = (z_right+z_left)/2;
          double phi   = (rphi_hi+rphi_lo)/2/sM.rmin;
          Position pos(0,0,pz);// x: along R pointing out. y: along rphi, 
          RotationZYX rot(phi,0,0);
          Transform3D tr(rot, pos);
          Tube c_tube(sM.rmin, sM.rmin + sM.uthickness, dz/2, -dphi/2, dphi/2);
          c_vol   = Volume(c_nam, c_tube, description.material(x_comp.materialStr()));
          pv_sens = m_vol.placeVolume(c_vol, tr);
          
          // ** inactive areas **
          // biasing
          dphi = BiasingWidth/sM.rmin;
          Tube f_tube1(sM.rmin, sM.rmin + sM.uthickness, dz/2, -dphi/2, dphi/2);
          Volume f_vol1(c_nam+"_biasing", f_tube1, description.material(sM.umaterial));
          pv = m_vol.placeVolume(f_vol1, Transform3D(RotationZYX(rphi_lo/sM.rmin/2,0,0),pos));
          f_vol1.setVisAttributes(description, sM.uvis);

          // read out
          dphi = ReadoutPadsWidth/sM.rmin;
          Tube f_tube2(sM.rmin, sM.rmin + sM.uthickness, dz/2, -dphi/2, dphi/2);
          Volume f_vol2(c_nam+"_readout", f_tube2, description.material(sM.umaterial));
          pv = m_vol.placeVolume(f_vol2, Transform3D(RotationZYX((rphi_hihi+rphi_hi)/sM.rmin/2,0,0),pos));
          f_vol2.setVisAttributes(description, sM.uvis);

          // backbone 
          dphi = BackboneWidth/sM.rmin;
          Tube f_tube3(sM.rmin, sM.rmin + sM.uthickness, BackboneLength/2, -dphi/2, dphi/2);
          Volume f_vol3(c_nam+"_backbone", f_tube3, description.material(sM.umaterial));
          pv = m_vol.placeVolume(f_vol3, Transform3D(RotationZYX(rphi_hihi/sM.rmin/2,0,0),Position(0,0,z_left-BackboneLength/2)));
          f_vol3.setVisAttributes(description, sM.uvis);
        } 
      }
      else{ // other component, place at given thickness, and module center
        double dphi      = comp_width / comp_rmin;
        Tube c_tube(comp_rmin, comp_rmin + comp_thickness, comp_length / 2, -dphi/2, dphi/2);
        c_vol = Volume(c_nam, c_tube, description.material(x_comp.materialStr()));
        pv_sens = m_vol.placeVolume(c_vol, Transform3D(RotationZYX(0,0,0),Position(0,0,0)));
      }
      c_vol.setVisAttributes(description, x_comp.visStr());
      
      if (x_comp.isSensitive()) {
        // // ***** SENSITIVE VOLUME
        pv_sens.addPhysVolID("sensor", sensor_number++);
        c_vol.setSensitiveDetector(sens);
        sensitives[sM.name].push_back(pv_sens);

        // -------- create a measurement plane for the tracking surface attached to the sensitive volume -----
        Vector3D u(-1., 0., 0.);
        Vector3D v(0., -1., 0.);
        Vector3D n(0., 0. , 1.);

        // Compute the inner (i.e. thickness until mid-sensitive-volume) and
        //             outer (from mid-sensitive-volume to top)
        // thicknesses that need to be assigned to the tracking surface
        double inner_thickness, outer_thickness;
        if (is_tile){ // always put the tile at inner most layer
          inner_thickness = 0 + sM.uthickness/2;
          outer_thickness = total_thickness - sM.uthickness/2;
        } else{
          inner_thickness = thickness_so_far + comp_thickness/2 ;
          outer_thickness = total_thickness - inner_thickness;
        }

        SurfaceType type(SurfaceType::Sensitive);
        VolPlane surf(c_vol, type, inner_thickness, outer_thickness, u, v, n); //,o ) ;
        volplane_surfaces[sM.name].push_back(surf);
      }
      if (!is_tile){
        comp_rmin        += comp_thickness;
        thickness_so_far += comp_thickness;
      }
    } //end of module component loop
      //place the frame
    // if (sM.has_rsu){
      // Volume rsu_vol(sM.uname, rsu_frame, description.material(sM.umaterial));
      // pv       = m_vol.placeVolume(rsu_vol, Position(0, 0, 0));
      // rsu_vol.setVisAttributes(description, sM.uvis);
    // }
  } //end of module loop
   
  // ********** LAYER
  // ***** RETRIEVE PARAMETERS
  for (xml_coll_t li(x_det, _U(layer)); li; ++li) {
    xml_comp_t x_layer   = li;
    int lay_id           = x_layer.id();
    string m_nam         = x_layer.moduleStr();
    Volume& module_vol   = volumes[m_nam];
    Placements& sensVols = sensitives[m_nam];

    xml_comp_t x_barrel  = x_layer.child(_U(barrel_envelope));
    double barrel_length = x_barrel.z_length();
    double barrel_z0     = getAttrOrDefault(x_barrel, _U(z0), 0.);
    // Calculate barrel ends
    double barrel_start = barrel_z0 - (barrel_length / 2);
    double barrel_end   = barrel_z0 + (barrel_length / 2);
    if (det_start > barrel_start || det_end < barrel_end) {
      printout(ERROR, "SVTBarrelTracker",
               "Layer #%d is outside the detector envelop\n"
               "det_start < barrel_start:  %.4f < %.4f\n "
               "det_end   > barrel_end  :  %.4f > %.4f\n  ",
               det_start, barrel_start, det_end, barrel_end);
      throw runtime_error("Logics error in building modules.");
    }

    // ***** LAYOUTS
    xml_comp_t x_layout = x_layer.child(_U(rphi_layout));
    double phi0         = x_layout.phi0();             // Starting phi of first module.
    xml_comp_t z_layout = x_layer.child(_U(z_layout)); // Get the <z_layout> element.
    double layer_zstart = getAttrOrDefault<double>(z_layout, _U(zstart), barrel_start);
    if (barrel_start > layer_zstart) {
      printout(ERROR, "SVTBarrelTracker", "Layer #%d: barrel_start>layer_zstart: %.4f > %.4f\n ",
               barrel_start, layer_zstart);
      throw runtime_error("Logics error in building modules.");
    }

    // ********** LAYER
    // ***** ENVELOPE
    string lay_nam = det_name + _toString(x_layer.id(), "_layer%d");
    Tube lay_tub(x_barrel.inner_r(), x_barrel.outer_r(), barrel_length / 2);
    Volume lay_vol(lay_nam, lay_tub, air); // Create the layer envelope volume.
    Position lay_pos(0, 0, barrel_z0);
    lay_vol.setVisAttributes(description.visAttributes(x_layer.visStr()));
    printout(DEBUG, "SVTBarrelTracker", "Layer \"%s\": rmin,max = %.2f,%.2f cm 1/2length = %.2f cm",
             lay_nam.c_str(), x_barrel.inner_r(), x_barrel.outer_r(), barrel_length / 2);

    DetElement lay_elt(sdet, lay_nam, lay_id);

    // the local coordinate systems of modules in dd4hep and acts differ
    // see http://acts.web.cern.ch/ACTS/latest/doc/group__DD4hepPlugins.html
    auto& layerParams =
        DD4hepDetectorHelper::ensureExtension<dd4hep::rec::VariantParameters>(lay_elt);
    for (xml_coll_t lmat(x_layer, _Unicode(layer_material)); lmat; ++lmat) {
      xml_comp_t x_layer_material = lmat;
      DD4hepDetectorHelper::xmlToProtoSurfaceMaterial(x_layer_material, layerParams,
                                                      "layer_material");
    }
    // ********** LOOP OVER THE SECTORS IN z and phi
    auto sM = staveModules[m_nam];
    // Determine "nphi" from RSU width and radius
    int nphi = int(2 * M_PI * sM.rmin / sM.width); // always round down to avoid overlap

    // ***** SECTOR POSITIONS ALONG Z
    // double modz_pos[1] = {-barrel_length / 2.0};
    int nz          = int((barrel_end - layer_zstart) / sM.length);
    double zstart   = layer_zstart; // for modules in z
    double z_incr   = sM.length;
    double module_z = zstart + z_incr / 2; // module center in z
    int nModules    = 0;

    for (int iz = 0; iz < nz; iz++) {
      double phi_incr = 2 * M_PI / nphi; // Phi increment for one module.
      // ***** LOOP OVER THE STAVES IN phi.
      int iphi;
      double phic;
      for (iphi = 0, phic = phi0; iphi < nphi; iphi++, phic += phi_incr, nModules++) {
        double x1, y1; // Coordinates of the centre of curvature of
        x1 = 0;        //rc * std::cos(phic);
        y1 = 0;        //rc * std::sin(phic);
        // int module_id       = 100 * iz + iphi;
        string module_name = _toString(nModules, "module%d");
        DetElement mod_elt(lay_elt, module_name, nModules);
        RotationZYX rot;
        rot = RotationZYX(phic, 0, 0);
        Transform3D tr(rot, Position(x1, y1, module_z));
        pv = lay_vol.placeVolume(module_vol, tr);
        pv.addPhysVolID("module", nModules);
        mod_elt.setPlacement(pv);

        // ***** SENSITIVE COMPONENT
        for (size_t ic = 0; ic < sensVols.size(); ++ic) {
          PlacedVolume sens_pv = sensVols[ic];
          DetElement comp_de(mod_elt, std::string("de_") + sens_pv.volume().name(), nModules);
          comp_de.setPlacement(sens_pv);
          auto& comp_de_params =
              DD4hepDetectorHelper::ensureExtension<dd4hep::rec::VariantParameters>(comp_de);
          comp_de_params.set<string>("axis_definitions", "XYZ");
          volSurfaceList(comp_de)->push_back(volplane_surfaces[m_nam][ic]);
        }
      }
      module_z += z_incr;
    }

    // ***** CREATE THE PhysicalVolume FOR THE LAYER.
    pv = assembly.placeVolume(lay_vol, lay_pos); // Place layer in mother
    pv.addPhysVolID("layer", lay_id);            // Set the layer ID.
    lay_elt.setAttributes(description, lay_vol, x_layer.regionStr(), x_layer.limitsStr(),
                          x_layer.visStr());
    lay_elt.setPlacement(pv);
  }
  sdet.setAttributes(description, assembly, x_det.regionStr(), x_det.limitsStr(), x_det.visStr());
  assembly.setVisAttributes(description.invisible());
  pv = description.pickMotherVolume(sdet).placeVolume(assembly);
  pv.addPhysVolID("system", det_id); // Set the subdetector system ID.
  sdet.setPlacement(pv);

  // #ifdef DEBUG_SVTBarrelTracker
  // // Reset initial print level before exiting
  // setPrintLevel(priorPrintLevel);
  // #endif

  return sdet;
}

//@}
// clang-format off
DECLARE_DETELEMENT(epic_CylinderSVTBarrel, create_SVTBarrelTracker)
