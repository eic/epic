//==========================================================================
//  A general implementation for homogeneous calorimeter
//  it supports three types of placements
//  1. Individual module placement with module dimensions and positions
//  2. Array placement with module dimensions and numbers of row and columns
//  3. Disk placement with module dimensions and (Rmin, Rmax), and (Phimin, Phimax)
//  4. Lines placement with module dimensions and (mirrorx, mirrory)
//     (NOTE: anchor point is the 0th block of the line instead of line center)
//--------------------------------------------------------------------------
//  Author: Chao Peng (ANL)
//  Date: 06/09/2021
//==========================================================================

#include "GeometryHelpers.h"
#include "DD4hep/DetFactoryHelper.h"
#include <XML/Helper.h>
#include <iostream>
#include <algorithm>
#include <tuple>
#include <math.h>
#include <fmt/core.h>

using namespace dd4hep;

// main
static Ref_t create_detector(Detector& desc, xml::Handle_t handle, SensitiveDetector sens)
{

  using namespace std;
  using namespace fmt;

    xml::DetElement detElem = handle;
    std::string detName = detElem.nameStr();
    int detID = detElem.id();
    DetElement det(detName, detID);
    sens.setType("calorimeter");

    auto glass_material = desc.material("SciGlass");
    auto crystal_material = desc.material("PbWO4");
    auto air_material = desc.material("Air");

    double ROut = desc.constantAsDouble("EcalEndcapN_rmax");
    double RIn_el = desc.constantAsDouble("EcalEndcapN_rmin1");
    double ionCutout_dphi = desc.constantAsDouble("EcalEndcapNIonCutout_dphi");
    double RIn = desc.constantAsDouble("EcalEndcapN_rmin2");
    double SizeZ = desc.constantAsDouble("EcalEndcapN_thickness");
    double thickness = desc.constantAsDouble("EcalEndcapN_thickness");
    double trans_radius = desc.constantAsDouble("EcalEndcapNCrystal_rmax");
    double Glass_z0 = desc.constantAsDouble("GlassModule_z0");
    double Glass_Width = desc.constantAsDouble("GlassModule_width");
    double Glass_thickness = desc.constantAsDouble("GlassModule_length");
    double Glass_Gap = desc.constantAsDouble("GlassModule_wrap");
    double glass_distance = desc.constantAsDouble("GlassModule_distance");

    double Crystal_Width = desc.constantAsDouble("CrystalModule_width");
    double Crystal_thickness = desc.constantAsDouble("CrystalModule_length");
    double Crystal_Gap = desc.constantAsDouble("CrystalModule_wrap");
    double crystal_distance = desc.constantAsDouble("CrystalModule_distance");
    double Crystal_z0 = desc.constantAsDouble("CrystalModule_z0");

    // RIn and ROut will define outer tube embedding the calorimeter
    // centers_rmin/out define the maximum radius of module centers
    // so that modules are not overlapping with mother tube volume
    const double glassHypotenuse = std::hypot(glass_distance, glass_distance)/2;
    const double crystalHypotenuse = glassHypotenuse/2;
    // Offset these values a bit so we don't miss edge-blocks
    const double glassCenters_rmax  = ROut - glassHypotenuse + 1 * mm;
    const double crystalCenters_rmin = RIn + crystalHypotenuse - 1 * mm;
    // Also limits of the inner crystal blocks fill
    const double cutout_tan = tan(ionCutout_dphi/2);
    const double cutout_rmin = RIn_el + crystalHypotenuse - 1 * mm;

    // Offset to align the modules at the zmin of the endcap,
    const double Crystal_offset = -0.5 * (Crystal_thickness - thickness);
    const double Glass_offset = -0.5*(Glass_thickness - thickness);

    // envelope
    // consists of an glass tube of the full thickness, and a crystal inner tube
    // for the crystal that allows us to get closet to the beampipe without
    // overlaps, and a partial electron tube that allows us to get closer to the
    // electron beampipe in the region where there is no ion beampipe
    Tube glass_tube(min(RIn + glassHypotenuse*2, trans_radius), ROut, SizeZ / 2.0, 0., 360.0 * deg);
    Tube crystal_tube(RIn, min(RIn + glassHypotenuse*2, trans_radius), Crystal_thickness/ 2.0, 0., 360.0 * deg);
    Tube electron_tube(RIn_el, RIn, Crystal_thickness / 2., ionCutout_dphi / 2., 360.0 * deg - ionCutout_dphi / 2);
    UnionSolid outer_envelope(glass_tube, crystal_tube, Position(0,0,Crystal_offset));
    UnionSolid envelope(outer_envelope, electron_tube, Position(0,0,Crystal_offset));
    Volume ecal_vol("negative_ecal", envelope, air_material);
    ecal_vol.setVisAttributes(desc.visAttributes("HybridEcalOuterVis"));

    // TODO why 1cm and not something else?
    double Glass_OuterR = ROut - 1 * cm ;
    double Glass_InnerR = trans_radius;

    // Geometry of modules
    Box glass_box("glass_box", Glass_Width * 0.5, Glass_Width * 0.5, Glass_thickness * 0.5);
    Volume glass_module("glass_module", glass_box, glass_material);
    glass_module.setVisAttributes(desc.visAttributes("EcalEndcapNModuleVis"));
    glass_module.setSensitiveDetector(sens);
    
    Box crystal_box("crystal_box",  Crystal_Width* 0.5, Crystal_Width * 0.5, Crystal_thickness * 0.5);
    Volume crystal_module("crystal_module", crystal_box, crystal_material);
    crystal_module.setVisAttributes(desc.visAttributes("EcalEndcapNModuleVis"));
    crystal_module.setSensitiveDetector(sens);

    // GLASS
    double diameter = 2 * Glass_OuterR;

    // Can we fit an even or odd amount of glass blocks within our rmax?
    // This determines the transition points between crystal and glass as we need the
    // outer crystal arangement to be in multiples of 2 (aligned with glass)
    const int towersInRow = std::ceil((diameter + Glass_Gap) /  (Glass_Width + Glass_Gap));
    const bool align_even = (towersInRow % 2);

    // Is it odd or even number of towersInRow
    double leftTowerPos, topTowerPos;
    if(towersInRow%2) {
      //             |
      //      [ ][ ][ ][ ][ ]
      //       ^     |
      int towersInHalfRow = std::ceil(towersInRow/2.0);
      topTowerPos = leftTowerPos = -towersInHalfRow * (Glass_Width + Glass_Gap);

    } else {
      //               |
      //      [ ][ ][ ][ ][ ][ ]
      //       ^      |
      int towersInHalfRow = towersInRow/2;
      topTowerPos = leftTowerPos = -(towersInHalfRow - 0.5) * (Glass_Width + Glass_Gap);
    }

    int moduleIndex = 0;

    int glass_module_index = 0;
    int cryst_module_index = 0;
    for(int rowIndex=0; rowIndex < towersInRow; rowIndex++) {
      for(int colIndex=0; colIndex < towersInRow; colIndex++) {
        const double glass_x = leftTowerPos + colIndex * (Glass_Width + Glass_Gap);
        const double glass_y = topTowerPos + rowIndex * (Glass_Width + Glass_Gap);
        const double r = std::hypot(glass_x, glass_y);
        // crystal if within the transition radius (as defined by the equivalent glass
        // block)
        if (r < trans_radius) {
          for (const auto dx : {-1, 1}) {
            for (const auto& dy : {-1, 1}) {
              const double crystal_x = glass_x + dx * crystal_distance / 2;
              const double crystal_y = glass_y + dy * crystal_distance / 2;
              const double crystal_r = std::hypot(crystal_x, crystal_y);
              // check if crystal in the main crystal ring?
              const bool mainRing = (crystal_r > crystalCenters_rmin);
              const bool innerRing = !mainRing && crystal_r > cutout_rmin;
              const bool ionCutout = crystal_x > 0 && fabs(crystal_y / crystal_x) < cutout_tan;
              if (mainRing || (innerRing && !ionCutout)) {
                auto placement =
                    ecal_vol.placeVolume(crystal_module, Position(crystal_x, crystal_y, Crystal_z0 + Crystal_offset));
                placement.addPhysVolID("sector", 1);
                placement.addPhysVolID("module", cryst_module_index++);
              }
            }
          }
        // Glass block if within the rmax
        } else if (r < glassCenters_rmax) {
          // glass module
          auto placement = ecal_vol.placeVolume(glass_module, Position(glass_x, glass_y, Glass_z0 + Glass_offset));
          placement.addPhysVolID("sector", 2);
          placement.addPhysVolID("module", glass_module_index++);
        }
      }
    }

    desc.add(Constant("EcalEndcapN_NModules_Sector1", std::to_string(cryst_module_index)));
    desc.add(Constant("EcalEndcapN_NModules_Sector2", std::to_string(glass_module_index)));
//    fmt::print("Total Glass modules: {}\n", towerIndex);
//    fmt::print("CE EMCAL GLASS END\n\n");

     // detector position and rotation
     auto pos = detElem.position();
     auto rot = detElem.rotation();
     Volume motherVol = desc.pickMotherVolume(det);
     Transform3D tr(RotationZYX(rot.z(), rot.y(), rot.x()), Position(pos.x(), pos.y(), pos.z()));
     PlacedVolume envPV = motherVol.placeVolume(ecal_vol, tr);
     envPV.addPhysVolID("system", detID);
     det.setPlacement(envPV);
     return det;

}

//@}
DECLARE_DETELEMENT(HybridCalorimeter, create_detector)

