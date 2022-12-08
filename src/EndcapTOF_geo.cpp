// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Nicolas Schmidt

/** \addtogroup Trackers Trackers
 * \brief Type: **Endcap Tracker with TOF**.
 * \author N. Schmidt
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
#include <map>
#include "DD4hepDetectorHelper.h"

using namespace dd4hep;
using namespace dd4hep::rec;
using namespace dd4hep::detail;

static Ref_t create_detector(Detector& description, xml_h e, SensitiveDetector sens)
{

  xml_det_t    x_det    = e;
  int          det_id   = x_det.id();
  std::string  det_name = x_det.nameStr();
  DetElement   ttl_detEl(det_name, det_id);
  Material     air = description.material("Air");
  PlacedVolume pv;

  // Set detector type flag
  dd4hep::xml::setDetectorTypeFlag(x_det, sdet);
  auto &params = DD4hepDetectorHelper::ensureExtension<dd4hep::rec::VariantParameters>(
      sdet);

  // Add the volume boundary material if configured
  for (xml_coll_t bmat(x_det, _Unicode(boundary_material)); bmat; ++bmat) {
    xml_comp_t x_boundary_material = bmat;
    DD4hepDetectorHelper::xmlToProtoSurfaceMaterial(x_boundary_material, params,
                                         "boundary_material");
  }

  xml_comp_t  x_mod         = x_det.child(_Unicode(module));
  std::string m_nam         = x_mod.nameStr();
  xml_comp_t  diskdimension = x_mod.dimensions();

  // load all information from the diskdimension definitions
  double disk_zPos              = getAttrOrDefault(diskdimension, _Unicode(zPos), 0.);
  double disk_rMin              = getAttrOrDefault(diskdimension, _Unicode(rMin), 0.);
  double disk_rMax              = getAttrOrDefault(diskdimension, _Unicode(rMax), 100.);
  double disk_xOffset           = getAttrOrDefault(diskdimension, _Unicode(xOffset), 0.);
  double disk_det_height        = getAttrOrDefault(diskdimension, _Unicode(det_height), 0.);
  double cooling_plate_height   = getAttrOrDefault(diskdimension, _Unicode(cooling_plate_height), 0.);
  double cooling_tube_thickness = getAttrOrDefault(diskdimension, _Unicode(wallthickness_coolingtube), 0.);
  double cooling_tube_diameter  = getAttrOrDefault(diskdimension, _Unicode(diameter_coolingtube), 0.);

  Assembly assembly(det_name);
  assembly.setVisAttributes(description.invisible());
  sens.setType("tracker");

  double disk_face_rel_z = cooling_tube_diameter / 2. + cooling_plate_height / 2.;
  // NOTE: Creation of aluminum disks with offset cutout for beampipe
  {
    // create a solid for the beampipe cutout
    Tube beampipe_cutout(0, disk_rMin, disk_det_height, 0, 2. * M_PI);
    // create solids for the cooling/mountin disks (two Al disks)
    Tube ttl_plate_1(0, disk_rMax, cooling_plate_height / 2., 0, 2. * M_PI);
    Tube ttl_plate_2(0, disk_rMax, cooling_plate_height / 2., 0, 2. * M_PI);

    // create the offset cutout in the disks
    Position         cutoutOffset(disk_xOffset, 0, 0);
    SubtractionSolid ttl_cooling_plate_1(ttl_plate_1, beampipe_cutout, cutoutOffset);
    SubtractionSolid ttl_cooling_plate_2(ttl_plate_2, beampipe_cutout, cutoutOffset);

    // set the material (Al)
    Material disk_mat = description.material(diskdimension.materialStr());

    // create the volumes and set visualization options
    Volume ttl_cooling_plate_volume_1("ttl_cooling_plate_1", ttl_cooling_plate_1, disk_mat);
    Volume ttl_cooling_plate_volume_2("ttl_cooling_plate_2", ttl_cooling_plate_2, disk_mat);
    ttl_cooling_plate_volume_1.setVisAttributes(description.visAttributes("TOFAluminum"));
    ttl_cooling_plate_volume_2.setVisAttributes(description.visAttributes("TOFAluminum"));

    // place disks in assembly
    pv = assembly.placeVolume(ttl_cooling_plate_volume_1, Position(0, 0, -disk_face_rel_z));
    pv = assembly.placeVolume(ttl_cooling_plate_volume_2, Position(0, 0, disk_face_rel_z));
  }
  disk_face_rel_z += cooling_plate_height / 2.;

  xml_comp_t sensparams = x_mod.child(_Unicode(sensorparameters));
  ;
  double sensor_ydim      = 21.2 * mm; // getAttrOrDefault(sensparams, _Unicode(width), 0.);
  double sensor_xdim      = 42.0 * mm; // getAttrOrDefault(sensparams, _Unicode(length), 0.);
  double sensor_thickness = 0.3 * mm;  // getAttrOrDefault(sensparams, _Unicode(thickness), 0.);
  double sensor_margin    = 0.15 * mm; // getAttrOrDefault(sensparams, _Unicode(margin), 0.);
  double baseplate_xdim   = getAttrOrDefault(sensparams, _Unicode(baseplate_length), 100.);
  double baseplate_ydim   = getAttrOrDefault(sensparams, _Unicode(baseplate_width), 0.);
  printout(DEBUG, "EndcapTOF", "sensor_xdim: %f\tsensor_ydim: %f", sensor_xdim, sensor_ydim);
  double sensor_hybrid_ydim = baseplate_ydim + baseplate_ydim / 2.;

  std::string strLayerName[10];
  Material    materialLayer[10];
  double      thicknessLayer[10];
  double      widthLayer[10];
  double      offsetLayer[10];
  double      lengthLayer[10];
  int         nLayers = 0;
  for (xml_coll_t l_iter(x_det, _U(layer)); l_iter; ++l_iter) {
    xml_comp_t x_layer = l_iter;
    if (x_layer.id() != 2)
      continue;

    nLayers = getAttrOrDefault(x_layer, _Unicode(numslices), 0.);

    int i_slice = 0;
    for (xml_coll_t s_iter(x_layer, _U(slice)); s_iter; ++s_iter, ++i_slice) {
      // If slices are only given a thickness attribute, they are radially concentric slices
      // If slices are given an inner_z attribute, they are longitudinal slices with equal rmin
      xml_comp_t x_slice      = s_iter;
      materialLayer[i_slice]  = description.material(x_slice.materialStr());
      strLayerName[i_slice]   = getAttrOrDefault<std::string>(x_slice, _U(name), "slice" + std::to_string(i_slice));
      thicknessLayer[i_slice] = x_slice.thickness();
      widthLayer[i_slice]     = getAttrOrDefault(x_slice, _Unicode(width), 0.);
      offsetLayer[i_slice]    = getAttrOrDefault(x_slice, _Unicode(offset), 0.);
      lengthLayer[i_slice]    = getAttrOrDefault(x_slice, _Unicode(length), 0.);
    }
  }

  // need to define two thicknesses here as the sensors have to be placed separately from the rest of the stack
  // this is done to properly assign sensor IDs to the placed volumes
  double thicknessDet1 = 0;
  double thicknessDet2 = 0;
  for (int ilay = 0; ilay < nLayers; ilay++) {
    // layer 5 would be the sensor
    if (ilay < 5)
      thicknessDet1 += thicknessLayer[ilay];
    if (ilay > 5)
      thicknessDet2 += thicknessLayer[ilay];
  }

  // ---------------------------------------
  // Create Service Hybrid Module Stack:
  // ---------------------------------------
  std::string strLayerName_SH[10];
  Material    materialLayer_SH[10];
  double      thicknessLayer_SH[10];
  double      widthLayer_SH[10];
  double      offsetLayer_SH[10];
  double      lengthLayer_SH[10];
  int         nLayers_SH = 0;
  for (xml_coll_t l_iter(x_det, _U(layer)); l_iter; ++l_iter) {
    xml_comp_t x_layer = l_iter;
    if (x_layer.id() != 3)
      continue;

    nLayers_SH = getAttrOrDefault(x_layer, _Unicode(numslices), 0.);

    int i_slice = 0;
    for (xml_coll_t s_iter(x_layer, _U(slice)); s_iter; ++s_iter, ++i_slice) {
      // If slices are only given a thickness attribute, they are radially concentric slices
      // If slices are given an inner_z attribute, they are longitudinal slices with equal rmin
      xml_comp_t x_slice         = s_iter;
      materialLayer_SH[i_slice]  = description.material(x_slice.materialStr());
      strLayerName_SH[i_slice]   = getAttrOrDefault<std::string>(x_slice, _U(name), "slice" + std::to_string(i_slice));
      thicknessLayer_SH[i_slice] = x_slice.thickness();
      widthLayer_SH[i_slice]     = getAttrOrDefault(x_slice, _Unicode(width), 0.);
      offsetLayer_SH[i_slice]    = getAttrOrDefault(x_slice, _Unicode(offset), 0.);
      lengthLayer_SH[i_slice]    = getAttrOrDefault(x_slice, _Unicode(length), 0.);
    }
  }

  int      layer_id = 0;
  Assembly layer_assembly("layer_assembly" + _toString(layer_id, "_%d"));

  DetElement layer_detEl(ttl_detEl, _toString(layer_id, "layer_%d"), layer_id);

  // NOTE: ACTS extension for the disk layer of the TTL
  // also defining the coordinate system that differs between ACTS and Geant4 (zyx vs xyz)
  auto &layerParams = DD4hepDetectorHelper::ensureExtension<dd4hep::rec::VariantParameters>(
                      layer_detEl);
  layerParams.set<double>("envelope_r_min", -80. * mm);
  layerParams.set<double>("envelope_r_max", 670. * mm);
  layerParams.set<double>("envelope_z_min", 10. * mm);
  layerParams.set<double>("envelope_z_max", 10. * mm);

  pv = assembly.placeVolume(layer_assembly);
  pv.addPhysVolID("layer", layer_id);

  layer_detEl.setPlacement(pv);

  // NOTE: create envelopes for the sensor supporting layers (glue, kapton, etc.)
  Box    box_sensor_stack1(baseplate_xdim / 2, baseplate_ydim / 2, thicknessDet1 / 2);
  Box    box_sensor_stack2(baseplate_xdim / 2, baseplate_ydim / 2, thicknessDet2 / 2);
  Volume log_sensor_stack1("log_sensor_stack1", box_sensor_stack1, air);
  Volume log_sensor_stack2("log_sensor_stack2", box_sensor_stack2, air);
  log_sensor_stack1.setVisAttributes(description.visAttributes("InvisibleWithDaughters"));
  log_sensor_stack2.setVisAttributes(description.visAttributes("InvisibleWithDaughters"));

  // individual layers will be placed in subsequent order and z_start is increased by the thickness of the previous
  // layer
  double z_start1 = -thicknessDet1 / 2;
  double z_start2 = -thicknessDet2 / 2;
  for (int ilay = 0; ilay < nLayers; ilay++) {
    const std::string layer_name_sns = "sensor_stack_" + strLayerName[ilay];

    Box    box_Module_Layer_Raw(lengthLayer[ilay] / 2, widthLayer[ilay] / 2, thicknessLayer[ilay] / 2);
    Solid  sol_Module_Layer_Raw = box_Module_Layer_Raw;
    Volume Log_Layer(layer_name_sns + "_Log", sol_Module_Layer_Raw, materialLayer[ilay]);

    if (ilay < 5) {
      pv = log_sensor_stack1.placeVolume(Log_Layer,
                                         Position(0, -offsetLayer[ilay], z_start1 + thicknessLayer[ilay] / 2));
      z_start1 += thicknessLayer[ilay];
    }
    if (ilay > 5) {
      pv = log_sensor_stack2.placeVolume(Log_Layer,
                                         Position(0, -offsetLayer[ilay], z_start2 + thicknessLayer[ilay] / 2));
      z_start2 += thicknessLayer[ilay];
    }
    Log_Layer.setVisAttributes(description.visAttributes("TOFLayers"));
  }

  double baseSH_width    = baseplate_ydim / 2;
  double thicknessDet_SH = 0;
  for (int ilay = 0; ilay < nLayers_SH; ilay++) {
    thicknessDet_SH += thicknessLayer_SH[ilay];
  }

  Box    sol_SH_stack(baseplate_xdim / 2, baseSH_width / 2, thicknessDet_SH / 2);
  Volume log_SH_stack("log_SH_stack", sol_SH_stack, air);
  log_SH_stack.setVisAttributes(description.visAttributes("TOFSensorAndReadoutLadder"));

  // place all layers into the SH stack
  double z_start_SH = -thicknessDet_SH / 2;
  for (int ilay = 0; ilay < nLayers_SH; ilay++) {
    const std::string layer_name       = "SH_stack_" + strLayerName_SH[ilay];
    const std::string layer_name_Solid = "sol_" + layer_name;

    Box    sol_Module_Layer_Raw(lengthLayer_SH[ilay] / 2, widthLayer_SH[ilay] / 2, thicknessLayer_SH[ilay] / 2);
    Volume Log_Layer(layer_name + "_Log", sol_Module_Layer_Raw, materialLayer_SH[ilay]);

    pv = log_SH_stack.placeVolume(Log_Layer,
                                  Position(0, -offsetLayer_SH[ilay], z_start_SH + thicknessLayer_SH[ilay] / 2));

    z_start_SH += thicknessLayer_SH[ilay];
    Log_Layer.setVisAttributes(description.visAttributes("TOFLayers"));
  }

  // NOTE: create sensor volume itself and assign a surface
  int         senslayer         = 5;
  std::string sensor_layer_name = det_name + "_sensor_layer_" + std::to_string(layer_id);
  Box         sensor_box(lengthLayer[senslayer] / 2, widthLayer[senslayer] / 2, thicknessLayer[senslayer] / 2);
  Volume      sens_vol(sensor_layer_name, sensor_box, materialLayer[senslayer]);
  sens_vol.setAttributes(description, x_det.regionStr(), x_det.limitsStr(), "TOFActiveMat");
  sens_vol.setSensitiveDetector(sens);
  // sens_vol.setVisAttributes(description.visAttributes("TOFAluminum"));

  // Definition of the sensitive surface of the sensor
  // the sensitive volume -----
  Vector3D u(0., 0., -1.);
  Vector3D v(-1., 0., 0.);
  Vector3D n(0., 1., 0.);

  // compute the inner and outer thicknesses that need to be assigned to the
  // tracking surface
  // depending on wether the support is above or below the sensor
  double inner_thickness = sensor_thickness / 2.0;
  double outer_thickness = disk_det_height - sensor_thickness / 2.0;

  SurfaceType type(SurfaceType::Sensitive);

  // if (isStripDetector) type.setProperty(SurfaceType::Measurement1D, true);

  VolPlane surf_front(sens_vol, type, inner_thickness, outer_thickness, u, v, n);
  VolPlane surf_back(sens_vol, type, outer_thickness, inner_thickness, u, v, n);

  // NOTE: create cooling tube and coolant to be placed with each sensor+SH element
  std::string cooling_tube_name = det_name + "_cooling_tube_" + std::to_string(layer_id);
  Tube        sol_cooling_tube(cooling_tube_diameter / 2.0 - cooling_tube_thickness, cooling_tube_diameter / 2.0,
                               baseplate_xdim / 2.0);
  Volume      log_cooling_tube(cooling_tube_name, sol_cooling_tube, description.material("Aluminum"));
  log_cooling_tube.setVisAttributes(description.visAttributes("TOFAluminum"));

  std::string cooling_liquid_name = det_name + "_cooling_liquid_" + std::to_string(layer_id);
  Tube        sol_cooling_liquid(0, cooling_tube_diameter / 2.0 - cooling_tube_thickness, baseplate_xdim / 2.0);
  Volume      log_cooling_liquid(cooling_liquid_name, sol_cooling_liquid, description.material("Water"));
  log_cooling_liquid.setVisAttributes(description.visAttributes("TOFWater"));

  int isensor = 0;
  // number of sensors+halfSH in y-radial direction
  int nSensY = (int)((disk_rMax - (sensor_hybrid_ydim / 2)) / sensor_hybrid_ydim);

  for (int ilay_y = nSensY; ilay_y >= -nSensY; ilay_y--) {
    int nSensX =
        (int)((2 * sqrt(pow(disk_rMax, 2) - pow((abs(ilay_y) * sensor_hybrid_ydim + sensor_hybrid_ydim / 2), 2))) /
              baseplate_xdim);

    if (nSensX == 0)
      continue;
    // we want an odd number of towers to be symmetrically centered around 0
    if (nSensX % 2 == 0)
      nSensX -= 1;

    for (int ilay_x = -(nSensX - 1) / 2; ilay_x <= (nSensX - 1) / 2; ++ilay_x) {
      if (ilay_x < 0) {
        if (ilay_y < 0) {
          if (sqrt(pow(ilay_x * baseplate_xdim + baseplate_xdim / 2 - disk_xOffset, 2) +
                   pow(ilay_y * sensor_hybrid_ydim + sensor_hybrid_ydim / 2.0, 2)) < disk_rMin)
            continue;
        } else {
          if (sqrt(pow(ilay_x * baseplate_xdim + baseplate_xdim / 2 - disk_xOffset, 2) +
                   pow(ilay_y * sensor_hybrid_ydim - sensor_hybrid_ydim / 2.0, 2)) < disk_rMin)
            continue;
        }
      } else {
        if (ilay_y < 0) {
          if (sqrt(pow(ilay_x * baseplate_xdim - baseplate_xdim / 2 - disk_xOffset, 2) +
                   pow(ilay_y * sensor_hybrid_ydim + sensor_hybrid_ydim / 2.0, 2)) < disk_rMin)
            continue;
        } else {
          if (sqrt(pow(ilay_x * baseplate_xdim - baseplate_xdim / 2 - disk_xOffset, 2) +
                   pow(ilay_y * sensor_hybrid_ydim - sensor_hybrid_ydim / 2.0, 2)) < disk_rMin)
            continue;
        }
      }
      // NOTE: placement of sensors on front side of disk
      std::string sens_name_front = det_name + "_sensor_" + std::to_string(ilay_x) + "_" + std::to_string(ilay_y) +
                                    "_front_" + std::to_string(isensor);
      RotationZYX rot_front(0, 0, 0);

      pv = layer_assembly.placeVolume(
          sens_vol,
          Transform3D(rot_front, Position((ilay_x)*baseplate_xdim,
                                          ilay_y % 2 == 0 ? (ilay_y)*sensor_hybrid_ydim - sensor_hybrid_ydim / 2 +
                                                                sensor_ydim / 2 + sensor_margin
                                                          : (ilay_y)*sensor_hybrid_ydim + sensor_hybrid_ydim / 2 -
                                                                sensor_ydim / 2 - sensor_margin,
                                          -disk_face_rel_z - thicknessDet1 - sensor_thickness / 2.0)));

      pv.addPhysVolID("module", layer_id).addPhysVolID("sensor", isensor);

      DetElement sensor_detEl_front(layer_detEl, sens_name_front, isensor);
      sensor_detEl_front.setPlacement(pv);
      volSurfaceList(sensor_detEl_front)->push_back(surf_front);

      auto &params = DD4hepDetectorHelper::ensureExtension<dd4hep::rec::VariantParameters>(sensor_detEl_front);
      params.set<string>("axis_definitions", "XZY");

      isensor++;

      // NOTE: placement of sensor supports on front side of disk
      {
        RotationZYX rot_front_flip(M_PI, M_PI, 0);
        RotationZYX rot_front_flip2(0, M_PI, 0);

        pv = layer_assembly.placeVolume(
            log_sensor_stack1,
            Transform3D(ilay_y % 2 == 0 ? rot_front_flip2 : rot_front_flip,
                        Position((ilay_x)*baseplate_xdim,
                                 ilay_y % 2 == 0
                                     ? (ilay_y)*sensor_hybrid_ydim - sensor_hybrid_ydim / 2 + baseplate_ydim / 2
                                     : (ilay_y)*sensor_hybrid_ydim + sensor_hybrid_ydim / 2 - baseplate_ydim / 2,
                                 -disk_face_rel_z - thicknessDet1 / 2.0)));

        pv = layer_assembly.placeVolume(
            log_sensor_stack2,
            Transform3D(ilay_y % 2 == 0 ? rot_front_flip2 : rot_front_flip,
                        Position((ilay_x)*baseplate_xdim,
                                 ilay_y % 2 == 0
                                     ? (ilay_y)*sensor_hybrid_ydim - sensor_hybrid_ydim / 2 + baseplate_ydim / 2
                                     : (ilay_y)*sensor_hybrid_ydim + sensor_hybrid_ydim / 2 - baseplate_ydim / 2,
                                 -disk_face_rel_z - thicknessDet1 - sensor_thickness - thicknessDet2 / 2.0)));
      }
      // NOTE: placement of service hybrids on front side
      pv = layer_assembly.placeVolume(
          log_SH_stack,
          Transform3D(rot_front,
                      Position((ilay_x)*baseplate_xdim,
                               ilay_y % 2 == 0
                                   ? (ilay_y)*sensor_hybrid_ydim + sensor_hybrid_ydim / 2 - baseplate_ydim / 4
                                   : (ilay_y)*sensor_hybrid_ydim - sensor_hybrid_ydim / 2 + baseplate_ydim / 4,
                               -disk_face_rel_z - thicknessDet_SH / 2)));

      // NOTE: placement of sensors on back side of disk
      std::string sens_name_back = det_name + "_sensor_" + std::to_string(ilay_x) + "_" + std::to_string(ilay_y) +
                                   "_back_" + std::to_string(isensor);
      RotationZYX rot_back(0, M_PI, 0);

      pv = layer_assembly.placeVolume(
          sens_vol,
          Transform3D(rot_back, Position((ilay_x)*baseplate_xdim,
                                         ilay_y % 2 == 0 ? (ilay_y)*sensor_hybrid_ydim + sensor_hybrid_ydim / 2 -
                                                               sensor_ydim / 2 - sensor_margin
                                                         : (ilay_y)*sensor_hybrid_ydim - sensor_hybrid_ydim / 2 +
                                                               sensor_ydim / 2 + sensor_margin,
                                         disk_face_rel_z + thicknessDet1 + sensor_thickness / 2.0)));

      pv.addPhysVolID("module", layer_id).addPhysVolID("sensor", isensor);

      DetElement sensor_detEl_back(layer_detEl, sens_name_back, isensor);
      sensor_detEl_back.setPlacement(pv);
      volSurfaceList(sensor_detEl_back)->push_back(surf_back);

      auto &params = DD4hepDetectorHelper::ensureExtension<dd4hep::rec::VariantParameters>(sensor_detEl_front);
      params.set<string>("axis_definitions", "XZY");

      isensor++;

      // NOTE: placement of sensor supports on back side of disk
      {
        RotationZYX rot_back_flip(0, 0, 0);
        RotationZYX rot_back_flip2(M_PI, 0, 0);

        pv = layer_assembly.placeVolume(
            log_sensor_stack1,
            Transform3D(ilay_y % 2 == 0 ? rot_back_flip2 : rot_back_flip,
                        Position((ilay_x)*baseplate_xdim,
                                 ilay_y % 2 == 0
                                     ? (ilay_y)*sensor_hybrid_ydim + sensor_hybrid_ydim / 2 - baseplate_ydim / 2
                                     : (ilay_y)*sensor_hybrid_ydim - sensor_hybrid_ydim / 2 + baseplate_ydim / 2,
                                 disk_face_rel_z + thicknessDet1 / 2.0)));

        pv = layer_assembly.placeVolume(
            log_sensor_stack2,
            Transform3D(ilay_y % 2 == 0 ? rot_back_flip2 : rot_back_flip,
                        Position((ilay_x)*baseplate_xdim,
                                 ilay_y % 2 == 0
                                     ? (ilay_y)*sensor_hybrid_ydim + sensor_hybrid_ydim / 2 - baseplate_ydim / 2
                                     : (ilay_y)*sensor_hybrid_ydim - sensor_hybrid_ydim / 2 + baseplate_ydim / 2,
                                 disk_face_rel_z + thicknessDet1 + sensor_thickness + thicknessDet2 / 2.0)));
      }

      // NOTE: placement of service hybrids on back side
      pv = layer_assembly.placeVolume(
          log_SH_stack,
          Transform3D(rot_back,
                      Position((ilay_x)*baseplate_xdim,
                               ilay_y % 2 == 0
                                   ? (ilay_y)*sensor_hybrid_ydim - sensor_hybrid_ydim / 2 + baseplate_ydim / 4
                                   : (ilay_y)*sensor_hybrid_ydim + sensor_hybrid_ydim / 2 - baseplate_ydim / 4,
                               disk_face_rel_z + thicknessDet_SH / 2.0)));

      // NOTE: place cooling lines and liquid
      RotationZYX rot_cooling(0, M_PI / 2.0, 0);

      pv = layer_assembly.placeVolume(
          log_cooling_liquid,
          Transform3D(rot_cooling,
                      Position((ilay_x)*baseplate_xdim, (ilay_y)*sensor_hybrid_ydim + sensor_hybrid_ydim / 2, 0.)));
      pv = layer_assembly.placeVolume(
          log_cooling_tube,
          Transform3D(rot_cooling,
                      Position((ilay_x)*baseplate_xdim, (ilay_y)*sensor_hybrid_ydim + sensor_hybrid_ydim / 2, 0.)));
    }
  }
  layer_assembly->GetShape()->ComputeBBox();

  layer_id++;

  pv = description.pickMotherVolume(ttl_detEl).placeVolume(assembly, Position(0, 0, disk_zPos));
  pv.addPhysVolID("system", det_id);
  ttl_detEl.setPlacement(pv);

  return ttl_detEl;
}

//@}
// clang-format off
DECLARE_DETELEMENT(epic_TOFEndcap, create_detector)
