/** \addtogroup Trackers Trackers
 * \brief Type: **BarrelTrackerWithFrame**.
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

#if defined(USE_ACTSDD4HEP)
#include "ActsDD4hep/ActsExtension.hpp"
#include "ActsDD4hep/ConvertMaterial.hpp"
#else
#include "Acts/Plugins/DD4hep/ActsExtension.hpp"
#include "Acts/Plugins/DD4hep/ConvertDD4hepMaterial.hpp"
#endif

using namespace dd4hep;
using namespace dd4hep::rec;
using namespace dd4hep::detail;

/** Endcap Trapezoidal Tracker.
 *
 * @author Whitney Armstrong
 *
 */
static Ref_t create_detector(Detector& description, xml_h e, SensitiveDetector sens)
{
  typedef std::vector<PlacedVolume> Placements;
  xml_det_t x_det = e;
  int det_id   = x_det.id();
  std::string det_name = x_det.nameStr();
  bool reflect  = x_det.reflect(false);
  DetElement sdet(det_name, det_id);
  Assembly assembly(det_name);
  Material air = description.material("Air");
  Volume motherVol = description.pickMotherVolume(sdet);

  int m_id = 0;
  std::map<std::string, Placements> sensitives;
  std::map<std::string, std::vector<VolPlane>> volplane_surfaces;
  PlacedVolume pv;

  // ACTS extension
  {
    Acts::ActsExtension* detWorldExt = new Acts::ActsExtension();
    detWorldExt->addType("endcap", "detector");
    // SJJ probably need to set the envelope here, as ACTS can't figure
    // that out for Assembly volumes. May also need binning to properly pick up
    // on the support material @TODO
    //
    // Add the volume boundary material if configured
    for (xml_coll_t bmat(x_det, _Unicode(boundary_material)); bmat; ++bmat) {
      xml_comp_t x_boundary_material = bmat;
      Acts::xmlToProtoSurfaceMaterial(x_boundary_material, *detWorldExt, "boundary_material");
    }
    sdet.addExtension<Acts::ActsExtension>(detWorldExt);
  }

  assembly.setVisAttributes(description.invisible());
  sens.setType("tracker");

  // this loop creates the forward and the backward disk
  for (xml_coll_t mi(x_det, _U(module)); mi; ++mi, ++m_id) {
    xml_comp_t x_mod = mi;
    std::string     m_nam = x_mod.nameStr();
    xml_comp_t diskdimension       = x_mod.dimensions();

    // load all information from the diskdimension definitions
    double     disk_zPos            = getAttrOrDefault(diskdimension, _Unicode(zPos), 0.);
    double     disk_rMin            = getAttrOrDefault(diskdimension, _Unicode(rMin), 0.);
    double     disk_rMax            = getAttrOrDefault(diskdimension, _Unicode(rMax), 100.);
    double     disk_xOffset         = getAttrOrDefault(diskdimension, _Unicode(xOffset), 0.);
    double     disk_det_height      = getAttrOrDefault(diskdimension, _Unicode(det_height), 0.);
    double     cooling_plate_height = getAttrOrDefault(diskdimension, _Unicode(cooling_plate_height), 0.);
    double     cooling_tube_thickness    = getAttrOrDefault(diskdimension, _Unicode(wallthickness_coolingtube), 0.);
    double     cooling_tube_diameter    = getAttrOrDefault(diskdimension, _Unicode(diameter_coolingtube), 0.);
    // for the electron-going direction we invert the z position
    if(reflect) disk_zPos = -disk_zPos;

    // some debug output
    // std::cout << "\tdisk_rMin = " << disk_rMin
    // << "\n\tdisk_rMax = " << disk_rMax
    // << "\n\tdisk_xOffset = " << disk_xOffset
    // << "\n\tdisk_zPos = " << disk_zPos
    // << "\n\tdisk_det_height = " << disk_det_height
    // << "\n\tdisk_det_height = " << disk_det_height
    // << "\n\tcooling_tube_thickness = " << cooling_tube_thickness
    // << "\n\tcooling_tube_diameter = " << cooling_tube_diameter
    // << std::endl;

    // create a solid for the beampipe cutout
    Tube beampipe_cutout(0, disk_rMin,disk_det_height, 0, 2. * M_PI);
    // create solids for the cooling/mountin disks (two Al disks)
    Tube ttl_plate_1(0, disk_rMax,cooling_plate_height / 2., 0, 2. * M_PI);
    Tube ttl_plate_2(0, disk_rMax,cooling_plate_height / 2., 0, 2. * M_PI);

    // create the offset cutout in the disks
    Position cutoutOffset(disk_xOffset, 0,0);
    SubtractionSolid ttl_cooling_plate_1(ttl_plate_1, beampipe_cutout,cutoutOffset);
    SubtractionSolid ttl_cooling_plate_2(ttl_plate_2, beampipe_cutout,cutoutOffset);

    // set the material (Al)
    Material disk_mat = description.material(diskdimension.materialStr());

    // create the volumes and set visualization options
    Volume ttl_cooling_plate_volume_1("ttl_cooling_plate_1", ttl_cooling_plate_1, disk_mat);
    Volume ttl_cooling_plate_volume_2("ttl_cooling_plate_2", ttl_cooling_plate_2, disk_mat);
    ttl_cooling_plate_volume_1.setVisAttributes(description.visAttributes("TOFAluminum"));
    ttl_cooling_plate_volume_2.setVisAttributes(description.visAttributes("TOFAluminum"));

    // place disks in assembly
    pv = assembly.placeVolume(ttl_cooling_plate_volume_1, Position(0, 0, disk_zPos-cooling_tube_diameter/2.-cooling_plate_height/2.));
    pv = assembly.placeVolume(ttl_cooling_plate_volume_2, Position(0, 0, disk_zPos+cooling_tube_diameter/2+cooling_plate_height/2.));

    // ---------------------------------------
    // Create Sensor Module Stack:
    // ---------------------------------------

    xml_comp_t sensparams   = x_mod.child(_Unicode(sensorparameters));;

    // load all information from the sensparams definitions
    // double sensor_width     = getAttrOrDefault(sensparams, _Unicode(width), 0.);
    // double sensor_length    = getAttrOrDefault(sensparams, _Unicode(length), 0.);
    double baseplate_length = getAttrOrDefault(sensparams, _Unicode(baseplate_length), 100.);
    double baseplate_width  = getAttrOrDefault(sensparams, _Unicode(baseplate_width), 0.);



    std::string strLayerName[10];
    Material materialLayer[10];
    double thicknessLayer[10];
    double widthLayer[10];
    double offsetLayer[10];
    double lengthLayer[10];
    int nLayers = 0;
    for (xml_coll_t l_iter(x_det, _U(layer)); l_iter; ++l_iter) {
      xml_comp_t x_layer = l_iter;
      if(x_layer.id()!=2) continue;

      nLayers = getAttrOrDefault(x_layer, _Unicode(numslices), 0.);

      int    i_slice     = 0;
      for (xml_coll_t s_iter(x_layer, _U(slice)); s_iter; ++s_iter, ++i_slice) {
        // If slices are only given a thickness attribute, they are radially concentric slices
        // If slices are given an inner_z attribute, they are longitudinal slices with equal rmin
        xml_comp_t x_slice      = s_iter;
        materialLayer[i_slice]  = description.material(x_slice.materialStr());
        strLayerName[i_slice]   = getAttrOrDefault<std::string>(x_slice, _U(name), "slice" + std::to_string(i_slice));
        thicknessLayer[i_slice] = x_slice.thickness();
        widthLayer[i_slice] = getAttrOrDefault(x_slice, _Unicode(width), 0.);
        offsetLayer[i_slice] = getAttrOrDefault(x_slice, _Unicode(offset), 0.);
        lengthLayer[i_slice] = getAttrOrDefault(x_slice, _Unicode(length), 0.);
      }
    }




    // need to define two thicknesses here as the sensors have to be placed separately from the rest of the stack
    // this is done to properly assign sensor IDs to the placed volumes
    double thicknessDet1 = 0;
    double thicknessDet2 = 0;
    for (int ilay = 0; ilay < nLayers; ilay++)
    {
      // layer 5 would be the sensor
      if(ilay<5 )thicknessDet1 += thicknessLayer[ilay];
      if(ilay>5 )thicknessDet2 += thicknessLayer[ilay];
    }

    // create the stack volumes in the same length as the baseplate
    double segmentlength = baseplate_length;  //(detlength - 10 * cm) / 6;//
    Box box_sensor_stack1(baseplate_length / 2, baseplate_width / 2, thicknessDet1 / 2);
    Box box_sensor_stack2(baseplate_length / 2, baseplate_width / 2, thicknessDet2 / 2);
    Volume log_sensor_stack1("log_sensor_stack1", box_sensor_stack1, air);
    Volume log_sensor_stack2("log_sensor_stack2", box_sensor_stack2, air);
    log_sensor_stack1.setVisAttributes(description.visAttributes("InvisibleWithDaughters"));
    log_sensor_stack2.setVisAttributes(description.visAttributes("InvisibleWithDaughters"));

    // individual layers will be placed in subsequent order and z_start is increased by the thickness of the previous layer
    double z_start1 = -thicknessDet1 / 2;
    double z_start2 = -thicknessDet2 / 2;
    for (int ilay = 0; ilay < nLayers; ilay++)
    {
      const std::string layer_name = "sensor_stack_" + strLayerName[ilay];
      const std::string layer_name_Solid = "sol_" + layer_name;

      Box box_Module_Layer_Raw(lengthLayer[ilay] / 2,widthLayer[ilay] / 2,thicknessLayer[ilay] / 2);
      Solid     sol_Module_Layer_Raw = box_Module_Layer_Raw;
      Volume Log_Layer(layer_name + "_Log", sol_Module_Layer_Raw, materialLayer[ilay]);

      if(ilay<5 ){
        pv = log_sensor_stack1.placeVolume(Log_Layer, Position(0, -offsetLayer[ilay], z_start1 + thicknessLayer[ilay] / 2));
        z_start1 += thicknessLayer[ilay];
      }
      if(ilay>5 ){
        pv = log_sensor_stack2.placeVolume(Log_Layer, Position(0, -offsetLayer[ilay], z_start2 + thicknessLayer[ilay] / 2));
        z_start2 += thicknessLayer[ilay];
      }
      Log_Layer.setVisAttributes(description.visAttributes("TOFLayers"));
    }
  
  
    // ---------------------------------------
    // Create Service Hybrid Module Stack:
    // ---------------------------------------
    std::string strLayerName_SH[10];
    Material materialLayer_SH[10];
    double thicknessLayer_SH[10];
    double widthLayer_SH[10];
    double offsetLayer_SH[10];
    double lengthLayer_SH[10];
    int nLayers_SH = 0;
    for (xml_coll_t l_iter(x_det, _U(layer)); l_iter; ++l_iter) {
      xml_comp_t x_layer = l_iter;
      if(x_layer.id()!=3) continue;

      nLayers_SH = getAttrOrDefault(x_layer, _Unicode(numslices), 0.);

      int    i_slice     = 0;
      for (xml_coll_t s_iter(x_layer, _U(slice)); s_iter; ++s_iter, ++i_slice) {
        // If slices are only given a thickness attribute, they are radially concentric slices
        // If slices are given an inner_z attribute, they are longitudinal slices with equal rmin
        xml_comp_t x_slice      = s_iter;
        materialLayer_SH[i_slice]  = description.material(x_slice.materialStr());
        strLayerName_SH[i_slice]   = getAttrOrDefault<std::string>(x_slice, _U(name), "slice" + std::to_string(i_slice));
        thicknessLayer_SH[i_slice] = x_slice.thickness();
        widthLayer_SH[i_slice] = getAttrOrDefault(x_slice, _Unicode(width), 0.);
        offsetLayer_SH[i_slice] = getAttrOrDefault(x_slice, _Unicode(offset), 0.);
        lengthLayer_SH[i_slice] = getAttrOrDefault(x_slice, _Unicode(length), 0.);
      }
    }


    double baseSH_width = baseplate_width / 2;
    double thicknessDet_SH = 0;
    for (int ilay = 0; ilay < nLayers_SH; ilay++)
    {
      thicknessDet_SH += thicknessLayer_SH[ilay];
    }

    Box sol_SH_stack(baseplate_length / 2,baseSH_width / 2, thicknessDet_SH / 2);
    Volume log_SH_stack("log_SH_stack", sol_SH_stack, air);
    log_SH_stack.setVisAttributes(description.visAttributes("TOFSensorAndReadoutLadder"));

    // place all layers into the SH stack
    double z_start_SH = -thicknessDet_SH / 2;
    for (int ilay = 0; ilay < nLayers_SH; ilay++)
    {
      const std::string layer_name = "SH_stack_" + strLayerName_SH[ilay];
      const std::string layer_name_Solid = "sol_" + layer_name;

      Box sol_Module_Layer_Raw(lengthLayer_SH[ilay] / 2,widthLayer_SH[ilay] / 2,thicknessLayer_SH[ilay] / 2);
      Volume Log_Layer(layer_name + "_Log", sol_Module_Layer_Raw, materialLayer_SH[ilay]);

      pv = log_SH_stack.placeVolume(Log_Layer, Position(0, -offsetLayer_SH[ilay], z_start_SH + thicknessLayer_SH[ilay] / 2));

      z_start_SH += thicknessLayer_SH[ilay];
      Log_Layer.setVisAttributes(description.visAttributes("TOFLayers"));
    }

    // define the width of the whole sensor+SH stack
    double fullsensor_width = baseSH_width+baseplate_width;

    // stacks are placed on both sides of the detector and require different offsets
    double offsetzFront = cooling_tube_diameter/2 + cooling_plate_height + thicknessDet_SH / 2;
    double offsetzBack = -cooling_tube_diameter/2 - cooling_plate_height - thicknessDet_SH / 2;

    //--------------------------------------------------------------------------------
    // Main loop for sensor and SH placement:
    //--------------------------------------------------------------------------------

    // number of towers in radial direction (on y axis)
    int rowYdir = (int) ( (disk_rMax-(fullsensor_width/2)) / fullsensor_width);
    int sensorCount = 0;
    for(int row=rowYdir;row>=-rowYdir;row--){
      // pythagoras -> get available length in circular mother volume for towers
      // divide given length by tower width -> get number of towers that can be placed
      int numSensorsRow = (int) ( ( 2* sqrt(pow(disk_rMax,2)-pow( (abs(row)*fullsensor_width + fullsensor_width/2) ,2)) ) / segmentlength );
      if(numSensorsRow==0) continue;
      // we want an odd number of towers to be symmetrically centered around 0
      if ( numSensorsRow % 2 == 0) numSensorsRow-=1;

      if( ( (abs(row)*fullsensor_width) -(fullsensor_width/2)) < disk_rMin ){
        // pythagoras -> get available length in circular mother volume for towers
        // divide given length by tower width -> get number of towers that can be placed
        int numSensorLeftAdd = ceil( (disk_xOffset -(segmentlength/2) - sqrt(pow(disk_rMin,2)-pow( (abs(row)*fullsensor_width)-(fullsensor_width/2) ,2)) ) / segmentlength );
        int numSensorRightAdd = ceil( (disk_xOffset -(segmentlength/2) + sqrt(pow(disk_rMin,2)-pow( (abs(row)*fullsensor_width)-(fullsensor_width/2) ,2)) ) / segmentlength );

        // we want an odd number of towers to be symmetrically centered around 0
        // create mother volume with space for numSensorsRow towers along x-axis
        int nSensLeft = ((numSensorsRow-1) /2 + numSensorLeftAdd);
        Box TTLDetRowLeftSolid((nSensLeft) * segmentlength / 2.0,fullsensor_width / 2.0,thicknessDet_SH / 2.0);
        Volume TTLDetRowLeftLogical_Front("TTLDetRowLeftLogical_Front" + std::to_string(row), TTLDetRowLeftSolid, air);
        TTLDetRowLeftLogical_Front.setVisAttributes(description.visAttributes("InvisibleWithDaughters"));

        //place log_sensor_and_readout_L in TTLDetRowLeftLogical_Front
        for(int sensor = 0; sensor < nSensLeft; sensor++) {
          Box box_sensor_and_readout_L(segmentlength / 2,fullsensor_width / 2,thicknessDet_SH / 2);
          Volume log_sensor_and_readout_L("log_sensor_and_readout_LF_" + std::to_string(row)+"_"+ std::to_string(sensorCount), box_sensor_and_readout_L, air);

          log_sensor_and_readout_L.setVisAttributes(description.visAttributes("InvisibleWithDaughters"));

          log_sensor_and_readout_L.placeVolume(log_sensor_stack1, Position(0, -fullsensor_width/2 + baseplate_width/2, -thicknessDet_SH/2+thicknessDet1/2));
          log_sensor_and_readout_L.placeVolume(log_sensor_stack2, Position(0, -fullsensor_width/2 + baseplate_width/2, -thicknessDet_SH/2+thicknessDet1+thicknessLayer[5]+thicknessDet2/2));
          log_sensor_and_readout_L.placeVolume(log_SH_stack, Position(0,  fullsensor_width/2 - baseplate_width/4, 0));


          // make active sensor solid and volume which is separately placed in the detector compared to the other layers
          int senslayer=5;
          Box sol_ACLGAD_Sens(lengthLayer[senslayer] / 2,widthLayer[senslayer] / 2,thicknessLayer[senslayer] / 2);
          Volume log_ACLGAD_Sens_L("log_ACLGAD_Sens_LF_" + std::to_string(row)+"_"+  std::to_string(sensorCount), sol_ACLGAD_Sens, materialLayer[senslayer]);
          log_ACLGAD_Sens_L.setVisAttributes(description.visAttributes("TOFActiveMat"));

          pv = log_sensor_and_readout_L.placeVolume(log_ACLGAD_Sens_L, Position(0, -fullsensor_width/2 + baseplate_width/2 - offsetLayer[senslayer], -thicknessDet_SH/2+thicknessDet1+thicknessLayer[5]/2));

          // make sensor sensitive, set ID, type and add to sensitive volumes
          pv.addPhysVolID("sensor", sensorCount);
          sens.setType("tracker");
          log_ACLGAD_Sens_L.setSensitiveDetector(sens);
          sensitives[m_nam].push_back(pv);

          // TODO: the block below is taken from TrapEndcapTracker_geo.cpp and should be checked if it is needed here
          // -------- create a measurement plane for the tracking surface attched to the sensitive volume -----
          Vector3D u(0., 0., -1.);
          Vector3D v(-1., 0., 0.);
          Vector3D n(0., 1., 0.);

          // compute the inner and outer thicknesses that need to be assigned to the tracking surface
          // depending on wether the support is above or below the sensor
          double inner_thickness = 0;
          double outer_thickness = thicknessLayer[senslayer];

          SurfaceType type(SurfaceType::Sensitive);
          VolPlane surf(log_ACLGAD_Sens_L, type, inner_thickness, outer_thickness, u, v, n); //,o ) ;
          volplane_surfaces[m_nam].push_back(surf);
          // TODO: end

          //place log_sensor_and_readout_L in TTLDetRowLeftLogical_Front
          Position pos_sensor((sensor - (nSensLeft - 1) / 2.0) * segmentlength, 0, 0);
          pv = TTLDetRowLeftLogical_Front.placeVolume(log_sensor_and_readout_L, pos_sensor);
          sensorCount++;

        }

        //place TTLDetRowLeftLogical_Front in assembly, roate every second row
        Transform3D  transfrow(RotationZYX(row%2==0 ? - M_PI : 0,0, 0), Translation3D( - ( (nSensLeft) * segmentlength / 2.0 + segmentlength / 2.0 - (numSensorLeftAdd* segmentlength)), (row*fullsensor_width), disk_zPos+offsetzFront));
        assembly.placeVolume(TTLDetRowLeftLogical_Front, transfrow);
        pv.addPhysVolID("layer", 0);
        pv.addPhysVolID("module", row);


        Volume TTLDetRowLeftLogical_Back("TTLDetRowLeftLogical_Back" + std::to_string(row), TTLDetRowLeftSolid, air);
        TTLDetRowLeftLogical_Back.setVisAttributes(description.visAttributes("InvisibleWithDaughters"));

        //place log_sensor_and_readout_L in TTLDetRowLeftLogical_Back
        for(int sensor = 0; sensor < nSensLeft; sensor++) {
          Box box_sensor_and_readout_L(segmentlength / 2,fullsensor_width / 2,thicknessDet_SH / 2);
          Volume log_sensor_and_readout_L("log_sensor_and_readout_LB_" + std::to_string(row)+"_"+ std::to_string(sensorCount), box_sensor_and_readout_L, air);

          log_sensor_and_readout_L.setVisAttributes(description.visAttributes("InvisibleWithDaughters"));

          log_sensor_and_readout_L.placeVolume(log_sensor_stack1, Position(0, -fullsensor_width/2 + baseplate_width/2, -thicknessDet_SH/2+thicknessDet1/2));
          log_sensor_and_readout_L.placeVolume(log_sensor_stack2, Position(0, -fullsensor_width/2 + baseplate_width/2, -thicknessDet_SH/2+thicknessDet1+thicknessLayer[5]+thicknessDet2/2));
          log_sensor_and_readout_L.placeVolume(log_SH_stack, Position(0,  fullsensor_width/2 - baseplate_width/4, 0));


          // make active sensor solid and volume
          int senslayer=5;
          Box sol_ACLGAD_Sens(lengthLayer[senslayer] / 2,widthLayer[senslayer] / 2,thicknessLayer[senslayer] / 2);
          Volume log_ACLGAD_Sens_L("log_ACLGAD_Sens_LB_" + std::to_string(row)+"_"+  std::to_string(sensorCount), sol_ACLGAD_Sens, materialLayer[senslayer]);
          log_ACLGAD_Sens_L.setVisAttributes(description.visAttributes("TOFActiveMat"));


          pv = log_sensor_and_readout_L.placeVolume(log_ACLGAD_Sens_L, Position(0, -fullsensor_width/2 + baseplate_width/2 -offsetLayer[senslayer], -thicknessDet_SH/2+thicknessDet1+thicknessLayer[5]/2));


          pv.addPhysVolID("sensor", sensorCount);
          sens.setType("tracker");
          log_ACLGAD_Sens_L.setSensitiveDetector(sens);
          sensitives[m_nam].push_back(pv);
          // pv.addPhysVolID("layer", 0);
          // pv.addPhysVolID("module", row);
          // -------- create a measurement plane for the tracking surface attched to the sensitive volume -----
          Vector3D u(0., 0., -1.);
          Vector3D v(-1., 0., 0.);
          Vector3D n(0., 1., 0.);

          // compute the inner and outer thicknesses that need to be assigned to the tracking surface
          // depending on wether the support is above or below the sensor
          double inner_thickness = 0;
          double outer_thickness = thicknessLayer[senslayer];

          SurfaceType type(SurfaceType::Sensitive);
          VolPlane surf(log_ACLGAD_Sens_L, type, inner_thickness, outer_thickness, u, v, n); //,o ) ;
          volplane_surfaces[m_nam].push_back(surf);

          //place log_sensor_and_readout_L in TTLDetRowLeftLogical_Back
          Position pos_sensor((sensor - (nSensLeft - 1) / 2.0) * segmentlength, 0, 0);
          pv = TTLDetRowLeftLogical_Back.placeVolume(log_sensor_and_readout_L, pos_sensor);
          sensorCount++;

        }

        //place TTLDetRowLeftLogical_Back on other side in assembly, roate every second row
        Transform3D  transfrow2(RotationZYX(row%2==0 ? - M_PI : 0,0,-M_PI), Translation3D( - ( (nSensLeft) * segmentlength / 2.0 + segmentlength / 2.0 - (numSensorLeftAdd* segmentlength)), (row*fullsensor_width), disk_zPos+offsetzBack));
        assembly.placeVolume(TTLDetRowLeftLogical_Back, transfrow2);
        pv.addPhysVolID("layer", 1);
        pv.addPhysVolID("module", row);
      

        // we want an odd number of towers to be symmetrically centered around 0
        // create mother volume with space for numSensorsRow towers along x-axis
        int nSensRight = ((numSensorsRow-1) /2 - numSensorRightAdd);
        Box TTLDetRowRightSolid( nSensRight * segmentlength / 2.0,fullsensor_width / 2.0,thicknessDet_SH / 2.0);
        Volume TTLDetRowRightLogical_Front("TTLDetRowRightLogical_Front" + std::to_string(row), TTLDetRowRightSolid, air);
        TTLDetRowRightLogical_Front.setVisAttributes(description.visAttributes("InvisibleWithDaughters"));

        //place log_sensor_and_readout_R in TTLDetRowRightLogical_Front
        for(int sensor = 0; sensor < nSensRight; sensor++) {
          Box box_sensor_and_readout_R(segmentlength / 2,fullsensor_width / 2,thicknessDet_SH / 2);
          Volume log_sensor_and_readout_R("log_sensor_and_readout_RF_" + std::to_string(row)+"_"+ std::to_string(sensorCount), box_sensor_and_readout_R, air);

          log_sensor_and_readout_R.setVisAttributes(description.visAttributes("InvisibleWithDaughters"));

          log_sensor_and_readout_R.placeVolume(log_sensor_stack1, Position(0, -fullsensor_width/2 + baseplate_width/2, -thicknessDet_SH/2+thicknessDet1/2));
          log_sensor_and_readout_R.placeVolume(log_sensor_stack2, Position(0, -fullsensor_width/2 + baseplate_width/2, -thicknessDet_SH/2+thicknessDet1+thicknessLayer[5]+thicknessDet2/2));
          log_sensor_and_readout_R.placeVolume(log_SH_stack, Position(0,  fullsensor_width/2 - baseplate_width/4, 0));


          // make active sensor solid and volume
          int senslayer=5;
          Box sol_ACLGAD_Sens(lengthLayer[senslayer] / 2,widthLayer[senslayer] / 2,thicknessLayer[senslayer] / 2);
          Volume log_ACLGAD_Sens_R("log_ACLGAD_Sens_RF_" + std::to_string(row)+"_"+ std::to_string(sensorCount), sol_ACLGAD_Sens, materialLayer[senslayer]);
          log_ACLGAD_Sens_R.setVisAttributes(description.visAttributes("TOFActiveMat"));

          pv = log_sensor_and_readout_R.placeVolume(log_ACLGAD_Sens_R, Position(0, -fullsensor_width/2 + baseplate_width/2 -offsetLayer[senslayer], -thicknessDet_SH/2+thicknessDet1+thicknessLayer[5]/2));
          pv.addPhysVolID("sensor", sensorCount);
          sens.setType("tracker");
          log_ACLGAD_Sens_R.setSensitiveDetector(sens);
          sensitives[m_nam].push_back(pv);
          // -------- create a measurement plane for the tracking surface attched to the sensitive volume -----
          Vector3D u(0., 0., -1.);
          Vector3D v(-1., 0., 0.);
          Vector3D n(0., 1., 0.);

          // compute the inner and outer thicknesses that need to be assigned to the tracking surface
          // depending on wether the support is above or below the sensor
          double inner_thickness = 0;
          double outer_thickness = thicknessLayer[senslayer];

          SurfaceType type(SurfaceType::Sensitive);
          VolPlane surf(log_ACLGAD_Sens_R, type, inner_thickness, outer_thickness, u, v, n); //,o ) ;
          volplane_surfaces[m_nam].push_back(surf);

          //place log_sensor_and_readout_R in TTLDetRowRightLogical_Front
          Position pos_sensor((sensor - (nSensRight - 1) / 2.0) * segmentlength, 0, 0);
          // PlacedVolume sensor_phys = 
          pv = TTLDetRowRightLogical_Front.placeVolume(log_sensor_and_readout_R, pos_sensor);
          sensorCount++;
        }

        //place TTLDetRowRightLogical_Front in assembly, roate every second row
        Transform3D  transfrowR(RotationZYX(row%2==0 ? - M_PI : 0,0, 0), Translation3D(( nSensRight * segmentlength / 2.0 + segmentlength / 2.0 + (numSensorRightAdd* segmentlength)), (row*fullsensor_width), disk_zPos+offsetzFront));
        pv = assembly.placeVolume(TTLDetRowRightLogical_Front, transfrowR);
        pv.addPhysVolID("layer", 0);
        pv.addPhysVolID("module", row);


        Volume TTLDetRowRightLogical_Back("TTLDetRowRightLogical_Back" + std::to_string(row), TTLDetRowRightSolid, air);
        TTLDetRowRightLogical_Back.setVisAttributes(description.visAttributes("InvisibleWithDaughters"));

        //place log_sensor_and_readout_R in TTLDetRowRightLogical_Back
        for(int sensor = 0; sensor < nSensRight; sensor++) {
          Box box_sensor_and_readout_R(segmentlength / 2,fullsensor_width / 2,thicknessDet_SH / 2);
          Volume log_sensor_and_readout_R("log_sensor_and_readout_RB_" + std::to_string(row)+"_"+ std::to_string(sensorCount), box_sensor_and_readout_R, air);

          log_sensor_and_readout_R.setVisAttributes(description.visAttributes("InvisibleWithDaughters"));

          log_sensor_and_readout_R.placeVolume(log_sensor_stack1, Position(0, -fullsensor_width/2 + baseplate_width/2, -thicknessDet_SH/2+thicknessDet1/2));
          log_sensor_and_readout_R.placeVolume(log_sensor_stack2, Position(0, -fullsensor_width/2 + baseplate_width/2, -thicknessDet_SH/2+thicknessDet1+thicknessLayer[5]+thicknessDet2/2));
          log_sensor_and_readout_R.placeVolume(log_SH_stack, Position(0,  fullsensor_width/2 - baseplate_width/4, 0));


          // make active sensor solid and volume
          int senslayer=5;
          Box sol_ACLGAD_Sens(lengthLayer[senslayer] / 2,widthLayer[senslayer] / 2,thicknessLayer[senslayer] / 2);
          Volume log_ACLGAD_Sens_R("log_ACLGAD_Sens_RB_" + std::to_string(row)+"_"+ std::to_string(sensorCount), sol_ACLGAD_Sens, materialLayer[senslayer]);
          log_ACLGAD_Sens_R.setVisAttributes(description.visAttributes("TOFActiveMat"));

          pv = log_sensor_and_readout_R.placeVolume(log_ACLGAD_Sens_R, Position(0, -fullsensor_width/2 + baseplate_width/2 -offsetLayer[senslayer], -thicknessDet_SH/2+thicknessDet1+thicknessLayer[5]/2));
          pv.addPhysVolID("sensor", sensorCount);
          sens.setType("tracker");
          log_ACLGAD_Sens_R.setSensitiveDetector(sens);
          sensitives[m_nam].push_back(pv);
          // -------- create a measurement plane for the tracking surface attched to the sensitive volume -----
          Vector3D u(0., 0., -1.);
          Vector3D v(-1., 0., 0.);
          Vector3D n(0., 1., 0.);

          // compute the inner and outer thicknesses that need to be assigned to the tracking surface
          // depending on wether the support is above or below the sensor
          double inner_thickness = 0;
          double outer_thickness = thicknessLayer[senslayer];

          SurfaceType type(SurfaceType::Sensitive);
          VolPlane surf(log_ACLGAD_Sens_R, type, inner_thickness, outer_thickness, u, v, n); //,o ) ;
          volplane_surfaces[m_nam].push_back(surf);

          //place log_sensor_and_readout_R in TTLDetRowRightLogical_Back
          Position pos_sensor((sensor - (nSensRight - 1) / 2.0) * segmentlength, 0, 0);
          // PlacedVolume sensor_phys = 
          pv = TTLDetRowRightLogical_Back.placeVolume(log_sensor_and_readout_R, pos_sensor);
          sensorCount++;
        }


        //place TTLDetRowRightLogical_Back on other side in assembly, roate every second row
        Transform3D  transfrowR2(RotationZYX(row%2==0 ? - M_PI : 0,0,-M_PI), Translation3D(( nSensRight * segmentlength / 2.0 + segmentlength / 2.0 + (numSensorRightAdd* segmentlength)), (row*fullsensor_width), disk_zPos+offsetzBack));
        pv = assembly.placeVolume(TTLDetRowRightLogical_Back, transfrowR2);
        pv.addPhysVolID("layer", 1);
        pv.addPhysVolID("module", row);



        Box sol_cutout_tube_left(1.3*(sqrt(pow(disk_rMax,2)-pow( (abs(row)*fullsensor_width) - (fullsensor_width/2.0) ,2))- sqrt(pow(disk_rMin,2)-pow( (abs(row)*fullsensor_width)-(fullsensor_width/2) ,2)) + disk_xOffset)/2,(cooling_tube_diameter - 2*cooling_tube_thickness) / 2,(cooling_tube_diameter - 2*cooling_tube_thickness) / 2);
        Box box_cooling_tube_left(0.86*(sqrt(pow(disk_rMax,2)-pow( (abs(row)*fullsensor_width) - (fullsensor_width/2.0) ,2))- sqrt(pow(disk_rMin,2)-pow( (abs(row)*fullsensor_width)-(fullsensor_width/2) ,2)) + disk_xOffset)/2,cooling_tube_diameter / 2,cooling_tube_diameter / 2);

        SubtractionSolid sol_cooling_tube_left(box_cooling_tube_left, sol_cutout_tube_left);

        Volume Log_cooling_tube_left("Log_cooling_tube_left" + std::to_string(row), sol_cooling_tube_left, description.material("Aluminum"));
        Log_cooling_tube_left.setVisAttributes(description.visAttributes("TOFAluminum"));


        Box sol_water_cooling_left(0.85*(sqrt(pow(disk_rMax,2)-pow( (abs(row)*fullsensor_width) - (fullsensor_width/2.0) ,2))- sqrt(pow(disk_rMin,2)-pow( (abs(row)*fullsensor_width)-(fullsensor_width/2) ,2)) + disk_xOffset)/2, 0.99*(cooling_tube_diameter - 2*cooling_tube_thickness) / 2, 0.99*(cooling_tube_diameter - 2*cooling_tube_thickness) / 2);
        Volume Log_water_cooling_left("Log_water_cooling_left" + std::to_string(row), sol_water_cooling_left, description.material("Water"));

        Log_water_cooling_left.setVisAttributes(description.visAttributes("TOFWater"));

        //place Log_water_cooling_left in assembly
        Transform3D  transfrowCL(RotationZYX(row%2==0 ? - M_PI : 0,0, 0), Translation3D(- ( (nSensLeft) * segmentlength / 2.0 + segmentlength / 2.0 - (numSensorLeftAdd* segmentlength)), row>0 ? (row*fullsensor_width) - (fullsensor_width/2.0) : (row*fullsensor_width) + (fullsensor_width/2.0), disk_zPos));
        if(row!=0){
          assembly.placeVolume(Log_cooling_tube_left, transfrowCL);
          assembly.placeVolume(Log_water_cooling_left, transfrowCL);
        }


        Box sol_cutout_tube_right(1.3*(sqrt(pow(disk_rMax,2)-pow( (abs(row)*fullsensor_width) - (fullsensor_width/2.0) ,2))- sqrt(pow(disk_rMin,2)-pow( (abs(row)*fullsensor_width)-(fullsensor_width/2) ,2)) - disk_xOffset)/2,(cooling_tube_diameter - 2*cooling_tube_thickness) / 2,(cooling_tube_diameter - 2*cooling_tube_thickness) / 2);
        Box box_cooling_tube_right(0.86*(sqrt(pow(disk_rMax,2)-pow( (abs(row)*fullsensor_width) - (fullsensor_width/2.0) ,2))- sqrt(pow(disk_rMin,2)-pow( (abs(row)*fullsensor_width)-(fullsensor_width/2) ,2)) - disk_xOffset)/2,cooling_tube_diameter / 2,cooling_tube_diameter / 2);

        SubtractionSolid sol_cooling_tube_right(box_cooling_tube_right, sol_cutout_tube_right);

        Volume Log_cooling_tube_right("Log_cooling_tube_right" + std::to_string(row), sol_cooling_tube_right, description.material("Aluminum"));
        Log_cooling_tube_right.setVisAttributes(description.visAttributes("TOFAluminum"));


        Box sol_water_cooling_right(0.85*(sqrt(pow(disk_rMax,2)-pow( (abs(row)*fullsensor_width) - (fullsensor_width/2.0) ,2))- sqrt(pow(disk_rMin,2)-pow( (abs(row)*fullsensor_width)-(fullsensor_width/2) ,2)) - disk_xOffset)/2, 0.99*(cooling_tube_diameter - 2*cooling_tube_thickness) / 2, 0.99*(cooling_tube_diameter - 2*cooling_tube_thickness) / 2);
        Volume Log_water_cooling_right("Log_water_cooling_right" + std::to_string(row), sol_water_cooling_right, description.material("Water"));

        Log_water_cooling_right.setVisAttributes(description.visAttributes("TOFWater"));

        //place Log_water_cooling_right in assembly
        Transform3D  transfrowCR(RotationZYX(row%2==0 ? - M_PI : 0,0, 0), Translation3D(( (nSensRight) * segmentlength / 2.0 + segmentlength / 2.0 + (numSensorRightAdd* segmentlength)), row>0 ? (row*fullsensor_width) - (fullsensor_width/2.0) : (row*fullsensor_width) + (fullsensor_width/2.0), disk_zPos));
        if(row!=0){
          assembly.placeVolume(Log_cooling_tube_right, transfrowCR);
          assembly.placeVolume(Log_water_cooling_right, transfrowCR);
        }
      } else {
          // create mother volume with space for numSensorsRow towers along x-axis
          Box TTLDetRowSolid(numSensorsRow * segmentlength / 2.0,fullsensor_width / 2.0,thicknessDet_SH / 2.0);
          Volume TTLDetRowLogical_Front("TTLDetRowLogical_Front" + std::to_string(row), TTLDetRowSolid, air);
          TTLDetRowLogical_Front.setVisAttributes(description.visAttributes("InvisibleWithDaughters"));

          //place log_sensor_and_readout_C in TTLDetRowLogical_Front
          for(int sensor = 0; sensor < numSensorsRow; sensor++) {
            Box box_sensor_and_readout_C(segmentlength / 2,fullsensor_width / 2,thicknessDet_SH / 2);
            Volume log_sensor_and_readout_C("log_sensor_and_readout_CF_" + std::to_string(row)+"_"+ std::to_string(sensorCount), box_sensor_and_readout_C, air);

            log_sensor_and_readout_C.setVisAttributes(description.visAttributes("InvisibleWithDaughters"));

            log_sensor_and_readout_C.placeVolume(log_sensor_stack1, Position(0, -fullsensor_width/2 + baseplate_width/2, -thicknessDet_SH/2+thicknessDet1/2));
            log_sensor_and_readout_C.placeVolume(log_sensor_stack2, Position(0, -fullsensor_width/2 + baseplate_width/2, -thicknessDet_SH/2+thicknessDet1+thicknessLayer[5]+thicknessDet2/2));
            log_sensor_and_readout_C.placeVolume(log_SH_stack, Position(0,  fullsensor_width/2 - baseplate_width/4, 0));

            // make active sensor solid and volume
            int senslayer=5;
            Box sol_ACLGAD_Sens(lengthLayer[senslayer] / 2,widthLayer[senslayer] / 2,thicknessLayer[senslayer] / 2);
            Volume log_ACLGAD_Sens_C("log_ACLGAD_Sens_CF_" + std::to_string(row)+"_"+ std::to_string(sensorCount), sol_ACLGAD_Sens, materialLayer[senslayer]);
            log_ACLGAD_Sens_C.setVisAttributes(description.visAttributes("TOFActiveMat"));

            pv = log_sensor_and_readout_C.placeVolume(log_ACLGAD_Sens_C, Position(0, -fullsensor_width/2 + baseplate_width/2 -offsetLayer[senslayer], -thicknessDet_SH/2+thicknessDet1+thicknessLayer[5]/2));

            pv.addPhysVolID("sensor", sensorCount);
            sens.setType("tracker");
            log_ACLGAD_Sens_C.setSensitiveDetector(sens);
            sensitives[m_nam].push_back(pv);

            // -------- create a measurement plane for the tracking surface attched to the sensitive volume -----
            Vector3D u(0., 0., -1.);
            Vector3D v(-1., 0., 0.);
            Vector3D n(0., 1., 0.);

            // compute the inner and outer thicknesses that need to be assigned to the tracking surface
            // depending on wether the support is above or below the sensor
            double inner_thickness = 0;
            double outer_thickness = thicknessLayer[senslayer];

            SurfaceType type(SurfaceType::Sensitive);
            VolPlane surf(log_ACLGAD_Sens_C, type, inner_thickness, outer_thickness, u, v, n); //,o ) ;
            volplane_surfaces[m_nam].push_back(surf);


            //place log_sensor_and_readout_C in TTLDetRowLogical_Front
            Position pos_sensor((sensor - (numSensorsRow - 1) / 2.0) * segmentlength, 0, 0);
            // PlacedVolume sensor_phys = 
            pv = TTLDetRowLogical_Front.placeVolume(log_sensor_and_readout_C, pos_sensor);
            sensorCount++;
          }

          //place TTLDetRowLogical_Front in assembly, roate every second row
          Transform3D  transfrow(RotationZYX(row%2==0 ? - M_PI : 0,0, 0), Translation3D(0, row * fullsensor_width, disk_zPos+offsetzFront));
          pv = assembly.placeVolume(TTLDetRowLogical_Front, transfrow);
          pv.addPhysVolID("layer", 0);
          pv.addPhysVolID("module", row);


          Volume TTLDetRowLogical_Back("TTLDetRowLogical_Back" + std::to_string(row), TTLDetRowSolid, air);
          TTLDetRowLogical_Back.setVisAttributes(description.visAttributes("InvisibleWithDaughters"));

          //place log_sensor_and_readout_C in TTLDetRowLogical_Back
          for(int sensor = 0; sensor < numSensorsRow; sensor++) {
            Box box_sensor_and_readout_C(segmentlength / 2,fullsensor_width / 2,thicknessDet_SH / 2);
            Volume log_sensor_and_readout_C("log_sensor_and_readout_CB_" + std::to_string(row)+"_"+ std::to_string(sensorCount), box_sensor_and_readout_C, air);

            log_sensor_and_readout_C.setVisAttributes(description.visAttributes("InvisibleWithDaughters"));

            log_sensor_and_readout_C.placeVolume(log_sensor_stack1, Position(0, -fullsensor_width/2 + baseplate_width/2, -thicknessDet_SH/2+thicknessDet1/2));
            log_sensor_and_readout_C.placeVolume(log_sensor_stack2, Position(0, -fullsensor_width/2 + baseplate_width/2, -thicknessDet_SH/2+thicknessDet1+thicknessLayer[5]+thicknessDet2/2));
            log_sensor_and_readout_C.placeVolume(log_SH_stack, Position(0,  fullsensor_width/2 - baseplate_width/4, 0));

            // make active sensor solid and volume
            int senslayer=5;
            Box sol_ACLGAD_Sens(lengthLayer[senslayer] / 2,widthLayer[senslayer] / 2,thicknessLayer[senslayer] / 2);
            Volume log_ACLGAD_Sens_C("log_ACLGAD_Sens_CB_" + std::to_string(row)+"_"+ std::to_string(sensorCount), sol_ACLGAD_Sens, materialLayer[senslayer]);
            log_ACLGAD_Sens_C.setVisAttributes(description.visAttributes("TOFActiveMat"));

            pv = log_sensor_and_readout_C.placeVolume(log_ACLGAD_Sens_C, Position(0, -fullsensor_width/2 + baseplate_width/2 -offsetLayer[senslayer], -thicknessDet_SH/2+thicknessDet1+thicknessLayer[5]/2));

            pv.addPhysVolID("sensor", sensorCount);
            sens.setType("tracker");
            log_ACLGAD_Sens_C.setSensitiveDetector(sens);
            sensitives[m_nam].push_back(pv);

            // -------- create a measurement plane for the tracking surface attched to the sensitive volume -----
            Vector3D u(0., 0., -1.);
            Vector3D v(-1., 0., 0.);
            Vector3D n(0., 1., 0.);

            // compute the inner and outer thicknesses that need to be assigned to the tracking surface
            // depending on wether the support is above or below the sensor
            double inner_thickness = 0;
            double outer_thickness = thicknessLayer[senslayer];

            SurfaceType type(SurfaceType::Sensitive);
            VolPlane surf(log_ACLGAD_Sens_C, type, inner_thickness, outer_thickness, u, v, n); //,o ) ;
            volplane_surfaces[m_nam].push_back(surf);


            //place log_sensor_and_readout_C in TTLDetRowLogical_Back
            Position pos_sensor((sensor - (numSensorsRow - 1) / 2.0) * segmentlength, 0, 0);
            // PlacedVolume sensor_phys = 
            pv = TTLDetRowLogical_Back.placeVolume(log_sensor_and_readout_C, pos_sensor);
            sensorCount++;
          }

          //place TTLDetRowLogical_Back on other side in assembly, roate every second row
          Transform3D  transfrow2(RotationZYX(row%2==0 ? - M_PI : 0,0,-M_PI), Translation3D(0, row * fullsensor_width, disk_zPos+offsetzBack));
          pv = assembly.placeVolume(TTLDetRowLogical_Back, transfrow2);
          pv.addPhysVolID("layer", 1);
          pv.addPhysVolID("module", row);


          Box sol_cutout_tube(2.5*sqrt(pow(disk_rMax,2)-pow( (abs(row)*fullsensor_width) - (fullsensor_width/2.0) + cooling_tube_diameter ,2))/2, (cooling_tube_diameter - 2*cooling_tube_thickness) / 2, (cooling_tube_diameter - 2*cooling_tube_thickness) / 2);
          Box box_cooling_tube(0.98*2*sqrt(pow(disk_rMax,2)-pow( (abs(row)*fullsensor_width) - (fullsensor_width/2.0) + cooling_tube_diameter ,2))/2,cooling_tube_diameter / 2, cooling_tube_diameter / 2);

          SubtractionSolid sol_cooling_tube(box_cooling_tube, sol_cutout_tube);

          Volume Log_cooling_tube("Log_cooling_tube" + std::to_string(row), sol_cooling_tube, description.material("Aluminum"));
          Log_cooling_tube.setVisAttributes(description.visAttributes("TOFAluminum"));


          Box sol_water_cooling(0.97*2*sqrt(pow(disk_rMax,2)-pow( (abs(row)*fullsensor_width) - (fullsensor_width/2.0) + cooling_tube_diameter ,2))/2,0.99*(cooling_tube_diameter - 2*cooling_tube_thickness) / 2,0.99*(cooling_tube_diameter - 2*cooling_tube_thickness) / 2);
          Volume Log_water_cooling("Log_water_cooling" + std::to_string(row), sol_water_cooling, description.material("Water"));

          Log_water_cooling.setVisAttributes(description.visAttributes("TOFWater"));

          //place Log_water_cooling in assembly
          Transform3D  transfrowC(RotationZYX(row%2==0 ? - M_PI : 0,0, 0), Translation3D(0, row>0 ? (row*fullsensor_width) - (fullsensor_width/2.0) : (row*fullsensor_width) + (fullsensor_width/2.0), disk_zPos));
          assembly.placeVolume(Log_cooling_tube, transfrowC);
          assembly.placeVolume(Log_water_cooling, transfrowC);
      }
    }
  }
  pv = motherVol.placeVolume(assembly, Position(0, 0, (reflect ? -1.0e-9 : 1.0e-9)));
  pv.addPhysVolID("system", det_id);
  sdet.setPlacement(pv);
  return sdet;
}

//@}
// clang-format off
DECLARE_DETELEMENT(epic_TOFEndcap, create_detector)
