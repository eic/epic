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

using namespace std;
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
  typedef vector<PlacedVolume> Placements;
  xml_det_t                    x_det    = e;
  // Material                     vacuum   = description.vacuum();
  int                          det_id   = x_det.id();
  string                       det_name = x_det.nameStr();
  bool                         reflect  = x_det.reflect(false);
  DetElement                   sdet(det_name, det_id);
  Assembly                     assembly(det_name);

  Material air = description.material("Air");
  // Volume      assembly    (det_name,Box(10000,10000,10000),vacuum);
  Volume                             motherVol = description.pickMotherVolume(sdet);
  int                                m_id = 0;//, c_id = 0, n_sensor = 0;
  map<string, Volume>                modules;
  map<string, Placements>            sensitives;
  map<string, std::vector<VolPlane>> volplane_surfaces;
  map<string, std::array<double, 2>> module_thicknesses;
  PlacedVolume                       pv;

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

  for (xml_coll_t mi(x_det, _U(module)); mi; ++mi, ++m_id) {
    xml_comp_t x_mod = mi;
    string     m_nam = x_mod.nameStr();
    xml_comp_t diskdimension       = x_mod.child(_Unicode(diskdimension));

    double     disk_zPos            = getAttrOrDefault(diskdimension, _Unicode(zPos), 0.);
    double     disk_rMin            = getAttrOrDefault(diskdimension, _Unicode(rMin), 0.);
    double     disk_rMax            = getAttrOrDefault(diskdimension, _Unicode(rMax), 100.);
    double     disk_xOffset         = getAttrOrDefault(diskdimension, _Unicode(xOffset), 0.);
    double     disk_det_height      = getAttrOrDefault(diskdimension, _Unicode(det_height), 0.);
    double     cooling_plate_height = getAttrOrDefault(diskdimension, _Unicode(cooling_plate_height), 0.);
    double     cooling_tube_thickness    = getAttrOrDefault(diskdimension, _Unicode(wallthickness_coolingtube), 0.);
    double     cooling_tube_diameter    = getAttrOrDefault(diskdimension, _Unicode(diameter_coolingtube), 0.);
    // cout << "\tdisk_rMin = " << disk_rMin
    // << "\n\tdisk_rMax = " << disk_rMax
    // << "\n\tdisk_xOffset = " << disk_xOffset
    // << "\n\tdisk_zPos = " << disk_zPos
    // << "\n\tdisk_det_height = " << disk_det_height
    // << "\n\tdisk_det_height = " << disk_det_height
    // << "\n\tcooling_tube_thickness = " << cooling_tube_thickness
    // << "\n\tcooling_tube_diameter = " << cooling_tube_diameter
    // << endl;

    Tube beampipe_cutout(0, disk_rMin,disk_det_height, 0, 2. * M_PI);
    // Tube ttl_disk(0, disk_rMax,disk_det_height / 2., 0, 2. * M_PI);
    Tube ttl_plate_1(0, disk_rMax,cooling_plate_height / 2., 0, 2. * M_PI);
    Tube ttl_plate_2(0, disk_rMax,cooling_plate_height / 2., 0, 2. * M_PI);

    Position cutoutOffset(disk_xOffset, 0,0);
    SubtractionSolid ttl_cooling_plate_1(ttl_plate_1, beampipe_cutout,cutoutOffset);
    SubtractionSolid ttl_cooling_plate_2(ttl_plate_2, beampipe_cutout,cutoutOffset);

    Solid     ttl_cooling_plate_solid_1 = ttl_cooling_plate_1;
    Solid     ttl_cooling_plate_solid_2 = ttl_cooling_plate_2;

    Material disk_mat = description.material(diskdimension.materialStr());


    Volume ttl_cooling_plate_volume_1("ttl_cooling_plate_1", ttl_cooling_plate_solid_1, disk_mat);
    Volume ttl_cooling_plate_volume_2("ttl_cooling_plate_2", ttl_cooling_plate_solid_2, disk_mat);
    ttl_cooling_plate_volume_1.setVisAttributes(description.visAttributes("TOFAluminum"));
    ttl_cooling_plate_volume_2.setVisAttributes(description.visAttributes("TOFAluminum"));


    pv = assembly.placeVolume(ttl_cooling_plate_volume_1, Position(0, 0, disk_zPos-cooling_tube_diameter/2.-cooling_plate_height/2.));
    pv = assembly.placeVolume(ttl_cooling_plate_volume_2, Position(0, 0, disk_zPos+cooling_tube_diameter/2+cooling_plate_height/2.));

    // Sensor Module:
    double sensor_width = 21.2 * mm;
    double sensor_length = 42.0 * mm;
    double baseplate_length = 43.1 * mm;
    double baseplate_width = 56.5 * mm / 2;

    const int nLayers = 8;
    std::string strLayerName[nLayers] = {
      "ThermalPad",
      "ALN",
      "LairdFilm",
      "ROC",
      "Solder",
      "Sensor",
      "Epoxy",
      "AIN"
    };
    Material materialLayer[nLayers] = {
      description.material("Graphite"), description.material("AluminumNitrate"), description.material("Graphite"), description.material("Plexiglass"), description.material("Tin"), description.material("Silicon"), description.material("Epoxy"), description.material("AluminumNitrate")
    };
    double thicknessLayer[nLayers] = {
        0.25 * mm, 0.79 * mm, 0.08 * mm, 0.25 * mm, 0.03 * mm, 0.3 * mm, 0.08 * mm, 0.51 * mm};
    double widthLayer[nLayers] = {
        baseplate_width - 0.3 * mm,
        baseplate_width,
        sensor_width + 1 * mm,
        sensor_width + 1 * mm,
        sensor_width - 0.2 * mm,
        sensor_width,
        sensor_width,
        baseplate_width - 4 * mm};
    double offsetLayer[nLayers] = {
        0.3 * mm / 2,
        0,
        (baseplate_width - widthLayer[2]) / 2 - 0.2 * mm,
        (baseplate_width - widthLayer[3]) / 2 - 0.2 * mm,
        (baseplate_width - widthLayer[4]) / 2 - 0.1 * mm,
        (baseplate_width - widthLayer[5] - 0.1 * mm) / 2,
        (baseplate_width - widthLayer[6] - 0.1 * mm) / 2,
        4 * mm / 2 - 0.2 * mm};
    double lengthLayer[nLayers] = {
        baseplate_length - 0.2 * mm,
        baseplate_length,
        sensor_length + 0.2 * mm,
        sensor_length + 0.2 * mm,
        sensor_length - 0.2 * mm,
        sensor_length,
        sensor_length,
        baseplate_length - 0.2 * mm};
    bool layerActive[nLayers] = {
        false, false, false, false, false, true, false, false};

    double thicknessDet = 0;
    for (int ilay = 0; ilay < nLayers; ilay++)
    {
      thicknessDet += thicknessLayer[ilay];
    }

    double segmentlength = baseplate_length;  //(detlength - 10 * cm) / 6;//
    Box box_sensor_stack(baseplate_length / 2, baseplate_width / 2, thicknessDet / 2);
    Solid     sol_sensor_stack = box_sensor_stack;
    Volume log_sensor_stack("log_sensor_stack", sol_sensor_stack, air);
    log_sensor_stack.setVisAttributes(description.visAttributes("TOFActiveMat"));


    double z_start = -thicknessDet / 2;
    for (int ilay = 0; ilay < nLayers; ilay++)
    {
      const std::string layer_name = "sensor_stack_" + strLayerName[ilay];
      const std::string layer_name_Solid = "sol_" + layer_name;

      Box box_Module_Layer_Raw(lengthLayer[ilay] / 2,widthLayer[ilay] / 2,thicknessLayer[ilay] / 2);
      Solid     sol_Module_Layer_Raw = box_Module_Layer_Raw;
      Volume Log_Layer(layer_name + "_Log", sol_Module_Layer_Raw, materialLayer[ilay]);

      pv = log_sensor_stack.placeVolume(Log_Layer, Position(0, -offsetLayer[ilay], z_start + thicknessLayer[ilay] / 2));
      z_start += thicknessLayer[ilay];
      Log_Layer.setVisAttributes(description.visAttributes("TOFLayers"));

      if (layerActive[ilay])
      {
        Log_Layer.setSensitiveDetector(sens);
        sens.setType("tracker");
        sensitives[m_nam].push_back(pv);
        pv.addPhysVolID("sensor", 1);
        // -------- create a measurement plane for the tracking surface attched to the sensitive volume -----
        Vector3D u(0., 0., -1.);
        Vector3D v(-1., 0., 0.);
        Vector3D n(0., 1., 0.);

        // compute the inner and outer thicknesses that need to be assigned to the tracking surface
        // depending on wether the support is above or below the sensor
        double inner_thickness = 0;
        double outer_thickness = thicknessLayer[ilay];

        SurfaceType type(SurfaceType::Sensitive);

        // if( isStripDetector )
        //  type.setProperty( SurfaceType::Measurement1D , true ) ;

        VolPlane surf(Log_Layer, type, inner_thickness, outer_thickness, u, v, n); //,o ) ;
        volplane_surfaces[m_nam].push_back(surf);

      //   //--------------------------------------------
      }
    }
  
  
    const int nLayers_SH = 4;

    double baseSH_width = baseplate_width / 2;
    std::string strLayerName_SH[nLayers_SH] = {
        "ThermalPad",
        "HighSpeedBoard",
        "ConnectorSpace",
        "Powerboard"};

    Material materialLayer_SH[nLayers_SH] = {
      description.material("Graphite"), description.material("Polystyrene"), air, description.material("Polystyrene")
    };
    double thicknessLayer_SH[nLayers_SH] = {
        0.25 * mm,
        1.00 * mm,
        1.50 * mm,
        3.10 * mm};
    double widthLayer_SH[nLayers_SH] = {
        baseSH_width - 0.2 * mm,
        baseSH_width - 0.2 * mm,
        baseSH_width - 0.35 * mm,
        baseSH_width};
    double offsetLayer_SH[nLayers_SH] = {
        0.2 * mm / 2,
        0.2 * mm / 2,
        0.35 * mm / 2,
        0};
    double lengthLayer_SH[nLayers_SH] = {
        baseplate_length,
        baseplate_length,
        baseplate_length,
        baseplate_length};

    double thicknessDet_SH = 0;
    for (int ilay = 0; ilay < nLayers_SH; ilay++)
    {
      thicknessDet_SH += thicknessLayer_SH[ilay];
    }

    Box sol_SH_stack(baseplate_length / 2,baseSH_width / 2, thicknessDet_SH / 2);
    Volume log_SH_stack("log_SH_stack", sol_SH_stack, air);
    log_SH_stack.setVisAttributes(description.visAttributes("TOFSensorAndReadoutLadder"));

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

    double fullsensor_width = baseSH_width+baseplate_width;

    Box box_sensor_and_readout(segmentlength / 2,fullsensor_width / 2,thicknessDet_SH / 2);
    Volume log_sensor_and_readout("log_sensor_and_readout", box_sensor_and_readout, air);

    log_sensor_and_readout.setVisAttributes(description.visAttributes("InvisibleWithDaughters"));

    log_sensor_and_readout.placeVolume(log_sensor_stack, Position(0, -fullsensor_width/2 + baseplate_width/2, -thicknessDet_SH/2+thicknessDet/2));
    log_sensor_and_readout.placeVolume(log_SH_stack, Position(0,  fullsensor_width/2 - baseplate_width/4, 0));

    // RegisterPhysicalVolume(new G4PVPlacement(0, G4ThreeVector(0, fullsensor_width/2 - baseplate_width/4, 0),
    //                     log_SH_stack, "ServiceHybridPlacedPhysical", log_sensor_and_readout, false, 0, overlapcheck_sector),false);


  double offsetzFront = cooling_tube_diameter/2 + cooling_plate_height + thicknessDet_SH / 2;
  double offsetzBack = -cooling_tube_diameter/2 - cooling_plate_height - thicknessDet_SH / 2;
  // number of towers in radial direction (on y axis)
  int rowYdir = (int) ( (disk_rMax-(fullsensor_width/2)) / fullsensor_width);
  for(int row=rowYdir;row>=-rowYdir;row--){
  // for(int row=0;row>=-rowYdir;row--){
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
        Volume TTLDetRowLeftLogical("TTLDetRowLeftLogical" + std::to_string(row), TTLDetRowLeftSolid, air);
        TTLDetRowLeftLogical.setVisAttributes(description.visAttributes("InvisibleWithDaughters"));

        //place log_sensor_and_readout in TTLDetRowLeftLogical
        for(int sensor = 0; sensor < nSensLeft; sensor++) {
            //place log_sensor_and_readout in TTLDetRowLeftLogical
            Position pos_sensor((sensor - (nSensLeft - 1) / 2.0) * segmentlength, 0, 0);
            // PlacedVolume sensor_phys = 
            TTLDetRowLeftLogical.placeVolume(log_sensor_and_readout, pos_sensor);
        }

        //place TTLDetRowLeftLogical in assembly, roate every second row
        Transform3D  transfrow(RotationZYX(row%2==0 ? - M_PI : 0,0, 0), Translation3D( - ( (nSensLeft) * segmentlength / 2.0 + segmentlength / 2.0 - (numSensorLeftAdd* segmentlength)), (row*fullsensor_width), disk_zPos+offsetzFront));
        assembly.placeVolume(TTLDetRowLeftLogical, transfrow);

        //place TTLDetRowLeftLogical on other side in assembly, roate every second row
        Transform3D  transfrow2(RotationZYX(row%2==0 ? - M_PI : 0,0,-M_PI), Translation3D( - ( (nSensLeft) * segmentlength / 2.0 + segmentlength / 2.0 - (numSensorLeftAdd* segmentlength)), (row*fullsensor_width), disk_zPos+offsetzBack));
        assembly.placeVolume(TTLDetRowLeftLogical, transfrow2);
      

        // we want an odd number of towers to be symmetrically centered around 0
        // create mother volume with space for numSensorsRow towers along x-axis
        int nSensRight = ((numSensorsRow-1) /2 - numSensorRightAdd);
        Box TTLDetRowRightSolid( nSensRight * segmentlength / 2.0,fullsensor_width / 2.0,thicknessDet_SH / 2.0);
        Volume TTLDetRowRightLogical("TTLDetRowRightLogical" + std::to_string(row), TTLDetRowLeftSolid, air);
        TTLDetRowRightLogical.setVisAttributes(description.visAttributes("InvisibleWithDaughters"));

        //place log_sensor_and_readout in TTLDetRowRightLogical
        for(int sensor = 0; sensor < nSensRight; sensor++) {
            //place log_sensor_and_readout in TTLDetRowRightLogical
            Position pos_sensor((sensor - (nSensRight - 1) / 2.0) * segmentlength, 0, 0);
            // PlacedVolume sensor_phys = 
            TTLDetRowRightLogical.placeVolume(log_sensor_and_readout, pos_sensor);
        }

        //place TTLDetRowRightLogical in assembly, roate every second row
        Transform3D  transfrowR(RotationZYX(row%2==0 ? - M_PI : 0,0, 0), Translation3D(( nSensRight * segmentlength / 2.0 + segmentlength / 2.0 + (numSensorRightAdd* segmentlength)), (row*fullsensor_width), disk_zPos+offsetzFront));
        assembly.placeVolume(TTLDetRowRightLogical, transfrowR);

        //place TTLDetRowRightLogical on other side in assembly, roate every second row
        Transform3D  transfrowR2(RotationZYX(row%2==0 ? - M_PI : 0,0,-M_PI), Translation3D(( nSensRight * segmentlength / 2.0 + segmentlength / 2.0 + (numSensorRightAdd* segmentlength)), (row*fullsensor_width), disk_zPos+offsetzBack));
        assembly.placeVolume(TTLDetRowRightLogical, transfrowR2);



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
        assembly.placeVolume(Log_cooling_tube_left, transfrowCL);
        assembly.placeVolume(Log_water_cooling_left, transfrowCL);



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
        assembly.placeVolume(Log_cooling_tube_right, transfrowCR);
        assembly.placeVolume(Log_water_cooling_right, transfrowCR);

    } else {
        // create mother volume with space for numSensorsRow towers along x-axis
        Box TTLDetRowSolid(numSensorsRow * segmentlength / 2.0,fullsensor_width / 2.0,thicknessDet_SH / 2.0);
        Volume TTLDetRowLogical("TTLDetRowLogical" + std::to_string(row), TTLDetRowSolid, air);
        TTLDetRowLogical.setVisAttributes(description.visAttributes("InvisibleWithDaughters"));

        //place log_sensor_and_readout in TTLDetRowLogical
        for(int sensor = 0; sensor < numSensorsRow; sensor++) {
            //place log_sensor_and_readout in TTLDetRowLogical
            Position pos_sensor((sensor - (numSensorsRow - 1) / 2.0) * segmentlength, 0, 0);
            // PlacedVolume sensor_phys = 
            TTLDetRowLogical.placeVolume(log_sensor_and_readout, pos_sensor);
        }

        //place TTLDetRowLogical in assembly, roate every second row
        Transform3D  transfrow(RotationZYX(row%2==0 ? - M_PI : 0,0, 0), Translation3D(0, row * fullsensor_width, disk_zPos+offsetzFront));
        assembly.placeVolume(TTLDetRowLogical, transfrow);

        //place TTLDetRowLogical on other side in assembly, roate every second row
        Transform3D  transfrow2(RotationZYX(row%2==0 ? - M_PI : 0,0,-M_PI), Translation3D(0, row * fullsensor_width, disk_zPos+offsetzBack));
        assembly.placeVolume(TTLDetRowLogical, transfrow2);



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
