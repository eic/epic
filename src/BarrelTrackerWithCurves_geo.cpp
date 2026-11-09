// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 Whitney Armstrong, Tuna Tasali, Sam Henry, Shujie Li

/** \addtogroup Trackers Trackers
 * \brief Type: **BarrelTrackerWithFrame**.
 * \author W. Armstrong
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
#include "TGDMLParse.h"

/*
This version is an adaptation of BarrelTrackerWithFrame_geo.cpp to handle the curved silicon surfaces of the outer barrels, and import
passive components from GDML files.

A curved component (identified by the word "Curved" in the name) is a segment of a cylinder. This is constructed by a series of flat segments
The number of segments is set in the xml file, together with the radius, anglular range, length, thickness, and the inner and outer thicknesses
used to set the sensitive surface. This approach is taken as it is not possibe to use the DD4HEP CylindricalGridPhiZ readout for a cylindrical surface
where the axis is displaced from the beam line. Therefore the curved surface is modelled as a series of flat segments.
Unlike in BarrelTrackerWithFrame, the position of components (related to the centre of the stave) must now be given in the xml files.

Passive components are modelled as a tessellated solid, imported from a .gdml file, with the filename given in the xml file.

8 September 2026. Sam Henry
*/

using namespace std;
using namespace dd4hep;
using namespace dd4hep::rec;
using namespace dd4hep::detail;

/** Barrel Tracker with space frame.
 *
 * - Optional "support" tag within the detector element.
 *
 * The shapes are created using createShape which can be one of many basic geomtries.
 * See the examples Check_shape_*.xml in
 * [dd4hep's examples/ClientTests/compact](https://github.com/AIDASoft/DD4hep/tree/master/examples/ClientTests/compact)
 * directory.
 *
 *
 * - Optional "frame" tag within the module element.
 *
 * \ingroup trackers
 *
 * \code
 * \endcode
 *
 *
 * @author Whitney Armstrong
 */
  struct SensorStackConfig {
    bool enabled = false;
    double copper = 0.0;
    double activeSi = 0.0;
    double inactiveSi = 0.0;

    double totalThickness() const {
        return copper + activeSi + inactiveSi;
    }
};

struct SensitiveSegment {
    PlacedVolume placement;
    VolPlane surface;
    std::string name;
};

struct SensitivePath {
    std::vector<PlacedVolume> placements;
    VolPlane surface;
    std::string name;
    string halfName;
};

struct BuiltAssembly {
    Assembly volume;
    std::vector<SensitivePath> sensitivePaths;
};

struct HalfRSUInstance {
    std::string name;
    PlacedVolume placement;
    std::vector<SensitivePath> sensitivePaths;
};

struct LASInstance {
    std::string name;
    PlacedVolume placement;
    std::vector<SensitivePath> sensitivePaths;
};

Transform3D getXMLTransform(xml_comp_t x_comp)
{
    xml_comp_t x_pos = x_comp.position(false);
    xml_comp_t x_rot = x_comp.rotation(false);
    Position pos(
        x_pos ? x_pos.x(0) : 0.0,
        x_pos ? x_pos.y(0) : 0.0,
        x_pos ? x_pos.z(0) : 0.0
    );
    RotationZYX rot(
        x_rot ? x_rot.z(0) : 0.0,
        x_rot ? x_rot.y(0) : 0.0,
        x_rot ? x_rot.x(0) : 0.0
    );
    return Transform3D(rot, pos);
}


void attachSensorSurface( DetElement& sensorDE,  const VolPlane& surface) {
    auto& params =  DD4hepDetectorHelper::ensureExtension<dd4hep::rec::VariantParameters>(sensorDE);
    params.set<std::string>("axis_definitions", "XYZ");
    volSurfaceList(sensorDE)->push_back(surface);
}


void attachLASDetElements( DetElement& moduleDE, const LASInstance& instance, int lasID) {
    // LAS is a daughter of the stave.
    DetElement lasDE( moduleDE, instance.name,  lasID );
    lasDE.setPlacement(instance.placement);

    // Each sensitive path is relative to the LAS.
    // Paths of length 2:
    //   [HalfRSU placement, sensitive segment placement]
    //
    // Paths of length 1:
    //   [direct sensitive segment placement]

    std::map<std::string, DetElement> halfElements;
    int halfID = 1;
    int directSensorID = 1;
    for (const auto& path : instance.sensitivePaths) {
        if (path.placements.size() == 1) {
            // Direct sensitive component inside LAS.
            DetElement sensorDE(lasDE,  "de_" + path.name,  directSensorID++  );
            sensorDE.setPlacement( path.placements[0]);
            attachSensorSurface( sensorDE, path.surface);
        } else if (path.placements.size() == 2) {
            PlacedVolume halfPV = path.placements[0];
            PlacedVolume sensorPV = path.placements[1];
            // Identify this HalfRSU instance by its physical placement within the LAS.
            const std::string& halfName = path.halfName;
            auto it = halfElements.find(halfName);
            if (it == halfElements.end()) {
                DetElement halfDE( lasDE,  halfName,  halfID++  );
                halfDE.setPlacement(halfPV);
                it = halfElements.emplace( halfName, halfDE).first;
            }
            DetElement& halfDE = it->second;
            DetElement sensorDE( halfDE, "de_" + path.name,  static_cast<int>( halfDE.children().size() + 1) );
            sensorDE.setPlacement(sensorPV);
            attachSensorSurface( sensorDE, path.surface);
        } else {
            throw std::runtime_error(
                "Unexpected LAS sensitive path depth");
        }
    }
}



void placeCurvedComponent(
    Detector& description,
    SensitiveDetector& sens,
    Volume& parent,
    xml_comp_t x_comp,
    const SensorStackConfig& stack,
    int& sensor_number,
    std::vector<PlacedVolume>& sensitive_placements,
    std::vector<VolPlane>& volplane_surface)
{
    const bool isSensitive = x_comp.isSensitive();
    const bool buildStack = stack.enabled && isSensitive;

    std::vector<std::string> layerNames;
    std::vector<std::string> layerMaterials;
    std::vector<double> layerThicknesses;
    std::vector<double> layerOffsets;
    std::vector<bool> layerSensitive;

    if (buildStack) {
        if (std::abs(x_comp.thickness() - stack.totalThickness()) > 1e-9 * mm) {  throw std::runtime_error("Curved component thickness differs from sensor stack"); }
        layerNames = {"InactiveSi", "ActiveSi", "Copper"};
        layerMaterials = {"Silicon", "Silicon", "Copper"};
        layerThicknesses = { stack.inactiveSi, stack.activeSi, stack.copper };
        layerOffsets = {stack.inactiveSi / 2 - stack.totalThickness() / 2,  stack.inactiveSi + stack.activeSi / 2 - stack.totalThickness() / 2,  stack.totalThickness() / 2 - stack.copper / 2};
        layerSensitive = {false, true, false};
    } else {
        layerNames = {"Component"};
        layerMaterials = {x_comp.materialStr()};
        layerThicknesses = {x_comp.thickness()};
        layerOffsets = {0.0};
        layerSensitive = {isSensitive};
    }

    const int nseg = x_comp.nsegments();
    if (nseg <= 0)
        throw std::runtime_error("Invalid number of curved segments");

    const double phi0 = x_comp.phi0();
    const double dphi = (x_comp.phi1() - phi0) / nseg;

    xml_comp_t x_pos = x_comp.position(false);
    xml_comp_t x_rot = x_comp.rotation(false);

    const Position offset(
        x_pos ? x_pos.x(0) : 0.0,
        x_pos ? x_pos.y(0) : 0.0,
        x_pos ? x_pos.z(0) : 0.0
    );

    const double rx = x_rot ? x_rot.x(0) : 0.0;
    const double ry = x_rot ? x_rot.y(0) : 0.0;
    const double rz = x_rot ? x_rot.z(0) : 0.0;

    for (size_t layer = 0; layer < layerNames.size(); ++layer) {
        const double radius = x_comp.radius() + layerOffsets[layer];
        const double thickness = layerThicknesses[layer];

        if (radius <= 0)
            throw std::runtime_error("Invalid curvature radius");

        for (int i = 0; i < nseg; ++i) {
            const double phiA = phi0 + i * dphi;
            const double phiB = phiA + dphi;
            const double phiMid = (phiA + phiB) / 2;
            const double xA = radius * std::sin(phiA);
            const double xB = radius * std::sin(phiB);
            const double zA = radius * std::cos(phiA);
            const double zB = radius * std::cos(phiB);
            const double width = std::hypot(xB - xA, zB - zA);
            const double midX = (xA + xB) / 2;
            const double midZ = (zA + zB) / 2;

            Box segmentShape( width / 2 * (radius - thickness / 2) / radius,  x_comp.length() / 2, thickness / 2 );
            const std::string name = x_comp.nameStr() + "_" + layerNames[layer] + "_" + std::to_string(i);
            Volume segment( name, segmentShape, description.material(layerMaterials[layer]) );
            segment.setRegion(description, x_comp.regionStr());
            segment.setLimitSet(description, x_comp.limitsStr());
            segment.setVisAttributes(description, x_comp.visStr());

            const Position position( midX + offset.x(),  offset.y(),  midZ + offset.z() );
            const RotationZYX rotation( rz, ry + phiMid, rx);
            PlacedVolume pv = parent.placeVolume(segment, Transform3D(rotation, position) );

            if (layerSensitive[layer]) {
                pv.addPhysVolID("sensor", sensor_number++);
                segment.setSensitiveDetector(sens);
		sensitive_placements.push_back(pv);
                const double inner = getAttrOrDefault<double>(x_comp, _Unicode(innerthickness), 0.5 * mm);
                const double outer = getAttrOrDefault<double>(x_comp, _Unicode(outerthickness), 0.5 * mm);
                SurfaceType type(SurfaceType::Sensitive);
                VolPlane surface(
                    segment, type, inner, outer,
                    Vector3D(-1, 0, 0),
                    Vector3D(0, -1, 0),
                    Vector3D(0, 0, 1)
                );
		volplane_surface.push_back(surface);
            }
        }
    }
}


BuiltAssembly buildHalfRSU(
    Detector& description,
    SensitiveDetector& sens,
    xml_comp_t x_half,
    const SensorStackConfig& stack,
    int& sensor_number) {

    Assembly half(x_half.nameStr());
    std::vector<PlacedVolume> placements;
    std::vector<VolPlane> surfaces;
    for (xml_coll_t ci(x_half, _Unicode(rsu_component)); ci; ++ci) {
        xml_comp_t x_comp = ci;
        placeCurvedComponent(description, sens, half, x_comp, stack, sensor_number, placements, surfaces );
    }
     if (placements.size() != surfaces.size()) { throw std::runtime_error( "HalfRSU sensitive placements/surfaces mismatch");  }

    std::vector<SensitivePath> paths;
    for (size_t i = 0; i < placements.size(); ++i) {
        paths.push_back({ {placements[i]},   surfaces[i],  placements[i].volume().name(), "" });
    }
    return {half, paths};
}

 
BuiltAssembly buildLAS(
    Detector& description,
    SensitiveDetector& sens,
    xml_comp_t x_las,
    const std::map<std::string, BuiltAssembly>& halfRSUs) {

    Assembly las(x_las.nameStr());
    std::vector<SensitivePath> paths;

    SensorStackConfig stack;
    stack.enabled = getAttrOrDefault<bool>( x_las, _Unicode(sensor_stack), false);
    stack.copper = getAttrOrDefault<double>( x_las, _Unicode(sensor_cu_thickness), 0.0);
    stack.activeSi = getAttrOrDefault<double>( x_las, _Unicode(sensor_active_si_thickness), 0.0);
    stack.inactiveSi = getAttrOrDefault<double>( x_las, _Unicode(sensor_inactive_si_thickness), 0.0);

    int halfIndex = 0;
    int sensor_number = 1;

    for (xml_coll_t ci(x_las, _Unicode(LAS_component));  ci; ++ci) {
        xml_comp_t x_comp = ci;
        const std::string ref = getAttrOrDefault<std::string>( x_comp, _Unicode(ref), "");
        if (!ref.empty()) {
            auto it = halfRSUs.find(ref);
            if (it == halfRSUs.end()) {
                throw std::runtime_error(
                    "Unknown HalfRSU reference: " + ref);
            }
            const BuiltAssembly& half = it->second;
            PlacedVolume halfPV = las.placeVolume( half.volume,  getXMLTransform(x_comp) );
            halfPV.addPhysVolID("halfrsu", halfIndex);
            for (const auto& halfPath : half.sensitivePaths) {
	    // Initialise all members, including VolPlane.
    		SensitivePath path{  halfPath.placements, halfPath.surface, x_comp.nameStr() + "_" + halfPath.name, x_comp.nameStr() };
            // Prepend the placement of this HalfRSU within the LAS.
    		path.placements.insert(path.placements.begin(), halfPV  );
                paths.push_back(std::move(path));
            }
            ++halfIndex;
        } else {
            // Curved component belonging directly to LAS
            std::vector<PlacedVolume> placements;
            std::vector<VolPlane> surfaces;
            placeCurvedComponent( description, sens, las, x_comp, stack, sensor_number, placements, surfaces );
            if (placements.size() != surfaces.size()) {
                throw std::runtime_error(
                    "LAS sensitive placement/surface mismatch");
            }
            for (size_t i = 0; i < placements.size(); ++i) {
                paths.push_back({ {placements[i]},  surfaces[i],  placements[i].volume().name(), ""  });
            }
        }
    }
    return {las, paths};
}

 
 

static Ref_t create_BarrelTrackerWithCurves(Detector& description, xml_h e,
                                            SensitiveDetector sens) {
  typedef vector<PlacedVolume> Placements;
  xml_det_t x_det = e;
  Material air    = description.air();
  int det_id      = x_det.id();
  string det_name = x_det.nameStr();
  DetElement sdet(det_name, det_id);

  map<string, Volume> volumes;
  map<string, Placements> sensitives;
  map<string, std::vector<VolPlane>> volplane_surfaces;
  map<string, std::array<double, 2>> module_thicknesses;
   map<std::string, std::vector<SensitivePath>> nested_sensitives;
  map<std::string, std::vector<HalfRSUInstance>> halfRSU_instances; // maybe obsolete
  map<std::string, std::vector<LASInstance>> las_instances;

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
  Assembly assembly(det_name);
  sens.setType("tracker");

  // Read HalfRSU definitions from XML
  std::map<std::string, BuiltAssembly> halfRSUs;
  for (xml_coll_t hi(x_det.parent(), _Unicode(HalfRSU)); hi; ++hi) {
      xml_comp_t x_half = hi;
      const std::string halfName = x_half.nameStr();
      SensorStackConfig stack;
      stack.enabled = getAttrOrDefault<bool>( x_half, _Unicode(sensor_stack), false);
      stack.copper = getAttrOrDefault<double>(x_half, _Unicode(sensor_cu_thickness), 0.0);
      stack.activeSi = getAttrOrDefault<double>( x_half, _Unicode(sensor_active_si_thickness), 0.0);
      stack.inactiveSi = getAttrOrDefault<double>( x_half, _Unicode(sensor_inactive_si_thickness), 0.0);
      if (stack.enabled &&  (stack.copper <= 0 ||  stack.activeSi <= 0 ||  stack.inactiveSi <= 0)) {
          throw std::runtime_error(  "Invalid sensor stack thicknesses in HalfRSU " +  halfName);
      }
      int sensor_number = 1;
      BuiltAssembly half = buildHalfRSU( description, sens, x_half, stack, sensor_number );
      halfRSUs.emplace(halfName, std::move(half));
  }

  // Construct LAS
  std::map<std::string, BuiltAssembly> lasAssemblies;
  for (xml_coll_t li(x_det, _Unicode(LAS)); li; ++li) { 
      xml_comp_t x_las = li;
      const std::string lasName = x_las.nameStr();
      BuiltAssembly las = buildLAS( description, sens, x_las, halfRSUs );
      auto inserted = lasAssemblies.emplace(lasName, std::move(las));
      if (!inserted.second) {
          throw std::runtime_error( "Duplicate LAS definition: " + lasName);
      } 
  }

  TGDMLParse* parser = new TGDMLParse();

  // loop over the modules
  for (xml_coll_t mi(x_det, _U(module)); mi; ++mi) {
    xml_comp_t x_mod = mi;
    string m_nam     = x_mod.nameStr();

    if (volumes.find(m_nam) != volumes.end()) {
      printout(ERROR, "BarrelTrackerWithFrame",
               string((string("Module with named ") + m_nam + string(" already exists."))).c_str());
      throw runtime_error("Logics error in building modules.");
    }

    const bool sensor_stack = getAttrOrDefault<bool>(x_mod, _Unicode(sensor_stack), false);
    const double copper_thickness =
        getAttrOrDefault<double>(x_mod, _Unicode(sensor_cu_thickness), 0.0);
    const double active_si_thickness =
        getAttrOrDefault<double>(x_mod, _Unicode(sensor_active_si_thickness), 0.0);
    const double inactive_si_thickness =
        getAttrOrDefault<double>(x_mod, _Unicode(sensor_inactive_si_thickness), 0.0);

    if (sensor_stack &&
        (copper_thickness <= 0 || active_si_thickness <= 0 || inactive_si_thickness <= 0)) {
      throw runtime_error("CurvedTrackerBarrel: sensor stack thicknesses must be positive");
    }

    SensorStackConfig stack;
    stack.enabled = sensor_stack;
    stack.copper = copper_thickness;
    stack.activeSi = active_si_thickness;
    stack.inactiveSi = inactive_si_thickness;

    int ncomponents   = 0;
    int sensor_number = 1;

    // the module assembly volume
    Assembly m_vol(m_nam);
    volumes[m_nam] = m_vol;
    m_vol.setVisAttributes(description.visAttributes(x_mod.visStr()));

    for (xml_coll_t mci(x_mod, _U(module_component)); mci; ++mci, ++ncomponents) {
      xml_comp_t x_comp  = mci;
      xml_comp_t x_pos   = x_comp.position(false);
      xml_comp_t x_rot   = x_comp.rotation(false);
      const string c_nam = _toString(ncomponents, "component%d");

	// Place LAS
      const std::string ref = getAttrOrDefault<std::string>( x_comp, _Unicode(ref), "");
      if (!ref.empty()) {
	 auto it = lasAssemblies.find(ref);
         if (it == lasAssemblies.end()) {  
             throw std::runtime_error( "Unknown LAS reference: " + ref);
         }
         const BuiltAssembly& las = it->second;
         const int lasIndex = static_cast<int>(las_instances[m_nam].size());
         PlacedVolume lasPV = m_vol.placeVolume( las.volume, getXMLTransform(x_comp)  );
         lasPV.addPhysVolID("las", lasIndex);
         LASInstance instance;
         instance.name = x_comp.nameStr();
         instance.placement = lasPV;
         instance.sensitivePaths = las.sensitivePaths;  // Paths are relative to the LAS assembly.
         las_instances[m_nam].push_back(std::move(instance)  );
	
      } else  if (x_comp.nameStr().find("Curved") != std::string::npos) {

		cout << "place curved component  " << getAttrOrDefault<std::string>(x_comp, _Unicode(name), " ") << endl;
       		placeCurvedComponent( description, sens, m_vol, x_comp, stack, sensor_number, sensitives[m_nam], volplane_surfaces[m_nam]);
          

      } else { // not a curved component

        // New code that constructs component from GDML files. Import the GDML file from the "file" attribute of the module_component
        std::string gdml_file = getAttrOrDefault<std::string>(x_comp, _Unicode(file), " ");

        Volume c_vol(c_nam);
        c_vol = parser->GDMLReadFile(gdml_file.c_str());



        //check the validity of the volume
        if (!c_vol.isValid()) {
          printout(WARNING, "BarrelTrackerOuter", "%s", gdml_file.c_str());
          printout(WARNING, "BarrelTrackerOuter", "c_vol invalid, GDML parser failed!");
          std::_Exit(EXIT_FAILURE);
        }
        c_vol.import();
        c_vol.setMaterial(description.material(x_comp.materialStr()));
        //risky bit, might quickly fill up memory if too many solids are imported
        TessellatedSolid c_sol(c_vol.solid());
        //note: c_sol gets casted automatically to parent class TGeoTessellated by the copy constructor. IDK why?
        c_sol->CloseShape(true, true, true); //otherwise you get an infinite bounding box
        c_sol->CheckClosure(
            true, true); //fix any flipped orientation in facets, the second 'true' is for verbose
        c_vol.setSolid(c_sol);

        c_vol.setRegion(description, x_comp.regionStr());
        c_vol.setLimitSet(description, x_comp.limitsStr());
        c_vol.setVisAttributes(description, x_comp.visStr());

        //place the volume
        if (x_pos && x_rot) {
          Position c_pos(x_pos.x(0), x_pos.y(0), x_pos.z(0));
          RotationZYX c_rot(x_rot.z(0), x_rot.y(0), x_rot.x(0) + M_PI / 2);
          pv = m_vol.placeVolume(c_vol, Transform3D(c_rot, c_pos));
        } else if (x_rot) {
          Position c_pos(0, 0, 0);
          pv = m_vol.placeVolume(
              c_vol,
              Transform3D(RotationZYX(x_rot.z(0), x_rot.y(0), x_rot.x(0) + M_PI / 2), c_pos));
        } else if (x_pos) {
          pv = m_vol.placeVolume(c_vol, Position(x_pos.x(0), x_pos.y(0), x_pos.z(0)));
        } else {
          RotationZYX c_rot(0, 0, +M_PI / 2);
          pv = m_vol.placeVolume(c_vol, Transform3D(c_rot, Position(0, 0, 0)));
        }
        c_vol.setVisAttributes(description, x_comp.visStr());
      }
    }
  }


  delete parser;

  // now build the layers
  for (xml_coll_t li(x_det, _U(layer)); li; ++li) {
    xml_comp_t x_layer  = li;
    xml_comp_t x_barrel = x_layer.child(_U(barrel_envelope));
    xml_comp_t x_layout = x_layer.child(_U(rphi_layout));
    xml_comp_t z_layout = x_layer.child(_U(z_layout)); // Get the <z_layout> element.
    int lay_id          = x_layer.id();
    string m_nam        = x_layer.moduleStr();
    string lay_nam      = det_name + _toString(x_layer.id(), "_layer%d");
    Tube lay_tub(x_barrel.inner_r(), x_barrel.outer_r(), x_barrel.z_length() / 2.0);
    Volume lay_vol(lay_nam, lay_tub, air); // Create the layer envelope volume.
    Position lay_pos(0, 0, getAttrOrDefault(x_barrel, _U(z0), 0.));
    lay_vol.setVisAttributes(description.visAttributes(x_layer.visStr()));

    double phi0     = x_layout.phi0();     // Starting phi of first module.
    double phi_tilt = x_layout.phi_tilt(); // Phi tilt of a module.
    double rc       = x_layout.rc();       // Radius of the module center.
    int nphi        = x_layout.nphi();     // Number of modfules in phi.
    double rphi_dr  = x_layout.dr();       // The delta radius of every other module.
    double phi_incr = (M_PI * 2) / nphi;   // Phi increment for one module.
    double phic     = phi0;                // Phi of the module center.
    double z0       = z_layout.z0();       // Z position of first module in phi.
    double nz       = z_layout.nz();       // Number of modules to place in z.
    double z_dr     = z_layout.dr();       // Radial displacement parameter, of every other module.

    Volume module_env = volumes[m_nam];
    DetElement lay_elt(sdet, lay_nam, lay_id);
    Placements& sensVols = sensitives[m_nam];

    // the local coordinate systems of modules in dd4hep and acts differ
    // see http://acts.web.cern.ch/ACTS/latest/doc/group__DD4hepPlugins.html
    auto& layerParams =
        DD4hepDetectorHelper::ensureExtension<dd4hep::rec::VariantParameters>(lay_elt);

    for (xml_coll_t lmat(x_layer, _Unicode(layer_material)); lmat; ++lmat) {
      xml_comp_t x_layer_material = lmat;
      DD4hepDetectorHelper::xmlToProtoSurfaceMaterial(x_layer_material, layerParams,
                                                      "layer_material");
    }

    // Z increment for module placement along Z axis.
    // Adjust for z0 at center of module rather than
    // the end of cylindrical envelope.
    double z_incr = nz > 1 ? (2.0 * z0) / (nz - 1) : 0.0;
    // Starting z for module placement along Z axis.
    double module_z = -z0;
    int module      = 1;

    // Loop over the number of modules in phi.
    for (int ii = 0; ii < nphi; ii++) {
      double dx = z_dr * std::cos(phic + phi_tilt); // Delta x of module position.
      double dy = z_dr * std::sin(phic + phi_tilt); // Delta y of module position.
      double x  = rc * std::cos(phic);              // Basic x module position.
      double y  = rc * std::sin(phic);              // Basic y module position.

      // Loop over the number of modules in z.
      for (int j = 0; j < nz; j++) {
        string module_name = _toString(module, "module%d");
        DetElement mod_elt(lay_elt, module_name, module);

        Transform3D tr(RotationZYX(0, ((M_PI / 2) - phic - phi_tilt), -M_PI / 2),
                       Position(x, y, module_z));

        pv = lay_vol.placeVolume(module_env, tr);
        pv.addPhysVolID("module", module);
        mod_elt.setPlacement(pv);
        
        int lasID = 1;
	for (const auto& instance : las_instances[m_nam]) {
	    attachLASDetElements( mod_elt,  instance,  lasID++  );
  	}
        
    	int halfID = 1;
	for (const auto& instance : halfRSU_instances[m_nam]) {
		DetElement halfDE(mod_elt, instance.name,  halfID++);
    		halfDE.setPlacement(instance.placement);
	        int segmentID = 1;
    		for (const auto& path : instance.sensitivePaths) {
			if (path.placements.size() != 1) {
            			throw std::runtime_error( "Expected direct sensitive segment inside HalfRSU");
        		}
			PlacedVolume sensPV = path.placements.front();
			DetElement sensorDE( halfDE,  "de_" + path.name,   segmentID++ );
			sensorDE.setPlacement(sensPV);
			auto& params2 =  DD4hepDetectorHelper::ensureExtension<dd4hep::rec::VariantParameters>(sensorDE);
			params2.set<std::string>("axis_definitions", "XYZ");
			volSurfaceList(sensorDE)->push_back(path.surface);
    		}
	}    
    
        for (size_t ic = 0; ic < sensVols.size(); ++ic) {
          PlacedVolume sens_pv = sensVols[ic];
          DetElement comp_de(mod_elt, std::string("de_") + sens_pv.volume().name(), module);
          comp_de.setPlacement(sens_pv);

          auto& comp_de_params =
              DD4hepDetectorHelper::ensureExtension<dd4hep::rec::VariantParameters>(comp_de);
          comp_de_params.set<string>("axis_definitions", "XYZ");
          // comp_de.setAttributes(description, sens_pv.volume(), x_layer.regionStr(), x_layer.limitsStr(),
          //                       xml_det_t(xmleles[m_nam]).visStr());
          //

          volSurfaceList(comp_de)->push_back(volplane_surfaces[m_nam][ic]);
        }

        /// Increase counters etc.
        module++;
        // Adjust the x and y coordinates of the module.
        x += dx;
        y += dy;
        // Flip sign of x and y adjustments.
        dx *= -1;
        dy *= -1;
        // Add z increment to get next z placement pos.
        module_z += z_incr;
      }
      phic += phi_incr; // Increment the phi placement of module.
      rc += rphi_dr;    // Increment the center radius according to dr parameter.
      rphi_dr *= -1;    // Flip sign of dr parameter.
      module_z = -z0;   // Reset the Z placement parameter for module.
    }
    // Create the PhysicalVolume for the layer.
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
  return sdet;
}

//@}
// clang-format off
DECLARE_DETELEMENT(epic_CurvedTrackerBarrel,   create_BarrelTrackerWithCurves)
