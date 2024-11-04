// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Whitney Armstrong

/** \addtogroup Trackers Trackers
 * \brief Type: **BarrelTrackerWithFrame**.
 * \author W. Armstrong
 *
 * \ingroup trackers
 *
 * @{
 */

#include "XML/Layering.h"
#include "XML/Utilities.h"
#include "TGDMLParseBiggerFiles.h"
//#include "CADPlugins.h"
#include "FileLoaderHelper.h"
#include "ImportCADHelper.h" //include utilities
#include <fstream>

using namespace std;
using namespace dd4hep;
using namespace dd4hep::rec;
using namespace dd4hep::detail;

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

//some global variables
TessellatedSolid::Vertex xhat(1., 0., 0.);
TessellatedSolid::Vertex yhat(0., 1., 0.);
TessellatedSolid::Vertex zhat(0., 0., 1.);

class TriangularFacet
{
  public:
  //some attributes of the facet we are working with
  bool initialized = false;
  bool extruded = false;
  vector<TessellatedSolid::Vertex> vertices;
  TessellatedSolid::Vertex centroid;
  //indices of the vertices
  int iv0 = -1;
  int iv1 = -1;
  int iv2 = -1;
  //TessellatedSolid solid = NULL; //pointer to the solid object, somehow does not like nullptr?
  double facet_area = 0;
  TessellatedSolid::Vertex normal;
  TessellatedSolid::Vertex extrusionVector;
  //now this might be the weirdest varible
  //+1 if the position of the facet and its normal are both pointing up from the yz plane
  //-1 if one is pointing up and the other is pointing down or vice-versa
  int concavity_score = 0;
  void compute_facet_properties(TessellatedSolid::Facet facet, TessellatedSolid c_sol){
    //get the vertex indices
    iv0 = facet.GetVertexIndex(0);
    iv1 = facet.GetVertexIndex(1);
    iv2 = facet.GetVertexIndex(2);
    //add the vertices to a vector
    vertices.push_back(c_sol->GetVertex(iv0));
    vertices.push_back(c_sol->GetVertex(iv1));
    vertices.push_back(c_sol->GetVertex(iv2));
    centroid = 0.33333*(vertices[0] + vertices[1] + vertices[2]);
    facet_area = compute_area();
    normal = compute_normal(); 

    //assign the concavity score:
    double normalscore = TessellatedSolid::Vertex::Dot(xhat, normal);
    //shift the centroid to a new coordinate system, whose origin is to the same side of all components of the stave
    //this is a geometric requirement when determining concavity this way
    TessellatedSolid::Vertex shift(0., -1000., 1.);
    centroid -= shift;
    double centroidscore = TessellatedSolid::Vertex::Dot(xhat, centroid);
    concavity_score = sgn(normalscore*centroidscore);
    //set initialization flag    
    initialized = true;
  }
  TessellatedSolid::Vertex compute_normal()
  {
    //figure out the castings here
    TessellatedSolid::Vertex vec1 = vertices.at(1) - vertices.at(2);
    TessellatedSolid::Vertex vec2 = vertices.at(0) - vertices.at(2);
    
    //take the cross product
    //note that cross is stored in the struct vertex so I needed to access it the static way
    TessellatedSolid::Vertex normal_return = TessellatedSolid::Vertex::Cross(vec1, vec2);
    if(!normal_return.IsNormalized()) {
      normal_return.Normalize();
    }
    return normal_return;
  }
  double compute_area()
  {
    //figure out the castings here
    TessellatedSolid::Vertex vec1 = vertices.at(1) - vertices.at(2);
    TessellatedSolid::Vertex vec2 = vertices.at(0) - vertices.at(2);
    
    //take the cross product
    //note that cross is stored in the struct vertex so I needed to access it the static way
    TessellatedSolid::Vertex normal_return = TessellatedSolid::Vertex::Cross(vec1, vec2);
    return(normal_return.Mag()/2);
  }
  TessellatedSolid create_TriangularPrism(double extrusion_length) {
    
    if(!initialized) {
      printout(ERROR, "BarrelTrackerOuter_standardized", "Try running struct TriangularPrism construct_from_facet_and_solid()");
      throw runtime_error("Vertices not initialized! Triangular prisim construction failed.");
      return NULL;
    }
    TessellatedSolid extruded_prism("prism", 6);
    if(vertices.size() > 3) {
      printout(ERROR, "BarrelTrackerOuter_standardized", "Trying to construct triangular prism with more or less than 3 vertices");
      throw runtime_error("Triangular prisim construction failed.");
    }
    if(facet_area < 0.05) return NULL;
    else {
      extrusionVector = extrusion_length * normal;
      vector<TessellatedSolid::Vertex> extruded_vertices;
      for(auto& element : vertices) 
      {
        extruded_vertices.push_back(element + extrusionVector);
      }
      //vector<TessellatedSolid::Vertex> all_vertices(vertices.size() + extruded_vertices.size());
      //merge(vertices.begin(), vertices.end(), extruded_vertices.begin(), extruded_vertices.end(), 
          //all_vertices.begin()); 
      //Now add the facets:
      //depending on the extrusion, change the orientation
      if(extrusion_length < 0)
      {
        //Top and bottom
        extruded_prism.addFacet(vertices.at(2), vertices.at(1), vertices.at(0));
        //note the reversal of order to keep the normals well
        extruded_prism.addFacet(extruded_vertices.at(0), extruded_vertices.at(1), extruded_vertices.at(2));
        //sides
        for(unsigned long i = 0; i < vertices.size(); i++) 
        {
          int next_i = (i + 1) % vertices.size();
          extruded_prism.addFacet(vertices.at(i), vertices.at(next_i), extruded_vertices.at(next_i), extruded_vertices.at(i));
          //extruded_prism.addFacet(extruded_vertices.at(i), extruded_vertices.at(next_i), vertices.at(next_i), vertices.at(i)); 
        }
      }
      else if(extrusion_length > 0)
      {
        //Top and bottom
        extruded_prism.addFacet(vertices.at(0), vertices.at(1), vertices.at(2));
        //note the reversal of order to keep the normals well
        extruded_prism.addFacet(extruded_vertices.at(2), extruded_vertices.at(1), extruded_vertices.at(0));
        //sides
        for(unsigned long i = 0; i < vertices.size(); i++) 
        {
          int next_i = (i + 1) % vertices.size();
          //extruded_prism.addFacet(vertices.at(i), vertices.at(next_i), extruded_vertices.at(next_i), extruded_vertices.at(i));
          extruded_prism.addFacet(extruded_vertices.at(i), extruded_vertices.at(next_i), vertices.at(next_i), vertices.at(i)); 
        }
      }
      else
      {
        printout(ERROR, "BarrelTrackerOuter_standardized", "Trying to construct triangular prism with zero extrusion_length!!!");
        throw runtime_error("Triangular prisim construction failed.");
      }
      extruded = true;
      return extruded_prism;
    }
  }
};

/** Barrel Tracker imported from GDML file
 *
 *
 * The shapes are created using createShape which can be one of many basic geomtries.
 * See the examples Check_shape_*.xml in
 * [dd4hep's examples/ClientTests/compact](https://github.com/AIDASoft/DD4hep/tree/master/examples/ClientTests/compact)
 * directory.
 *
 *
 * 
 *
 * \ingroup trackers
 *
 * \code
 * \endcode
 *
 *
 * @author Whitney Armstrong, Tuna Tasali
 */

//a conversion function between Vector3D and TessellatedSolid::Vertex
Vector3D vertex_to_vector3D(TessellatedSolid::Vertex vertex) 
{
  Vector3D vec(vertex[0], vertex[1], vertex[2]);
  return(vec);
}

static Ref_t create_BarrelTrackerOuterStandardized(Detector& description, xml_h e, SensitiveDetector sens) {

  typedef vector<PlacedVolume> Placements;
  xml_det_t x_det = e;
  Material air    = description.air();
  //Material vacuum = description.vacuum();
  //Material silicon = description.material();
  int det_id      = x_det.id();
  string det_name = x_det.nameStr();
  DetElement sdet(det_name, det_id);

  map<string, Volume> volumes;
  vector<string> module_names;
  map<string, Placements> sensitives;
  map<string, std::vector<VolPlane>> volplane_surfaces;
  map<string, std::array<double, 2>> module_thicknesses;

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

  // dd4hep::xml::Dimension dimensions(x_det.dimensions());
  // Tube topVolumeShape(dimensions.rmin(), dimensions.rmax(), dimensions.length() * 0.5);
  // Volume assembly(det_name,topVolumeShape,air);
  Assembly assembly(det_name);

  //you can pick between these two options
  sens.setType("tracker");
  //sens.setType("calorimeter");
  //
  //

  // Loop over the suports
  for (xml_coll_t su(x_det, _U(support)); su; ++su) {
    xml_comp_t x_support     = su;
    double support_thickness = getAttrOrDefault(x_support, _U(thickness), 2.0 * mm);
    double support_length    = getAttrOrDefault(x_support, _U(length), 2.0 * mm);
    double support_rmin      = getAttrOrDefault(x_support, _U(rmin), 2.0 * mm);
    double support_zstart    = getAttrOrDefault(x_support, _U(zstart), 2.0 * mm);
    std::string support_name =
        getAttrOrDefault<std::string>(x_support, _Unicode(name), "support_tube");
    std::string support_vis = getAttrOrDefault<std::string>(x_support, _Unicode(vis), "AnlRed");
    xml_dim_t pos(x_support.child(_U(position), false));
    xml_dim_t rot(x_support.child(_U(rotation), false));
    Solid support_solid;
    if (x_support.hasChild(_U(shape))) {
      xml_comp_t shape(x_support.child(_U(shape)));
      string shape_type = shape.typeStr();
      support_solid     = xml::createShape(description, shape_type, shape);
    } else {
      support_solid = Tube(support_rmin, support_rmin + support_thickness, support_length / 2);
    }
    Transform3D tr =
        Transform3D(Rotation3D(), Position(0, 0, (support_zstart + support_length / 2)));
    if (pos.ptr() && rot.ptr()) {
      Rotation3D rot3D(RotationZYX(rot.z(0), rot.y(0), rot.x(0)));
      Position pos3D(pos.x(0), pos.y(0), pos.z(0));
      tr = Transform3D(rot3D, pos3D);
    } else if (pos.ptr()) {
      tr = Transform3D(Rotation3D(), Position(pos.x(0), pos.y(0), pos.z(0)));
    } else if (rot.ptr()) {
      Rotation3D rot3D(RotationZYX(rot.z(0), rot.y(0), rot.x(0)));
      tr = Transform3D(rot3D, Position());
    }
    Material support_mat = description.material(x_support.materialStr());
    Volume support_vol(support_name, support_solid, support_mat);
    support_vol.setVisAttributes(description.visAttributes(support_vis));
    pv = assembly.placeVolume(support_vol, tr);
    // pv = assembly.placeVolume(support_vol, Position(0, 0, support_zstart + support_length / 2));
  }



  // loop over the modules
  // to parse GDML files

  //needed for some reports
  double sensitive_area = 0;
  double total_surface_area = 0;

  TGDMLParseBiggerFiles* parser = new TGDMLParseBiggerFiles();
  for (xml_coll_t mi(x_det, _U(module)); mi; ++mi) {
    xml_comp_t x_mod = mi;
    string m_nam     = x_mod.nameStr();
    module_names.push_back(m_nam);
    if (volumes.find(m_nam) != volumes.end()) {
      printout(ERROR, "BarrelTrackerOuter",
               string((string("Module with named ") + m_nam + string(" already exists."))).c_str());
      throw runtime_error("Logics error in building modules.");
    }

    int ncomponents        = 0;
    int sensor_number      = 1;
    double total_thickness = 0;

    // Compute module total thickness from components
    xml_coll_t ci(x_mod, _U(module_component));
    for (ci.reset(), total_thickness = 0.0; ci; ++ci) {
      double thic = getAttrOrDefault(ci, _U(thickness), 0 * mm);
      total_thickness += thic;
    }
    Assembly m_vol(m_nam);
    volumes[m_nam] = m_vol;
    m_vol.setVisAttributes(description.visAttributes(x_mod.visStr()));

    double thickness_so_far = 0.0;
    double thickness_sum    = -total_thickness / 2.0;
    for (xml_coll_t mci(x_mod, _U(module_component)); mci; ++mci, ++ncomponents) {
      xml_comp_t x_comp  = mci;
      xml_comp_t x_pos   = x_comp.position(false);
      xml_comp_t x_rot   = x_comp.rotation(false);
      const string c_nam = _toString(ncomponents, "component%d");
      const string c_nam_mesh = _toString(ncomponents, "component%d_MESH");

      //new code that consturcts them from GDML files
      //import the GDML file from the "file" attribute of the module_component
      std::string gdml_file =
          getAttrOrDefault<std::string>(x_comp, _Unicode(file), " ");
      //printout(WARNING, "BarrelTrackerOuter", gdml_file);
      std::string given_name = getAttrOrDefault<std::string>(x_comp, _Unicode(name), " ");

      Volume c_vol(c_nam);
      //printout(WARNING, "BarrelTrackerOuter", "Parsing a large GDML file may lead to segfault or heap overflow.");
      c_vol = parser->GDMLReadFile(gdml_file.c_str());
      //check the validity of the volume
      if (!c_vol.isValid()) {
        printout(WARNING, "BarrelTrackerOuter", "%s", gdml_file.c_str());
        printout(WARNING, "BarrelTrackerOuter", "c_vol invalid, GDML parser failed!");
        std::_Exit(EXIT_FAILURE);
      }
      c_vol.import();
      c_vol.setMaterial(description.material(x_comp.materialStr()));
      printout(WARNING, "BarrelTrackerOuter", "%s", x_comp.materialStr().c_str());
      //risky bit, might quickly fill up memory if too many solids are imported
      TessellatedSolid c_sol(c_vol.solid());
      //note: c_sol gets casted automatically to parent class TGeoTessellated by the copy constructor. IDK why?
      c_sol->CloseShape(true, true, true); //otherwise you get an infinite bounding box
      c_sol->CheckClosure(true, true); //fix any flipped orientation in facets, the second 'true' is for verbose
      c_vol.setSolid(c_sol);
      c_vol.setRegion(description, x_comp.regionStr());
      c_vol.setLimitSet(description, x_comp.limitsStr());

      //an important offset variable that allows radial offsetting to prevent overlaps:
      const double radial_offset =
          getAttrOrDefault<double>(x_comp, _Unicode(offset), 0);
      
      // Utility variable for the relative z-offset based off the previous components
      const double zoff = thickness_sum + radial_offset / 2.0;

      // now, the code branches in the following way:
      // if the volume is sensitive, build the tesselated solid into thickened pixels of extruded polyogons
      // and place those under m_vol
      // else, just place c_vol under m_vol
      if (x_comp.isSensitive()) {
        //loop over the facets of c_vol to define volumes and set them as sensitive
        //note these facets won't be actual pixels, but rather just bits of the mesh
        //some variables to monitor if any sensitive area is lost
        double extrusion_length = x_comp.thickness(); //thickness will control the extrusion length
        //note that the negative sign has the same effect of flipping the normal of the inner side
        vector<TriangularFacet> accepted_facets;
        //calculate the concavity
        //double component_direction = 0;
        double costheta_threshold = 0.88;
        double facet_area_threshold = 0.05;
        for(int facet_index = 0; facet_index < c_sol->GetNfacets(); facet_index++)  {
          if (!(c_sol->GetFacet(facet_index).GetNvert() == 3)) throw runtime_error("BarrelTrackerOuterStandardized: Non triangular facets not supported. Please revise your mesh.");
          //construct a TriangularFacet from the facet
          TriangularFacet tri_facet;
          //printout(WARNING, "BarrelTrackerOuterSENSSSSSSSSSSSSS", gdml_file);
          tri_facet.compute_facet_properties(c_sol->GetFacet(facet_index), c_sol);
          double facet_area = tri_facet.facet_area;
          total_surface_area += facet_area;
          //determine if the facet is a good facet
          //note the standard in importing the cad models I defined when importing them:
            //stave long axis: z
            //stave thickness: y
            //stave width: x
            //extract some parameters
          double costheta = TessellatedSolid::Vertex::Dot(yhat, tri_facet.normal);
          //component_direction += sgn(costheta) * tri_facet.facet_area;
          bool is_good_facet = abs(costheta) >= costheta_threshold && 
                            facet_area > facet_area_threshold &&
                            tri_facet.concavity_score == -1; //always extrude the concave side
          if(!is_good_facet) continue;
          else accepted_facets.push_back(tri_facet);
        }
        //now work with the accepted prisms
        for(auto& tri_facet: accepted_facets) {
          sensitive_area += tri_facet.facet_area; //log the area
          //printout(WARNING, "BarrelTrackingOuter", "%f", TessellatedSolid::Vertex::Dot(yhat, tri_facet.normal));
          //now define a facet solid and place it under sc_vol
          TessellatedSolid extruded_facet = tri_facet.create_TriangularPrism(extrusion_length); //extrude the facet
          if(!extruded_facet) throw runtime_error("BarrelTrackerOuterStandardized: Error when extruding facet!");
          extruded_facet->CloseShape(true, true, true); //otherwise you get an infinite bounding box
          extruded_facet->CheckClosure(true, true); //fix any flipped orientation in facets, the second 'true' is for verbose
          Volume sc_vol_facet("facet" + to_string(sensor_number));
          sc_vol_facet.setSolid(extruded_facet); //note: the dereferenced pointer is also a pointer
          sc_vol_facet.setMaterial(description.material(x_comp.materialStr()));
          sc_vol_facet.setRegion(description, x_comp.regionStr());
          sc_vol_facet.setLimitSet(description, x_comp.limitsStr());
          //now place the volume
          RotationZYX c_rot(0, 0, -M_PI/2);
          pv = m_vol.placeVolume(sc_vol_facet, Transform3D(c_rot, Position(0, 0, zoff)));
          sc_vol_facet.setVisAttributes(description, x_comp.visStr());
          pv.addPhysVolID("sensor", sensor_number);
          sensor_number = sensor_number + 1;
          sc_vol_facet.setSensitiveDetector(sens);
          sensitives[m_nam].push_back(pv);

          //SURFACE WORK
          //module_thicknesses[m_nam] = {extrusion_length,
                                      //0};
          module_thicknesses[m_nam] = {thickness_so_far + x_comp.thickness() / 2.0,
                                        total_thickness - thickness_so_far - x_comp.thickness() / 2.0};                   
          // -------- create a measurement plane for the tracking surface attched to the sensitive volume -----
          vector<TessellatedSolid::Vertex> vertices = tri_facet.vertices;
          Vector3D u = vertex_to_vector3D(vertices.at(1) - vertices.at(2));
          Vector3D v = vertex_to_vector3D(vertices.at(0) - vertices.at(2));
          Vector3D n = vertex_to_vector3D(tri_facet.normal);
          Vector3D o = vertex_to_vector3D(0.333333333 * (vertices.at(0) + vertices.at(1) + vertices.at(2)));// + extrusion_length*n;
          // compute the inner and outer thicknesses that need to be assigned to the tracking surface
          // depending on whether the support is above or below the sensor
          double inner_thickness = module_thicknesses[m_nam][0];
          double outer_thickness = module_thicknesses[m_nam][1];
          SurfaceType type(SurfaceType::Sensitive);
          VolPlane surf(c_vol, type, inner_thickness, outer_thickness, u, v, n, o);
          volplane_surfaces[m_nam].push_back(surf);

        }
        //printout(WARNING, "BarrelTrackingOuter", "Please check the following values, they should be very close: ");
        //printout(WARNING, "BarrelTrackingOuter", "Total Area of Facets / 2 : %f", total_surface_area/2);
        //printout(WARNING, "BarrelTrackingOuter", "Total Sensitive Area: %f", sensitive_area/2);
      }
      else { // not a sensitive volume
        //place the volume
        if (x_pos && x_rot) {
          Position c_pos(x_pos.x(0), x_pos.y(0), x_pos.z(0) + zoff);
          RotationZYX c_rot(x_rot.z(0), x_rot.y(0), x_rot.x(0));
          pv = m_vol.placeVolume(c_vol, Transform3D(c_rot, c_pos));
        } else if (x_rot) {
          Position c_pos(0, 0, zoff);
          pv = m_vol.placeVolume(c_vol,
                                Transform3D(RotationZYX(x_rot.z(0), x_rot.y(0), x_rot.x(0)), c_pos));
        } else if (x_pos) {
          pv = m_vol.placeVolume(c_vol, Position(x_pos.x(0), x_pos.y(0), x_pos.z(0) + zoff));
        } else {
          //the c_rot is a temporary adjustment I added
          RotationZYX c_rot(0, 0, -M_PI/2);
          pv = m_vol.placeVolume(c_vol, Transform3D(c_rot, Position(0, 0, zoff)));
        
        }
        c_vol.setVisAttributes(description, x_comp.visStr());
        
      }
      //double thic = getAttrOrDefault<double>(x_comp, _U(thickness), 0 * mm);
      thickness_sum += radial_offset;
      thickness_so_far += radial_offset;
      // apply relative offsets in z-position used to stack components side-by-side
      if (x_pos) {
        thickness_sum += x_pos.z(0);
        thickness_so_far += x_pos.z(0);
      }
    }
    
  }
  delete parser;
  

  //dump reports on the components:
  ofstream area_report("Sens_Area_Report.txt");
  if(!area_report.is_open()) throw std::runtime_error("Unable to open or create the file: Sens_Area_Report.txt");
  area_report << "Sensitive area loss report from extrusion..." << std::endl;
  area_report << "Check if these values are close to each other (former will be slightly bigger):" << std::endl;
  area_report << "Total Surface Area / 2: " << total_surface_area/2 << std::endl;
  area_report << "Sensitive area: " << sensitive_area << std::endl;

  // now build the layers: changes added to allow multimodule layers
  for (xml_coll_t li(x_det, _U(layer)); li; ++li) {
    xml_comp_t x_layer  = li;
    xml_comp_t x_barrel = x_layer.child(_U(barrel_envelope));
    xml_comp_t x_layout = x_layer.child(_U(rphi_layout));
    xml_comp_t z_layout = x_layer.child(_U(z_layout)); // Get the <z_layout> element.
    int lay_id          = x_layer.id();
    string lay_nam      = det_name + _toString(x_layer.id(), "_layer%d");
    Tube lay_tub(x_barrel.inner_r(), x_barrel.outer_r(), x_barrel.z_length() / 2.0);
    Volume lay_vol(lay_nam, lay_tub, air); // Create the layer envelope volume.
    Position lay_pos(0, 0, getAttrOrDefault(x_barrel, _U(z0), 0.));
    lay_vol.setVisAttributes(description.visAttributes(x_layer.visStr()));

    //code that I used previously in lieu of an envelope, for visualization purposes
    //Assembly lay_vol(lay_nam);
    //Position lay_pos(0, 0, getAttrOrDefault(x_barrel, _U(z0), 0.));
    DetElement lay_elt(sdet, lay_nam, lay_id);
    // the local coordinate systems of modules in dd4hep and acts differ
    // see http://acts.web.cern.ch/ACTS/latest/doc/group__DD4hepPlugins.html
    auto& layerParams =
        DD4hepDetectorHelper::ensureExtension<dd4hep::rec::VariantParameters>(lay_elt);
    double phi0     = x_layout.phi0();     // Starting phi of first module.
    double phi_tilt = x_layout.phi_tilt(); // Phi tilt of a module.
    double rc       = x_layout.rc();       // Radius of the module center.
    int nphi        = x_layout.nphi();     // Number of modules in phi.
    double rphi_dr  = x_layout.dr();       // The delta radius of every other module.
    double phi_incr = (M_PI * 2) / nphi;   // Phi increment for one module.
    double phic     = phi0;                // Phi of the module center.
    double z0       = z_layout.z0();       // Z position of first module in phi.
    double nz       = z_layout.nz();       // Number of modules to place in z.
    double z_dr     = z_layout.dr();       // Radial displacement parameter, of every other module.
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
    for(auto& m_nam : module_names) {

      Volume module_env = volumes[m_nam];
      // Loop over the number of modules in phi.
      for (int ii = 0; ii < nphi; ii++) {
        double dx = z_dr * std::cos(phic + phi_tilt); // Delta x of module position.
        double dy = z_dr * std::sin(phic + phi_tilt); // Delta y of module position.
        double x  = rc * std::cos(phic);              // Basic x module position.
        double y  = rc * std::sin(phic);              // Basic y module position.

        // Loop over the number of modules in z. note that for a system working with stave designs nz = 1 always
        for (int j = 0; j < nz; j++) {
          string module_name = _toString(module, "module%d");
          DetElement mod_elt(lay_elt, module_name, module);

          Transform3D tr(RotationZYX(0, ((M_PI / 2) - phic - phi_tilt), -M_PI/2), /*RotationZ(phic)*/
                          Position(x, y, module_z));

          pv = lay_vol.placeVolume(module_env, tr);
          pv.addPhysVolID("module", module);
          mod_elt.setPlacement(pv);
          if(sensitives.count(m_nam)) 
          {
            Placements& sensVols = sensitives[m_nam];
            for (size_t ic = 0; ic < sensVols.size(); ++ic) {
            PlacedVolume sens_pv = sensVols[ic];
            DetElement comp_de(mod_elt, std::string("de_") + sens_pv.volume().name(), module);
            comp_de.setPlacement(sens_pv);

            auto& comp_de_params =
                DD4hepDetectorHelper::ensureExtension<dd4hep::rec::VariantParameters>(comp_de);
            comp_de_params.set<string>("axis_definitions", "XYZ");
            //comp_de.setAttributes(description, sens_pv.volume(), x_layer.regionStr(), x_layer.limitsStr(),
                                  //xml_det_t(xmleles[m_nam]).visStr());
            //

            volSurfaceList(comp_de)->push_back(volplane_surfaces[m_nam][ic]);
            }
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
    }
    // Create the PhysicalVolume for the layer.
    pv = assembly.placeVolume(lay_vol, lay_pos); // Place layer in mother
    pv.addPhysVolID("layer", lay_id);            // Set the layer ID.
    lay_elt.setAttributes(description, lay_vol, x_layer.regionStr(), x_layer.limitsStr(),
                          x_layer.visStr());
    //lay_vol.setVisAttributes(description.invisible());
    lay_elt.setPlacement(pv); 
  }
  
  //finally, place the world
  sdet.setAttributes(description, assembly, x_det.regionStr(), x_det.limitsStr(), x_det.visStr());
  assembly.setVisAttributes(description.invisible());
  pv = description.pickMotherVolume(sdet).placeVolume(assembly);
  pv.addPhysVolID("system", det_id); // Set the subdetector system ID.
  sdet.setPlacement(pv);
  //printout(WARNING, "BarrelTrackerOuter", "DetElement instance \"sdet\" might be corrupted if the GDML design file is too big.");	
  
  return sdet;
}

//@}
// clang-format off
//Macros to access the XML files
DECLARE_DETELEMENT(epic_SiliconBarrelStandardized,    create_BarrelTrackerOuterStandardized)
