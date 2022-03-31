/** \addtogroup PID
 * \brief Type: **Barrel bar detector (think DIRC) with longitudinal frame**.
 * \author S. Joosten
 *
 * \ingroup PID
 *
 *
 * \code
 * \endcode
 *
 * @{
 */
#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/Printout.h"
#include "DD4hep/Shapes.h"
#include "DDRec/DetectorData.h"
#include "DDRec/Surface.h"
#include "XML/Layering.h"

using namespace std;
using namespace dd4hep;

/** Barrel Bar detector with optional framek
 *
 * - Optional "frame" tag within the module element (frame eats part of the module width
 *   and is longitudinal at both sides)
 * - Detector is setup as a "tracker" so we can use the hits to perform fast MC
 *   for e.g. the DIRC
 *
 */
static Ref_t create_BarrelBarDetectorWithSideFrame(Detector& description, xml_h e, SensitiveDetector sens)
{
  typedef vector<PlacedVolume>            Placements;
  xml_det_t                               x_det    = e;
  Material                                air      = description.air();
  int                                     det_id   = x_det.id();
  string                                  det_name = x_det.nameStr();
  DetElement                              sdet(det_name, det_id);
  map<string, Volume>                     volumes;
  map<string, Placements>                 sensitives;
  map<string, std::vector<rec::VolPlane>> volplane_surfaces;
  PlacedVolume                            pv;
  dd4hep::xml::Dimension                  dimensions(x_det.dimensions());
  xml_dim_t                               dirc_pos = x_det.position();

  map<string, std::array<double, 2>> module_thicknesses;

  Tube   topVolumeShape(dimensions.rmin(), dimensions.rmax(), dimensions.length() * 0.5);
  Volume assembly(det_name, topVolumeShape, air);

  sens.setType("tracker");

  // loop over the modules
  for (xml_coll_t mi(x_det, _U(module)); mi; ++mi) {
    xml_comp_t x_mod = mi;
    string     m_nam = x_mod.nameStr();

    if (volumes.find(m_nam) != volumes.end()) {
      printout(ERROR, "BarrelBarDetectorWithSideFrame",
               string((string("Module with named ") + m_nam + string(" already exists."))).c_str());
      throw runtime_error("Logics error in building modules.");
    }

    int    ncomponents     = 0;
    int    sensor_number   = 1;
    double total_thickness = 0;

    // Compute module total thickness from components
    xml_coll_t ci(x_mod, _U(module_component));
    for (ci.reset(), total_thickness = 0.0; ci; ++ci) {
      total_thickness += xml_comp_t(ci).thickness();
    }
    // the module assembly volume
    Assembly m_vol(m_nam);
    volumes[m_nam] = m_vol;
    m_vol.setVisAttributes(description, x_mod.visStr());

    // Optional module frame.
    // frame is 2 longitudinal bars at the side of the module
    //
    // |___| <-- example module cross section (x-y plane), frame is flush with the
    //           bottom of the module and protrudes on the top if needed
    //

    // Get frame width, as it impacts the main module for being built. We
    // construct the actual frame structure later (once we know the module width)
    double frame_width = 0;
    if (x_mod.hasChild("frame")) {
      xml_comp_t m_frame = x_mod.child(_U(frame));
      frame_width        = m_frame.width();
    }

    double thickness_so_far    = 0.0;
    double thickness_sum       = -total_thickness / 2.0;
    double max_component_width = 0;
    for (xml_coll_t ci(x_mod, _U(module_component)); ci; ++ci, ++ncomponents) {
      xml_comp_t x_comp = ci;
      xml_comp_t x_pos  = x_comp.position(false);
      string     c_nam  = _toString(ncomponents, "component%d");

      double box_width    = x_comp.width() - 2 * frame_width;
      max_component_width = fmax(max_component_width, box_width);

      Box    c_box{box_width / 2, x_comp.length() / 2, x_comp.thickness() / 2};
      Volume c_vol{c_nam, c_box, description.material(x_comp.materialStr())};

      pv = m_vol.placeVolume(c_vol, Position(0, 0, thickness_sum + x_comp.thickness() / 2.0));

      c_vol.setRegion(description, x_comp.regionStr());
      c_vol.setLimitSet(description, x_comp.limitsStr());
      c_vol.setVisAttributes(description, x_comp.visStr());
      if (x_comp.isSensitive()) {
        c_vol.setSensitiveDetector(sens);
        sensitives[m_nam].push_back(pv);
        module_thicknesses[m_nam] = {thickness_so_far + x_comp.thickness() / 2.0,
                                     total_thickness - thickness_so_far - x_comp.thickness() / 2.0};
        // -------- create a measurement plane for the tracking surface attched to the sensitive volume -----
        rec::Vector3D u(0., 1., 0.);
        rec::Vector3D v(0., 0., 1.);
        rec::Vector3D n(1., 0., 0.);

        // compute the inner and outer thicknesses that need to be assigned to the tracking surface
        // depending on wether the support is above or below the sensor
        double inner_thickness = module_thicknesses[m_nam][0];
        double outer_thickness = module_thicknesses[m_nam][1];

        rec::SurfaceType type(rec::SurfaceType::Sensitive);

        rec::VolPlane surf(c_vol, type, inner_thickness, outer_thickness, u, v, n); //,o ) ;
        volplane_surfaces[m_nam].push_back(surf);
      }
      thickness_sum += x_comp.thickness();
      thickness_so_far += x_comp.thickness();
    }
    // Now add-on the frame
    if (x_mod.hasChild("frame")) {
      xml_comp_t m_frame         = x_mod.child(_U(frame));
      double     frame_thickness = getAttrOrDefault<double>(m_frame, _U(thickness), total_thickness);

      Box lframe_box{m_frame.width() / 2., m_frame.length() / 2., frame_thickness / 2.};
      Box rframe_box{m_frame.width() / 2., m_frame.length() / 2., frame_thickness / 2.};

      // Keep track of frame with so we can adjust the module bars appropriately

      Volume lframe_vol{"left_frame", lframe_box, description.material(m_frame.materialStr())};
      Volume rframe_vol{"right_frame", rframe_box, description.material(m_frame.materialStr())};
      lframe_vol.setVisAttributes(description, m_frame.visStr());
      rframe_vol.setVisAttributes(description, m_frame.visStr());

      m_vol.placeVolume(lframe_vol, Position(frame_width / 2. + max_component_width / 2, 0.,
                                             frame_thickness / 2. - total_thickness / 2.0));
      m_vol.placeVolume(rframe_vol, Position(-frame_width / 2. - max_component_width / 2, 0.,
                                             frame_thickness / 2. - total_thickness / 2.0));
    }
  }

  // now build the layers
  for (xml_coll_t li(x_det, _U(layer)); li; ++li) {
    xml_comp_t x_layer  = li;
    xml_comp_t x_barrel = x_layer.child(_U(barrel_envelope));
    xml_comp_t x_layout = x_layer.child(_U(rphi_layout));
    xml_comp_t z_layout = x_layer.child(_U(z_layout)); // Get the <z_layout> element.
    int        lay_id   = x_layer.id();
    string     m_nam    = x_layer.moduleStr();
    string     lay_nam  = _toString(x_layer.id(), "layer%d");
    Tube       lay_tub(x_barrel.inner_r(), x_barrel.outer_r(), x_barrel.z_length() / 2.0);
    Volume     lay_vol(lay_nam, lay_tub, air); // Create the layer envelope volume.
    lay_vol.setVisAttributes(description, x_layer.visStr());

    double phi0     = x_layout.phi0();     // Starting phi of first module.
    double phi_tilt = x_layout.phi_tilt(); // Phi tilt of a module.
    double rc       = x_layout.rc();       // Radius of the module center.
    int    nphi     = x_layout.nphi();     // Number of modules in phi.
    double rphi_dr  = x_layout.dr();       // The delta radius of every other module.
    double phi_incr = (M_PI * 2) / nphi;   // Phi increment for one module.
    double phic     = phi0;                // Phi of the module center.
    double z0       = z_layout.z0();       // Z position of first module in phi.
    double nz       = z_layout.nz();       // Number of modules to place in z.
    double z_dr     = z_layout.dr();       // Radial displacement parameter, of every other module.

    Volume      module_env = volumes[m_nam];
    DetElement  lay_elt(sdet, _toString(x_layer.id(), "layer%d"), lay_id);
    Placements& sensVols = sensitives[m_nam];

    // Z increment for module placement along Z axis.
    // Adjust for z0 at center of module rather than
    // the end of cylindrical envelope.
    double z_incr = nz > 1 ? (2.0 * z0) / (nz - 1) : 0.0;
    // Starting z for module placement along Z axis.
    double module_z = -z0;
    int    module   = 1;

    // Loop over the number of modules in phi.
    for (int ii = 0; ii < nphi; ii++) {
      double dx = z_dr * std::cos(phic + phi_tilt); // Delta x of module position.
      double dy = z_dr * std::sin(phic + phi_tilt); // Delta y of module position.
      double x  = rc * std::cos(phic);              // Basic x module position.
      double y  = rc * std::sin(phic);              // Basic y module position.

      // Loop over the number of modules in z.
      for (int j = 0; j < nz; j++) {
        string     module_name = _toString(module, "module%d");
        DetElement mod_elt(lay_elt, module_name, module);

        Transform3D tr(RotationZYX(0, ((M_PI / 2) - phic - phi_tilt), -M_PI / 2), Position(x, y, module_z));

        pv = lay_vol.placeVolume(module_env, tr);
        pv.addPhysVolID("module", module);
        pv.addPhysVolID("section", j);
        mod_elt.setPlacement(pv);
        for (size_t ic = 0; ic < sensVols.size(); ++ic) {
          PlacedVolume sens_pv = sensVols[ic];
          DetElement   comp_de(mod_elt, std::string("de_") + sens_pv.volume().name(), module);
          comp_de.setPlacement(sens_pv);
          rec::volSurfaceList(comp_de)->push_back(volplane_surfaces[m_nam][ic]);
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
    pv = assembly.placeVolume(lay_vol); // Place layer in mother
    pv.addPhysVolID("layer", lay_id);   // Set the layer ID.
    lay_elt.setAttributes(description, lay_vol, x_layer.regionStr(), x_layer.limitsStr(), x_layer.visStr());
    lay_elt.setPlacement(pv);
  }
  sdet.setAttributes(description, assembly, x_det.regionStr(), x_det.limitsStr(), x_det.visStr());
  assembly.setVisAttributes(description.invisible());
  pv = description.pickMotherVolume(sdet).placeVolume(assembly, Position(0, 0, dirc_pos.z()));
  pv.addPhysVolID("system", det_id); // Set the subdetector system ID.
  sdet.setPlacement(pv);
  return sdet;
}

//@}
// clang-format off
DECLARE_DETELEMENT(BarrelBarDetectorWithSideFrame, create_BarrelBarDetectorWithSideFrame)
DECLARE_DETELEMENT(athena_FakeDIRC, create_BarrelBarDetectorWithSideFrame)
