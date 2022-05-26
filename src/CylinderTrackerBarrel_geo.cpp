//==========================================================================
// Specialized generic detector constructor
//==========================================================================

#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/Printout.h"

#if defined(USE_ACTSDD4HEP)
#include "ActsDD4hep/ActsExtension.hpp"
#include "ActsDD4hep/ConvertMaterial.hpp"
#else
#include "Acts/Plugins/DD4hep/ActsExtension.hpp"
#include "Acts/Plugins/DD4hep/ConvertDD4hepMaterial.hpp"
#endif

using namespace std;
using namespace dd4hep;
using namespace dd4hep::detail;

/** A barrel tracker with a module that is curved (not flat).
 *
 * \ingroup tracking
 */
static Ref_t CylinderTrackerBarrel_create_detector(Detector& description, xml_h e, SensitiveDetector sens)
{
  typedef vector<PlacedVolume> Placements;
  xml_det_t                    x_det = e;
  Material                     air   = description.air();

  int        det_id   = x_det.id();
  string     det_name = x_det.nameStr();
  DetElement sdet(det_name, det_id);

  //Acts::ActsExtension* barrelExtension = new Acts::ActsExtension();
  //barrelExtension->addType("barrel", "detector");
  //sdet.addExtension<Acts::ActsExtension>(barrelExtension);

  Assembly                assembly(det_name);
  map<string, Volume>     mod_volumes;
  map<string, Placements> sensitives;
  PlacedVolume            pv;

  sens.setType("tracker");
  int n_modules = 0;
  for (xml_coll_t mi(x_det, _U(module)); mi; ++mi) {
    n_modules++;
    xml_comp_t x_mod = mi;
    xml_comp_t m_env = x_mod.child(_U(module_envelope));
    string     m_nam = x_mod.nameStr();

    Assembly module_assembly(_toString(n_modules, "mod_assembly_%d"));
    auto     module_rmin      = m_env.rmin();
    auto     module_thickness = m_env.thickness();
    auto     module_length    = m_env.length();
    auto     module_phi       = getAttrOrDefault(m_env, _Unicode(phi), 90.0);

    Volume m_vol(
        m_nam,
        Tube(module_rmin, module_rmin + module_thickness, module_length / 2, -module_phi / 2.0, module_phi / 2.0), air);
    int    ncomponents = 0, sensor_number = 1;
    module_assembly.placeVolume(m_vol, Position(-module_rmin, 0, 0));
    mod_volumes[m_nam] = module_assembly;
    m_vol.setVisAttributes(description.visAttributes(x_mod.visStr()));

    auto comp_rmin = module_rmin;
    for (xml_coll_t ci(x_mod, _U(module_component)); ci; ++ci, ++ncomponents) {
      xml_comp_t x_comp = ci;
      xml_comp_t x_pos  = x_comp.position(false);
      xml_comp_t x_rot  = x_comp.rotation(false);
      string     c_nam  = _toString(ncomponents, "component%d");

      auto comp_thickness = x_comp.thickness();
      comp_rmin           = getAttrOrDefault(x_comp, _Unicode(rmin), comp_rmin);
      auto comp_phi       = getAttrOrDefault(x_comp, _Unicode(phi), module_phi);
      auto comp_phi0      = getAttrOrDefault(x_comp, _Unicode(phi0), 0.0);
      auto comp_length    = getAttrOrDefault(x_comp, _Unicode(length), module_length);

      Tube         c_tube(comp_rmin, comp_rmin + comp_thickness, comp_length / 2, -comp_phi / 2.0 + comp_phi0,
                  comp_phi / 2.0 + comp_phi0);
      Volume       c_vol(c_nam, c_tube, description.material(x_comp.materialStr()));
      PlacedVolume c_pv;

      if (x_pos && x_rot) {
        Position    c_pos(x_pos.x(0), x_pos.y(0), x_pos.z(0));
        RotationZYX c_rot(x_rot.z(0), x_rot.y(0), x_rot.x(0));
        c_pv = m_vol.placeVolume(c_vol, Transform3D(c_rot, c_pos));
      } else if (x_rot) {
        c_pv = m_vol.placeVolume(c_vol, RotationZYX(x_rot.z(0), x_rot.y(0), x_rot.x(0)));
      } else if (x_pos) {
        c_pv = m_vol.placeVolume(c_vol, Position(x_pos.x(0), x_pos.y(0), x_pos.z(0)));
      } else {
        c_pv = m_vol.placeVolume(c_vol);
      }
      c_vol.setRegion(description, x_comp.regionStr());
      c_vol.setLimitSet(description, x_comp.limitsStr());
      c_vol.setVisAttributes(description, x_comp.visStr());
      if (x_comp.isSensitive()) {
        c_pv.addPhysVolID(_U(sensor), sensor_number++);
        c_vol.setSensitiveDetector(sens);
        sensitives[m_nam].push_back(c_pv);
      }
      comp_rmin = comp_rmin + comp_thickness;
    }
  }
  for (xml_coll_t li(x_det, _U(layer)); li; ++li) {
    xml_comp_t x_layer  = li;
    xml_comp_t x_barrel = x_layer.child(_U(barrel_envelope));
    xml_comp_t x_layout = x_layer.child(_U(rphi_layout));
    xml_comp_t z_layout = x_layer.child(_U(z_layout)); // Get the <z_layout> element.
    int        lay_id   = x_layer.id();
    string     m_nam    = x_layer.moduleStr();
    string     lay_nam  = _toString(x_layer.id(), "layer%d");
    Tube       lay_tub(x_barrel.inner_r(), x_barrel.outer_r(), x_barrel.z_length() / 2);
    Volume     lay_vol(lay_nam, lay_tub, air); // Create the layer envelope volume.
    lay_vol.setVisAttributes(description.visAttributes(x_layer.visStr()));
    double     phi0     = x_layout.phi0();     // Starting phi of first module.
    double     phi_tilt = x_layout.phi_tilt(); // Phi tilt of a module.
    double     rc       = x_layout.rc();       // Radius of the module center.
    int        nphi     = x_layout.nphi();     // Number of modules in phi.
    double     rphi_dr  = x_layout.dr();       // The delta radius of every other module.
    double     phi_incr = (M_PI * 2) / nphi;   // Phi increment for one module.
    double     phic     = phi0;                // Phi of the module center.
    double     z0       = z_layout.z0();       // Z position of first module in phi.
    double     nz       = z_layout.nz();       // Number of modules to place in z.
    double     z_dr     = z_layout.dr();       // Radial displacement parameter, of every other module.
    Volume     m_env    = mod_volumes[m_nam];
    DetElement lay_elt(sdet, _toString(x_layer.id(), "layer%d"), lay_id);

    ///Acts::ActsExtension* layerExtension = new Acts::ActsExtension();
    ///layerExtension->addType("sensitive cylinder", "layer");
    /////// layerExtension->addValue(10. * Acts::UnitConstants::mm, "r", "envelope");
    ///lay_elt.addExtension<Acts::ActsExtension>(layerExtension);

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
        // Module PhysicalVolume.
        //         Transform3D
        //         tr(RotationZYX(0,-((M_PI/2)-phic-phi_tilt),M_PI/2),Position(x,y,module_z));
        // NOTE (Nikiforos, 26/08 Rotations needed to be fixed so that component1 (silicon) is on the
        // outside
        Transform3D tr(RotationZYX(phic - phi_tilt, 0, 0), Position(x, y, module_z));

        pv = lay_vol.placeVolume(m_env, tr);
        pv.addPhysVolID("module", module);
        mod_elt.setPlacement(pv);
        for (size_t ic = 0; ic < sensVols.size(); ++ic) {
          PlacedVolume sens_pv = sensVols[ic];
          DetElement   comp_elt(mod_elt, sens_pv.volume().name(), module);
          comp_elt.setPlacement(sens_pv);
          //Acts::ActsExtension* moduleExtension = new Acts::ActsExtension();
          //comp_elt.addExtension<Acts::ActsExtension>(moduleExtension);
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
  // assembly.setVisAttributes(description.invisible());
  pv = description.pickMotherVolume(sdet).placeVolume(assembly);
  pv.addPhysVolID("system", det_id); // Set the subdetector system ID.
  pv.addPhysVolID("barrel", 1);      // Flag this as a barrel subdetector.
  sdet.setPlacement(pv);
  return sdet;
}

// clang-format off
DECLARE_DETELEMENT(athena_CylinderTrackerBarrel, CylinderTrackerBarrel_create_detector)
DECLARE_DETELEMENT(athena_MMTrackerBarrel,       CylinderTrackerBarrel_create_detector)
DECLARE_DETELEMENT(athena_RWellTrackerBarrel,    CylinderTrackerBarrel_create_detector)
DECLARE_DETELEMENT(athena_CylinderVertexBarrel,  CylinderTrackerBarrel_create_detector)

