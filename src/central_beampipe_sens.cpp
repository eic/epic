// SPDX-License-Identifier: LGPL-3.0-or-later

#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/Printout.h"
#include "DD4hep/Shapes.h"
#include "DDRec/DetectorData.h"
#include "DDRec/Surface.h"
#include "XML/Layering.h"
#include <XML/Helper.h>

using namespace std;
using namespace dd4hep;

static Ref_t create_detector(Detector& description, xml_h e, SensitiveDetector sens) {
  sens.setType("tracker");

  xml_det_t x_det       = e;
  string det_name       = x_det.nameStr();
  string mat_name = dd4hep::getAttrOrDefault<string>(x_det, _U(material), "Vacuum");

  // Helper function to create polycone pairs (wall, coating, and vacuum)
  auto zplane_to_polycones = [](xml::Component& x_pipe) {
    std::vector<double> zero, z;
    std::vector<double> rmax_wall, rmax_coating, rmax_vacuum;
    std::vector<double> rmin_wall, rmin_coating, rmin_vacuum;
    // thickness
    auto wall_thickness    = getAttrOrDefault(x_pipe, _Unicode(wall_thickness), 1 * mm);
    auto coating_thickness = getAttrOrDefault(x_pipe, _Unicode(coating_thickness), 30 * um);

    for (xml_coll_t x_zplane_i(x_pipe, _Unicode(zplane)); x_zplane_i; ++x_zplane_i) {
      xml_comp_t x_zplane = x_zplane_i;
      // z position
      z.push_back(x_zplane.attr<double>(_Unicode(z)));
      // outer radius
      rmax_wall.push_back(x_zplane.attr<double>(_Unicode(ID)) / 2.0 + wall_thickness);
      rmax_coating.push_back(x_zplane.attr<double>(_Unicode(ID)) / 2.0);
      rmax_vacuum.push_back(x_zplane.attr<double>(_Unicode(ID)) / 2.0 - coating_thickness);
      // inner radius
      rmin_wall.push_back(x_zplane.attr<double>(_Unicode(ID)) / 2.0);
      rmin_coating.push_back(x_zplane.attr<double>(_Unicode(ID)) / 2.0 - coating_thickness);
      rmin_vacuum.push_back(0);
    }
    return std::tuple<Polycone, Polycone, Polycone>({0, 2.0 * M_PI, rmin_wall, rmax_wall, z},
                                                    {0, 2.0 * M_PI, rmin_coating, rmax_coating, z},
                                                    {0, 2.0 * M_PI, rmin_vacuum, rmax_vacuum, z});
  };

  // Sensitive pipe
  xml::Component sens_pipe_c = x_det.child(_Unicode(sens_pipe));
  auto sens_pipe_polycones = zplane_to_polycones(sens_pipe_c);
  Volume vol(det_name + "_vol",std::get<0>(sens_pipe_polycones), description.material(mat_name));
  vol.setVisAttributes(description.visAttributes(x_det.visStr()));
  vol.setSensitiveDetector(sens);

  DetElement det(det_name, x_det.id());
  Volume motherVol = description.pickMotherVolume(det);
  PlacedVolume phv = motherVol.placeVolume(vol, 0);
  phv.addPhysVolID("system", x_det.id());

  det.setPlacement(phv);

  return det;
}

DECLARE_DETELEMENT(BeamPipeSens, create_detector)
