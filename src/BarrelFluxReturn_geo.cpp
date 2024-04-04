// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2024 Wouter Deconinck

#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/Printout.h"
#include "XML/Utilities.h"

using namespace dd4hep;

static dd4hep::Ref_t create_detector(dd4hep::Detector& description, xml_h e, [[maybe_unused]] dd4hep::SensitiveDetector sens)
{
  xml_det_t x_det  = e;
  xml_comp_t x_dim = x_det.dimensions();

  // Create element
  dd4hep::DetElement sdet(x_det.nameStr(), x_det.id());

  // Create envelope
  dd4hep::Material air = description.air();
  Tube env_solid(x_dim.rmin(), x_dim.rmax(), x_dim.z() / 2.0);
  Volume env_vol(x_det.nameStr() + "_env", env_solid, air);

  // Create volume map
  std::map<std::string, dd4hep::Volume> volumes_by_name;
  for (xml_coll_t shape(e, _U(shape)); shape; ++shape) {
    xml_comp_t x_shape = shape;
    Solid solid = xml::createShape(description, x_shape.typeStr(), shape);
    Material mat = description.material(x_shape.materialStr());
    volumes_by_name[x_shape.nameStr()] = Volume(x_shape.nameStr(), solid, mat);
  }

  // Replicate volumes
  for (xml_coll_t repl(e, _U(replicate)); repl; ++repl) {
    xml_comp_t x_repl = repl;
    Volume& vol = volumes_by_name[x_repl.attr<std::string>(_U(shape))];
    Transform3D tf = xml::createTransformation(x_repl);
    for (int i = 0; i < x_repl.count(); ++i) {
      double phi = x_repl.phi0() + i * x_repl.attr<double>(_Unicode(dphi));
      env_vol.placeVolume(vol, Transform3D(RotationZ(phi)) * tf);
    }
  }

  // Get position and place volume
  PlacedVolume pv = description.pickMotherVolume(sdet).placeVolume(env_vol);
  sdet.setPlacement(pv);
  return sdet;
}

DECLARE_DETELEMENT(epic_BarrelFluxReturn, create_detector)
