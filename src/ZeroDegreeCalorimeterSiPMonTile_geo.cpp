// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2023 Sebouh J. Paul

//==========================================================================
//  Implementation of SiPM-on-tile Zero-Degree Calorimeter
//--------------------------------------------------------------------------
//  Author: Sebouh J. Paul (UCR)
//==========================================================================

#include "DD4hep/DetFactoryHelper.h"
#include <XML/Helper.h>
#include <XML/Layering.h>

using namespace dd4hep;

static Ref_t createDetector(Detector& desc, xml_h handle, SensitiveDetector sens)
{
  xml_det_t     detElem = handle;
  std::string   detName = detElem.nameStr();
  int           detID   = detElem.id();

  xml_dim_t  dim        = detElem.dimensions();
  double     width      = dim.x(); // Size along x-axis
  double     height     = dim.y(); // Size along y-axis
  double     length     = dim.z(); // Size along z-axis

  xml_dim_t  pos        = detElem.position(); // Position in global coordinates
  xml_dim_t rot         = detElem.rotation();

  Material   air        = desc.material("Air");

  // Defining envelope
  Box envelope(width / 2.0, height / 2.0, length / 2.0);

  // Defining envelope volume
  Volume envelopeVol(detName, envelope, air);
  // Setting envelope attributes
  envelopeVol.setAttributes(
    desc,
    detElem.regionStr(),
    detElem.limitsStr(),
    detElem.visStr()
  );

  PlacedVolume pv;

  double z_distance_traversed = 0.;

  int layer_num = 1;

  // Looping through all the different layer sections
  for(xml_coll_t c(detElem,_U(layer)); c; ++c)
  {
    xml_comp_t x_layer = c;
    int repeat = x_layer.repeat();
    double layer_thickness = x_layer.thickness();

    // Looping through the number of repeated layers in each section
    for(int i = 0; i < repeat; i++)
    {
      std::string layer_name = detName + _toString(layer_num, "_layer%d");

      Box layer(width / 2., height / 2., layer_thickness / 2.);

      Volume layer_vol(layer_name, layer, air);

      int slice_num = 1;
      double slice_z = -layer_thickness / 2.; // Keeps track of slices' z locations in each layer

      // Looping over each layer's slices
      for(xml_coll_t l(x_layer,_U(slice)); l; ++l)
      {
        xml_comp_t x_slice = l;
        double slice_thickness = x_slice.thickness();
        std::string slice_name = layer_name + _toString(slice_num, "slice%d");
        Material slice_mat = desc.material(x_slice.materialStr());
        slice_z += slice_thickness/2.; // Going to slice halfway point

        Box slice(width/2., height/2., slice_thickness/2.);

        Volume slice_vol (slice_name, slice, slice_mat);

        // Setting appropriate slices as sensitive
        if(x_slice.isSensitive())
        {
          sens.setType("calorimeter");
          slice_vol.setSensitiveDetector(sens);
        }

        // Setting slice attributes
        slice_vol.setAttributes(
          desc,
          x_slice.regionStr(),
          x_slice.limitsStr(),
          x_slice.visStr()
        );

        // Placing slice within layer
        pv = layer_vol.placeVolume(
          slice_vol,
          Transform3D(
            RotationZYX(0, 0, 0),
            Position(
              0.,
              0.,
              slice_z
            )
          )
        );
        pv.addPhysVolID("slice", slice_num);
        slice_z += slice_thickness/2.;
        z_distance_traversed += slice_thickness;
        ++slice_num;
      }

      // Setting layer attributes
      layer_vol.setAttributes(
        desc,
        x_layer.regionStr(),
        x_layer.limitsStr(),
        x_layer.visStr()
      );
      // Placing each layer inside the envelope volume
      // -length/2. is front of detector in global coordinate system
      // + (z_distance_traversed - layer_thickness) goes to the front of each layer
      // + layer_thickness/2. places layer in correct spot
      // Example: After placement of slices in first layer, z_distance_traversed = layer_thickness
      //          Subtracting layer_thickness goes back to the front of the first slice (Now, z = -length/2)
      //          Adding layer_thickness/2. goes to half the first layer thickness (proper place to put layer)
      //          Each loop over repeat will increases z_distance_traversed by layer_thickness
      pv = envelopeVol.placeVolume(
        layer_vol,
        Transform3D(
          RotationZYX(0, 0, 0),
          Position(
            0.,
            0.,
            -length/2.  + (z_distance_traversed - layer_thickness) + layer_thickness/2.
          )
        )
      );
      pv.addPhysVolID("layer", layer_num);
      layer_num++;
    }
  }

  DetElement   det(detName, detID);
  Volume motherVol = desc.pickMotherVolume(det);

  // Placing ZDC in world volume
  auto tr = Transform3D(RotationZYX(rot.z(), rot.y(), rot.x()),Position(pos.x(), pos.y(), pos.z()));
  PlacedVolume phv = motherVol.placeVolume(envelopeVol, tr);
  phv.addPhysVolID("system", detID);
  det.setPlacement(phv);

  return det;
}
DECLARE_DETELEMENT(ZeroDegreeCalorimeterSiPMonTile, createDetector)
