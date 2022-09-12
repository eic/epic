//==========================================================================
//  Implementation of forward calorimeter with the insert shape cut out
//--------------------------------------------------------------------------
//  Author: Ryan Milton (UCR)
//==========================================================================

#include "DD4hep/DetFactoryHelper.h"
#include <XML/Helper.h>
#include <XML/Layering.h>

using namespace dd4hep;

static Ref_t createDetector(Detector& desc, xml_h handle, SensitiveDetector sens)
{
  xml_det_t   detElem = handle;
  std::string detName = detElem.nameStr();
  int         detID   = detElem.id();

  xml_dim_t dim    = detElem.dimensions();
  double    rmin   = dim.rmin(); // Dummy variable. Set to 0 since cutting out insert
  double    rmax   = dim.rmax(); // Max radius of endcap
  double    length = dim.z();    // Size along z-axis

  xml_dim_t pos = detElem.position();

  Material air = desc.material("Air");

  // Getting insert dimensions
  const xml::Component& insert_xml       = detElem.child(_Unicode(insert));
  xml_dim_t             insert_dim       = insert_xml.dimensions();
  xml_dim_t             insert_local_pos = insert_xml.position();

  // Defining envelope
  Tube envelope(rmin, rmax, length / 2.0);

  // Removing insert shape from envelope
  Box              insert(insert_dim.x() / 2., insert_dim.y() / 2., length / 2.);
  SubtractionSolid envelope_with_inserthole(envelope, insert, Position(insert_local_pos.x(), insert_local_pos.y(), 0.));
  Volume           envelopeVol(detName, envelope_with_inserthole, air);
  // Setting envelope attributes
  envelopeVol.setAttributes(desc, detElem.regionStr(), detElem.limitsStr(), detElem.visStr());
  PlacedVolume pv;

  int    layer_num = 1;
  double layer_z   = -length / 2.; // Keeps track of layers' local z locations
  // Looping through all the different layer sections (Steel/Sc, W/Sc)
  for (xml_coll_t c(detElem, _U(layer)); c; ++c) {
    xml_comp_t x_layer         = c;
    int        repeat          = x_layer.repeat();
    double     layer_thickness = x_layer.thickness();

    // Looping through the number of repeated layers in each section
    for (int i = 0; i < repeat; i++) {
      layer_z += layer_thickness / 2.; // Going to halfway point in layer

      std::string layer_name = detName + _toString(layer_num, "_layer%d");
      Tube        layer(rmin, rmax, layer_thickness / 2.);

      // Removing insert shape from each layer
      Box              layer_insert(insert_dim.x() / 2., insert_dim.y() / 2., layer_thickness / 2.);
      SubtractionSolid layer_with_inserthole(layer, layer_insert,
                                             Position(insert_local_pos.x(), insert_local_pos.y(), 0.));
      Volume           layer_vol(layer_name, layer_with_inserthole, air);

      int    slice_num = 1;
      double slice_z   = -layer_thickness / 2.; // Keeps track of slices' z locations in each layer

      // Looping over each layer's slices
      for (xml_coll_t l(x_layer, _U(slice)); l; ++l) {
        xml_comp_t  x_slice         = l;
        double      slice_thickness = x_slice.thickness();
        std::string slice_name      = layer_name + _toString(slice_num, "slice%d");
        Material    slice_mat       = desc.material(x_slice.materialStr());
        slice_z += slice_thickness / 2.; // Going to slice halfway point

        Tube slice(rmin, rmax, slice_thickness / 2.);

        // Removing insert shape from each slice
        Box              slice_insert(insert_dim.x() / 2., insert_dim.y() / 2., slice_thickness / 2.0);
        SubtractionSolid slice_with_inserthole(slice, slice_insert,
                                               Position(insert_local_pos.x(), insert_local_pos.y(), 0.));
        Volume           slice_vol(slice_name, slice_with_inserthole, slice_mat);

        // Setting appropriate slices as sensitive
        if (x_slice.isSensitive()) {
          sens.setType("calorimeter");
          slice_vol.setSensitiveDetector(sens);
        }

        // Setting slice attributes
        slice_vol.setAttributes(desc, x_slice.regionStr(), x_slice.limitsStr(), x_slice.visStr());

        // Placing slice within layer
        pv = layer_vol.placeVolume(slice_vol, Transform3D(RotationZYX(0, 0, 0), Position(0., 0., slice_z)));
        pv.addPhysVolID("slice", slice_num);
        slice_z += slice_thickness / 2.; // Going to end of slice
        ++slice_num;
      }

      // Setting layer attributes
      layer_vol.setAttributes(desc, x_layer.regionStr(), x_layer.limitsStr(), x_layer.visStr());

      // Placing each layer inside the envelope volume
      // layer_z starts as -length/2 (front of endcap)
      // layer_z is increased by layer_thickness/2 at the start of the layer loop
      // This moves to the middle of the layer (location of placement)
      // It's then increased by layer_thickness/2 again to move to end of layer
      pv = envelopeVol.placeVolume(layer_vol, Transform3D(RotationZYX(0, 0, 0), Position(0., 0., layer_z)));
      pv.addPhysVolID("layer", layer_num);
      layer_num++;
      layer_z += layer_thickness / 2.;
    }
  }

  DetElement det(detName, detID);
  Volume     motherVol = desc.pickMotherVolume(det);

  // Placing endcap in world volume
  auto tr = Transform3D(Position(pos.x(), pos.y(), pos.z() + length / 2.));

  PlacedVolume phv = motherVol.placeVolume(envelopeVol, tr);
  phv.addPhysVolID("system", detID);
  det.setPlacement(phv);

  return det;
}
DECLARE_DETELEMENT(epic_EndcapCalorimeterWithInsert, createDetector)
