// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2023 Friederike Bock

//==========================================================================
//  Implementation of longitudinally separated forward calorimeter
//--------------------------------------------------------------------------
//  Author: Friederike Bock (ORNL)
//==========================================================================

#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/Printout.h"
#include "DD4hep/Shapes.h"
#include "DDRec/DetectorData.h"
#include "DDRec/Surface.h"
#include <XML/Helper.h>
#include "XML/Layering.h"
#include "XML/Utilities.h"
using namespace dd4hep;

struct moduleParamsStrct{
  moduleParamsStrct(): mod_BIwidth(0.), mod_BIheight(0.), mod_SWThick(0.),   mod_TWThick(0.), mod_FWThick (0.),
                      mod_BWThick(0.), mod_width(0.), mod_height(0.),
                      mod_notchDepth(0.), mod_notchHeight(0.), mod_foilThick(0.),  mod_pcbLength(0.), mod_pcbThick(0.), mod_pcbWidth(0.), mod_visStr(""), mod_regStr(""), mod_limStr("")
                      {}
  moduleParamsStrct(   double BIwidth, double BIheight, double SWThick, double TWThick, double FWThick, double BWThick, double width, double height,
                       double notchDepth, double notchHeight, double foilThick,
                       double pcbLegth, double pcbThick, double pcbWidth,
                       std::string visStr, std::string regStr, std::string limStr){
      mod_BIwidth       = BIwidth;
      mod_BIheight      = BIheight;
      mod_SWThick       = SWThick;
      mod_TWThick       = TWThick;
      mod_FWThick       = FWThick;
      mod_BWThick       = BWThick;
      mod_width         = width;
      mod_height        = height;
      mod_notchDepth      = notchDepth;
      mod_notchHeight     = notchHeight;
      mod_foilThick       = foilThick;
      mod_pcbLength       = pcbLegth;
      mod_pcbThick        = pcbThick;
      mod_pcbWidth        = pcbWidth;
      mod_visStr          = visStr;
      mod_regStr          = regStr;
      mod_limStr          = limStr;
  }
  double      mod_BIwidth   = 0.;
  double      mod_BIheight  = 0.;
  double      mod_SWThick   = 0.;
  double      mod_TWThick   = 0.;
  double      mod_FWThick   = 0.;
  double      mod_BWThick   = 0.;
  double      mod_width     = 0.;
  double      mod_height    = 0.;
  double      mod_notchDepth  = 0.;
  double      mod_notchHeight  = 0.;
  double      mod_foilThick    = 0.;
  double      mod_pcbLength    = 0.;
  double      mod_pcbThick     = 0.;
  double      mod_pcbWidth     = 0.;
  std::string mod_visStr      = "";
  std::string mod_regStr      = "";
  std::string mod_limStr      = "";
} ;

struct sliceParamsStrct{
  sliceParamsStrct(): layer_ID(0), slice_ID(0), slice_partID(0), slice_thick(0.), slice_offset(0.), slice_readoutLayer(0), slice_matStr(""), slice_visStr(""), slice_regStr(""), slice_limStr("")
                      {}
  sliceParamsStrct(
                      int l_ID, int sl_ID, int sl_partID, double sl_thick, double sl_off, int l_rl, std::string sl_matStr, std::string sl_visStr, std::string sl_regStr, std::string sl_limStr ){
      layer_ID             = l_ID;
      slice_ID            = sl_ID;
      slice_partID        = sl_partID;
      slice_thick         = sl_thick;
      slice_offset        = sl_off;
      slice_readoutLayer  = l_rl;
      slice_matStr        = sl_matStr;
      slice_visStr        = sl_visStr;
      slice_regStr        = sl_regStr;
      slice_limStr        = sl_limStr;

  }
  int         layer_ID      = 0;
  int         slice_ID      = 0;
  int         slice_partID  = 0;
  double      slice_thick   = 0.;
  double      slice_offset  = 0.;
  int         slice_readoutLayer  = 0;
  std::string slice_matStr  = "";
  std::string slice_visStr  = "";
  std::string slice_regStr  = "";
  std::string slice_limStr  = "";
};

//************************************************************************************************************
//************************** Assembly for absorber plates  ***************************************************
//************************************************************************************************************
Volume createAbsorberPlate(Detector& desc,
                                   std::string basename,
                                   double h_mod,
                                   double w_mod,
                                   double t_mod_tp,
                                   double t_mod_sp,
                                   double t_slice,
                                   double w_notch,
                                   double h_notch,
                                   Material slice_mat,
                                   std::string region,
                                   std::string limit,
                                   std::string vis,
                                   bool renderComp
){

  double w_plate  = (w_mod/2-t_mod_sp)*2;
  double l_A = -w_plate/2;
  double l_B = -(w_plate/2-w_notch);
  double r_A = w_plate/2;
                                        // 0      1     2     3     4
  const std::vector<double> xCoord      = { l_A,  r_A,  r_A,  l_A,   l_A,
                                        // 5     6     7
                                           l_B,   l_B, l_A
                                          };
                                        // 0      1     2     3      4

  double topA = h_mod/2-t_mod_tp;
  double topB = h_notch/2;
  double botA = -(h_mod/2-t_mod_tp);
  double botB = -(h_notch/2);
                                          // 0     1       2      3       4
  const std::vector<double> yCoord      = { topA,  topA,   botA,  botA,   botB,
                                        // 5       6       7      8       9
                                          botB,   topB,   topB
                                        };

  const std::vector<double> zStep       = {-t_slice/2, t_slice/2};
  const std::vector<double> zStepX      = {0., 0.};
  const std::vector<double> zStepY      = {0., 0.};
  const std::vector<double> zStepScale  = {1., 1.};

  ExtrudedPolygon absplate = ExtrudedPolygon( xCoord, yCoord, zStep, zStepX, zStepY, zStepScale);

  Volume      absplate_vol(basename, absplate, slice_mat);
  // Setting slice attributes
  if (renderComp){
    absplate_vol.setAttributes(desc, region, limit, vis);
  } else {
    absplate_vol.setAttributes(desc, region, limit, "InvisibleNoDaughters");
  }


  return absplate_vol;

}

//************************************************************************************************************
//************************** Filler plate i.e. air & kapton & PCB & ESR
//************************************************************************************************************
Volume createFillerPlate( Detector& desc,
                                   std::string basename,
                                   double h_mod,
                                   double w_mod,
                                   double t_mod_tp,
                                   double t_mod_sp,
                                   double t_slice,
                                   double w_notch,
                                   Material slice_mat,
                                   std::string region,
                                   std::string limit,
                                   std::string vis,
                                   bool renderComp
){
  double w_plate     = w_mod-2*t_mod_sp-w_notch;
  double h_plate     = h_mod-2*t_mod_tp;

  Box         filler( w_plate / 2., h_plate / 2., t_slice / 2.);
  Volume      filler_vol(basename, filler, slice_mat);
  // Setting slice attributes
  if (renderComp){
    filler_vol.setAttributes(desc, region, limit, vis);
  } else {
    filler_vol.setAttributes(desc, region, limit, "InvisibleNoDaughters");
  }

  return filler_vol;
}

//************************************************************************************************************
//************************** single scintillator plate for tower *********************************************
//************************************************************************************************************
Volume createScintillatorTower( Detector& desc,
                                   std::string basename,
                                   double w_tow,
                                   double h_tow,
                                   double t_slice,
                                   Material slice_mat,
                                   std::string region,
                                   std::string limit,
                                   std::string vis,
                                  SensitiveDetector sens,
                                  bool renderComp
){

  Box         scintplate( w_tow / 2., h_tow / 2., t_slice / 2.);
  Volume      slice_vol(basename, scintplate, slice_mat);
    // Setting appropriate slices as sensitive
  sens.setType("calorimeter");
  slice_vol.setSensitiveDetector(sens);
  // Setting slice attributes
  if (renderComp){
    slice_vol.setAttributes(desc, region, limit, vis);
  } else {
    slice_vol.setAttributes(desc, region, limit, "InvisibleNoDaughters");
  }
  return slice_vol;

}

//************************************************************************************************************
//************************** create scintillator plate with separations for 8M *******************************
//************************************************************************************************************
Assembly createScintillatorPlateEightM( Detector& desc,
                                        std::string basename,
//                                         int modID,
                                        int layerID,
                                        double h_mod,
                                        double w_mod,
                                        double t_mod_tp,
                                        double t_mod_sp,
                                        double t_slice,
                                        double w_notch,
                                        double t_foil,
                                        Material slice_mat,
                                        int roLayer,
                                        std::string region,
                                        std::string limit,
                                        std::string vis,
                                        SensitiveDetector sens,
                                        bool renderComp
){
  // Tower placement in 8M module
  //======================================================================
  //||              ||              ||                ||                ||
  //||      0       ||      1       ||        2       ||        3       ||
  //||              ||              ||                ||                ||
  //======================================================================
  //||              ||              ||                ||                ||
  //||      4       ||      5       ||        6       ||        7       ||
  //||              ||              ||                ||                ||
  //======================================================================
  Assembly modScintAssembly(basename);
  double w_plate  = w_mod-w_notch-2*t_mod_sp-2*t_foil;
  double h_plate  = h_mod-2*t_mod_tp-2*t_foil;
  double w_tow    = (w_plate-6*t_foil)/4;
  double h_tow    = (h_plate-2*t_foil)/2;

  // placement volumes
  PlacedVolume pvm;

  // foil separations                     // 0              1                   2                   3                   4
  const std::vector<double> xCoordTi = {  -(w_plate/2.),    -(w_tow+3*t_foil),  -(w_tow+3*t_foil),  -(w_tow+1*t_foil),  -(w_tow+1*t_foil),
                                          // 5              6                   7                   8                   9
                                          -t_foil,          -t_foil,            t_foil,             t_foil,             w_tow+1*t_foil,
                                          // 10             11                  12                  13                  14
                                          w_tow+1*t_foil,   w_tow+3*t_foil,     w_tow+3*t_foil,     w_plate/2.,         w_plate/2.,
                                          // 15             16                  17                  18                  19
                                          w_tow+3*t_foil,   w_tow+3*t_foil,     w_tow+1*t_foil,     w_tow+1*t_foil,     t_foil,
                                          // 20             21                  22                  23                  24
                                          t_foil,           -t_foil,            -t_foil,            -(w_tow+1*t_foil),  -(w_tow+1*t_foil),
                                          // 25             26                  27
                                          -(w_tow+3*t_foil),-(w_tow+3*t_foil),  -(w_plate/2.)
                                       };
                                          // 0              1                   2                   3                   4
  const std::vector<double> yCoordTi = {  t_foil,           t_foil,             (h_plate/2.),       (h_plate/2.),       t_foil,
                                          // 5              6                   7                   8                   9
                                          t_foil,           (h_plate/2.),       (h_plate/2.),       t_foil,             t_foil,
                                          // 10             11                  12                  13                  14
                                          (h_plate/2.),     (h_plate/2.),       t_foil,             t_foil,             -t_foil,
                                          // 15             16                  17                  18                  19
                                          -t_foil,          -(h_plate/2.),      -(h_plate/2.),      -t_foil,            -t_foil,
                                          // 20             21                  22                  23                  24
                                          -(h_plate/2.),    -(h_plate/2.),      -t_foil,            -t_foil,            -(h_plate/2.),
                                          // 25             26                  27
                                          -(h_plate/2.),    -t_foil,            -t_foil
                                        };

  const std::vector<double> zStepTi       = {-t_slice/2, t_slice/2};
  const std::vector<double> zStepXTi      = {0., 0.};
  const std::vector<double> zStepYTi      = {0., 0.};
  const std::vector<double> zStepScaleTi  = {1., 1.};

  ExtrudedPolygon foilgrid = ExtrudedPolygon( xCoordTi, yCoordTi, zStepTi, zStepXTi, zStepYTi, zStepScaleTi);
  Box         foil_t( (w_plate+2*t_foil) / 2., t_foil / 2., t_slice / 2.);
  Box         foil_s( t_foil / 2., h_plate / 2., t_slice / 2.);
  Volume      foilgrid_vol(basename+"_ESRFoil_"+_toString(layerID, "_layer_%d"), foilgrid, slice_mat);
  Volume      foil_t_vol(basename+"_ESRFoilT_"+_toString(layerID, "_layer_%d"), foil_t, slice_mat);
  Volume      foil_b_vol(basename+"_ESRFoilB_"+_toString(layerID, "_layer_%d"), foil_t, slice_mat);
  Volume      foil_l_vol(basename+"_ESRFoilL_"+_toString(layerID, "_layer_%d"), foil_s, slice_mat);
  Volume      foil_r_vol(basename+"_ESRFoilR_"+_toString(layerID, "_layer_%d"), foil_s, slice_mat);
  // Setting slice attributes
  if (renderComp){
    foilgrid_vol.setAttributes(desc, region, limit, "LFHCALLayerSepVis");
    foil_t_vol.setAttributes(desc, region, limit, "LFHCALLayerSepVis");
    foil_b_vol.setAttributes(desc, region, limit, "LFHCALLayerSepVis");
    foil_l_vol.setAttributes(desc, region, limit, "LFHCALLayerSepVis");
    foil_r_vol.setAttributes(desc, region, limit, "LFHCALLayerSepVis");
  } else {
    foilgrid_vol.setAttributes(desc, region, limit, "InvisibleNoDaughters");
    foil_t_vol.setAttributes(desc, region, limit, "InvisibleNoDaughters");
    foil_b_vol.setAttributes(desc, region, limit, "InvisibleNoDaughters");
    foil_l_vol.setAttributes(desc, region, limit, "InvisibleNoDaughters");
    foil_r_vol.setAttributes(desc, region, limit, "InvisibleNoDaughters");
  }
  pvm = modScintAssembly.placeVolume(foilgrid_vol, Position(0, 0, 0 ));
  pvm = modScintAssembly.placeVolume(foil_t_vol, Position(0, 1.5*t_foil+h_tow, 0 ));
  pvm = modScintAssembly.placeVolume(foil_b_vol, Position(0, -(1.5*t_foil+h_tow), 0 ));
  pvm = modScintAssembly.placeVolume(foil_l_vol, Position(-(3.5*t_foil+2*w_tow), 0, 0 ));
  pvm = modScintAssembly.placeVolume(foil_r_vol, Position((3.5*t_foil+2*w_tow), 0, 0 ));

  // 8M module placement of scintillator for tower
  double rotZ[8] = {0,  0,  0,  0,  0,  0,  0,  0};
  double rotY[8] = {0,  0,  0,  0,  0,  0,  0,  0};
  double rotX[8] = {0,  0,  0,  0,  0,  0,  0,  0};
  double posX[8] = {(w_tow*1.5+3*t_foil),     (w_tow*0.5+t_foil),     -(w_tow*0.5+t_foil),    -(w_tow*1.5+3*t_foil),
                    (w_tow*1.5+3*t_foil),     (w_tow*0.5+t_foil),     -(w_tow*0.5+t_foil),    -(w_tow*1.5+3*t_foil)};
  double posY[8] = {0.5*(h_tow)+t_foil,       0.5*(h_tow)+t_foil,     0.5*(h_tow)+t_foil,     0.5*(h_tow)+t_foil,
                    -(0.5*(h_tow)+t_foil),    -(0.5*(h_tow)+t_foil),  -(0.5*(h_tow)+t_foil),  -(0.5*(h_tow)+t_foil)};
  double posZ[8] = {0,                        0,                      0,                      0,
                    0,                        0,                      0,                      0};
  int towerx  = 0;
  int towery  = 0;

  // loop over all towers within same module
  for (int i = 0; i < 8; i++){
    // printout(DEBUG, "LFHCAL_geo", basename + _toString(i, "_tower_%d") + "\t" + _toString(modID) + "\t" + _toString(i) + "\t" + _toString(layerID));
    Volume modScintTowerAss = createScintillatorTower( desc,  basename+ _toString(i, "_tower_%d"),
                                                            w_tow, h_tow, t_slice,
                                                            slice_mat, region, limit, vis, sens, renderComp);
    pvm = modScintAssembly.placeVolume(modScintTowerAss, Transform3D(RotationZYX(rotZ[i], rotY[i], rotX[i]), Position(posX[i], posY[i], posZ[i] )));
    towerx = i%4;
    towery = 0;
    if (i > 3) towery = 1;
    pvm.addPhysVolID("towerx", towerx).addPhysVolID("towery", towery).addPhysVolID("layerz", layerID).addPhysVolID("passive", 0).addPhysVolID("rlayerz", roLayer);
  }
  return modScintAssembly;
}

//************************************************************************************************************
//************************** create scintillator plate with separations for 4M *******************************
//************************************************************************************************************
Assembly createScintillatorPlateFourM( Detector& desc,
                                        std::string basename,
//                                         int modID,
                                        int layerID,
                                        double h_mod,
                                        double w_mod,
                                        double t_mod_tp,
                                        double t_mod_sp,
                                        double t_slice,
                                        double w_notch,
                                        double t_foil,
                                        Material slice_mat,
                                        int roLayer,
                                        std::string region,
                                        std::string limit,
                                        std::string vis,
                                        SensitiveDetector sens,
                                        bool renderComp
){
  // Tower placement in 4M module
  //--------------------------------
  //|              ||              |
  //|      0       ||      1       |
  //|              ||              |
  //|==============================|
  //|              ||              |
  //|      2       ||      3       |
  //|              ||              |
  //--------------------------------
  Assembly modScintAssembly(basename);

  double w_plate  = w_mod-w_notch-2*t_mod_sp-2*t_foil;
  double h_plate  = h_mod-2*t_mod_tp-2*t_foil;
  double w_tow    = (w_plate-2*t_foil)/2;
  double h_tow    = (h_plate-2*t_foil)/2;

  // placement volumes
  PlacedVolume pvm;

  // foil separations
                                          // 0            1             2             3               4
  const std::vector<double> xCoordTi = {  -(w_plate/2.),  -t_foil,      -t_foil,      t_foil,         t_foil,
                                          // 5            6             7             8               9
                                          w_plate/2.,     w_plate/2.,   t_foil,       t_foil,         -t_foil,
                                          // 10           11
                                          -t_foil,        -(w_plate/2.)
                                       };
                                          // 0            1             2             3               4
  const std::vector<double> yCoordTi = {  t_foil,         t_foil,       (h_plate/2.), (h_plate/2.),   t_foil,
                                          // 5            6             7             8               9
                                          t_foil,         -t_foil,      -t_foil,      -(h_plate/2.),  -(h_plate/2.),
                                          // 10           11
                                          -t_foil,        -t_foil,
                                        };

  const std::vector<double> zStepTi       = {-t_slice/2, t_slice/2};
  const std::vector<double> zStepXTi      = {0., 0.};
  const std::vector<double> zStepYTi      = {0., 0.};
  const std::vector<double> zStepScaleTi  = {1., 1.};

  ExtrudedPolygon foilgrid = ExtrudedPolygon( xCoordTi, yCoordTi, zStepTi, zStepXTi, zStepYTi, zStepScaleTi);
  Box         foil_t( (w_plate+2*t_foil) / 2., t_foil / 2., t_slice / 2.);
  Box         foil_s( t_foil / 2., h_plate / 2., t_slice / 2.);
  Volume      foilgrid_vol(basename+"_ESRFoil_"+_toString(layerID, "_layer_%d"), foilgrid, slice_mat);
  Volume      foil_t_vol(basename+"_ESRFoilT_"+_toString(layerID, "_layer_%d"), foil_t, slice_mat);
  Volume      foil_b_vol(basename+"_ESRFoilB_"+_toString(layerID, "_layer_%d"), foil_t, slice_mat);
  Volume      foil_l_vol(basename+"_ESRFoilL_"+_toString(layerID, "_layer_%d"), foil_s, slice_mat);
  Volume      foil_r_vol(basename+"_ESRFoilR_"+_toString(layerID, "_layer_%d"), foil_s, slice_mat);
  // Setting slice attributes
  if (renderComp){
    foilgrid_vol.setAttributes(desc, region, limit, "LFHCALLayerSepVis");
    foil_t_vol.setAttributes(desc, region, limit, "LFHCALLayerSepVis");
    foil_b_vol.setAttributes(desc, region, limit, "LFHCALLayerSepVis");
    foil_l_vol.setAttributes(desc, region, limit, "LFHCALLayerSepVis");
    foil_r_vol.setAttributes(desc, region, limit, "LFHCALLayerSepVis");
  } else {
    foilgrid_vol.setAttributes(desc, region, limit, "InvisibleNoDaughters");
    foil_t_vol.setAttributes(desc, region, limit, "InvisibleNoDaughters");
    foil_b_vol.setAttributes(desc, region, limit, "InvisibleNoDaughters");
    foil_l_vol.setAttributes(desc, region, limit, "InvisibleNoDaughters");
    foil_r_vol.setAttributes(desc, region, limit, "InvisibleNoDaughters");
  }
  pvm = modScintAssembly.placeVolume(foilgrid_vol, Position(0, 0, 0 ));
  pvm = modScintAssembly.placeVolume(foil_t_vol, Position(0, 1.5*t_foil+h_tow, 0 ));
  pvm = modScintAssembly.placeVolume(foil_b_vol, Position(0, -(1.5*t_foil+h_tow), 0 ));
  pvm = modScintAssembly.placeVolume(foil_l_vol, Position(-(1.5*t_foil+w_tow), 0, 0 ));
  pvm = modScintAssembly.placeVolume(foil_r_vol, Position((1.5*t_foil+w_tow), 0, 0 ));


  // 4M module placement of scintillator for tower
  double rotZ[4] = {0,    0, 0,     0   };
  double rotY[4] = {0,    0, 0,     0   };
  double rotX[4] = {0,    0, 0,     0   };
  double posX[4] = {(w_tow*0.5+t_foil),     -(w_tow*0.5+t_foil),
                    (w_tow*0.5+t_foil),     -(w_tow*0.5+t_foil)};
  double posY[4] = {0.5*(h_tow)+t_foil,    0.5*(h_tow)+t_foil,
                    -(0.5*(h_tow)+t_foil), -(0.5*(h_tow)+t_foil)};
  double posZ[4] = {0,                              0,
                    0,                              0};
  int towerx  = 0;
  int towery  = 0;
  // loop over all towers within same module

  for (int i = 0; i < 4; i++){
    // printout(DEBUG, "LFHCAL_geo", basename + _toString(i, "_tower_%d") + "\t" + _toString(modID) + "\t" + _toString(i) + "\t" + _toString(layerID));
    Volume modScintTowerAss = createScintillatorTower( desc,  basename+ _toString(i, "_tower_%d"),
                                                             w_tow, h_tow, t_slice,
                                                            slice_mat, region, limit, vis, sens, renderComp);
    pvm = modScintAssembly.placeVolume(modScintTowerAss, Transform3D(RotationZYX(rotZ[i], rotY[i], rotX[i]), Position(posX[i], posY[i], posZ[i] )));
    towerx = i%2;
    towery = 0;
    if (i > 1) towery = 1;
    pvm.addPhysVolID("towerx", towerx).addPhysVolID("towery", towery).addPhysVolID("layerz", layerID).addPhysVolID("passive", 0).addPhysVolID("rlayerz", roLayer);
  }
  return modScintAssembly;
}

//************************************************************************************************************
//************************** create 8M module assembly  ******************************************************
//************************************************************************************************************
Volume createEightMModule ( Detector& desc,
                              moduleParamsStrct mod_params,
                              std::vector<sliceParamsStrct> sl_params,
//                               int modID,
                              double length,
                              SensitiveDetector sens,
                              bool renderComp,
                              bool allSen
){
  std::string baseName = "LFHCAL_8M";

  // assembly definition
  Box         modBox( mod_params.mod_width / 2., mod_params.mod_height / 2., length / 2.);
  Volume  vol_mod(baseName,modBox,desc.material("Air"));
  vol_mod.setVisAttributes(desc.visAttributes(mod_params.mod_visStr.data()));

  // placement operator
  PlacedVolume pvm;
  // ********************************************************************************
  // Casing definition
  // ********************************************************************************
  // geom definition 8M module casing
  Box         modFrontPlate( mod_params.mod_width / 2., mod_params.mod_height / 2., mod_params.mod_FWThick / 2.);
  Box         modSidePlateL( mod_params.mod_SWThick / 2., mod_params.mod_height / 2., (length-mod_params.mod_FWThick-mod_params.mod_BWThick) / 2.);
  Box         modSidePlateR( mod_params.mod_SWThick / 2., mod_params.mod_height / 2., (length-mod_params.mod_FWThick-mod_params.mod_BWThick) / 2.);
  Box         modTopPlate( (mod_params.mod_width-2*mod_params.mod_SWThick) / 2., mod_params.mod_TWThick / 2., (length-mod_params.mod_FWThick-mod_params.mod_BWThick) / 2.);
  Box         modBottomPlate( (mod_params.mod_width-2*mod_params.mod_SWThick) / 2., mod_params.mod_TWThick / 2., (length-mod_params.mod_FWThick-mod_params.mod_BWThick) / 2.);
  Box         modBackCutOut( mod_params.mod_BIwidth / 2., mod_params.mod_BIheight / 2., mod_params.mod_BWThick / 2.);
  Box         modBackPlateFull( mod_params.mod_width / 2., mod_params.mod_height / 2., mod_params.mod_BWThick / 2.);
  SubtractionSolid modBackPlate(modBackPlateFull, modBackCutOut);

  // volume definition 8M module casing
  Volume  vol_modFrontPlate(baseName+"_FrontPlate",modFrontPlate,desc.material("Steel235"));
  Volume  vol_modBackPlate(baseName+"_BackPlate",modBackPlate,desc.material("Steel235"));
  Volume  vol_modSidePlateL(baseName+"_LeftSidePlate",modSidePlateL,desc.material("Steel235"));
  Volume  vol_modSidePlateR(baseName+"_RightSidePlate",modSidePlateR,desc.material("Steel235"));
  Volume  vol_modTopPlate(baseName+"_TopPlate",modTopPlate,desc.material("Steel235"));
  Volume  vol_modBottomPlate(baseName+"_BottomPlate",modBottomPlate,desc.material("Steel235"));

  if (allSen){
     sens.setType("calorimeter");
     vol_modFrontPlate.setSensitiveDetector(sens);
     vol_modBackPlate.setSensitiveDetector(sens);
     vol_modSidePlateL.setSensitiveDetector(sens);
     vol_modSidePlateR.setSensitiveDetector(sens);
     vol_modTopPlate.setSensitiveDetector(sens);
     vol_modBottomPlate.setSensitiveDetector(sens);
  }

  if (renderComp){
    vol_modFrontPlate.setAttributes(desc, mod_params.mod_regStr, mod_params.mod_limStr, mod_params.mod_visStr);
    vol_modBackPlate.setAttributes(desc, mod_params.mod_regStr, mod_params.mod_limStr, mod_params.mod_visStr);
    vol_modSidePlateL.setAttributes(desc, mod_params.mod_regStr, mod_params.mod_limStr, mod_params.mod_visStr);
    vol_modSidePlateR.setAttributes(desc, mod_params.mod_regStr, mod_params.mod_limStr, mod_params.mod_visStr);
    vol_modTopPlate.setAttributes(desc, mod_params.mod_regStr, mod_params.mod_limStr, mod_params.mod_visStr);
    vol_modBottomPlate.setAttributes(desc, mod_params.mod_regStr, mod_params.mod_limStr, mod_params.mod_visStr);
  } else {
    vol_modFrontPlate.setAttributes(desc, mod_params.mod_regStr, mod_params.mod_limStr, "InvisibleNoDaughters");
    vol_modBackPlate.setAttributes(desc, mod_params.mod_regStr, mod_params.mod_limStr, "InvisibleNoDaughters");
    vol_modSidePlateL.setAttributes(desc, mod_params.mod_regStr, mod_params.mod_limStr, "InvisibleNoDaughters");
    vol_modSidePlateR.setAttributes(desc, mod_params.mod_regStr, mod_params.mod_limStr, "InvisibleNoDaughters");
    vol_modTopPlate.setAttributes(desc, mod_params.mod_regStr, mod_params.mod_limStr, "InvisibleNoDaughters");
    vol_modBottomPlate.setAttributes(desc, mod_params.mod_regStr, mod_params.mod_limStr, "InvisibleNoDaughters");
  }

  // ********************************************************************************
  // long PCB
  // ********************************************************************************
  Box         modPCB( mod_params.mod_pcbThick / 2., mod_params.mod_pcbWidth / 2., (mod_params.mod_pcbLength) / 2.);
  Volume  vol_modPCB(baseName+"_PCB",modPCB,desc.material("Fr4"));
  if (renderComp){
    vol_modPCB.setAttributes(desc, mod_params.mod_regStr, mod_params.mod_limStr, "LFHCALModPCB");
  } else {
    vol_modPCB.setAttributes(desc, mod_params.mod_regStr, mod_params.mod_limStr, "InvisibleNoDaughters");
  }

  int    layer_num = 0;
  double slice_z   = -length/2+mod_params.mod_FWThick; // Keeps track of layers' local z locations
  // Looping through the number of repeated layers & slices in each section
  for (int i = 0; i < (int)sl_params.size(); i++){
    slice_z += sl_params[i].slice_offset + sl_params[i].slice_thick / 2.; // Going to halfway point in layer
    layer_num = sl_params[i].layer_ID;
    //*************************************************
    // absorber plates
    //*************************************************
    Material slice_mat = desc.material(sl_params[i].slice_matStr);
    if (sl_params[i].slice_partID == 1 ){
      Volume modAbsAssembly = createAbsorberPlate( desc,
                                                          baseName+"_Abs"+_toString(sl_params[i].layer_ID, "_layer_%d"),
                                                          mod_params.mod_height, mod_params.mod_width, mod_params.mod_TWThick, mod_params.mod_SWThick,
                                                          sl_params[i].slice_thick, mod_params.mod_notchDepth,
                                                          mod_params.mod_notchHeight,
                                                          slice_mat, sl_params[i].slice_regStr, sl_params[i].slice_limStr, sl_params[i].slice_visStr, renderComp);
      // Placing slice within layer
      if (allSen) modAbsAssembly.setSensitiveDetector(sens);
      pvm = vol_mod.placeVolume(modAbsAssembly, Transform3D(RotationZYX(0, 0, 0), Position(0., 0., slice_z)));
      if (allSen) pvm.addPhysVolID("towerx", 0).addPhysVolID("towery", 0).addPhysVolID("rlayerz", sl_params[i].slice_readoutLayer).addPhysVolID("layerz", layer_num).addPhysVolID("passive", 1);
    //*************************************************
    // air & kapton & PCB & ESR
    //*************************************************
    } else if (sl_params[i].slice_partID == 2 ){
      Volume modFillAssembly =  createFillerPlate( desc,
                                                   baseName+"_Fill"+_toString(sl_params[i].layer_ID, "_layer_%d")+_toString(sl_params[i].slice_ID, "slice_%d"),
                                                   mod_params.mod_height, mod_params.mod_width, mod_params.mod_TWThick, mod_params.mod_SWThick,
                                                   sl_params[i].slice_thick, mod_params.mod_notchDepth,
                                                   slice_mat, sl_params[i].slice_regStr, sl_params[i].slice_limStr, sl_params[i].slice_visStr, renderComp);
      // Placing slice within layer
      if (allSen) modFillAssembly.setSensitiveDetector(sens);
      pvm = vol_mod.placeVolume(modFillAssembly, Transform3D(RotationZYX(0, 0, 0), Position((mod_params.mod_notchDepth)/2., 0., slice_z)));
      if (allSen) pvm.addPhysVolID("towerx", 1).addPhysVolID("towery", 0).addPhysVolID("rlayerz", sl_params[i].slice_readoutLayer).addPhysVolID("layerz", layer_num).addPhysVolID("passive", 1);
    //*************************************************
    // scintillator
    //*************************************************
    } else {
      Assembly modScintAssembly =  createScintillatorPlateEightM(  desc,
                                                                  baseName+"_ScintAssembly"+_toString(sl_params[i].layer_ID, "_layer_%d"),
                                                                  layer_num,
                                                                  mod_params.mod_height, mod_params.mod_width, mod_params.mod_TWThick, mod_params.mod_SWThick,
                                                                  sl_params[i].slice_thick, mod_params.mod_notchDepth, mod_params.mod_foilThick,
                                                                  slice_mat, sl_params[i].slice_readoutLayer ,sl_params[i].slice_regStr, sl_params[i].slice_limStr, sl_params[i].slice_visStr, sens, renderComp);
      // Placing slice within layer
      pvm = vol_mod.placeVolume(modScintAssembly, Transform3D(RotationZYX(0, 0, 0), Position((mod_params.mod_notchDepth)/2., 0, slice_z)));
    }
    slice_z += sl_params[i].slice_thick / 2.;
  }

  // placement 8M module casing
  pvm = vol_mod.placeVolume(vol_modFrontPlate, Position(0, 0, -( length-mod_params.mod_FWThick) / 2. ));
  if (allSen) pvm.addPhysVolID("towerx", 2).addPhysVolID("towery", 0).addPhysVolID("layerz", 0).addPhysVolID("passive", 1);
  pvm = vol_mod.placeVolume(vol_modBackPlate, Position(0, 0, ( length-mod_params.mod_BWThick) / 2. ));
  if (allSen) pvm.addPhysVolID("towerx", 2).addPhysVolID("towery", 0).addPhysVolID("layerz", layer_num).addPhysVolID("passive", 1);
  pvm = vol_mod.placeVolume(vol_modSidePlateL, Position(-(mod_params.mod_width-mod_params.mod_SWThick)/2., 0,  (mod_params.mod_FWThick-mod_params.mod_BWThick)/2));
  if (allSen) pvm.addPhysVolID("towerx", 3).addPhysVolID("towery", 0).addPhysVolID("layerz", 0).addPhysVolID("passive", 1);
  pvm = vol_mod.placeVolume(vol_modSidePlateR, Position((mod_params.mod_width-mod_params.mod_SWThick)/2., 0,(mod_params.mod_FWThick-mod_params.mod_BWThick)/2));
  if (allSen) pvm.addPhysVolID("towerx", 0).addPhysVolID("towery", 1).addPhysVolID("layerz", 0).addPhysVolID("passive", 1);
  pvm = vol_mod.placeVolume(vol_modTopPlate, Position(0, (mod_params.mod_height-mod_params.mod_TWThick)/2., (mod_params.mod_FWThick-mod_params.mod_BWThick)/2));
  if (allSen) pvm.addPhysVolID("towerx", 1).addPhysVolID("towery", 1).addPhysVolID("layerz", 0).addPhysVolID("passive", 1);
  pvm = vol_mod.placeVolume(vol_modBottomPlate, Position(0, -(mod_params.mod_height-mod_params.mod_TWThick)/2., (mod_params.mod_FWThick-mod_params.mod_BWThick)/2));
  if (allSen) pvm.addPhysVolID("towerx", 2).addPhysVolID("towery", 1).addPhysVolID("layerz", 0).addPhysVolID("passive", 1);

  double lengthA      = length-mod_params.mod_FWThick-mod_params.mod_BWThick;
  double z_offSetPCB  = (mod_params.mod_FWThick-mod_params.mod_BWThick)/2-(lengthA-mod_params.mod_pcbLength)/2.;

  pvm = vol_mod.placeVolume(vol_modPCB, Position(-(mod_params.mod_width-2*mod_params.mod_SWThick-mod_params.mod_notchDepth)/2., 0, z_offSetPCB));

  return vol_mod;
}


//************************************************************************************************************
//************************** create 8M module assembly  ******************************************************
//************************************************************************************************************
Volume createFourMModule ( Detector& desc,
                              moduleParamsStrct mod_params,
                              std::vector<sliceParamsStrct> sl_params,
//                               int modID,
                              double length,
                              SensitiveDetector sens,
                              bool renderComp,
                              bool allSen
){

  std::string baseName = "LFHCAL_4M";

  // assembly definition
  Box         modBox( mod_params.mod_width / 2., mod_params.mod_height / 2., length / 2.);
  Volume  vol_mod(baseName,modBox,desc.material("Air"));
  printout(DEBUG, "LFHCAL_geo", "visualization string module: " + mod_params.mod_visStr);
  vol_mod.setVisAttributes(desc.visAttributes(mod_params.mod_visStr.data()));

  // placement operator
  PlacedVolume pvm;
  // ********************************************************************************
  // Casing definition
  // ********************************************************************************
  // geom definition 8M module casing
  Box         modFrontPlate( mod_params.mod_width / 2., mod_params.mod_height / 2., mod_params.mod_FWThick / 2.);
  Box         modSidePlateL( mod_params.mod_SWThick / 2., mod_params.mod_height / 2., (length-mod_params.mod_FWThick-mod_params.mod_BWThick) / 2.);
  Box         modSidePlateR( mod_params.mod_SWThick / 2., mod_params.mod_height / 2., (length-mod_params.mod_FWThick-mod_params.mod_BWThick) / 2.);
  Box         modTopPlate( (mod_params.mod_width-2*mod_params.mod_SWThick) / 2., mod_params.mod_TWThick / 2., (length-mod_params.mod_FWThick-mod_params.mod_BWThick) / 2.);
  Box         modBottomPlate( (mod_params.mod_width-2*mod_params.mod_SWThick) / 2., mod_params.mod_TWThick / 2., (length-mod_params.mod_FWThick-mod_params.mod_BWThick) / 2.);
  Box         modBackCutOut( mod_params.mod_BIwidth / 2., mod_params.mod_BIheight / 2., mod_params.mod_BWThick / 2.);
  Box         modBackPlateFull( mod_params.mod_width / 2., mod_params.mod_height / 2., mod_params.mod_BWThick / 2.);
  SubtractionSolid modBackPlate(modBackPlateFull, modBackCutOut);

  // volume definition 8M module casing
  Volume  vol_modFrontPlate(baseName+"_FrontPlate",modFrontPlate,desc.material("Steel235"));
  Volume  vol_modBackPlate(baseName+"_BackPlate",modBackPlate,desc.material("Steel235"));
  Volume  vol_modSidePlateL(baseName+"_LeftSidePlate",modSidePlateL,desc.material("Steel235"));
  Volume  vol_modSidePlateR(baseName+"_RightSidePlate",modSidePlateR,desc.material("Steel235"));
  Volume  vol_modTopPlate(baseName+"_TopPlate",modTopPlate,desc.material("Steel235"));
  Volume  vol_modBottomPlate(baseName+"_BottomPlate",modBottomPlate,desc.material("Steel235"));

  if (allSen){
     sens.setType("calorimeter");
     vol_modFrontPlate.setSensitiveDetector(sens);
     vol_modBackPlate.setSensitiveDetector(sens);
     vol_modSidePlateL.setSensitiveDetector(sens);
     vol_modSidePlateR.setSensitiveDetector(sens);
     vol_modTopPlate.setSensitiveDetector(sens);
     vol_modBottomPlate.setSensitiveDetector(sens);
  }


  if (renderComp){
    vol_modFrontPlate.setAttributes(desc, mod_params.mod_regStr, mod_params.mod_limStr, mod_params.mod_visStr);
    vol_modBackPlate.setAttributes(desc, mod_params.mod_regStr, mod_params.mod_limStr, mod_params.mod_visStr);
    vol_modSidePlateL.setAttributes(desc, mod_params.mod_regStr, mod_params.mod_limStr, mod_params.mod_visStr);
    vol_modSidePlateR.setAttributes(desc, mod_params.mod_regStr, mod_params.mod_limStr, mod_params.mod_visStr);
    vol_modTopPlate.setAttributes(desc, mod_params.mod_regStr, mod_params.mod_limStr, mod_params.mod_visStr);
    vol_modBottomPlate.setAttributes(desc, mod_params.mod_regStr, mod_params.mod_limStr, mod_params.mod_visStr);
  } else {
    vol_modFrontPlate.setAttributes(desc, mod_params.mod_regStr, mod_params.mod_limStr, "InvisibleNoDaughters");
    vol_modBackPlate.setAttributes(desc, mod_params.mod_regStr, mod_params.mod_limStr, "InvisibleNoDaughters");
    vol_modSidePlateL.setAttributes(desc, mod_params.mod_regStr, mod_params.mod_limStr, "InvisibleNoDaughters");
    vol_modSidePlateR.setAttributes(desc, mod_params.mod_regStr, mod_params.mod_limStr, "InvisibleNoDaughters");
    vol_modTopPlate.setAttributes(desc, mod_params.mod_regStr, mod_params.mod_limStr, "InvisibleNoDaughters");
    vol_modBottomPlate.setAttributes(desc, mod_params.mod_regStr, mod_params.mod_limStr, "InvisibleNoDaughters");
  }

  // ********************************************************************************
  // long PCB
  // ********************************************************************************
  Box         modPCB( mod_params.mod_pcbThick / 2., mod_params.mod_pcbWidth / 2., (mod_params.mod_pcbLength) / 2.);
  Volume  vol_modPCB(baseName+"_PCB",modPCB,desc.material("Fr4"));
  if (renderComp){
    vol_modPCB.setAttributes(desc, mod_params.mod_regStr, mod_params.mod_limStr, "LFHCALModPCB");
  } else {
    vol_modPCB.setAttributes(desc, mod_params.mod_regStr, mod_params.mod_limStr, "InvisibleNoDaughters");
  }

  int    layer_num = 0;
  double slice_z   = -length/2+mod_params.mod_FWThick; // Keeps track of layers' local z locations

  // Looping through the number of repeated layers & slices in each section
  for (int i = 0; i < (int)sl_params.size(); i++){
    slice_z += sl_params[i].slice_offset + sl_params[i].slice_thick/2. ; // Going to halfway point in layer
    layer_num = sl_params[i].layer_ID;
    //*************************************************
    // absorber plates
    //*************************************************
    Material slice_mat = desc.material(sl_params[i].slice_matStr);
    if (sl_params[i].slice_partID == 1 ){
      Volume modAbsAssembly = createAbsorberPlate( desc, baseName+"_Abs"+_toString(sl_params[i].layer_ID, "_layer_%d"),
                                                              mod_params.mod_height, mod_params.mod_width, mod_params.mod_TWThick, mod_params.mod_SWThick,
                                                              sl_params[i].slice_thick, mod_params.mod_notchDepth, mod_params.mod_notchHeight,
                                                              slice_mat, sl_params[i].slice_regStr, sl_params[i].slice_limStr, sl_params[i].slice_visStr, renderComp);
      // Placing slice within layer
      if (allSen) modAbsAssembly.setSensitiveDetector(sens);
      pvm = vol_mod.placeVolume(modAbsAssembly, Transform3D(RotationZYX(0, 0, 0), Position(0., 0., slice_z)));
      if (allSen) pvm.addPhysVolID("towerx", 0).addPhysVolID("towery", 0).addPhysVolID("rlayerz", sl_params[i].slice_readoutLayer).addPhysVolID("layerz", layer_num).addPhysVolID("passive", 1);
    //*************************************************
    // air & kapton & PCB & ESR
    //*************************************************
    } else if (sl_params[i].slice_partID == 2 ){
      Volume modFillAssembly =  createFillerPlate( desc,
                                              baseName+"_Fill"+_toString(sl_params[i].layer_ID, "_layer_%d")+_toString(sl_params[i].slice_ID, "slice_%d"),
                                              mod_params.mod_height, mod_params.mod_width, mod_params.mod_TWThick, mod_params.mod_SWThick,
                                              sl_params[i].slice_thick, mod_params.mod_notchDepth,
                                              slice_mat, sl_params[i].slice_regStr, sl_params[i].slice_limStr, sl_params[i].slice_visStr, renderComp);
      // Placing slice within layer
      if (allSen) modFillAssembly.setSensitiveDetector(sens);
      pvm = vol_mod.placeVolume(modFillAssembly, Transform3D(RotationZYX(0, 0, 0), Position((mod_params.mod_notchDepth)/2.,0. , slice_z)));
      if (allSen) pvm.addPhysVolID("towerx", 1).addPhysVolID("towery", 0).addPhysVolID("rlayerz", sl_params[i].slice_readoutLayer).addPhysVolID("layerz", layer_num).addPhysVolID("passive", 1);
    //*************************************************
    // scintillator
    //*************************************************
    } else {
      Assembly modScintAssembly =  createScintillatorPlateFourM( desc,baseName+"_ScintAssembly"+_toString(sl_params[i].layer_ID, "_layer_%d"),
                                                                 layer_num, mod_params.mod_height, mod_params.mod_width, mod_params.mod_TWThick, mod_params.mod_SWThick,
                                                                  sl_params[i].slice_thick, mod_params.mod_notchDepth, mod_params.mod_foilThick,
                                                                  slice_mat, sl_params[i].slice_readoutLayer, sl_params[i].slice_regStr, sl_params[i].slice_limStr, sl_params[i].slice_visStr, sens, renderComp);
      // Placing slice within layer
      pvm = vol_mod.placeVolume(modScintAssembly, Transform3D(RotationZYX(0, 0, 0), Position((mod_params.mod_notchDepth)/2., 0, slice_z)));
    }
    slice_z += sl_params[i].slice_thick/2.;
  }

  // placement 4M module casing
  pvm = vol_mod.placeVolume(vol_modFrontPlate, Position(0, 0, -( length-mod_params.mod_FWThick) / 2. ));
  if (allSen) pvm.addPhysVolID("towerx", 2).addPhysVolID("towery", 0).addPhysVolID("layerz", 0).addPhysVolID("passive", 1);
  pvm = vol_mod.placeVolume(vol_modBackPlate, Position(0, 0, ( length-mod_params.mod_BWThick) / 2. ));
  if (allSen) pvm.addPhysVolID("towerx", 2).addPhysVolID("towery", 0).addPhysVolID("layerz", layer_num).addPhysVolID("passive", 1);
  pvm = vol_mod.placeVolume(vol_modSidePlateL, Position(-(mod_params.mod_width-mod_params.mod_SWThick)/2., 0,  (mod_params.mod_FWThick-mod_params.mod_BWThick)/2));
  if (allSen) pvm.addPhysVolID("towerx", 3).addPhysVolID("towery", 0).addPhysVolID("layerz", 0).addPhysVolID("passive", 1);
  pvm = vol_mod.placeVolume(vol_modSidePlateR, Position((mod_params.mod_width-mod_params.mod_SWThick)/2., 0,(mod_params.mod_FWThick-mod_params.mod_BWThick)/2));
  if (allSen) pvm.addPhysVolID("towerx", 0).addPhysVolID("towery", 1).addPhysVolID("layerz", 0).addPhysVolID("passive", 1);
  pvm = vol_mod.placeVolume(vol_modTopPlate, Position(0, (mod_params.mod_height-mod_params.mod_TWThick)/2., (mod_params.mod_FWThick-mod_params.mod_BWThick)/2));
  if (allSen) pvm.addPhysVolID("towerx", 1).addPhysVolID("towery", 1).addPhysVolID("layerz", 0).addPhysVolID("passive", 1);
  pvm = vol_mod.placeVolume(vol_modBottomPlate, Position(0, -(mod_params.mod_height-mod_params.mod_TWThick)/2., (mod_params.mod_FWThick-mod_params.mod_BWThick)/2));
  if (allSen) pvm.addPhysVolID("towerx", 2).addPhysVolID("towery", 1).addPhysVolID("layerz", 0).addPhysVolID("passive", 1);

  double lengthA      = length-mod_params.mod_FWThick-mod_params.mod_BWThick;
  double z_offSetPCB  = (mod_params.mod_FWThick-mod_params.mod_BWThick)/2-(lengthA-mod_params.mod_pcbLength)/2.;

  pvm = vol_mod.placeVolume(vol_modPCB, Position(-(mod_params.mod_width-2*mod_params.mod_SWThick-mod_params.mod_notchDepth)/2., 0, z_offSetPCB));
  return vol_mod;
}

//********************************************************************************************
//*                                                                                          *
//*                              Create detector                                             *
//==============================  MAIN FUNCTION  =============================================
//*                                                                                          *
//********************************************************************************************
static Ref_t createDetector(Detector& desc, xml_h handle, SensitiveDetector sens)
{
  // global detector variables
  xml_det_t   detElem = handle;
  std::string detName = detElem.nameStr();
  int         detID   = detElem.id();

  // general detector dimensions
  xml_dim_t dim    = detElem.dimensions();
  double    length = dim.z();    // Size along z-axis

  // general detector position
  xml_dim_t pos = detElem.position();
  printout(DEBUG, "LFHCAL_geo", "global LFHCal position " + _toString(pos.x()) + "\t" + _toString(pos.y()) + "\t" + _toString(pos.z()));

  // envelope volume
  xml_comp_t x_env = detElem.child(_Unicode(envelope));
  Tube rmaxtube(0, dim.rmax(), dim.z() / 2);
  Box beampipe(dim.x() / 2, dim.y() / 2, dim.z() / 2);
  Solid env = SubtractionSolid(rmaxtube, beampipe, Position(dim.x0(),0,0));
  Volume env_vol(detName + "_env", env, desc.material(x_env.materialStr()));

  bool renderComponents = getAttrOrDefault(detElem, _Unicode(renderComponents), 0.);
  bool allSensitive     = getAttrOrDefault(detElem, _Unicode(allSensitive), 0.);
  if (renderComponents) {
    printout(DEBUG, "LFHCAL_geo", "enabled visualization");
  } else {
    printout(DEBUG, "LFHCAL_geo", "switchted off visualization");
  }

  // 8M module specific loading
  xml_comp_t  eightM_xml        = detElem.child(_Unicode(eightmodule));
  xml_dim_t eightMmod_dim        = eightM_xml.dimensions();
  moduleParamsStrct eightM_params(getAttrOrDefault(eightMmod_dim, _Unicode(widthBackInner), 0.),
                                  getAttrOrDefault(eightMmod_dim, _Unicode(heightBackInner), 0.),
                                  getAttrOrDefault(eightMmod_dim, _Unicode(widthSideWall), 0.),
                                  getAttrOrDefault(eightMmod_dim, _Unicode(widthTopWall), 0.),
                                  getAttrOrDefault(eightMmod_dim, _Unicode(thicknessFrontWall), 0.),
                                  getAttrOrDefault(eightMmod_dim, _Unicode(thicknessBackWall), 0.),
                                  getAttrOrDefault(eightMmod_dim, _Unicode(width), 0.),
                                  getAttrOrDefault(eightMmod_dim, _Unicode(height), 0.),
                                  getAttrOrDefault(eightMmod_dim, _Unicode(notchDepth), 0.),
                                  getAttrOrDefault(eightMmod_dim, _Unicode(notchHeight), 0.),
                                  getAttrOrDefault(eightMmod_dim, _Unicode(foilThick), 0.),
                                  getAttrOrDefault(eightMmod_dim, _Unicode(pcbLength), 0.),
                                  getAttrOrDefault(eightMmod_dim, _Unicode(pcbThick), 0.),
                                  getAttrOrDefault(eightMmod_dim, _Unicode(pcbWidth), 0.),
                                  eightM_xml.visStr(),
                                  eightM_xml.regionStr(),
                                  eightM_xml.limitsStr()
                                 );

  // 4M module specific loading
  xml_comp_t  fourM_xml        = detElem.child(_Unicode(fourmodule));
  xml_dim_t fourMmod_dim        = fourM_xml.dimensions();
  moduleParamsStrct fourM_params( getAttrOrDefault(fourMmod_dim, _Unicode(widthBackInner), 0.),
                                  getAttrOrDefault(fourMmod_dim, _Unicode(heightBackInner), 0.),
                                  getAttrOrDefault(fourMmod_dim, _Unicode(widthSideWall), 0.),
                                  getAttrOrDefault(fourMmod_dim, _Unicode(widthTopWall), 0.),
                                  getAttrOrDefault(fourMmod_dim, _Unicode(thicknessFrontWall), 0.),
                                  getAttrOrDefault(fourMmod_dim, _Unicode(thicknessBackWall), 0.),
                                  getAttrOrDefault(fourMmod_dim, _Unicode(width), 0.),
                                  getAttrOrDefault(fourMmod_dim, _Unicode(height), 0.),
                                  getAttrOrDefault(fourMmod_dim, _Unicode(notchDepth), 0.),
                                  getAttrOrDefault(fourMmod_dim, _Unicode(notchHeight), 0.),
                                  getAttrOrDefault(fourMmod_dim, _Unicode(foilThick), 0.),
                                  getAttrOrDefault(fourMmod_dim, _Unicode(pcbLength), 0.),
                                  getAttrOrDefault(fourMmod_dim, _Unicode(pcbThick), 0.),
                                  getAttrOrDefault(fourMmod_dim, _Unicode(pcbWidth), 0.),
                                  fourM_xml.visStr(),
                                  fourM_xml.regionStr(),
                                  fourM_xml.limitsStr());

  std::vector<sliceParamsStrct> slice_Params;
  int layer_num   = 0;
  int readLayerC  = 0;
  for (xml_coll_t c(detElem, _U(layer)); c; ++c) {
    xml_comp_t x_layer         = c;
    int        repeat          = x_layer.repeat();
    int        readlayer       = getAttrOrDefault(x_layer, _Unicode(readoutlayer), 0.);
    if (readLayerC != readlayer){
      readLayerC  = readlayer;
      layer_num   = 0;
    }
    // Looping through the number of repeated layers in each section
    for (int i = 0; i < repeat; i++) {
      int    slice_num = 1;

      // Looping over each layer's slices
      for (xml_coll_t l(x_layer, _U(slice)); l; ++l) {
        xml_comp_t  x_slice         = l;
        sliceParamsStrct slice_param( layer_num,  slice_num, getAttrOrDefault(l, _Unicode(type), 0.), x_slice.thickness(), getAttrOrDefault(l, _Unicode(offset), 0.), readlayer,
                                      x_slice.materialStr(), x_slice.visStr(), x_slice.regionStr(), x_slice.limitsStr());
        slice_Params.push_back(slice_param);
        ++slice_num;
      }
      layer_num++;
    }
  }

  // create mother volume
  DetElement det(detName, detID);
  Assembly assembly(detName);
  PlacedVolume phv;

  int moduleIDx = -1;
  int moduleIDy = -1;


  struct position {
    double x,y,z;
  };

  std::vector<position> pos8M;

  xml_coll_t eightMPos(detElem, _Unicode(eightmodulepositions));
  for (xml_coll_t position_i(eightMPos, _U(position)); position_i; ++position_i){
    xml_comp_t position_comp = position_i;
    if (! getAttrOrDefault(position_comp, _Unicode(if), true)) {
      printout(DEBUG, "LFHCAL_geo", "skipping x = %.1f cm, y = %.1f cm", position_comp.x(), position_comp.y());
      continue;
    }
    pos8M.push_back({position_comp.x(), position_comp.y(), position_comp.z()});
  }

  // create 8M modules
  Volume  eightMassembly = createEightMModule ( desc, eightM_params, slice_Params, length, sens, renderComponents, allSensitive);
  for (int e = 0; e < (int)pos8M.size(); e++){
    if(e%20 == 0 ) printout(DEBUG, "LFHCAL_geo", "LFHCAL placing 8M module: " + _toString(e) + "/" +  _toString((int)pos8M.size()) + "\t" + _toString(pos8M[e].x) + "\t" + _toString(pos8M[e].y) + "\t" + _toString(pos8M[e].z));
      if(moduleIDx<0 || moduleIDy<0){
      printout(DEBUG, "LFHCAL_geo", "LFHCAL WRONG ID FOR 8M module: " + _toString(e) + "/" + _toString((int)pos8M.size()) + "\t" + _toString(moduleIDx) + "\t"
                + _toString(moduleIDy));
    }
    moduleIDx             = ((pos8M[e].x + 270) / 10);
    moduleIDy             = ((pos8M[e].y + 265) / 10);

    // Placing modules in world volume
    auto tr8M = Transform3D(Position(-pos8M[e].x-0.5*eightM_params.mod_width, -pos8M[e].y, pos8M[e].z));
    phv = assembly.placeVolume(eightMassembly, tr8M);
    phv.addPhysVolID("moduleIDx", moduleIDx).addPhysVolID("moduleIDy", moduleIDy).addPhysVolID("moduletype", 0);
  }

  std::vector<position> pos4M;

  xml_coll_t fourMPos(detElem, _Unicode(fourmodulepositions));
  for (xml_coll_t position_i(fourMPos, _U(position)); position_i; ++position_i){
    xml_comp_t position_comp = position_i;
    if (! getAttrOrDefault(position_comp, _Unicode(if), true)) {
      printout(DEBUG, "LFHCAL_geo", "skipping x = %.1f cm, y = %.1f cm", position_comp.x(), position_comp.y());
      continue;
    }
    pos4M.push_back({position_comp.x(), position_comp.y(), position_comp.z()});
  }

  // create 4M modules
  Volume  fourMassembly = createFourMModule ( desc, fourM_params, slice_Params,  length, sens, renderComponents, allSensitive);
  for (int f = 0; f < (int)pos4M.size(); f++){
    if(f%20 == 0 ) printout(DEBUG, "LFHCAL_geo", "LFHCAL placing 4M module: " + _toString(f) + "/" + _toString((int)pos4M.size()) + "\t" + _toString(pos4M[f].x) + "\t" + _toString(pos4M[f].y) + "\t" + _toString(pos4M[f].z));

    moduleIDx             = ((pos4M[f].x + 265) / 10);
    moduleIDy             = ((pos4M[f].y + 265) / 10);
    if(moduleIDx<0 || moduleIDy<0){
      printout(DEBUG, "LFHCAL_geo", "LFHCAL WRONG ID FOR 4M module: " + _toString(f) + "/" + _toString((int)pos4M.size()) + "\t" + _toString(moduleIDx) + "\t"
                + _toString(moduleIDy));
    }
    auto tr4M = Transform3D(Position(-pos4M[f].x-0.5*fourM_params.mod_width, -pos4M[f].y, pos4M[f].z));
    phv = assembly.placeVolume(fourMassembly, tr4M);
    phv.addPhysVolID("moduleIDx", moduleIDx).addPhysVolID("moduleIDy", moduleIDy).addPhysVolID("moduletype", 1);
  }

  Volume motherVol = desc.pickMotherVolume(det);
  phv = env_vol.placeVolume(assembly);
  phv = motherVol.placeVolume(env_vol, Transform3D(Position(pos.x(), pos.y(), pos.z() + length / 2.)));
  phv.addPhysVolID("system", detID);
  det.setPlacement(phv);

  return det;
}
DECLARE_DETELEMENT(epic_LFHCAL, createDetector)
