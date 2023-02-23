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
                      mod_BWThick(0.), mod_width(0.), mod_height(0.), mod_notchWidthAbsA(0.), mod_notchWidthAbsB(0.), mod_notchWidthAbsC (0.),
                      mod_notchWidthScin(0.), mod_notchDepth(0.), mod_sepDepth(0.),  mod_visStr(""), mod_regStr(""), mod_limStr("")
                      {}
  moduleParamsStrct(   double BIwidth, double BIheight, double SWThick, double TWThick, double FWThick, double BWThick, double width, double height, 
                       double notchWidthAbsA, double notchWidthAbsB, double notchWidthAbsC, double notchWidthScin, double notchDepth, double sepDepth, 
                       std::string visStr, std::string regStr, std::string limStr){
      mod_BIwidth       = BIwidth;
      mod_BIheight      = BIheight;
      mod_SWThick       = SWThick;  
      mod_TWThick       = TWThick;
      mod_FWThick       = FWThick;
      mod_BWThick       = BWThick;
      mod_width         = width;
      mod_height        = height;
      mod_notchWidthAbsA  = notchWidthAbsA;
      mod_notchWidthAbsB  = notchWidthAbsB;
      mod_notchWidthAbsC  = notchWidthAbsC;
      mod_notchWidthScin  = notchWidthScin;
      mod_notchDepth      = notchDepth;
      mod_sepDepth        = sepDepth;
      mod_visStr          = visStr;
      mod_regStr          = regStr;
      mod_limStr          = limStr;
  }
  double      mod_BIwidth;
  double      mod_BIheight;
  double      mod_SWThick;
  double      mod_TWThick;
  double      mod_FWThick;
  double      mod_BWThick;
  double      mod_width;
  double      mod_height;
  double      mod_notchWidthAbsA;
  double      mod_notchWidthAbsB;
  double      mod_notchWidthAbsC;
  double      mod_notchWidthScin;
  double      mod_notchDepth;
  double      mod_sepDepth;
  std::string mod_visStr;
  std::string mod_regStr;
  std::string mod_limStr;
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
  int         layer_ID;
  int         slice_ID;
  int         slice_partID;
  double      slice_thick;
  double      slice_offset;
  int         slice_readoutLayer;
  std::string slice_matStr;
  std::string slice_visStr;
  std::string slice_regStr;
  std::string slice_limStr;
};  

//************************************************************************************************************
//************************** Assembly for absorber plates for 8M modules *************************************
//************************************************************************************************************
Volume createAbsorberPlateEightM(Detector& desc,
                                   std::string basename, 
                                   double h_mod,
                                   double w_mod,
                                   double t_mod_tp,
                                   double t_mod_sp,
                                   double t_slice,
                                   double h_notch,
                                   double w_notchA,
                                   double w_notchB,
                                   double w_notchC,
                                   Material slice_mat,
                                   std::string region,
                                   std::string limit,
                                   std::string vis,
                                   bool renderComp
){
  
  double w_plate  = (w_mod/2-t_mod_sp)*2;
  double l_A = -w_plate/2;
  double l_B = -(w_plate/2-w_notchC);
  double l_C = -(w_plate/4+w_notchA/2);
  double l_D = -(w_plate/4-w_notchA/2);
  double l_E = -w_notchB/2;
  double r_A = w_plate/2;
  double r_B = w_plate/2-w_notchC;
  double r_C = w_plate/4+w_notchA/2;
  double r_D = w_plate/4-w_notchA/2;
  double r_E = w_notchB/2; 
                                        // 0      1     2     3      4
  const std::vector<double> xCoord      = { l_A,  l_B,  l_B,  l_C,   l_C,
                                        // 5      6     7     8      9 
                                            l_D,  l_D,  l_E,  l_E,   r_E,
                                        // 10     11    12    13     14
                                            r_E,  r_D,  r_D,  r_C,   r_C,
                                        // 15     16    17    18     19
                                            r_B,  r_B,  r_A,  r_A,   r_B,
                                        // 20     21    22    23     24
                                            r_B,  r_C,  r_C,  r_D,   r_D,
                                        // 25    26     27    28     29
                                            r_E,  r_E,  l_E,  l_E,   l_D,
                                        // 30    31     32    33     34
                                            l_D, l_C,   l_C,  l_B,   l_B,
                                        //35
                                            l_A
                                        };
                                        
  double topA = h_mod/2-t_mod_tp;
  double topB = h_mod/2-t_mod_tp-h_notch;
  double botA = -(h_mod/2-t_mod_tp);
  double botB = -(h_mod/2-t_mod_tp-h_notch);
                                        // 0       1       2      3       4
  const std::vector<double> yCoord      = { topA,  topA,   topB,  topB,   topA,
                                        // 5       6       7      8       9 
                                          topA,   topB,   topB,   topA,   topA, 
                                        // 10      11     12      13      14
                                          topB,   topB,   topA,   topA,   topB,
                                        // 15     16      17      18      19
                                          topB,   topA,   topA,   botA,   botA,
                                        // 20     21      22      23      24
                                          botB,   botB,   botA,   botA,   botB,
                                         // 25    26      27      28      29
                                          botB,   botA,   botA,   botB,   botB,
                                         // 30    31      32      33      34
                                          botA,   botA,   botB,   botB,   botA,
                                         //35
                                          botA,
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
    absplate_vol.setRegion(desc,region);
    absplate_vol.setLimitSet(desc,limit);
  }

  
  return absplate_vol;
  
}

//************************************************************************************************************
//************************** Assembly for absorber plates for 4M modules *************************************
//************************************************************************************************************
Volume createAbsorberPlateFourM(Detector& desc,
                                   std::string basename, 
                                   double h_mod,
                                   double w_mod,
                                   double t_mod_tp,
                                   double t_mod_sp,
                                   double t_slice,
                                   double h_notch,
                                   double w_notchA,
                                   double w_notchC,
                                   Material slice_mat,
                                   std::string region,
                                   std::string limit,
                                   std::string vis,
                                   bool renderComp
){
                                          // 0                 1                         2                          3            4
  const std::vector<double> xCoord      = { -(w_mod/2-t_mod_sp), -(w_mod/2-t_mod_sp-w_notchC),  -(w_mod/2-t_mod_sp-w_notchC),  -w_notchA/2,  -w_notchA/2,  
    
    
                                        // 5                  6                         7                            8                              9 
                                            w_notchA/2,    w_notchA/2,                w_mod/2-t_mod_sp-w_notchC,  w_mod/2-t_mod_sp-w_notchC,        w_mod/2-t_mod_sp, 
                                        // 10                   11                          12                              13           14
                                            w_mod/2-t_mod_sp, w_mod/2-t_mod_sp-w_notchC,    w_mod/2-t_mod_sp-w_notchC, w_notchA/2, w_notchA/2,    
                                        // 15          16          17                         18                            19
                                         -w_notchA/2, -w_notchA/2,  -(w_mod/2-t_mod_sp-w_notchC), -(w_mod/2-t_mod_sp-w_notchC),-(w_mod/2-t_mod_sp)
                                        };
                                        // 0                 1                         2                      3                           4
  const std::vector<double> yCoord      = {h_mod/2-t_mod_tp, h_mod/2-t_mod_tp,  h_mod/2-t_mod_tp-h_notch,    h_mod/2-t_mod_tp-h_notch,   h_mod/2-t_mod_tp, 
                                        // 5                  6                         7                            8                    9 
                                          h_mod/2-t_mod_tp,  h_mod/2-t_mod_tp-h_notch, h_mod/2-t_mod_tp-h_notch,    h_mod/2-t_mod_tp,    h_mod/2-t_mod_tp,
                                        // 10                   11                      12                              13                          14
                                          -(h_mod/2-t_mod_tp),  -(h_mod/2-t_mod_tp),    -(h_mod/2-t_mod_tp-h_notch),    -(h_mod/2-t_mod_tp-h_notch), -(h_mod/2-t_mod_tp),
                                        // 15                   16                          17                              18                        19
                                         -(h_mod/2-t_mod_tp),   -(h_mod/2-t_mod_tp-h_notch), -(h_mod/2-t_mod_tp-h_notch),   -(h_mod/2-t_mod_tp),    -(h_mod/2-t_mod_tp)
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
    absplate_vol.setRegion(desc,region);
    absplate_vol.setLimitSet(desc,limit);
  }
  
  return absplate_vol;
}


//************************************************************************************************************
//************************** Filler plate i.e. Tyvek/Air for 8M module ***************************************
//************************************************************************************************************
Volume createFillerPlateEightM( Detector& desc,
                                   std::string basename, 
                                   double h_mod,
                                   double w_mod,
                                   double t_mod_tp,
                                   double t_mod_sp,
                                   double t_slice,
                                   double h_notch,
                                   double w_notch,
                                   Material slice_mat,
                                   std::string region,
                                   std::string limit,
                                   std::string vis,
                                   bool renderComp
){
  double w_plate     = w_mod-2*t_mod_sp;
  double h_plate     = h_mod-2*t_mod_tp;
                                            // 0                 1                         2                          3                       4
  const std::vector<double> xCoord      = { -w_plate/2,   -(w_plate/4+w_notch/2),       -(w_plate/4+w_notch/2),      -(w_plate/4-w_notch/2), -(w_plate/4-w_notch/2), 
                                        // 5                    6                         7                            8                      9 
                                          w_plate/4-w_notch/2,  w_plate/4-w_notch/2,    w_plate/4+w_notch/2,          w_plate/4+w_notch/2,    w_plate/2,
                                        // 10                   11                    12                        13                  14            
                                        w_plate/2,          w_plate/4+w_notch/2,      w_plate/4+w_notch/2,   w_plate/4-w_notch/2,  w_plate/4-w_notch/2,
                                        // 15                   16                        17                      18                      19
                                        -(w_plate/4-w_notch/2), -(w_plate/4-w_notch/2),   -(w_plate/4+w_notch/2), -(w_plate/4+w_notch/2), -w_plate/2
                                        };  
                                        // 0                     1                  2              3             4
  const std::vector<double> yCoord      = { 
                                        h_plate/2-h_notch,      h_plate/2-h_notch,  h_plate/2,   h_plate/2,     h_plate/2-h_notch,  
                                        // 5                     6                 7                8                    9 
                                         h_plate/2-h_notch,      h_plate/2,        h_plate/2,    h_plate/2-h_notch,    h_plate/2-h_notch,
                                        // 10                   11                    12                        13                  14            
                                        -(h_plate/2-h_notch),   -(h_plate/2-h_notch),  -h_plate/2,         -h_plate/2,  -(h_plate/2-h_notch),
                                        // 15                   16                   17              18                      19
                                         -(h_plate/2-h_notch),  -h_plate/2,      -h_plate/2,      -(h_plate/2-h_notch),        -(h_plate/2-h_notch)
                                        };
  
  const std::vector<double> zStep       = {-t_slice/2, t_slice/2};
  const std::vector<double> zStepX      = {0., 0.};
  const std::vector<double> zStepY      = {0., 0.};
  const std::vector<double> zStepScale  = {1., 1.};
  
  ExtrudedPolygon filler = ExtrudedPolygon( xCoord, yCoord, zStep, zStepX, zStepY, zStepScale);
    
  Volume      filler_vol(basename, filler, slice_mat);
  // Setting slice attributes
  if (renderComp){
    filler_vol.setAttributes(desc, region, limit, vis);
  } else {
    filler_vol.setRegion(desc,region);
    filler_vol.setLimitSet(desc,limit);
  }
  
  return filler_vol;
}

//************************************************************************************************************
//************************** Filler plate i.e. Tyvek/Air for 4M module ***************************************
//************************************************************************************************************
Volume createFillerPlateFourM( Detector& desc,
                                   std::string basename, 
                                   double h_mod,
                                   double w_mod,
                                   double t_mod_tp,
                                   double t_mod_sp,
                                   double t_slice,
                                   double h_notch,
                                   double w_notch,
                                   Material slice_mat,
                                   std::string region,
                                   std::string limit,
                                   std::string vis,
                                   bool renderComp
){
                                              // 0                 1                         2                          3            4
  const std::vector<double> xCoord      = { -(w_mod/2-t_mod_sp),   -w_notch/2,  -w_notch/2,  w_notch/2,    w_notch/2,      
    
    
                                        // 5                  6                         7                            8                              9 
                                        w_mod/2-t_mod_sp,   w_mod/2-t_mod_sp,  w_notch/2, w_notch/2,    -w_notch/2, 
                                        // 10                   11             
                                        -w_notch/2,        -(w_mod/2-t_mod_sp)
                                        };
                                        // 0                           1                         2                   3                4
  const std::vector<double> yCoord      = {h_mod/2-t_mod_tp-h_notch,  h_mod/2-t_mod_tp-h_notch,   h_mod/2-t_mod_tp, h_mod/2-t_mod_tp, h_mod/2-t_mod_tp-h_notch,
                                        // 5                          6                            7                            8                          9 
                                            h_mod/2-t_mod_tp-h_notch, -(h_mod/2-t_mod_tp-h_notch), -(h_mod/2-t_mod_tp-h_notch), -(h_mod/2-t_mod_tp),  -(h_mod/2-t_mod_tp), 
                                        // 10                   11                     
                                        -(h_mod/2-t_mod_tp-h_notch),    -(h_mod/2-t_mod_tp-h_notch)
                                        };
  
  const std::vector<double> zStep       = {-t_slice/2, t_slice/2};
  const std::vector<double> zStepX      = {0., 0.};
  const std::vector<double> zStepY      = {0., 0.};
  const std::vector<double> zStepScale  = {1., 1.};
  
  ExtrudedPolygon filler = ExtrudedPolygon( xCoord, yCoord, zStep, zStepX, zStepY, zStepScale);
    
  Volume      filler_vol(basename, filler, slice_mat);
  // Setting slice attributes
  if (renderComp){
    filler_vol.setAttributes(desc, region, limit, vis);
  } else {
    filler_vol.setRegion(desc,region);
    filler_vol.setLimitSet(desc,limit);
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
                                   double h_notch,
                                   double w_notch,
                                   Material slice_mat,
                                   std::string region,
                                   std::string limit,
                                   std::string vis,
                                  SensitiveDetector sens,
                                  bool renderComp
){
                                        // 0                1,              2,                3,          4
  const std::vector<double> xCoord      = { -(w_tow/2),   w_tow/2-w_notch,  w_tow/2-w_notch,  w_tow/2,    w_tow/2,      
                                        // 5                 
                                       -(w_tow/2)
                                        };
                                        // 0                 1                  2        3        4
  const std::vector<double> yCoord      = {h_tow/2,  h_tow/2,   h_tow/2+h_notch, h_tow/2+h_notch, -h_tow/2,
                                        // 5                         
                                          -h_tow/2,
                                        };
  
  const std::vector<double> zStep       = {-t_slice/2, t_slice/2};
  const std::vector<double> zStepX      = {0., 0.};
  const std::vector<double> zStepY      = {0., 0.};
  const std::vector<double> zStepScale  = {1., 1.};
  
  ExtrudedPolygon scintplate = ExtrudedPolygon( xCoord, yCoord, zStep, zStepX, zStepY, zStepScale);
    
  Volume      slice_vol(basename, scintplate, slice_mat);
    // Setting appropriate slices as sensitive
  sens.setType("calorimeter");
  slice_vol.setSensitiveDetector(sens);  
  // Setting slice attributes
  if (renderComp){
    slice_vol.setAttributes(desc, region, limit, vis);
  } else {
    slice_vol.setRegion(desc,region);
    slice_vol.setLimitSet(desc,limit);
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
                                        double h_notch,
                                        double w_notch,
                                        double t_sep,
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
  double w_tow   = (w_mod-2*t_mod_sp-3*t_sep)/4;
  double h_tow  = (h_mod-2*t_mod_tp-t_sep-2*h_notch)/2;

  // placement volumes
  PlacedVolume pvm;

  // titanium-dioxide separations                                        // 0                        1                     2                    3                     4
  const std::vector<double> xCoordTi      = {-(w_mod/2-t_mod_sp),  -(w_tow+1.5*t_sep),  -(w_tow+1.5*t_sep),    -(w_tow+0.5*t_sep),  -(w_tow+0.5*t_sep),
                                        // 5          6          7            8           9
                                        -t_sep/2,    -t_sep/2,   t_sep/2,    t_sep/2,  w_tow+0.5*t_sep,
                                       // 10              11               12               13                  14
                                       w_tow+0.5*t_sep,   w_tow+1.5*t_sep, w_tow+1.5*t_sep, w_mod/2-t_mod_sp,   w_mod/2-t_mod_sp,
                                        // 15             16                  17              18             19
                                        w_tow+1.5*t_sep, w_tow+1.5*t_sep, w_tow+0.5*t_sep,  w_tow+0.5*t_sep,   t_sep/2, 
                                        // 20             21              22              23            24
                                          t_sep/2,         -t_sep/2,     -t_sep/2,   -(w_tow+0.5*t_sep),  -(w_tow+0.5*t_sep),
                                        // 25               26              27
                                        -(w_tow+1.5*t_sep), -(w_tow+1.5*t_sep), -(w_mod/2-t_mod_sp)
                                        };
                                        // 0           1            2                         3                     4
  const std::vector<double> yCoordTi      = { t_sep/2.,  t_sep/2.,   (h_tow+h_notch+t_sep/2), (h_tow+h_notch+t_sep/2), t_sep/2.,
                                        // 5            6                    7                         8            9                         
                                          t_sep/2.,   (h_tow+t_sep/2),      (h_tow+t_sep/2),   t_sep/2.,    t_sep/2.,
                                        // 10              11                               12              13            14
                                         (h_tow+h_notch+t_sep/2), (h_tow+h_notch+t_sep/2),   t_sep/2.,    t_sep/2.,    -t_sep/2.,
                                         // 15         16                        17                       18            19
                                        -t_sep/2.,   -(h_tow+h_notch+t_sep/2), -(h_tow+h_notch+t_sep/2),  -t_sep/2.,    -t_sep/2., 
                                        // 20               21                       22              23            24
                                        -(h_tow+t_sep/2), -(h_tow+t_sep/2), -t_sep/2.,    -t_sep/2.,    -(h_tow+h_notch+t_sep/2),
                                        // 25             26              27
                                        -(h_tow+h_notch+t_sep/2), -t_sep/2.,    -t_sep/2.
                                        };
                                        
  const std::vector<double> zStepTi       = {-t_slice/2, t_slice/2};
  const std::vector<double> zStepXTi      = {0., 0.};
  const std::vector<double> zStepYTi      = {0., 0.};
  const std::vector<double> zStepScaleTi  = {1., 1.};
  
  ExtrudedPolygon tiOgrid = ExtrudedPolygon( xCoordTi, yCoordTi, zStepTi, zStepXTi, zStepYTi, zStepScaleTi);
    
  Volume      ti0grid_vol(basename+"_Ti02Epoxy_"+_toString(layerID, "_layer_%d"), tiOgrid, slice_mat);
  // Setting slice attributes
  if (renderComp){
    ti0grid_vol.setAttributes(desc, region, limit, "LFHCALLayerTiOVis");
  } else {
    ti0grid_vol.setRegion(desc,region);
    ti0grid_vol.setLimitSet(desc,limit);
  }
  pvm = modScintAssembly.placeVolume(ti0grid_vol, Position(0, 0, 0 ));
  
  // 8M module placement of scintillator for tower
  double rotZ[8] = {0,    0, 0,     0,  0,    0,    0,    0};
  double rotY[8] = {M_PI, 0, M_PI,  0,  M_PI, 0,    M_PI, 0};
  double rotX[8] = {0,    0, 0,     0,  M_PI, M_PI, M_PI, M_PI};
  double posX[8] = {(w_tow*1.5+1.5*t_sep),     (w_tow*0.5+0.5*t_sep),     -(w_tow*0.5+0.5*t_sep),    -(w_tow*1.5+1.5*t_sep),
                    (w_tow*1.5+1.5*t_sep),     (w_tow*0.5+0.5*t_sep),     -(w_tow*0.5+0.5*t_sep),    -(w_tow*1.5+1.5*t_sep)};
  double posY[8] = {0.5*(h_tow)+0.5*t_sep,    0.5*(h_tow)+0.5*t_sep,    0.5*(h_tow)+0.5*t_sep,    0.5*(h_tow)+0.5*t_sep,
                    -(0.5*(h_tow)+0.5*t_sep), -(0.5*(h_tow)+0.5*t_sep), -(0.5*(h_tow)+0.5*t_sep), -(0.5*(h_tow)+0.5*t_sep)};
  double posZ[8] = {0,                                        0,                                        0,                                        0,
                    0,                                        0,                                        0,                                        0};
  int towerx  = 0;
  int towery  = 0;

  // loop over all towers within same module
  for (int i = 0; i < 8; i++){
//     std::cout << basename << _toString(i, "_tower_%d") << "\t"<< modID << "\t" << i << "\t"<< layerID << std::endl;
    Volume modScintTowerAss = createScintillatorTower( desc,  basename+ _toString(i, "_tower_%d"),  
                                                            w_tow, h_tow, t_slice,
                                                            h_notch, (w_notch-t_sep)/2, 
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
                                        double h_notch,
                                        double w_notch,
                                        double t_sep,
                                        Material slice_mat,
                                        int roLayer,
                                        std::string region,
                                        std::string limit,
                                        std::string vis,
                                        SensitiveDetector sens,
                                        bool renderComp 
){
  // Tower placement in 4M module
  //==================================
  //||              ||              ||
  //||      0       ||      1       ||
  //||              ||              ||
  //==================================
  //||              ||              ||
  //||      2       ||      3       ||
  //||              ||              ||
  //==================================
  Assembly modScintAssembly(basename);
  double w_tow   = (w_mod-2*t_mod_sp-t_sep)/2;
  double h_tow  = (h_mod-2*t_mod_tp-t_sep-2*h_notch)/2;

  // placement volumes
  PlacedVolume pvm;
  
  // titanium-dioxide separations
  const std::vector<double> xCoordTi      = {-(w_mod/2-t_mod_sp),  -(t_sep/2.),        -(t_sep/2.),      t_sep/2.,  t_sep/2.,
                                        // 5                6                7                8           9
                                        w_mod/2-t_mod_sp, w_mod/2-t_mod_sp,  t_sep/2.,         t_sep/2.,   -(t_sep/2.),
                                       // 10                11
                                        -(t_sep/2.),        -(w_mod/2-t_mod_sp)
                                        };
                                        // 0           1            2                         3                     4
  const std::vector<double> yCoordTi      = { t_sep/2.,  t_sep/2.,   (h_tow+h_notch+t_sep/2), (h_tow+h_notch+t_sep/2), t_sep/2.,
                                        // 5            6           7             8                        9                         
                                          t_sep/2.,    -t_sep/2.,   -t_sep/2.,    -(h_tow+h_notch+t_sep/2), -(h_tow+h_notch+t_sep/2),
                                        // 10               11
                                         -t_sep/2.,   -t_sep/2.,
                                        };
  
  const std::vector<double> zStepTi       = {-t_slice/2, t_slice/2};
  const std::vector<double> zStepXTi      = {0., 0.};
  const std::vector<double> zStepYTi      = {0., 0.};
  const std::vector<double> zStepScaleTi  = {1., 1.};
  
  ExtrudedPolygon tiOgrid = ExtrudedPolygon( xCoordTi, yCoordTi, zStepTi, zStepXTi, zStepYTi, zStepScaleTi);
    
  Volume      ti0grid_vol(basename+"_Ti02Epoxy_"+_toString(layerID, "_layer_%d"), tiOgrid, slice_mat);
  // Setting slice attributes
  if (renderComp){
    ti0grid_vol.setAttributes(desc, region, limit, "LFHCALLayerTiOVis");
  } else {
   ti0grid_vol.setRegion(desc,region);
   ti0grid_vol.setLimitSet(desc,limit);
  }
  pvm = modScintAssembly.placeVolume(ti0grid_vol, Position(0, 0, 0 ));
  
  // 4M module placement of scintillator for tower
  double rotZ[4] = {0,    0, 0,     0   };
  double rotY[4] = {M_PI, 0, M_PI,  0   };
  double rotX[4] = {0,    0, M_PI,  M_PI};
  double posX[4] = {(w_tow*0.5+0.5*t_sep),     -(w_tow*0.5+0.5*t_sep),  
                    (w_tow*0.5+0.5*t_sep),     -(w_tow*0.5+0.5*t_sep)};
  double posY[4] = {0.5*(h_tow)+0.5*t_sep,    0.5*(h_tow)+0.5*t_sep,
                    -(0.5*(h_tow)+0.5*t_sep), -(0.5*(h_tow)+0.5*t_sep)};
  double posZ[4] = {0,                              0,                  
                    0,                              0};
  int towerx  = 0;
  int towery  = 0;
  // loop over all towers within same module
                    
  for (int i = 0; i < 4; i++){
//     std::cout << basename << _toString(i, "_tower_%d") << "\t"<< modID << "\t" << i << "\t"<< layerID << std::endl;
    Volume modScintTowerAss = createScintillatorTower( desc,  basename+ _toString(i, "_tower_%d"),  
                                                            w_tow, h_tow, t_slice,
                                                            h_notch, (w_notch-t_sep)/2., 
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
  vol_mod.setVisAttributes(mod_params.mod_visStr);

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
    vol_modFrontPlate.setRegion(desc,mod_params.mod_regStr);
    vol_modFrontPlate.setLimitSet(desc,mod_params.mod_limStr);
    vol_modBackPlate.setRegion(desc,mod_params.mod_regStr);
    vol_modBackPlate.setLimitSet(desc,mod_params.mod_limStr);
    vol_modSidePlateL.setRegion(desc,mod_params.mod_regStr);
    vol_modSidePlateL.setLimitSet(desc,mod_params.mod_limStr);
    vol_modSidePlateR.setRegion(desc,mod_params.mod_regStr);
    vol_modSidePlateR.setLimitSet(desc,mod_params.mod_limStr);
    vol_modTopPlate.setRegion(desc,mod_params.mod_regStr);
    vol_modTopPlate.setLimitSet(desc,mod_params.mod_limStr);
    vol_modBottomPlate.setRegion(desc,mod_params.mod_regStr);
    vol_modBottomPlate.setLimitSet(desc,mod_params.mod_limStr);
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
      Volume modAbsAssembly = createAbsorberPlateEightM( desc, 
                                                          baseName+"_Abs"+_toString(sl_params[i].layer_ID, "_layer_%d"),
                                                          mod_params.mod_height, mod_params.mod_width, mod_params.mod_TWThick, mod_params.mod_SWThick,
                                                          sl_params[i].slice_thick, mod_params.mod_notchDepth, mod_params.mod_notchWidthAbsA, mod_params.mod_notchWidthAbsB, mod_params.mod_notchWidthAbsC,
                                                          slice_mat, sl_params[i].slice_regStr, sl_params[i].slice_limStr, sl_params[i].slice_visStr, renderComp);
      // Placing slice within layer
      if (allSen) modAbsAssembly.setSensitiveDetector(sens);  
      pvm = vol_mod.placeVolume(modAbsAssembly, Transform3D(RotationZYX(0, 0, 0), Position(0., 0., slice_z)));
      if (allSen) pvm.addPhysVolID("towerx", 0).addPhysVolID("towery", 0).addPhysVolID("rlayerz", sl_params[i].slice_readoutLayer).addPhysVolID("layerz", layer_num).addPhysVolID("passive", 1);
    //*************************************************  
    // air & tyvek
    //*************************************************  
    } else if (sl_params[i].slice_partID == 2 ){
      Volume modFillAssembly =  createFillerPlateEightM( desc, 
                                                          baseName+"_Fill"+_toString(sl_params[i].layer_ID, "_layer_%d")+_toString(sl_params[i].slice_ID, "slice_%d"),
                                                          mod_params.mod_height, mod_params.mod_width, mod_params.mod_TWThick, mod_params.mod_SWThick,
                                                          sl_params[i].slice_thick, mod_params.mod_notchDepth, mod_params.mod_notchWidthScin,
                                                          slice_mat, sl_params[i].slice_regStr, sl_params[i].slice_limStr, sl_params[i].slice_visStr, renderComp);
      // Placing slice within layer
      if (allSen) modFillAssembly.setSensitiveDetector(sens);
      pvm = vol_mod.placeVolume(modFillAssembly, Transform3D(RotationZYX(0, 0, 0), Position(0., 0., slice_z)));
      if (allSen) pvm.addPhysVolID("towerx", 1).addPhysVolID("towery", 0).addPhysVolID("rlayerz", sl_params[i].slice_readoutLayer).addPhysVolID("layerz", layer_num).addPhysVolID("passive", 1);
    //*************************************************
    // scintillator
    //*************************************************
    } else {
      Assembly modScintAssembly =  createScintillatorPlateEightM(  desc,
                                                                  baseName+"_ScintAssembly"+_toString(sl_params[i].layer_ID, "_layer_%d"),
//                                                                   modID, 
                                                                  layer_num, 
                                                                  mod_params.mod_height, mod_params.mod_width, mod_params.mod_TWThick, mod_params.mod_SWThick,
                                                                  sl_params[i].slice_thick, mod_params.mod_notchDepth, mod_params.mod_notchWidthScin, mod_params.mod_sepDepth, 
                                                                  slice_mat, sl_params[i].slice_readoutLayer ,sl_params[i].slice_regStr, sl_params[i].slice_limStr, sl_params[i].slice_visStr, sens, renderComp);
      // Placing slice within layer
      pvm = vol_mod.placeVolume(modScintAssembly, Transform3D(RotationZYX(0, 0, 0), Position(0., 0., slice_z)));
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
  vol_mod.setVisAttributes(mod_params.mod_visStr);

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
    vol_modFrontPlate.setRegion(desc,mod_params.mod_regStr);
    vol_modFrontPlate.setLimitSet(desc,mod_params.mod_limStr);
    vol_modBackPlate.setRegion(desc,mod_params.mod_regStr);
    vol_modBackPlate.setLimitSet(desc,mod_params.mod_limStr);
    vol_modSidePlateL.setRegion(desc,mod_params.mod_regStr);
    vol_modSidePlateL.setLimitSet(desc,mod_params.mod_limStr);
    vol_modSidePlateR.setRegion(desc,mod_params.mod_regStr);
    vol_modSidePlateR.setLimitSet(desc,mod_params.mod_limStr);
    vol_modTopPlate.setRegion(desc,mod_params.mod_regStr);
    vol_modTopPlate.setLimitSet(desc,mod_params.mod_limStr);
    vol_modBottomPlate.setRegion(desc,mod_params.mod_regStr);
    vol_modBottomPlate.setLimitSet(desc,mod_params.mod_limStr);
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
      Volume modAbsAssembly = createAbsorberPlateFourM( desc, baseName+"_Abs"+_toString(sl_params[i].layer_ID, "_layer_%d"),
                                                              mod_params.mod_height, mod_params.mod_width, mod_params.mod_TWThick, mod_params.mod_SWThick,
                                                              sl_params[i].slice_thick, mod_params.mod_notchDepth, mod_params.mod_notchWidthAbsA, mod_params.mod_notchWidthAbsC,
                                                              slice_mat, sl_params[i].slice_regStr, sl_params[i].slice_limStr, sl_params[i].slice_visStr, renderComp);
      // Placing slice within layer
      if (allSen) modAbsAssembly.setSensitiveDetector(sens);  
      pvm = vol_mod.placeVolume(modAbsAssembly, Transform3D(RotationZYX(0, 0, 0), Position(0., 0., slice_z)));
      if (allSen) pvm.addPhysVolID("towerx", 0).addPhysVolID("towery", 0).addPhysVolID("rlayerz", sl_params[i].slice_readoutLayer).addPhysVolID("layerz", layer_num).addPhysVolID("passive", 1);
    //*************************************************  
    // air & tyvek
    //*************************************************  
    } else if (sl_params[i].slice_partID == 2 ){
      Volume modFillAssembly =  createFillerPlateFourM( desc, baseName+"_Fill"+_toString(sl_params[i].layer_ID, "_layer_%d")+_toString(sl_params[i].slice_ID, "slice_%d"),
                                                              mod_params.mod_height, mod_params.mod_width, mod_params.mod_TWThick, mod_params.mod_SWThick,
                                                              sl_params[i].slice_thick, mod_params.mod_notchDepth, mod_params.mod_notchWidthScin,
                                                              slice_mat, sl_params[i].slice_regStr, sl_params[i].slice_limStr, sl_params[i].slice_visStr, renderComp);
      // Placing slice within layer
      if (allSen) modFillAssembly.setSensitiveDetector(sens);  
      pvm = vol_mod.placeVolume(modFillAssembly, Transform3D(RotationZYX(0, 0, 0), Position(0., 0., slice_z)));
      if (allSen) pvm.addPhysVolID("towerx", 1).addPhysVolID("towery", 0).addPhysVolID("rlayerz", sl_params[i].slice_readoutLayer).addPhysVolID("layerz", layer_num).addPhysVolID("passive", 1);
    //*************************************************
    // scintillator
    //*************************************************
    } else {
      Assembly modScintAssembly =  createScintillatorPlateFourM( desc,baseName+"_ScintAssembly"+_toString(sl_params[i].layer_ID, "_layer_%d"),
//                                                                   modID,  
                                                                 layer_num, 
                                                                  mod_params.mod_height, mod_params.mod_width, mod_params.mod_TWThick, mod_params.mod_SWThick,
                                                                  sl_params[i].slice_thick, mod_params.mod_notchDepth, mod_params.mod_notchWidthScin, mod_params.mod_sepDepth, 
                                                                  slice_mat, sl_params[i].slice_readoutLayer, sl_params[i].slice_regStr, sl_params[i].slice_limStr, sl_params[i].slice_visStr, sens, renderComp);
      // Placing slice within layer
      pvm = vol_mod.placeVolume(modScintAssembly, Transform3D(RotationZYX(0, 0, 0), Position(0., 0., slice_z)));
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
  xml_dim_t pos = detElem.position();

  std::cout<< "global LFHCal position" << pos.x() << "\t" << pos.y() << "\t" << pos.z()  << std::endl;
  
  
  bool renderComponents = getAttrOrDefault(detElem, _Unicode(renderComponents), 0.);
  bool allSensitive     = getAttrOrDefault(detElem, _Unicode(allSensitive), 0.);
  if (renderComponents) {
    std::cout << "enabled visualization" << std::endl;
  } else {
    std::cout << "switchted off visualization" << std::endl;
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
                                  getAttrOrDefault(eightMmod_dim, _Unicode(notchWidthAbsA), 0.),
                                  getAttrOrDefault(eightMmod_dim, _Unicode(notchWidthAbsB), 0.),
                                  getAttrOrDefault(eightMmod_dim, _Unicode(notchWidthAbsC), 0.),
                                  getAttrOrDefault(eightMmod_dim, _Unicode(notchWidthScin), 0.),
                                  getAttrOrDefault(eightMmod_dim, _Unicode(notchDepth), 0.),
                                  getAttrOrDefault(eightMmod_dim, _Unicode(sepDepth), 0.), 
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
                                  getAttrOrDefault(fourMmod_dim, _Unicode(notchWidthAbsA), 0.),
                                  getAttrOrDefault(fourMmod_dim, _Unicode(notchWidthAbsB), 0.),
                                  getAttrOrDefault(fourMmod_dim, _Unicode(notchWidthAbsC), 0.),
                                  getAttrOrDefault(fourMmod_dim, _Unicode(notchWidthScin), 0.),
                                  getAttrOrDefault(fourMmod_dim, _Unicode(notchDepth), 0.),
                                  getAttrOrDefault(fourMmod_dim, _Unicode(sepDepth), 0.), 
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

  
  int    moduleID   = 0;
  std::vector<double> xpos8M;
  std::vector<double> ypos8M;
  std::vector<double> zpos8M;
  
  for(xml_coll_t i(handle, _Unicode(eightmodulepositions)); i; ++i){
    xml_comp_t  x_mtrx = i;

    std::string   mtrx_name       = getAttrOrDefault<std::string>(x_mtrx, _Unicode(name), " ");
    std::string   mtrx_values     = getAttrOrDefault<std::string>(x_mtrx, _Unicode(values), " ");

    std::vector<double> *aptr = NULL;

    if(mtrx_name == "xpos")
      aptr = &xpos8M;
    else if(mtrx_name == "ypos")
      aptr = &ypos8M;
    else if(mtrx_name == "zpos")
      aptr = &zpos8M;
    else{
      printout(WARNING, "LFHCAL", "unknown <eightmodulepositions> data!");
      continue;
    }

    std::string delimiter = " ";
    size_t posC = 0;
    std::string token;
    while ((posC = mtrx_values.find(delimiter)) != std::string::npos) {
      token = mtrx_values.substr(0, posC);
      aptr->push_back(atof(token.c_str()));
      mtrx_values.erase(0, posC + delimiter.length());
    }
    aptr->push_back(atof(mtrx_values.c_str()));
  }

  // create 8M modules
  Volume  eightMassembly = createEightMModule ( desc, eightM_params, slice_Params, length, sens, renderComponents, allSensitive);
  if (xpos8M.size() != ypos8M.size() || xpos8M.size() != zpos8M.size()){
    std::cout << xpos8M.size() << "\t" <<ypos8M.size() <<  "\t" << zpos8M.size() << std::endl;
    std::cout <<  "idiot you can't count" << std::endl;
  } else {
    for (int e = 0; e < (int)xpos8M.size(); e++){
      if(e%20 == 0 ) std::cout <<  "LFHCAL placing 8M module: " << e << "/"<< (int)xpos8M.size() << "\t" << xpos8M[e] << "\t" << ypos8M[e] << "\t" << zpos8M[e]<< std::endl;
        if(moduleIDx<0 || moduleIDy<0){
        std::cout << "LFHCAL WRONG ID FOR 8M module: " << e << "/" << (int)xpos8M.size() << "\t" << moduleIDx << "\t"
                  << moduleIDy << std::endl;
      }
      moduleIDx             = ((xpos8M[e] + 270) / 10);
      moduleIDy             = ((ypos8M[e] + 265) / 10);
      
      // Placing modules in world volume
      auto tr8M = Transform3D(Position(pos.x()-xpos8M[e]*dd4hep::cm-0.5*eightM_params.mod_width, pos.y()-ypos8M[e]*dd4hep::cm, pos.z() +zpos8M[e]*dd4hep::cm + length / 2.));
      phv = assembly.placeVolume(eightMassembly, tr8M);
      phv.addPhysVolID("moduleIDx", moduleIDx).addPhysVolID("moduleIDy", moduleIDy).addPhysVolID("moduletype", 0);
      moduleID++;
    }
  }

  std::vector<double> xpos4M;
  std::vector<double> ypos4M;
  std::vector<double> zpos4M;

  for(xml_coll_t i(handle, _Unicode(fourmodulepositions)); i; ++i){
    xml_comp_t  x_mtrx = i;

    std::string   mtrx_name       = getAttrOrDefault<std::string>(x_mtrx, _Unicode(name), " ");
    std::string   mtrx_values     = getAttrOrDefault<std::string>(x_mtrx, _Unicode(values), " ");
    std::vector<double> *aptr = NULL;

    if(mtrx_name == "xpos")
      aptr = &xpos4M;
    else if(mtrx_name == "ypos")
      aptr = &ypos4M;
    else if(mtrx_name == "zpos")
      aptr = &zpos4M;
    else{
      printout(WARNING, "LFHCAL", "unknown <fourmodulepositions> data!");
      continue;
    }

    std::string delimiter = " ";
    size_t posC = 0;
    std::string token;
    while ((posC = mtrx_values.find(delimiter)) != std::string::npos) {
      token = mtrx_values.substr(0, posC);
      aptr->push_back(atof(token.c_str()));
      mtrx_values.erase(0, posC + delimiter.length());
    }
    aptr->push_back(atof(mtrx_values.c_str()));
  }

  // create 4M modules
  Volume  fourMassembly = createFourMModule ( desc, fourM_params, slice_Params,  length, sens, renderComponents, allSensitive);
  if (xpos4M.size() != ypos4M.size() || xpos4M.size() != zpos4M.size()){
    std::cout << xpos4M.size() << "\t" <<ypos4M.size() <<  "\t" << zpos4M.size() << std::endl;
    std::cout <<  "idiot you can't count" << std::endl;
  } else {
    for (int f = 0; f < (int)xpos4M.size(); f++){
      if(f%20 == 0 )std::cout <<  "LFHCAL placing 4M module: " << f << "/"<< (int)xpos4M.size() << "\t" << xpos4M[f] << "\t" << ypos4M[f] << "\t" << zpos4M[f]<< std::endl;
      
      moduleIDx             = ((xpos4M[f] + 265) / 10);
      moduleIDy             = ((ypos4M[f] + 265) / 10);
      if(moduleIDx<0 || moduleIDy<0){
        std::cout << "LFHCAL WRONG ID FOR 4M module: " << f << "/" << (int)xpos4M.size() << "\t" << moduleIDx << "\t"
                  << moduleIDy << std::endl;
      }
      auto tr4M = Transform3D(Position(pos.x()-xpos4M[f]*dd4hep::cm-0.5*fourM_params.mod_width, pos.y()-ypos4M[f]*dd4hep::cm, pos.z() +zpos4M[f]*dd4hep::cm + length / 2.));
      phv = assembly.placeVolume(fourMassembly, tr4M);
      
      phv.addPhysVolID("moduleIDx", moduleIDx).addPhysVolID("moduleIDy", moduleIDy).addPhysVolID("moduletype", 1);
  
      moduleID++;
    }
  }
  
  Volume     motherVol = desc.pickMotherVolume(det);
  Transform3D  tr        = Translation3D(0., 0., 0.) * RotationZYX(0.,0.,0.);
  PlacedVolume envPV     = motherVol.placeVolume(assembly, tr);
  envPV.addPhysVolID("system", detID);
  det.setPlacement(envPV);
  
  return det;
}
DECLARE_DETELEMENT(epic_LFHCAL, createDetector)
