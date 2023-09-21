#include "DD4hep/DetFactoryHelper.h"
#include "TVector3.h"
#include "XML/Layering.h"
using namespace std;
using namespace dd4hep;
using namespace dd4hep::detail;

const int nOuterRows       = 2;
const int nStarEmcRows     = 12;
const int nInnerRowsGroup1 = 4; // innermost // 4 rows
const int nInnerRowsGroup2 = 2; // central inner //2 rows
const int nInnerRowsGroup3 = 3; // star emc size continuation //3 rows
// const double innerHcalEtaBins[11] =
// {3.48733, 3.33696, 3.19285, 3.05481, 2.92264, 2.79614, 2.67511, 2.55935, 2.44868, 2.34289};
// 4 mostinner rows 15 deg tiles   +   2 rows x 10 deg tiles   +   4 rows x 6 deg tiles

const unsigned nHcalSectors                            = 12;
const double   starEmcEtaBins[nStarEmcRows + 1]        = {2.0,    1.9008, 1.8065, 1.7168, 1.6317, 1.5507, 1.4738,
                                                          1.4007, 1.3312, 1.2651, 1.2023, 1.1427, 1.086};
const double   outerHcalEtaBins[nOuterRows + 1]        = {1.38221, 1.3219, 1.26371};
const double   innerHcalEtaBins1[nInnerRowsGroup1 + 1] = {3.48733, 3.33696, 3.19285, 3.05481, 2.92264};
const double   innerHcalEtaBins2[nInnerRowsGroup2 + 1] = {2.92264, 2.79614, 2.67511};
const double   innerHcalEtaBins3[nInnerRowsGroup3 + 1] = {2.67511, 2.55935, 2.44868, 2.34289};

const double starEmcZPositions[24] = {270.19,  271.695, 273.15,  274.555, 275.96,  277.365, 282.363, 283.768,
                                      285.173, 286.578, 287.983, 289.388, 290.793, 292.198, 293.603, 295.008,
                                      296.413, 297.818, 299.223, 300.628, 302.033, 303.438, 304.843, 306.158};
const int    tileMap[11]           = {0, 2, 4, 5, 6, 9, 11, 13, 16, 18, 20};
double       getR(const double& z, const double& eta) { return z / (sinh(eta)); }

void buildSubsector(Detector& description, xml_h e, SensitiveDetector sens, vector<PlacedVolume>& sensitives,
                    const int& layerNumber, const double& globalZ, const xml_comp_t x_slice, Volume& sectorVol,
                    const int& nSubSectors, const int& nRows, const double* etaBins, const int& sectornum,
                    const double& oldEmcGap = 0, bool isStarEmc = false)
{
  xml_det_t x_det = e;
  Material  air   = description.air();
  Layering  layering(x_det);
  double    slice_thickness = x_slice.thickness();

  double dphiSector   = 2. * M_PI / (nHcalSectors); // 30  deg
  double dphi         = dphiSector / nSubSectors;
  int    oldTileLayer = tileMap[layerNumber - 1];        // layer of the old tile corresponding to the current layer
  double zOldPos      = starEmcZPositions[oldTileLayer]; // center z position of the old tile
  double zPos         = globalZ;
  double rmax         = isStarEmc ? getR(zOldPos, starEmcEtaBins[nRows])
                                  : getR(zPos, etaBins[nRows]) - oldEmcGap; // closest tile distance to the beam pipe
  double rmin         = isStarEmc ? getR(zOldPos, starEmcEtaBins[0])
                                  : getR(zPos, etaBins[0]) - oldEmcGap; // furthest tile distance to the beam pipe
  if (oldEmcGap <= 0 || nRows == nInnerRowsGroup3) // tricky check if this is not the innermost group or second
                                                   // innermost group
    rmax *= (2 - cos(dphi / 2));
  // to account for difference between trapezoid and tube segment

  for (int iSubSector = 0; iSubSector < nSubSectors; iSubSector++) {
    int subSector_volID = iSubSector;
    if (nSubSectors == 5 && nRows == nStarEmcRows)
      subSector_volID = iSubSector + 10;
    if (nSubSectors == 5 && nRows == nInnerRowsGroup3)
      subSector_volID = iSubSector + 15;
    if (nSubSectors == 3 && nRows == nInnerRowsGroup2)
      subSector_volID = iSubSector + 20;
    if (nSubSectors == 2 && nRows == nInnerRowsGroup1)
      subSector_volID = iSubSector + 23;
    string subSectorName = Form("subsector_%i", subSector_volID);
    Volume subSectorVol(subSectorName, Tube(rmin, rmax, slice_thickness / 2, 0, dphi), air);

    for (int iRow = 0; iRow < nRows; iRow++) {
      // recreate tile size for the given layer:
      double rBottom = isStarEmc
                           ? getR(zOldPos, starEmcEtaBins[iRow])
                           : getR(zPos, etaBins[iRow]) - oldEmcGap; // current tile closest distance to the beam pipe
      double rTop    = isStarEmc
                           ? getR(zOldPos, starEmcEtaBins[iRow + 1])
                           : getR(zPos, etaBins[iRow + 1]) - oldEmcGap; // current tile closest distance to the beam pipe

      if (!isStarEmc && iRow == (nRows - 1) &&
          oldEmcGap > 0) // last row in 2 most inner groups - a temporary solution for avoiding overlaps
      {
        if (nRows == nInnerRowsGroup2) {
          rTop = rTop * cos(dphi / 2) * cos(dphi / 3);
        } else if (nRows == nInnerRowsGroup1) {
          rTop = rTop * 1.002 * cos(dphi / 2) * cos(dphi / 3);
        }
      }

      // cosine theorem for the angle between the two radii:
      //(2 dx1)^2 = r^2+r^2-2*r*r*cos(phi) = 2r^2(1-cos(phi));
      //   dx1=r*sqrt(2(1-cos(phi)))/2;
      double    dx1 = rBottom * sqrt(2 * (1 - cos(dphi))) / 2;
      double    dx2 = rTop * sqrt(2 * (1 - cos(dphi))) / 2;
      double    dz  = (rTop - rBottom) / 2;
      double    dy  = slice_thickness / 2;
      Trapezoid tile(dx1, dx2, dy, dy, dz);
      Material  s_mat      = description.material(x_slice.materialStr());
      int       tile_volID = iRow;
      if (nSubSectors == 5 && nRows == nStarEmcRows)
        tile_volID = iRow + 2;
      if (nSubSectors == 5 && nRows == nInnerRowsGroup3)
        tile_volID = iRow + 14;
      if (nSubSectors == 3 && nRows == nInnerRowsGroup2)
        tile_volID = iRow + 17;
      if (nSubSectors == 2 && nRows == nInnerRowsGroup1)
        tile_volID = iRow + 19;
      string tileName =
          isStarEmc ? _toString(100 * tile_volID + 10000 * subSector_volID + sectornum, "tileOld%i")
                    : Form("tile%i_%i_%i_%i", (int)(dphi / M_PI * 180 + 0.5), tile_volID, subSector_volID, sectornum);
      Volume tile_vol(tileName, tile, s_mat);
      tile_vol.setVisAttributes(description.visAttributes(x_slice.visStr()));
      double       distance  = (rTop + rBottom) / 2;
      PlacedVolume tile_phys = subSectorVol.placeVolume(
          tile_vol, Transform3D(RotationZYX(0, M_PI / 2 + dphi / 2, M_PI / 2),
                                Position(distance * cos(dphi / 2), distance * sin(dphi / 2), 0)));
      tile_phys.addPhysVolID("tile", tile_volID);
      sens.setType("calorimeter");
      tile_vol.setSensitiveDetector(sens);
      sensitives.push_back(tile_phys);
    }

    PlacedVolume subSectorVol_phys = sectorVol.placeVolume(subSectorVol, RotationZYX((iSubSector)*dphi, 0, 0));
    subSectorVol_phys.addPhysVolID("subsector", subSector_volID);
  }
}
static Ref_t create_detector(Detector& description, xml_h e, SensitiveDetector sens)
{
  xml_det_t x_det    = e;
  xml_dim_t dim      = x_det.dimensions();
  int       det_id   = x_det.id();
  bool      reflect  = x_det.reflect(true);
  string    det_name = x_det.nameStr();
  Material  air      = description.air();
  // int            numsides = dim.numsides();
  xml::Component pos  = x_det.position();
  double         rmin = dim.rmin();
  double         rmax = dim.rmax();
  double         zmin = dim.zmin();

  Layering layering(x_det);
  double   totalThickness = layering.totalThickness();

  Volume endcapVol("endcap", Tube(rmin, rmax + 2 * cm, totalThickness), air); // rmax + 2 * cm -to account for disk gap
  DetElement endcap("endcap", det_id);

  double layerZ = -totalThickness / 2;
  endcapVol.setAttributes(description, x_det.regionStr(), x_det.limitsStr(), x_det.visStr());

  int    iLayer  = 0;
  double globalZ = zmin;
  // Looping through all the different layer sections
  for (xml_coll_t xc(x_det, _U(layer)); xc; ++xc) {
    xml_comp_t x_layer         = xc;
    int        repeat          = x_layer.repeat();
    double     layer_thickness = x_layer.thickness();

    // Looping through the number of repeated Rows in each section
    for (int i = 0; i < repeat; i++) {
      iLayer++;
      string               l_name   = _toString(iLayer - 1, "layer%d");
      double               disksGap = getAttrOrDefault(x_layer, _Unicode(disksGap), 1 * cm);
      Volume               layer_volume("layer", Tube(rmin, rmax + disksGap, layer_thickness), air);
      vector<PlacedVolume> sensitives;
      int                  slice_number = 1;
      double               sliceZ       = -layer_thickness / 2;

      for (xml_coll_t xs(x_layer, _U(slice)); xs; ++xs) {
        xml_comp_t x_slice         = xs;
        string     slice_name      = _toString(slice_number, "slice%d");
        double     slice_thickness = x_slice.thickness();
        globalZ += slice_thickness / 2;

        Material s_mat = description.material(x_slice.materialStr());
        sliceZ += slice_thickness / 2;

        Volume slice_volume("slice", Tube(rmin, rmax + disksGap, slice_thickness / 2), air);

        if (x_slice.isSensitive()) {
          ConeSegment sectorShape(slice_thickness / 2, rmin, rmax, rmin, rmax, 0, 2 * M_PI / (nHcalSectors));
          ///////Volume      sectorVol{"sector", sectorShape, air};
          // void buildSubsector(Detector& description, xml_h e, SensitiveDetector sens, vector<PlacedVolume>&
          // sensitives,
          //     const double& slice_thickness, const int& layerNumber, const xml_comp_t x_slice, Volume& sectorVol,
          //     const int& nSubSectors, const int& nRows, const double* etaBins, bool isStarEmc = false)
          for (unsigned iSector = 0; iSector < nHcalSectors; iSector++) {
            string SectorName = Form("sector_%i", iSector);
            Volume sectorVol{SectorName, sectorShape, air};
            int    nSubsectors;
            nSubsectors = 10; // 3 degree tiles
            buildSubsector(description, e, sens, sensitives, iLayer, globalZ, x_slice, sectorVol, nSubsectors,
                           nOuterRows, outerHcalEtaBins, iSector, -2.18);
            nSubsectors = 5; // 6 degree tiles
            buildSubsector(description, e, sens, sensitives, iLayer, globalZ, x_slice, sectorVol, nSubsectors,
                           nStarEmcRows, starEmcEtaBins, iSector, 0., true);
            nSubsectors = 5; // 6 degree tiles
            buildSubsector(description, e, sens, sensitives, iLayer, globalZ, x_slice, sectorVol, nSubsectors,
                           nInnerRowsGroup3, innerHcalEtaBins3, iSector, 2.18);
            nSubsectors = 3; // 10 degree tiles
            buildSubsector(description, e, sens, sensitives, iLayer, globalZ, x_slice, sectorVol, nSubsectors,
                           nInnerRowsGroup2, innerHcalEtaBins2, iSector, 2.18);
            nSubsectors = 2; // 15 degree tiles
            buildSubsector(description, e, sens, sensitives, iLayer, globalZ, x_slice, sectorVol, nSubsectors,
                           nInnerRowsGroup1, innerHcalEtaBins1, iSector, 2.18);
            PlacedVolume sectorVol_phys = slice_volume.placeVolume(
                sectorVol, Transform3D(RotationZYX(M_PI / 2 + iSector * 2 * M_PI / (nHcalSectors), 0, 0),
                                       Position(((iSector + 1) > nHcalSectors / 2 ? disksGap : -disksGap), 0, 0)));
            sectorVol_phys.addPhysVolID("sector", iSector);
          }
          PlacedVolume s_phv = layer_volume.placeVolume(slice_volume, Position(0, 0, sliceZ));
          s_phv.addPhysVolID("slice", slice_number);
        } else {
          Volume halfdisk("halfdisk", Tube(rmin, rmax, slice_thickness / 2, M_PI / 2, M_PI * 3 / 2), s_mat);
          halfdisk.setVisAttributes(description.visAttributes(x_slice.visStr()));
          PlacedVolume s_phv1 = slice_volume.placeVolume(halfdisk, Position(-disksGap / 2, 0, 0));
          s_phv1.addPhysVolID("halfdisk", 0);
          PlacedVolume s_phv2 =
              slice_volume.placeVolume(halfdisk, Transform3D(RotationZYX(M_PI, 0, 0), Position(+disksGap / 2, 0, 0)));
          s_phv2.addPhysVolID("halfdisk", 1);
          PlacedVolume s_phv = layer_volume.placeVolume(slice_volume, Position(0, 0, sliceZ));
          s_phv.addPhysVolID("slice", slice_number);
        }
        globalZ += slice_thickness / 2;
        sliceZ += slice_thickness / 2;
        slice_number++;
      }

      layer_volume.setVisAttributes(description.visAttributes(x_layer.visStr()));

      string phys_lay = Form("layer%i", iLayer);
      layerZ += layer_thickness / 2;
      DetElement   layer_elt(endcap, phys_lay, iLayer);
      PlacedVolume pv = endcapVol.placeVolume(layer_volume, Position(0, 0, layerZ));
      pv.addPhysVolID("layer", iLayer);
      layer_elt.setPlacement(pv);
      for (size_t ic = 0; ic < sensitives.size(); ++ic) {
        PlacedVolume sens_pv = sensitives[ic];
        DetElement   comp_elt(layer_elt, sens_pv.volume().name(), iLayer);
        comp_elt.setPlacement(sens_pv);
      }
      layerZ += layer_thickness / 2;
    }
  }
  double       z_pos = zmin + totalThickness / 2;
  PlacedVolume pv;
  // Reflect it.
  Assembly   assembly(det_name);
  DetElement endcapAssyDE(det_name, det_id);
  Volume     motherVol = description.pickMotherVolume(endcapAssyDE);
  if (reflect) {
    pv = assembly.placeVolume(endcapVol, Transform3D(RotationZYX(0, M_PI, 0), Position(0, 0, -z_pos)));
    pv.addPhysVolID("barrel", 2);
    Ref_t(endcap)->SetName((det_name + "_backward").c_str());
    endcap.setPlacement(pv);
  } else {
    pv = assembly.placeVolume(endcapVol, Transform3D(RotationZYX(0, 0, 0), Position(0, 0, z_pos)));
    pv.addPhysVolID("barrel", 1);
    Ref_t(endcap)->SetName((det_name + "_forward").c_str());
    endcap.setPlacement(pv);
  }
  endcapAssyDE.add(endcap);
  pv = motherVol.placeVolume(assembly, Position(pos.x(), pos.y(), pos.z()));
  pv.addPhysVolID("system", det_id);
  endcapAssyDE.setPlacement(pv);
  return endcapAssyDE;
}

// clang-format off
DECLARE_DETELEMENT(epic_PolyhedraEndcapCalorimeter2, create_detector)
DECLARE_DETELEMENT(epic_PolyhedraEndcapCalorimeter, create_detector)
