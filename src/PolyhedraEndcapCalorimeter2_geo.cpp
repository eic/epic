#include "DD4hep/DetFactoryHelper.h"
#include "TVector3.h"
#include "XML/Layering.h"
using namespace std;
using namespace dd4hep;
using namespace dd4hep::detail;

const unsigned nHcalSectors    = 12;
const unsigned nHcalSubSectors = 5;

const double starEmcEtaBins[13] = {2.0,    1.9008, 1.8065, 1.7168, 1.6317, 1.5507, 1.4738,
                                   1.4007, 1.3312, 1.2651, 1.2023, 1.1427, 1.086};

const double newHcalEtaBins[7] = {3.13294, 2.99414, 2.86202, 2.73624, 2.61643, 2.50226, 2.39341};

const double starEmcZPositions[24] = {270.19,  271.695, 273.15,  274.555, 275.96,  277.365, 282.363, 283.768,
                                      285.173, 286.578, 287.983, 289.388, 290.793, 292.198, 293.603, 295.008,
                                      296.413, 297.818, 299.223, 300.628, 302.033, 303.438, 304.843, 306.158};

const int nNewLayers          = 11;
const int tileMap[nNewLayers] = {0, 2, 4, 5, 6, 8, 10, 12, 15, 17, 19};
double    getEta(const double& z, const double& r) { return asinh(z / r); }
double    getR(const double& z, const double& eta) { return z / (sinh(eta)); }

void buildSubsector(Detector& description, xml_h e, SensitiveDetector sens, vector<PlacedVolume>& sensitives,
                    const double& s_thick, const int& layerNumber, const xml_comp_t x_slice, Volume& sectorVol,
                    const int& nSubSectors, const int& nEtaBins, const double* etaBins, bool isStarEmc = false)
{
  xml_det_t x_det        = e;
  xml_dim_t dim          = x_det.dimensions();
  Material  air          = description.air();
  double    zmin         = dim.zmin();
  double    dphiSector   = 2 * M_PI / (nHcalSectors);
  double    dphi         = dphiSector / nSubSectors;
  int       oldTileLayer = tileMap[layerNumber - 1];        // layer of the old tile corresponding to the current layer
  double    zOldPos      = starEmcZPositions[oldTileLayer]; // center z position of the old tile
  double    zPos         = zmin + (layerNumber - 1) * s_thick + s_thick / 2;
  double    oldEmcGap    = 2;                               // cm

  double rmax = isStarEmc ? getR(zOldPos, starEmcEtaBins[nEtaBins])
                          : getR(zPos, etaBins[nEtaBins]) - oldEmcGap; // closest tile distance to the beam pipe
  double rmin = isStarEmc ? getR(zOldPos, starEmcEtaBins[0])
                          : getR(zPos, etaBins[0]) - oldEmcGap;        // furthest tile distance to the beam pipe

  Volume subSectorVol("subsector", Tube(rmin, rmax, s_thick / 2, 0, dphi), air);

  for (int iEtaBin = 0; iEtaBin < nEtaBins; iEtaBin++) {
    // recreate tile size for the given layer:
    double rBottom = isStarEmc
                         ? getR(zOldPos, starEmcEtaBins[iEtaBin])
                         : getR(zPos, etaBins[iEtaBin]) - oldEmcGap; // current tile closest distance to the beam pipe
    double rTop    = isStarEmc
                         ? getR(zOldPos, starEmcEtaBins[iEtaBin + 1])
                         : getR(zPos, etaBins[iEtaBin + 1]) - oldEmcGap; // current tile closest distance to the beam pipe
    // cosine theorem for the angle between the two radii:
    //(2 dx1)^2 = r^2+r^2-2*r*r*cos(phi) = 2r^2(1-cos(phi));
    //   dx1=r*sqrt(2(1-cos(phi)))/2;
    double    dx1 = rBottom * sqrt(2 * (1 - cos(dphi))) / 2;
    double    dx2 = rTop * sqrt(2 * (1 - cos(dphi))) / 2;
    double    dz  = (rTop - rBottom) / 2;
    double    dy  = s_thick / 2;
    Trapezoid tile(dx1, dx2, dy, dy, dz);
    Material  s_mat = description.material(x_slice.materialStr());
    string    tileName =
        isStarEmc ? _toString(iEtaBin, "tileOld%i") : Form("tile%i_%i", (int)(dphi / M_PI * 180 + 0.5), iEtaBin);
    Volume       tile_vol(tileName, tile, s_mat);
    PlacedVolume tile_phys = subSectorVol.placeVolume(tile_vol, Position(0, 0, (rTop + rBottom) / 2));
    tile_phys.addPhysVolID("etabin", iEtaBin);

    if (x_slice.isSensitive()) {
      sens.setType("calorimeter");
      tile_vol.setSensitiveDetector(sens);
      sensitives.push_back(tile_phys);
    }
  }

  for (int iSubSector = 0; iSubSector < nSubSectors; iSubSector++) {
    PlacedVolume subSectorVol_phys =
        sectorVol.placeVolume(subSectorVol, RotationZYX(0, (iSubSector + 0.5L) * dphi, M_PI / 2));
    subSectorVol_phys.addPhysVolID("subsector", iSubSector);
  }
}

static Ref_t create_detector(Detector& description, xml_h e, SensitiveDetector sens)
{
  xml_det_t      x_det    = e;
  xml_dim_t      dim      = x_det.dimensions();
  int            det_id   = x_det.id();
  bool           reflect  = x_det.reflect(true);
  string         det_name = x_det.nameStr();
  Material       air      = description.air();
  int            numsides = dim.numsides();
  xml::Component pos      = x_det.position();
  double         rmin     = dim.rmin();
  double         rmax     = dim.rmax();
  double         zmin     = dim.zmin();

  Layering   layering(x_det);
  double     totalThickness = layering.totalThickness();
  Volume     endcapVol("endcap", PolyhedraRegular(numsides, rmin, rmax, totalThickness), air);
  DetElement endcap("endcap", det_id);

  double layerZ = -totalThickness / 2;
  endcapVol.setAttributes(description, x_det.regionStr(), x_det.limitsStr(), x_det.visStr());

  xml_coll_t xf(x_det, _U(layer));
  xml_comp_t x_layer_first = xf;
  int        l_repeat      = x_layer_first.repeat();

  if (l_repeat <= 0) {
    throw std::runtime_error(x_det.nameStr() + "> Invalid repeat value");
  }
  for (int iLayer = 1; iLayer <= l_repeat; ++iLayer) {

    for (xml_coll_t xc(x_det, _U(layer)); xc; ++xc) {
      xml_comp_t           x_layer = xc;
<<<<<<< HEAD
      double               l_thick = layering.layer(iLayer - 1)->thickness();
      string               l_name  = _toString(iLayer-1, "layer%d");
=======
      double               l_thick = layering.layer(l_num - 1)->thickness();
      string               l_name  = _toString(layerType, "layer%d");
>>>>>>> cd1974dca74022f5cd744909da0711ae3c9149fc
      Volume               l_vol(l_name, PolyhedraRegular(numsides, rmin, rmax, l_thick), air);
      vector<PlacedVolume> sensitives;
      int                  s_num    = 1;
      double               sliceZ   = -l_thick / 2;
      double               disksGap = getAttrOrDefault(x_layer, _Unicode(disksGap), 1 * cm);

      for (xml_coll_t xs(x_layer, _U(slice)); xs; ++xs) {
        xml_comp_t x_slice = xs;
        string     s_name  = _toString(s_num, "slice%d");
        double     s_thick = x_slice.thickness();
        Material   s_mat   = description.material(x_slice.materialStr());
        sliceZ += s_thick / 2;

        Volume slice("slice", Tube(rmin, rmax, s_thick / 2), air);

        if (x_slice.visStr() == "HcalAbsorberVis") {
          Volume       absorber("absorber", Tube(rmin, rmax, s_thick / 2, M_PI / 2, M_PI * 3 / 2), s_mat);
          PlacedVolume s_phv1 = slice.placeVolume(absorber, Position(-disksGap / 2, 0, 0));
          s_phv1.addPhysVolID("absorberDisk", 0);
          PlacedVolume s_phv2 =
              slice.placeVolume(absorber, Transform3D(RotationZYX(M_PI, 0, 0), Position(+disksGap / 2, 0, 0)));
          s_phv2.addPhysVolID("absorberDisk", 1);

          PlacedVolume s_phv = l_vol.placeVolume(slice, Position(0, 0, sliceZ));
          s_phv.addPhysVolID("slice", s_num);
        }

        else {
          ConeSegment sectorShape(s_thick / 2, rmin, rmax, rmin, rmax, 0, 2 * M_PI / (nHcalSectors));
          Volume      sectorVol{"sector", sectorShape, air};
<<<<<<< HEAD
          // void buildSubsector(Detector& description, xml_h e, SensitiveDetector sens, vector<PlacedVolume>&
          // sensitives,
          //     const double& s_thick, const int& layerNumber, const xml_comp_t x_slice, Volume& sectorVol,
          //     const int& nSubSectors, const int& nEtaBins, const double* etaBins, bool isStarEmc = false)

          buildSubsector(description, e, sens, sensitives, s_thick, iLayer, x_slice, sectorVol, 5, 12, starEmcEtaBins,
                         true);               // star Emc implementation for new Z position with 12 eta bins
          buildSubsector(description, e, sens, sensitives, s_thick, iLayer, x_slice, sectorVol, 5, 3,
                         &newHcalEtaBins[3]); // 3 eta bins of inner sector with the same substructure as Star Emc with
                                              // 6 (30/5)degree tiles
          buildSubsector(description, e, sens, sensitives, s_thick, iLayer, x_slice, sectorVol, 3, 2,
                         &newHcalEtaBins[1]); // 2 eta bins of inner sector with 10 (30/3) degree tiles
          buildSubsector(description, e, sens, sensitives, s_thick, iLayer, x_slice, sectorVol, 2, 1,
=======
          // buildSubsector( ... const int& nSubSectors, const int& nEtaBins, const double* etaBins, const double&
          // s_thick, const int& layerNumber, const double & disksGap, bool isStarEmc = false)\

          // void buildSubsector(Detector& description, xml_h e, SensitiveDetector sens, vector<PlacedVolume>&
          // sensitives,
          //     const double& s_thick, const int& layerNumber, const xml_comp_t x_slice, Volume& sectorVol,
          //     const int& nSubSectors, const int& nEtaBins, const double* etaBins, bool isStarEmc = false)

          buildSubsector(description, e, sens, sensitives, s_thick, l_num, x_slice, sectorVol, 5, 12, starEmcEtaBins,
                         true);               // star Emc implementation for new Z position with 12 eta bins
          buildSubsector(description, e, sens, sensitives, s_thick, l_num, x_slice, sectorVol, 5, 3,
                         &newHcalEtaBins[3]); // 3 eta bins of inner sector with the same substructure as Star Emc with 6
                                              // (30/5)degree tiles
          buildSubsector(description, e, sens, sensitives, s_thick, l_num, x_slice, sectorVol, 3, 2,
                         &newHcalEtaBins[1]); // 2 eta bins of inner sector with 10 (30/3) degree tiles
          buildSubsector(description, e, sens, sensitives, s_thick, l_num, x_slice, sectorVol, 2, 1,
>>>>>>> cd1974dca74022f5cd744909da0711ae3c9149fc
                         &newHcalEtaBins[0]); // 1 eta bins of innermost sector with 15 (30/2) degree tiles

          for (unsigned iSector = 0; iSector < nHcalSectors; iSector++) {
            PlacedVolume sectorVol_phys = slice.placeVolume(
                sectorVol, Transform3D(RotationZYX(M_PI / 2, iSector * 2 * M_PI / (nHcalSectors), M_PI / 2),
                                       Position(((iSector + 1) > nHcalSectors / 2 ? -disksGap : disksGap), 0, 0)));
            sectorVol_phys.addPhysVolID("sector", iSector);
          }
          PlacedVolume s_phv = l_vol.placeVolume(slice, Position(0, 0, sliceZ));
          s_phv.addPhysVolID("slice", s_num);
        }
        sliceZ += s_thick / 2;
        s_num++;
      }

      l_vol.setVisAttributes(description.visAttributes(x_layer.visStr()));

<<<<<<< HEAD
      string phys_lay = Form("layer%i", iLayer);
=======
      string phys_lay = Form("layer%i", l_num);
>>>>>>> cd1974dca74022f5cd744909da0711ae3c9149fc
      layerZ += l_thick / 2;
      DetElement   layer_elt(endcap, phys_lay, iLayer);
      PlacedVolume pv = endcapVol.placeVolume(l_vol, Position(0, 0, layerZ));
      pv.addPhysVolID("layer", iLayer);
      layer_elt.setPlacement(pv);
      for (size_t ic = 0; ic < sensitives.size(); ++ic) {
        PlacedVolume sens_pv = sensitives[ic];
        DetElement   comp_elt(layer_elt, sens_pv.volume().name(), iLayer);
        comp_elt.setPlacement(sens_pv);
      }
      layerZ += l_thick / 2;
    }
  }

  double       z_pos = zmin + totalThickness / 2;
  PlacedVolume pv;
  // Reflect it.
  Assembly   assembly(det_name);
  DetElement endcapAssyDE(det_name, det_id);
  Volume     motherVol = description.pickMotherVolume(endcapAssyDE);
  if (reflect) {
    pv = assembly.placeVolume(endcapVol, Transform3D(RotationZYX(M_PI / numsides, M_PI, 0), Position(0, 0, -z_pos)));
    pv.addPhysVolID("barrel", 2);
    Ref_t(endcap)->SetName((det_name + "_backward").c_str());
    endcap.setPlacement(pv);
  } else {
    pv = assembly.placeVolume(endcapVol, Transform3D(RotationZYX(M_PI / numsides, 0, 0), Position(0, 0, z_pos)));
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
