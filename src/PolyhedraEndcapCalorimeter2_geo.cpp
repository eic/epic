
#include "DD4hep/DetFactoryHelper.h"
#include "TVector3.h"
#include "XML/Layering.h"
using namespace std;
using namespace dd4hep;
using namespace dd4hep::detail;

const double   phiShift        = 0 * M_PI / 180.0; // first sector (1 o'clock) phi value [index 0]
const unsigned nEtaBins        = 12;
const unsigned nHcalSectors    = 12;
const unsigned nHcalSubSectors = 5;

const unsigned nDefaultLayers     = 24;
const double   defaultEtaBins[13] = {2.0,    1.9008, 1.8065, 1.7168, 1.6317, 1.5507, 1.4738,
                                     1.4007, 1.3312, 1.2651, 1.2023, 1.1427, 1.086};

const double newEtaBins[28] = {3.9736842, 3.82929, 3.65542, 3.48947, 3.33106, 3.17987, 3.03556,
                               2.89782,   2.76635, 2.64086, 2.52109, 2.40676, 2.29764, 2.19349,
                               2.09299,   1.99721, 1.90582, 1.81882, 1.7357,  1.65646, 1.58079,
                               1.50851,   1.43941, 1.3734,  1.3104,  1.25011, 1.19095, 1.1158023};

const double defaultZpos[nDefaultLayers] = {270.19,  271.695, 273.15,  274.555, 275.96,  277.365, 282.363, 283.768,
                                            285.173, 286.578, 287.983, 289.388, 290.793, 292.198, 293.603, 295.008,
                                            296.413, 297.818, 299.223, 300.628, 302.033, 303.438, 304.843, 306.158};

const int nNewLayers          = 11;
const int       tileMap[nNewLayers] = {0, 2, 4, 5, 6, 8, 10, 13, 15, 17, 19};
double    getEta(const double& z, const double& r) { return asinh(z / r); }
double    getR(const double& z, const double& eta) { return z / (sinh(eta)); }
double    getNewEta(const double& zOld, const double& etaOld, double zNew) { return getEta(zNew, getR(zOld, etaOld)); }

inline double getEtaMean(unsigned eta) { return 0.5 * (defaultEtaBins[eta] + defaultEtaBins[eta + 1]); }
inline double getEtaHalfWidth(unsigned eta) { return 0.5 * fabs(defaultEtaBins[eta] - defaultEtaBins[eta + 1]); }
// adjust angle so it falls into [-pi,pi] interval
static inline double AdjustAngle(double alpha)
{
  while (alpha < -TMath::Pi())
    alpha += TMath::TwoPi();
  while (alpha > TMath::Pi())
    alpha -= TMath::TwoPi();
  return alpha;
}
inline double getPhiMean(unsigned sector)
{
  const double dPhi = TMath::TwoPi() / nHcalSectors;
  return AdjustAngle((sector + 0.5L) * dPhi + phiShift);
}

inline double getPhiHalfWidth()
{
  const double dPhi = TMath::TwoPi() / nHcalSectors;
  return (double)(0.5L / nHcalSubSectors * dPhi);
}
inline double getPhiMean(unsigned sector, unsigned subSector)
{
  const double dPhi = TMath::TwoPi() / nHcalSectors;
  return AdjustAngle((double(sector) + (subSector + 0.5L) / nHcalSubSectors) * dPhi + phiShift);
}

TVector3 getTowerCenter(const unsigned sector, const unsigned subSector, const unsigned etaBin, const double z)
{
  double phi = 0.0;
  double eta = -1.0;
  double rho = 0.0;

  phi = getPhiMean(sector, subSector);
  eta = getEtaMean(etaBin);

  rho = getR(z, eta);
  // create vector pointing toward the center of the tower
  return TVector3(rho * cos(phi), rho * sin(phi), z);
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

  int    l_num     = 1;
  int    layerType = 0;
  double layerZ    = -totalThickness / 2;

  endcapVol.setAttributes(description, x_det.regionStr(), x_det.limitsStr(), x_det.visStr());

  for (xml_coll_t xc(x_det, _U(layer)); xc; ++xc) {
    // std::cout << "l_num = " << l_num << "\n";
    // std::cout << "xc = " << xc << "\n";
    xml_comp_t x_layer = xc;
    double     l_thick = layering.layer(l_num - 1)->thickness();
    // std::cout << "xc = " << xc << "\n";
    string               l_name   = _toString(layerType, "layer%d");
    int                  l_repeat = x_layer.repeat();
    Volume               l_vol(l_name, PolyhedraRegular(numsides, rmin, rmax, l_thick), air);
    vector<PlacedVolume> sensitives;

    int    s_num  = 1;
    double sliceZ = -l_thick / 2;

    double disksGap = getAttrOrDefault(x_layer, _Unicode(disksGap), 0 * cm);

    // std::cout << "disksGap = " << disksGap << "\n";
    for (xml_coll_t xs(x_layer, _U(slice)); xs; ++xs) {
      xml_comp_t x_slice = xs;
      string     s_name  = _toString(s_num, "slice%d");
      double     s_thick = x_slice.thickness();
      Material   s_mat   = description.material(x_slice.materialStr());
      sliceZ += s_thick / 2;
      int         oldTileLayer = tileMap[l_num - 1]; // layer of the old tile corresponding to the current layer
      double      zOldPos      = defaultZpos[oldTileLayer];         // center z position of the old tile

       Volume absorber("aborber",Tube(rmin, rmax, s_thick / 2, M_PI / 2, M_PI * 3 / 2),s_mat);
      if (x_slice.visStr() == "HcalAbsorberVis") {

        PlacedVolume s_phv1 = l_vol.placeVolume(absorber, Position(-disksGap / 2, 0, sliceZ));
        s_phv1.addPhysVolID("absorberDisk", 0);
        PlacedVolume s_phv2 =
            l_vol.placeVolume(absorber, Transform3D(RotationZYX(M_PI, 0, 0), Position(+disksGap / 2, 0, sliceZ)));
        s_phv2.addPhysVolID("absorberDisk", 1);
      }

      else {
        double dphi = 2 * getPhiHalfWidth();
        Volume subSectorVol("subsector", Tube(rmin, rmax, s_thick / 2, 0, 2 * M_PI / (nHcalSectors * nHcalSubSectors)),
                            air);
        for (unsigned iEtaBin = 0; iEtaBin < nEtaBins; iEtaBin++) {
          // recreate old tile size for the given layer:
        
          double rBottom = getR(zOldPos, defaultEtaBins[iEtaBin]);  // closest current tile distance to the beam pipe
          double rTop = getR(zOldPos, defaultEtaBins[iEtaBin + 1]); // furthest current tile distance to the beam pipe
          // cosine theorem for the angle between the two radii:
          //(2 dx1)^2 = r^2+r^2-2*r*r*cos(phi) = 2r^2(1-cos(phi));
          //   dx1=r*sqrt(2(1-cos(phi)))/2;
          double       dx1 = rBottom * sqrt(2 * (1 - cos(dphi))) / 2;
          double       dx2 = rTop * sqrt(2 * (1 - cos(dphi))) / 2;
          double       dz  = (rTop - rBottom) / 2;
          double       dy  = s_thick / 2;
          Trapezoid    tile(dx1, dx2, dy, dy, dz);
          Volume       tile_vol(_toString(iEtaBin, "tile%i"), tile, s_mat);
          PlacedVolume tile_phys = subSectorVol.placeVolume(tile_vol,  Position(0, 0, (rTop + rBottom) / 2));

          tile_phys.addPhysVolID("etabin", iEtaBin);

          if (x_slice.isSensitive()) {
            sens.setType("calorimeter");
            tile_vol.setSensitiveDetector(sens);
            sensitives.push_back(tile_phys);
          }
        }

        
        
        double      rTop         = getR(zOldPos, defaultEtaBins[12]); // furthest current tile distance to the beam pipe
        ConeSegment sectorShape(s_thick / 2, rmin, rTop, rmin, rTop, 0, 2 * M_PI / nHcalSectors);
        Volume      sectorVol{"sector", sectorShape, air};

        for (unsigned iSubSector = 0; iSubSector < nHcalSubSectors; iSubSector++) {
          PlacedVolume subSectorVol_phys = sectorVol.placeVolume(
              subSectorVol, RotationZYX(0, (iSubSector+0.5L) * 2 * M_PI / (nHcalSectors * nHcalSubSectors), M_PI/2));
          subSectorVol_phys.addPhysVolID("subsector", iSubSector);
        }

        for (unsigned iSector = 0; iSector < nHcalSectors; iSector++) {
          PlacedVolume sectorVol_phys = l_vol.placeVolume(
              sectorVol, Transform3D(RotationZYX(M_PI/2, iSector * 2 * M_PI / (nHcalSectors), M_PI / 2),
                                     Position((iSector > nHcalSectors / 2 ? -disksGap : disksGap), 0, sliceZ)));
          sectorVol_phys.addPhysVolID("sector", iSector);
        }
      }
      sliceZ += s_thick / 2;
      s_num++;
    }

    l_vol.setVisAttributes(description.visAttributes(x_layer.visStr()));
    if (l_repeat <= 0)
      throw std::runtime_error(x_det.nameStr() + "> Invalid repeat value");
    for (int j = 0; j < l_repeat; ++j) {
      string phys_lay = _toString(l_num, "layer%d");
      layerZ += l_thick / 2;
      DetElement   layer_elt(endcap, phys_lay, l_num);
      PlacedVolume pv = endcapVol.placeVolume(l_vol, Position(0, 0, layerZ));
      pv.addPhysVolID("layer", l_num);
      layer_elt.setPlacement(pv);
      for (size_t ic = 0; ic < sensitives.size(); ++ic) {
        PlacedVolume sens_pv = sensitives[ic];
        DetElement   comp_elt(layer_elt, sens_pv.volume().name(), l_num);
        comp_elt.setPlacement(sens_pv);
      }
      layerZ += l_thick / 2;
      ++l_num;
    }
    ++layerType;
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
