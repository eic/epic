
#include <TFile.h>

#include "DD4hep/DetFactoryHelper.h"

#include <XML/Helper.h>

#define _WITH_OPTICS_

// Mimic ePIC pfRICH vessel length;
#define _FIDUCIAL_VOLUME_LENGTH_SIMPLE_    (500.0*mm)
// Choose ePIC dRICH location;
#define _FIDUCIAL_VOLUME_OFFSET_SIMPLE_   (1550.0*mm)
#define _FAKE_HRPPD_PITCH_                  (3.25*mm)
// This one does not really matter;
#define _FIDUCIAL_VOLUME_WIDTH_SIMPLE_    (705*_FAKE_HRPPD_PITCH_)

#define _VESSEL_FRONT_SIDE_THICKNESS_        (5.0*mm)
#define _SENSOR_AREA_LENGTH_                 (5.0*cm)

#define _BUILDING_BLOCK_CLEARANCE_           (1.0*mm)

#define _AEROGEL_THICKNESS_                  (2.5*cm)
#define _ACRYLIC_THICKNESS_                  (3.0*mm)

using namespace dd4hep;

#ifdef _WITH_OPTICS_
#include "IRT/CherenkovDetectorCollection.h"
#endif

#define _FAKE_HRPPD_WINDOW_THICKNESS_            (5.0*mm)
#define _FAKE_HRPPD_TILE_SIZE_          (_FIDUCIAL_VOLUME_WIDTH_SIMPLE_ - 3*_FAKE_HRPPD_PITCH_)
#define _FAKE_HRPPD_ACTIVE_AREA_SIZE_   (_FIDUCIAL_VOLUME_WIDTH_SIMPLE_ - 5*_FAKE_HRPPD_PITCH_)
#define _FAKE_HRPPD_CONTAINER_VOLUME_HEIGHT_    (32.0*mm)

#define _FAKE_SENSOR_PLANE_GEOMETRIC_EFFICIENCY_   (1.00)
#define _FAKE_SAFETY_FACTOR_                       (0.70)

static Ref_t createDetector(Detector& desc, xml::Handle_t handle, SensitiveDetector sens)
{
  xml::DetElement detElem       = handle;
  std::string detName           = detElem.nameStr();
  int detID                     = detElem.id();
  DetElement det(detName, detID);

  // Whatever this means;
  sens.setType("tracker");

  // Start optical configuration if needed;
#ifdef _WITH_OPTICS_
  auto output_file = TFile::Open("qrich-optics.root", "RECREATE");
  auto geometry = new CherenkovDetectorCollection();
  auto cdet = geometry->AddNewDetector(detName.c_str());
#endif
  
  // Vessel material;
  auto vesselMaterial     = desc.material(detElem.attr<std::string>(_Unicode(vessel)));

  // Gas volume;
  auto gasvolMaterialName = detElem.attr<std::string>(_Unicode(gas));
  auto gasvolMaterial     = desc.material(gasvolMaterialName);
  
  // Aerogel;
  auto aerogelMaterialName = detElem.attr<std::string>(_Unicode(aerogel));
  auto aerogelMaterial     = desc.material(aerogelMaterialName);

  // Acrylic filter;
  auto acrylicMaterialName = detElem.attr<std::string>(_Unicode(acrylic));
  auto acrylicMaterial     = desc.material(acrylicMaterialName);
  
  // HRPPD window;
  auto windowMaterialName  = detElem.attr<std::string>(_Unicode(window));
  auto windowMaterial      = desc.material(windowMaterialName);
  
  // Photocathode;
  auto pcMaterialName      = detElem.attr<std::string>(_Unicode(photocathode));
  auto pcMaterial          = desc.material(pcMaterialName);

  auto resinMat       = desc.material("Epoxy");
  
  // - sensor readout
  auto readoutName = detElem.attr<std::string>(_Unicode(readout));

  // readout coder <-> unique sensor ID
  std::vector<std::string> sensorIDfields = {"pdu"};//, "sipm", "sector"};
  const auto& readoutCoder = *desc.readout(readoutName).idSpec().decoder();
  uint64_t cellMask = 0;
  for (const auto& idField : sensorIDfields)
    cellMask |= readoutCoder[idField].mask();
  // Do not mind to store it twice;
  desc.add(Constant("QRICH_cell_mask", std::to_string(cellMask)));
#ifdef _WITH_OPTICS_
  cdet->SetReadoutCellMask(cellMask);
#endif
  // create a unique sensor ID from a sensor's PlacedVolume::volIDs
  auto encodeSensorID = [&readoutCoder](auto ids) {
    uint64_t enc = 0;
    for (const auto& [idField, idValue] : ids)
      enc |= uint64_t(idValue) << readoutCoder[idField].offset();
    return enc;
  };

  const double fvLength = _FIDUCIAL_VOLUME_LENGTH_SIMPLE_, fvWidth  = _FIDUCIAL_VOLUME_WIDTH_SIMPLE_;
  const double fvOffset = _FIDUCIAL_VOLUME_OFFSET_SIMPLE_;
  
  // Fiducial volume / vessel;
  auto vesselSolid = Box(fvWidth/2, fvWidth/2, fvLength/2);
  Volume vesselVolume(detName, vesselSolid, vesselMaterial);

  // Place vessel into the mother volume;
  Volume motherVol = desc.pickMotherVolume(det);
  PlacedVolume vesselPV = motherVol.placeVolume(vesselVolume,  Position(0, 0, fvOffset));
  vesselPV.addPhysVolID("system", detID);
  det.setPlacement(vesselPV);
    
  // Gas volume, aerogel, acrylic filter;
  {
    const double gvWidth = fvWidth - 1*mm;
    const double gvLength = fvLength - _VESSEL_FRONT_SIDE_THICKNESS_ - _SENSOR_AREA_LENGTH_;
    auto gasvolSolid = Box(gvWidth/2, gvWidth/2, gvLength/2);
    Volume gasvolVol(detName + "_gas", gasvolSolid, gasvolMaterial);
    
    // Place gas volume into the vessel;
    double gvOffset = -(_SENSOR_AREA_LENGTH_ - _VESSEL_FRONT_SIDE_THICKNESS_)/2;
    PlacedVolume gasvolPV = vesselVolume.placeVolume(gasvolVol, Position(0, 0, gvOffset));
    DetElement gasvolDE(det, "gasvol_de", 0);
    gasvolDE.setPlacement(gasvolPV);

#ifdef _WITH_OPTICS_
    {
      // FIXME: Z-location does not really matter here, right?;
      auto boundary = new FlatSurface(TVector3(0,0,0), /*sign**/TVector3(1,0,0), TVector3(0,-1,0));
      
      auto radiator = geometry->SetContainerVolume(cdet, "GasVolume", 0, (G4LogicalVolume*)(0x0), 0, boundary);
      // FIXME: hardcoded;
      radiator->SetReferenceRefractiveIndex(1.00033);
      radiator->SetAlternativeMaterialName(gasvolMaterialName.c_str());
    }
#endif
    
    {
      double gzOffset = -gvLength/2 + _BUILDING_BLOCK_CLEARANCE_;
    
      // Aerogel volume;
      {
	const double agWidth = gvWidth - 1*mm, agThick = _AEROGEL_THICKNESS_;
	Box aerogelSolid(agWidth/2, agWidth/2, agThick/2);
	Volume aerogelVol(detName + "_aerogel", aerogelSolid, aerogelMaterial);
	
	// Place aerogel into gas volume; NB: 'gzOffset' used in a FlatSurface() description below; 
	gzOffset += agThick/2;
	gasvolVol.placeVolume(aerogelVol, Position(0, 0, gzOffset));
	
#ifdef _WITH_OPTICS_
	{
	  TVector3 nx(1,0,0), ny(0,-1,0);
	
	  auto surface = new FlatSurface((1/mm)*TVector3(0,0,fvOffset + gvOffset + gzOffset), nx, ny);
	  
	  auto radiator = geometry->AddFlatRadiator(cdet, "Aerogel", CherenkovDetector::Upstream, 
						      0, (G4LogicalVolume*)(0x1), 0, surface, agThick/mm);
	  radiator->SetAlternativeMaterialName(aerogelMaterialName.c_str());
	  // FIXME: what is it good for?;
	  geometry->AddRadiatorLogicalVolume(radiator, (G4LogicalVolume*)(0x1));//ag_log);
	}
#endif
    
	gzOffset += agThick/2 + _BUILDING_BLOCK_CLEARANCE_;
      }
      
      // Acrylic filter volume;
      {
	const double acWidth = gvWidth - 1*mm, acThick = _ACRYLIC_THICKNESS_;
	Box acrylicSolid(acWidth/2, acWidth/2, acThick/2);
	Volume acrylicVol(detName + "_filter", acrylicSolid, acrylicMaterial);
	
	// Place acrylic filter into gas volume; NB: 'gzOffset' used in a FlatSurface() description below;
	gzOffset += acThick/2;
	gasvolVol.placeVolume(acrylicVol, Position(0, 0, gzOffset));
	
#ifdef _WITH_OPTICS_
	{
	  TVector3 nx(1,0,0), ny(0,-1,0);
	
	  auto surface = new FlatSurface((1/mm)*TVector3(0,0,fvOffset + gvOffset + gzOffset), nx, ny);
	  
	  auto radiator = geometry->AddFlatRadiator(cdet, "Acrylic", CherenkovDetector::Upstream, 
						    0, (G4LogicalVolume*)(0x2), 0, surface, acThick/mm);
	  // FIXME: hardcoded;
	  radiator->SetReferenceRefractiveIndex(1.5017);
	  radiator->SetAlternativeMaterialName(acrylicMaterialName.c_str());
	}
#endif
      }
    }
  }

  // A fake HRPPD container volume, window, photocathode;
  {
    // Mimic BuildFakePhotonDetectorMatrix() used in the standalone code;
    const double fvzOffset = fvOffset, wzOffset = fvLength/2 - _SENSOR_AREA_LENGTH_;
    const double pdthick = 0.01*mm, zpdc = wzOffset + _FAKE_HRPPD_WINDOW_THICKNESS_ + pdthick/2;
    auto pd_box  = Box(/*"PhotoDetector",*/ _FAKE_HRPPD_ACTIVE_AREA_SIZE_/2,
		       _FAKE_HRPPD_ACTIVE_AREA_SIZE_/2, pdthick/2);
    Volume pcVol(detName + "_photodetector", pd_box, pcMaterial);
    pcVol.setSensitiveDetector(sens);
    
#ifdef _WITH_OPTICS_
    //auto pd = new CherenkovPhotonDetector(pd_box, m_Bialkali);
    // [0,0]: have neither access to G4VSolid nor to G4Material; IRT code does not care; fine;
    auto pd = new CherenkovPhotonDetector(0, 0);
    
    // FIXME: '0' stands for the unknown (and irrelevant) G4LogicalVolume;
    geometry->AddPhotonDetector(cdet, 0, pd);
#endif
    
    // Full size quartz window; FIXME: Sapphire, here and in all other places;
    auto wnd_box = Box(/*"QuartzWindow",*/ _FAKE_HRPPD_TILE_SIZE_/2, _FAKE_HRPPD_TILE_SIZE_/2,
		       _FAKE_HRPPD_WINDOW_THICKNESS_/2);
    Volume wnd_log(detName + "_wnd", wnd_box, windowMaterial);
    
    //+auto hrppd_log = BuildFakeHRPPD(wnd_log, pd_box, pd);
#ifdef _WITH_OPTICS_
    {	
      TVector3 nx(1,0,0), ny(0,-1,0);
      
      auto surface = 
	new FlatSurface((1/mm)*TVector3(0,0,fvzOffset + wzOffset + _FAKE_HRPPD_WINDOW_THICKNESS_/2), nx, ny);
      
      auto radiator = geometry->AddFlatRadiator(cdet, "QuartzWindow", CherenkovDetector::Downstream, 
						0, (G4LogicalVolume*)(0x3), 0, /*wnd_log, m_FusedSilica,*/ surface,
						_FAKE_HRPPD_WINDOW_THICKNESS_/mm);
      radiator->SetReferenceRefractiveIndex(1.46);//1.00029);
      radiator->SetAlternativeMaterialName(windowMaterialName.c_str());
    }	
#endif
    
    //?m_Geometry->AddPhotonDetector(cdet, pd->GetLogicalVolume(), pd);
    {
      auto hrppd_box =
	Box(_FAKE_HRPPD_TILE_SIZE_/2, _FAKE_HRPPD_TILE_SIZE_/2, _FAKE_HRPPD_CONTAINER_VOLUME_HEIGHT_/2);
      // FIXME: think (OLD: well, let it be vacuum inside);
      Volume hrppd_log(detName + "_hrppd", hrppd_box, gasvolMaterial);//vacuumMat);
      
#ifdef _WITH_OPTICS_
      // FIXME: really needed in EICrecon?;
      pd->SetCopyIdentifierLevel(1);
      //?pd->DefineLogicalVolume();
      // Cannot access GEANT shapes in the reconstruction code -> store this value;
      pd->SetActiveAreaSize(_FAKE_HRPPD_ACTIVE_AREA_SIZE_/mm);
      
      pd->SetGeometricEfficiency(_FAKE_SENSOR_PLANE_GEOMETRIC_EFFICIENCY_ * _FAKE_SAFETY_FACTOR_);
#endif
      
#if _LATER_
      auto qd_box  =
	Box(_FAKE_HRPPD_ACTIVE_AREA_SIZE_/2, _FAKE_HRPPD_ACTIVE_AREA_SIZE_/2, pdthick/2);
      Volume qd_log(detName + "_fake_photodetector", qd_box, bialkaliMat);
#endif
      
      {
	double accu = -_FAKE_HRPPD_CONTAINER_VOLUME_HEIGHT_/2;
	
	// Window layer;
	hrppd_log.placeVolume(wnd_log, Position(0, 0, accu + _FAKE_HRPPD_WINDOW_THICKNESS_/2));
	accu += _FAKE_HRPPD_WINDOW_THICKNESS_;
	
	// Fake photodector layer (photocathode);
	//new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, accu + pdthick/2), qd_log, "FakePhotoDetector", 
	//		    hrppd_log, false, 0);
	//+hrppd_log.placeVolume(qd_log, Position(0, 0, accu + pdthick/2));
	//+accu += pdthick;
	
	// Photodector layer (photocathode);
	auto pcPlacement = Position(0, 0, accu + pdthick/2);
	auto pcPV = hrppd_log.placeVolume(pcVol, pcPlacement);
	accu += pdthick;
	{
	  double pcThickness = pdthick;
	  double resinSide = _FAKE_HRPPD_ACTIVE_AREA_SIZE_ + 1*mm, resinThickness = pcThickness + 1*mm;
	  Box resinSolid(resinSide / 2., resinSide / 2., resinThickness / 2.);
	  
	  auto resinPlacement = Position(0, 0, accu + resinThickness/2);
	  Volume resinVol("_resin_", resinSolid/*Embedded*/, resinMat);
	  hrppd_log.placeVolume(resinVol, resinPlacement);
	}
	{
	  // sensor readout // NOTE: follow `sensorIDfields`
	  pcPV.addPhysVolID("pdu", 0);
	  
	  // sensor DetElement
	  auto sensorID = encodeSensorID(pcPV.volIDs());
	  printf("@S@: %lu\n", sensorID); 
	  //std::string sensorIDname =
	  //secName + "_pdu" + std::to_string(ipdu) + "_sipm" + std::to_string(isipm);
	  DetElement pcDE(det, "sensor_de_"/* + sensorIDname*/, sensorID);
	  pcDE.setPlacement(pcPV);
	}
      }
            
      vesselVolume.placeVolume(hrppd_log, Position(0, 0, wzOffset + _FAKE_HRPPD_CONTAINER_VOLUME_HEIGHT_/2));
    }
    
#ifdef _WITH_OPTICS_
    {
      // Photocathode surface;
      auto surface = new FlatSurface((1/mm)*TVector3(0.0, 0.0, fvzOffset + zpdc),
				     TVector3(1,0,0), TVector3(0,-1,0));
      
      {
	// Mimic det->CreatePhotonDetectorInstance();
	unsigned sector = 0, icopy = 0;
	auto irt = pd->AllocateIRT(sector, icopy);
	
	// Aerogel and acrylic;
	if (cdet->m_OpticalBoundaries[CherenkovDetector::Upstream].find(sector) != 
	    cdet->m_OpticalBoundaries[CherenkovDetector::Upstream].end())
	  for(auto boundary: cdet->m_OpticalBoundaries[CherenkovDetector::Upstream][sector])
	    irt->AddOpticalBoundary(boundary);
	
	// Fused silica windows;
	if (cdet->m_OpticalBoundaries[CherenkovDetector::Downstream].find(sector) != 
	    cdet->m_OpticalBoundaries[CherenkovDetector::Downstream].end())
	  for(auto boundary: cdet->m_OpticalBoundaries[CherenkovDetector::Downstream][sector])
	    irt->AddOpticalBoundary(boundary);
	
	// Terminate the optical path;
	pd->AddItselfToOpticalBoundaries(irt, surface);
      }
    } 
#endif
    
    //for(auto radiator: cdet->Radiators())
    //radiator.second->SetReferenceRefractiveIndex(radiator.second->GetMaterial()->RefractiveIndex(eV*_MAGIC_CFF_/_LAMBDA_NOMINAL_));
  }
    
  // Write out optical configuration;
#ifdef _WITH_OPTICS_
  geometry->Write();
  output_file->Close();
#endif
  
  return det;
}

DECLARE_DETELEMENT(epic_QRICH, createDetector)
