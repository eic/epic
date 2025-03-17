
#include <TFile.h>

#include "DD4hep/DetFactoryHelper.h"

#include <XML/Helper.h>

#define _WITH_OPTICS_

//#include "simple.default.h"
// Mimic ePIC pfRICH vessel length;
#define _FIDUCIAL_VOLUME_LENGTH_SIMPLE_    (500.0*mm)
// Choose ePIC dRICH location;
#define _FIDUCIAL_VOLUME_OFFSET_SIMPLE_   (1550.0*mm)
#define _FAKE_HRPPD_PITCH_                  (3.25*mm)
// This one does not really matter;
#define _FIDUCIAL_VOLUME_WIDTH_SIMPLE_    (705*_FAKE_HRPPD_PITCH_)//2000.0*mm)

#define _VESSEL_FRONT_SIDE_THICKNESS_        (5.0*mm)
#define _SENSOR_AREA_LENGTH_                 (5.0*cm)

#define _BUILDING_BLOCK_CLEARANCE_           (1.0*mm)

#define _AEROGEL_THICKNESS_                  (2.5*cm)
#define _ACRYLIC_THICKNESS_                  (3.0*mm)

#define _MAGIC_CFF_                          (1239.8)

using namespace dd4hep;

#ifdef _WITH_OPTICS_
#include "IRT/CherenkovDetectorCollection.h"
#endif

#define _FAKE_HRPPD_WINDOW_THICKNESS_            (5.0*mm)
#define _FAKE_HRPPD_TILE_SIZE_          (_FIDUCIAL_VOLUME_WIDTH_SIMPLE_ - 3*_FAKE_HRPPD_PITCH_)//10*mm)
#define _FAKE_HRPPD_ACTIVE_AREA_SIZE_   (_FIDUCIAL_VOLUME_WIDTH_SIMPLE_ - 5*_FAKE_HRPPD_PITCH_)//20*mm)
#define _FAKE_HRPPD_CONTAINER_VOLUME_HEIGHT_    (32.0*mm)

//#define _FAKE_QE_DOWNSCALING_FACTOR_          (30.0/37.0)
#define _FAKE_SENSOR_PLANE_GEOMETRIC_EFFICIENCY_   (1.00)
#define _FAKE_SAFETY_FACTOR_                       (0.70)

//static CherenkovDetectorCollection *geoptr = new CherenkovDetectorCollection();

//#include <SubtractionSolid.h>


static Ref_t createDetector(Detector& desc, xml::Handle_t handle, SensitiveDetector sens)
{
  printf("@@@ Ku-start!\n");

  // Start optical configuration;
#ifdef _WITH_OPTICS_
  auto output_file = TFile::Open("qrich-optics.root", "RECREATE");
  auto geometry = new CherenkovDetectorCollection();
  auto cdet = geometry->AddNewDetector("QRICH");
  //det->SetReadoutCellMask(~0x0);
  cdet->SetReadoutCellMask(0xFFFFFFFFFFFFFFFF);
#endif
  
  xml::DetElement detElem       = handle;
  std::string detName           = detElem.nameStr();
  int detID                     = detElem.id();
  OpticalSurfaceManager surfMgr = desc.surfaceManager();
  DetElement det(detName, detID);
  sens.setType("tracker");

  //auto vesselMat    = desc.material("AirOptical");//detElem.attr<std::string>(_Unicode(material)));
  //auto gasvolMat    = desc.material("AirOptical");
  auto vesselMat    = desc.material(detElem.attr<std::string>(_Unicode(material)));
  auto gasvolMat    = desc.material(detElem.attr<std::string>(_Unicode(gas)));
  auto resinMat       = desc.material("Epoxy");
  
  auto pssElem      = detElem.child(_Unicode(sensors)).child(_Unicode(pss));
  auto pssMat       = desc.material(pssElem.attr<std::string>(_Unicode(material)));
  //auto pssSurf      = surfMgr.opticalSurface(pssElem.attr<std::string>(_Unicode(surface)));
  auto pssSurf      = surfMgr.opticalSurface("SensorSurface_DRICH");
  
  auto aerogelMat   = desc.material("Aerogel_QRICH");
  auto acrylicMat   = desc.material("Acrylic_DRICH");
  auto quartzMat    = desc.material("QuartzOptical");
  //auto vacuumMat    = desc.material("AirOptical");//VacuumOptical");
  //auto bialkaliMat  = desc.material("Bialkali");
  
  // - sensor readout
  auto readoutName = detElem.attr<std::string>(_Unicode(readout));

  // readout coder <-> unique sensor ID
  /* - `sensorIDfields` is a list of readout fields used to specify a unique sensor ID
   * - `cellMask` is defined such that a hit's `cellID & cellMask` is the corresponding sensor's unique ID
   */
#if 1//_TODAY_
  std::vector<std::string> sensorIDfields = {"pdu"};//, "sipm", "sector"};
  const auto& readoutCoder                = *desc.readout(readoutName).idSpec().decoder();
  // create a unique sensor ID from a sensor's PlacedVolume::volIDs
  auto encodeSensorID = [&readoutCoder](auto ids) {
    uint64_t enc = 0;
    for (const auto& [idField, idValue] : ids)
      enc |= uint64_t(idValue) << readoutCoder[idField].offset();
    return enc;
  };
#endif

  // Fiducial volume / vessel;
  double fvLength = _FIDUCIAL_VOLUME_LENGTH_SIMPLE_;
  double fvWidth  = _FIDUCIAL_VOLUME_WIDTH_SIMPLE_;
  auto vesselSolid = Box(fvWidth/2, fvWidth/2, fvLength/2);
  Volume vesselVol(detName, vesselSolid, vesselMat);

  // Place vessel into the mother volume;
  Volume motherVol      = desc.pickMotherVolume(det);
  double fvOffset = _FIDUCIAL_VOLUME_OFFSET_SIMPLE_;
  auto vesselPos = Position(0, 0, fvOffset);
  PlacedVolume vesselPV = motherVol.placeVolume(vesselVol, vesselPos);
  vesselPV.addPhysVolID("system", detID);
  det.setPlacement(vesselPV);
    
  // Gas volume;
  {
    double gvWidth = fvWidth - 1*mm;
    double gvLength = fvLength - _VESSEL_FRONT_SIDE_THICKNESS_ - _SENSOR_AREA_LENGTH_;
    auto gasvolSolid = Box(gvWidth/2, gvWidth/2, gvLength/2);
    Volume gasvolVol(detName + "_gas", gasvolSolid, gasvolMat);
    
    // Place gas volume into the vessel;
    double gvOffset = -(_SENSOR_AREA_LENGTH_ - _VESSEL_FRONT_SIDE_THICKNESS_)/2;
    PlacedVolume gasvolPV = vesselVol.placeVolume(gasvolVol, Position(0, 0, gvOffset));
    DetElement gasvolDE(det, "gasvol_de", 0);
    gasvolDE.setPlacement(gasvolPV);

#ifdef _WITH_OPTICS_
    {
      // FIXME: Z-location does not really matter here, right?;
      auto boundary = new FlatSurface(TVector3(0,0,0), /*sign**/TVector3(1,0,0), TVector3(0,-1,0));
      
      geometry->SetContainerVolume(cdet, "GasVolume", 0, (G4LogicalVolume*)(0x0), 0, boundary);
    }
#endif
    
    {
      double gzOffset = -gvLength/2 + _BUILDING_BLOCK_CLEARANCE_;
    
      // Aerogel volume;
      {
	double agWidth = gvWidth - 1*mm, agThick = _AEROGEL_THICKNESS_;
	Box aerogelSolid(agWidth/2, agWidth/2, agThick/2);
	//Volume aerogelVol(detName + "_aerogel_" + secName, aerogelSolid, aerogelMat);
	Volume aerogelVol(detName + "_aerogel", aerogelSolid, aerogelMat);
	
	// Place aerogel into gas volume;
	gzOffset += agThick/2;
	gasvolVol.placeVolume(aerogelVol, Position(0, 0, gzOffset));
	
#ifdef _WITH_OPTICS_
	{
	  TVector3 nx(1,0,0), ny(0,-1,0);
	
	  auto surface = new FlatSurface((1/mm)*TVector3(0,0,fvOffset + gvOffset + gzOffset), nx, ny);
	  
	  auto radiator = geometry->AddFlatRadiator(cdet, "Aerogel", CherenkovDetector::Upstream, 
						      0, (G4LogicalVolume*)(0x1), 0, surface, agThick/mm);
	  geometry->AddRadiatorLogicalVolume(radiator, (G4LogicalVolume*)(0x1));//ag_log);
	}
#endif
    
	gzOffset += agThick/2 + _BUILDING_BLOCK_CLEARANCE_;
      }
      
      // Acrylic filter volume;
#if _TODAY_
      {
	double acWidth = gvWidth - 1*mm, acThick = _ACRYLIC_THICKNESS_;
	Box acrylicSolid(acWidth/2, acWidth/2, acThick/2);
	Volume acrylicVol(detName + "_filter", acrylicSolid, acrylicMat);
	
	// Place acrylic filter into gas volume;
	gzOffset += acThick/2;
	gasvolVol.placeVolume(acrylicVol, Position(0, 0, gzOffset));
	
#if 0//def _WITH_OPTICS_
	{
	  TVector3 nx(1,0,0), ny(0,-1,0);
	
	  auto surface = new FlatSurface((1/mm)*TVector3(0,0,fvOffset + gvOffset + gzOffset), nx, ny);
	  
	  /*auto radiator =*/ geometry->AddFlatRadiator(cdet, "Acrylic", CherenkovDetector::Upstream, 
							0, (G4LogicalVolume*)(0x2), 0, surface, acThick/mm);
	}
#endif
      }
#endif
    }

    {
      // Mimic BuildFakePhotonDetectorMatrix() used in the standalone code;
      double fvzOffset = fvOffset, wzOffset = fvLength/2 - _SENSOR_AREA_LENGTH_;
      
      double pdthick = 0.01*mm, zpdc = wzOffset + _FAKE_HRPPD_WINDOW_THICKNESS_ + pdthick/2;
      //G4Box *pd_box  = new G4Box("PhotoDetector", _FAKE_HRPPD_ACTIVE_AREA_SIZE_/2,
      //			     _FAKE_HRPPD_ACTIVE_AREA_SIZE_/2, pdthick/2);
      auto pd_box  = Box(/*"PhotoDetector",*/ _FAKE_HRPPD_ACTIVE_AREA_SIZE_/2,
			 _FAKE_HRPPD_ACTIVE_AREA_SIZE_/2, pdthick/2);
      Volume pssVol(detName + "_photodetector", pd_box, pssMat);//bialkaliMat);
      //gasvolVol.placeVolume(pd_log, Position(0, 0, 0));
      pssVol.setSensitiveDetector(sens);
  
#ifdef _WITH_OPTICS_
      //auto pd = new CherenkovPhotonDetector(pd_box, m_Bialkali);
      // [0,0]: have neither access to G4VSolid nor to G4Material; IRT code does not care; fine;
      auto pd = new CherenkovPhotonDetector(0, 0);
  
      // FIXME: '0' stands for the unknown (and irrelevant) G4LogicalVolume;
      //+geometry->AddPhotonDetector(detector, 0, pd);
      geometry->AddPhotonDetector(cdet, 0, pd);
#endif
	
      // Full size quartz window; FIXME: Sapphire, here and in all other places;
      //auto wnd_box = new G4Box("QuartzWindow", _FAKE_HRPPD_TILE_SIZE_/2, _FAKE_HRPPD_TILE_SIZE_/2,
      //			   _FAKE_HRPPD_WINDOW_THICKNESS_/2);
      auto wnd_box = Box(/*"QuartzWindow",*/ _FAKE_HRPPD_TILE_SIZE_/2, _FAKE_HRPPD_TILE_SIZE_/2,
			 _FAKE_HRPPD_WINDOW_THICKNESS_/2);

      //auto wnd_log =
      //new G4LogicalVolume(wnd_box, /*_HRPPD_WINDOW_MATERIAL_*/m_FusedSilica,  "QuartzWindow", 0, 0, 0);
      Volume wnd_log(detName + "_wnd", wnd_box, quartzMat);///*_HRPPD_WINDOW_MATERIAL_*/m_FusedSilica,  "QuartzWindow", 0, 0, 0);
      //Volume gasvolVol(detName + "_gas", gasvolSolid, gasvolMat);
      
      //+auto hrppd_log = BuildFakeHRPPD(wnd_log, pd_box, pd);
      //#if _LATER_
#ifdef _WITH_OPTICS_
      {	
	TVector3 nx(1,0,0), ny(0,-1,0);
	
	auto surface = 
	  new FlatSurface((1/mm)*TVector3(0,0,fvzOffset + wzOffset + _FAKE_HRPPD_WINDOW_THICKNESS_/2), nx, ny);
	
	geometry->AddFlatRadiator(cdet, "QuartzWindow", CherenkovDetector::Downstream, 
				  0, (G4LogicalVolume*)(0x3), 0, /*wnd_log, m_FusedSilica,*/ surface, _FAKE_HRPPD_WINDOW_THICKNESS_/mm);
      }	
#endif
      
      //?m_Geometry->AddPhotonDetector(cdet, pd->GetLogicalVolume(), pd);
      {
	//double pdthick = 1*cm;//0.01*mm;
	
	auto hrppd_box =
	  Box(_FAKE_HRPPD_TILE_SIZE_/2, _FAKE_HRPPD_TILE_SIZE_/2, _FAKE_HRPPD_CONTAINER_VOLUME_HEIGHT_/2);
	// Well, let it be vacuum inside;
	Volume hrppd_log(detName + "_hrppd", hrppd_box, gasvolMat);//vacuumMat);
	
#ifdef _WITH_OPTICS_
	pd->SetCopyIdentifierLevel(1);
	//?pd->DefineLogicalVolume();
	//pd->SetColor(G4Colour(1, 0, 0, 1.0));
	// Cannot access GEANT shapes in the reconstruction code -> store this value;
	pd->SetActiveAreaSize(_FAKE_HRPPD_ACTIVE_AREA_SIZE_);
#endif
	
#if _LATER_
	auto qd_box  =
	  Box(_FAKE_HRPPD_ACTIVE_AREA_SIZE_/2, _FAKE_HRPPD_ACTIVE_AREA_SIZE_/2, pdthick/2);
	Volume qd_log(detName + "_fake_photodetector", qd_box, bialkaliMat);
#endif
	
	{
	  double accu = -_FAKE_HRPPD_CONTAINER_VOLUME_HEIGHT_/2;
	  
#if 1//_LATER_
	  // Window layer;
	  //+auto wnd_phys = 
	  //new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, accu + _FAKE_HRPPD_WINDOW_THICKNESS_/2), wnd_log, 
	  //		    "QuartzWindow", hrppd_log, false, 0);
	  //+hrppd_log.placeVolume(wnd_log, Position(0, 0, accu + _FAKE_HRPPD_WINDOW_THICKNESS_/2));
	  //+accu += _FAKE_HRPPD_WINDOW_THICKNESS_;
	  
	  //#if _TODAY_
	  // Fake photodector layer (photocathode);
	  //new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, accu + pdthick/2), qd_log, "FakePhotoDetector", 
	  //		    hrppd_log, false, 0);
	  //+hrppd_log.placeVolume(qd_log, Position(0, 0, accu + pdthick/2));
	  //+accu += pdthick;
#endif
	  
	  // Photodector layer (photocathode);
	  //new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, accu + pdthick/2), pd->GetLogicalVolume(), "PhotoDetector", 
	  //		    hrppd_log, false, 0);
	  auto pssPlacement = Position(0, 0, accu + pdthick/2);
	  auto pssPV = hrppd_log.placeVolume(pssVol, pssPlacement);//Position(0, 0, accu + pdthick/2));
	  //gasvolVol.placeVolume(pd_log, Position(0, 0, 0));
	  accu += pdthick;
	  //#endif
	  {
	    double pssThickness = pdthick;
	    double resinSide = _FAKE_HRPPD_ACTIVE_AREA_SIZE_ + 1*mm, resinThickness = pssThickness + 1*mm;
	    Box resinSolid(resinSide / 2., resinSide / 2., resinThickness / 2.);

	    // embed pss solid in resin solid, by subtracting `pssSolid` from `resinSolid`
#if 1
	    //SubtractionSolid resinSolidEmbedded(
	    //					resinSolid, pssVol,//Solid,
	    //					Transform3D(Translation3D(0., 0., (resinThickness - pssThickness) / 2.)));
	    auto resinPlacement = Position(0, 0, accu + resinThickness/2);
	    Volume resinVol("_resin_", resinSolid/*Embedded*/, resinMat);
	    hrppd_log.placeVolume(resinVol, resinPlacement);
#endif
	  }
	  {
	    // sensor readout // NOTE: follow `sensorIDfields`
#if _TODAY_
	    pssPV.addPhysVolID("sector", isec)
	      .addPhysVolID("pdu", ipdu)
	      .addPhysVolID("sipm", isipm);
#endif
	    pssPV.addPhysVolID("pdu", 0);//0);//isec)
	    
	    // sensor DetElement
#if 1//_TODAY_
	    auto sensorID = encodeSensorID(pssPV.volIDs());
	    //std::string sensorIDname =
	    //secName + "_pdu" + std::to_string(ipdu) + "_sipm" + std::to_string(isipm);
	    DetElement pssDE(det, "sensor_de_"/* + sensorIDname*/, sensorID);
	    pssDE.setPlacement(pssPV);
#endif
	    {
	      SkinSurface pssSkin(desc, pssDE, "sensor_optical_surface_"/* + sensorIDname*/, pssSurf,
				  pssVol);
	      pssSkin.isValid();
	    }
	  }
	}

#if _THINK_
	// Is this stuff really needed here?;
	{                      
	  const unsigned qeEntries = 26;
	  
	  // Create HRPPD QE table; use LAPPD #126 from Alexey's March 2022 LAPPD Workshop presentation;
	  double WL[qeEntries] = { 160,  180,  200,  220,  240,  260,  280,  300,  320,  340,  360,  380,  400,  
				   420,  440,  460,  480,  500,  520,  540,  560,  580,  600,  620,  640,  660};
	  double QE[qeEntries] = {0.25, 0.26, 0.27, 0.30, 0.32, 0.35, 0.36, 0.36, 0.36, 0.36, 0.37, 0.35, 0.30, 
				  0.27, 0.24, 0.20, 0.18, 0.15, 0.13, 0.11, 0.10, 0.09, 0.08, 0.07, 0.05, 0.05};  
	  
	  double qemax = 0.0, qePhotonEnergy[qeEntries], qeData[qeEntries];
	  for(int iq=0; iq<qeEntries; iq++) {
	    qePhotonEnergy[iq] = eV * _MAGIC_CFF_ / (WL[qeEntries - iq - 1] + 0.0);
	    qeData        [iq] =                     QE[qeEntries - iq - 1] * _FAKE_QE_DOWNSCALING_FACTOR_;
	    
	    if (qeData[iq] > qemax) qemax = qeData[iq];
	  } //for iq
	  
	  pd->SetQE(eV * _MAGIC_CFF_ / WL[qeEntries-1], eV * _MAGIC_CFF_ / WL[0], 
		    // NB: last argument: want a built-in selection of unused photons, which follow the QE(lambda);
		    // see CherenkovSteppingAction::UserSteppingAction() for a usage case;
		    new G4DataInterpolation(qePhotonEnergy, qeData, qeEntries, 0.0, 0.0), qemax ? 1.0/qemax : 1.0);
	}
#endif
	pd->SetGeometricEfficiency(_FAKE_SENSOR_PLANE_GEOMETRIC_EFFICIENCY_ * _FAKE_SAFETY_FACTOR_);
   
	//new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, wzOffset + _FAKE_HRPPD_CONTAINER_VOLUME_HEIGHT_/2), 
	//		hrppd_log, "HRPPD", dbox->m_fiducial_volume_phys->GetLogicalVolume(), false, 0);
	gasvolVol.placeVolume(hrppd_log, Position(0, 0, wzOffset + _FAKE_HRPPD_CONTAINER_VOLUME_HEIGHT_/2));
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
    
#if _OLD_
    {
      //int isec = 0;
      std::string secName = "sec";// + std::to_string(isec);
      //int ipdu = 0;
      
      Box pssSolid(198*cm, 198*cm, 1*cm);
      // pss volume;
      Volume pssVol(detName + "_pss_" + secName, pssSolid, pssMat);
    
      // place PDU assembly
      //+gasvolVol.placeVolume(pduAssembly, pduAssemblyPlacement);
      auto pssPlacement   = Position(0, 0, 0);
      gasvolVol.placeVolume(pssVol, pssPlacement);
      // sensitivity
      pssVol.setSensitiveDetector(sens);
    }
#endif
  }
  
#if _TODAY_
  auto pduAssemblyPlacement = Position(0, 0, 0);
  Assembly pduAssembly(detName + "_pdu_" + secName);
  //int isipm                 = 0;
  Assembly sensorAssembly(detName + "_sensor_" + secName);
  // placement transformations
  // - placement of objects in `sensorAssembly`
#endif
#if _TODAY_
  auto sensorAssemblyPlacement = Position(0, 0, 0);
  //+auto pssPV =
  sensorAssembly.placeVolume(pssVol, pssPlacement);
  pduAssembly.placeVolume(sensorAssembly, sensorAssemblyPlacement);
#endif
  
  // sensor readout // NOTE: follow `sensorIDfields`
#if _TODAY_
  pssPV.addPhysVolID("sector", isec)
    .addPhysVolID("pdu", ipdu)
    .addPhysVolID("sipm", isipm);
#endif
  
  // sensor DetElement
#if _TODAY_
  auto sensorID = encodeSensorID(pssPV.volIDs());
  std::string sensorIDname =
    secName + "_pdu" + std::to_string(ipdu) + "_sipm" + std::to_string(isipm);
  DetElement pssDE(det, "sensor_de_" + sensorIDname, sensorID);
  pssDE.setPlacement(pssPV);
#endif
  
  // sensor surface properties
#if _TODAY_
  SkinSurface pssSkin(desc, pssDE, "sensor_optical_surface_" + sensorIDname, pssSurf,
  		      pssVol);
  pssSkin.isValid();
#endif 

  // Write out optical configuration;
#ifdef _WITH_OPTICS_
  geometry->Write();//"qrich.root");
  output_file->Close();
#endif
  
  printf("@@@ Ku-end!\n");
  
  return det;
}

DECLARE_DETELEMENT(epic_QRICH, createDetector)
