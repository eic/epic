
#include <TFile.h>

#include "DD4hep/DetFactoryHelper.h"

#include <XML/Helper.h>

// Mimic ePIC dRICH vessel length;
#define _FIDUCIAL_VOLUME_LENGTH_             (1400.0*mm)
// Choose ePIC dRICH location in the forward endcap;
#define _FIDUCIAL_VOLUME_OFFSET_             (3200.0*mm)

#define _MIRROR_RADIUS_                      (2750.0*mm)

// Use HRPPD-like photosensor flat wall;
#define _HRPPD_MATRIX_DIM_                          (19)
// FIXME: may want to use qrich.xml where 32 and 3.25mm are also defined;
#define _HRPPD_PITCH_                          (3.25*mm)
#define _HRPPD_ACTIVE_AREA_SIZE_      (32*_HRPPD_PITCH_)
#define _HRPPD_TILE_SIZE_                     (120.0*mm)
#define _HRPPD_INSTALLATION_GAP_                (3.0*mm)
#define _HRPPD_INSTALLATION_PITCH_ (_HRPPD_TILE_SIZE_ + _HRPPD_INSTALLATION_GAP_)
#define _FIDUCIAL_VOLUME_WIDTH_ (_HRPPD_TILE_SIZE_*_HRPPD_MATRIX_DIM_ + \
				 _HRPPD_INSTALLATION_GAP_*(_HRPPD_MATRIX_DIM_+1) + 50.0*mm)

#define _SENSOR_AREA_LENGTH_                   (50.0*mm)

#define _BUILDING_BLOCK_CLEARANCE_              (1.0*mm)

#define _AEROGEL_THICKNESS_                    (40.0*mm)
// Roughly match eta=3.0 at a distance of 3200-1400/2=250mm;
#define _AEROGEL_VOFFSET_                     (250.0*mm)
#define _ACRYLIC_THICKNESS_                     (3.0*mm)

#define _HRPPD_WINDOW_THICKNESS_                (5.0*mm)
#define _HRPPD_CONTAINER_VOLUME_HEIGHT_        (32.0*mm)

#define _SENSOR_PLANE_GEOMETRIC_EFFICIENCY_       (1.00)
//#define _SAFETY_FACTOR_                           (1.00)

using namespace dd4hep;

#ifdef _WITH_IRT_OPTICS_
#include "IRT/CherenkovDetectorCollection.h"
#include "IRT/SphericalSurface.h"
#endif

static Ref_t createDetector(Detector& desc, xml::Handle_t handle, SensitiveDetector sens)
{
  xml::DetElement detElem       = handle;
  std::string detName           = detElem.nameStr();
  int detID                     = detElem.id();
  DetElement det(detName, detID);

  OpticalSurfaceManager surfMgr = desc.surfaceManager();
  
  // Whatever this means;
  sens.setType("tracker");

  // Start optical configuration if needed;
#ifdef _WITH_IRT_OPTICS_
  // Output file with optics description;
  auto fname = detElem.attr<std::string>(_Unicode(optics));
  
  auto output_file = TFile::Open(fname.c_str(), "RECREATE");
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
  //auto acrylicMaterialName = detElem.attr<std::string>(_Unicode(acrylic));
  //auto acrylicMaterial     = desc.material(acrylicMaterialName);
  
  // HRPPD window;
  auto windowMaterialName  = detElem.attr<std::string>(_Unicode(window));
  auto windowMaterial      = desc.material(windowMaterialName);
  
  // Photocathode;
  auto pcMaterialName      = detElem.attr<std::string>(_Unicode(photocathode));
  auto pcMaterial          = desc.material(pcMaterialName);

  // A fake non-optical layer behind the photocathode;
  auto absorberMaterial    = desc.material(detElem.attr<std::string>(_Unicode(absorber)));

  // Mirror surface description;
  auto mirrorSurf          = surfMgr.opticalSurface(detElem.attr<std::string>(_Unicode(mirror)));
  
  // Sensor readout
  auto readoutName         = detElem.attr<std::string>(_Unicode(readout));

  // Readout coder <-> unique sensor ID
  std::vector<std::string> sensorIDfields = {"hrppd"};
  const auto& readoutCoder = *desc.readout(readoutName).idSpec().decoder();
  uint64_t cellMask = 0;
  for (const auto& idField : sensorIDfields)
    cellMask |= readoutCoder[idField].mask();
  desc.add(Constant("FRICH_cell_mask", std::to_string(cellMask)));
#ifdef _WITH_IRT_OPTICS_
  // Do not mind to store it twice;
  cdet->SetReadoutCellMask(cellMask);
#endif
  // create a unique sensor ID from a sensor's PlacedVolume::volIDs
  auto encodeSensorID = [&readoutCoder](auto ids) {
    uint64_t enc = 0;
    for (const auto& [idField, idValue] : ids)
      enc |= uint64_t(idValue) << readoutCoder[idField].offset();
    return enc;
  };

  const double fvLength = _FIDUCIAL_VOLUME_LENGTH_, fvWidth  = _FIDUCIAL_VOLUME_WIDTH_;
  const double fvOffset = _FIDUCIAL_VOLUME_OFFSET_;
  
  // Fiducial volume / vessel;
  auto vesselSolid = Box(fvWidth/2, fvWidth/2, fvLength/2);
  Volume vesselVolume(detName, vesselSolid, vesselMaterial);

  // Place vessel into the mother volume;
  Volume motherVol = desc.pickMotherVolume(det);
  PlacedVolume vesselPV = motherVol.placeVolume(vesselVolume,  Position(0, 0, fvOffset));
#ifdef _WITH_IRT_OPTICS_
  const double sign = 1.0;//flip ? -1.0 : 1.0;
#endif
  vesselPV.addPhysVolID("system", detID);
  det.setPlacement(vesselPV);

#ifdef _WITH_IRT_OPTICS_
  SphericalSurface *msurface = 0;
  OpticalBoundary *mboundary = 0;
#endif
    
  // Gas volume, aerogel, acrylic filter, mirrors;
  {
    const double gvWidth = fvWidth - 1*mm;
    const double gvLength = fvLength - _SENSOR_AREA_LENGTH_;
    auto gasvolSolid = Box(gvWidth/2, gvWidth/2, gvLength/2);
    Volume gasvolVol(detName + "_gas", gasvolSolid, gasvolMaterial);
    
    // Place gas volume into the vessel;
    double gvOffset = _SENSOR_AREA_LENGTH_/2;
    PlacedVolume gasvolPV = vesselVolume.placeVolume(gasvolVol, Position(0, 0, gvOffset));
    DetElement gasvolDE(det, "gasvol_de", 0);
    gasvolDE.setPlacement(gasvolPV);

#ifdef _WITH_IRT_OPTICS_
    {
      // FIXME: Z-location does not really matter here, right?;
      auto boundary = new FlatSurface(TVector3(0,0,0), sign*TVector3(1,0,0), TVector3(0,-1,0));
      
      auto radiator = geometry->SetContainerVolume(cdet, "GasVolume", 0, (G4LogicalVolume*)(0x0), 0, boundary);
      radiator->SetAlternativeMaterialName(gasvolMaterialName.c_str());
      // FIXME: what is it good for in ePIC IRT 2.0 implementation?;
      geometry->AddRadiatorLogicalVolume(radiator, (G4LogicalVolume*)(0x0));
    }
#endif
    
    {
      double gzOffset = -gvLength/2 + _BUILDING_BLOCK_CLEARANCE_;
    
      // Aerogel volume;
      {
	// FIXME: for now just tweak this size such that reflected light at {eta=2.5,phi=90}
	// does not pass through this aerogel layer twice; 
	const double agWidth = 5*cm, agThick = _AEROGEL_THICKNESS_;
	Box aerogelSolid(agWidth/2, agWidth/2, agThick/2);
	Volume aerogelVol(detName + "_aerogel", aerogelSolid, aerogelMaterial);
	
	// Place aerogel into gas volume; NB: 'gzOffset' used in a FlatSurface() description below; 
	gzOffset += agThick/2;
	gasvolVol.placeVolume(aerogelVol, Position(0, _AEROGEL_VOFFSET_, gzOffset));
	
#ifdef _WITH_IRT_OPTICS_
	{
	  TVector3 nx(1*sign,0,0), ny(0,-1,0);
	
	  auto surface = new FlatSurface(sign*(1/mm)*TVector3(0,0,fvOffset + gvOffset + gzOffset), nx, ny);
	  
	  auto radiator = geometry->AddFlatRadiator(cdet, "Aerogel", CherenkovDetector::Upstream, 
						      0, (G4LogicalVolume*)(0x1), 0, surface, agThick/mm);
	  radiator->SetAlternativeMaterialName(aerogelMaterialName.c_str());
	  // FIXME: what is it good for in ePIC IRT 2.0 implementation?;
	  geometry->AddRadiatorLogicalVolume(radiator, (G4LogicalVolume*)(0x1));
	}
#endif
    
	gzOffset += agThick/2 + _BUILDING_BLOCK_CLEARANCE_;
      }
      
      // Acrylic filter volume;
#if _LATER_
      {
	const double acWidth = 5*cm, acThick = _ACRYLIC_THICKNESS_;
	Box acrylicSolid(acWidth/2, acWidth/2, acThick/2);
	Volume acrylicVol(detName + "_filter", acrylicSolid, acrylicMaterial);
	
	// Place acrylic filter into gas volume; NB: 'gzOffset' used in a FlatSurface() description below;
	gzOffset += acThick/2;
	gasvolVol.placeVolume(acrylicVol, Position(0, _AEROGEL_VOFFSET_, gzOffset));
	
#ifdef _WITH_IRT_OPTICS_
	{
	  TVector3 nx(1*sign,0,0), ny(0,-1,0);
	
	  auto surface = new FlatSurface(sign*(1/mm)*TVector3(0,0,fvOffset + gvOffset + gzOffset), nx, ny);
	  
	  auto radiator = geometry->AddFlatRadiator(cdet, "Acrylic", CherenkovDetector::Upstream, 
						    0, (G4LogicalVolume*)(0x2), 0, surface, acThick/mm);
	  radiator->SetAlternativeMaterialName(acrylicMaterialName.c_str());
	  // FIXME: what is it good for in ePIC IRT 2.0 implementation?;
	  geometry->AddRadiatorLogicalVolume(radiator, (G4LogicalVolume*)(0x2));
	}
#endif
      }
#endif
    }

    {
      const double mifvWidth = gvWidth - 1*mm, mifvLength = 500*mm, miThickness = 1*mm;
      double miRadius = _MIRROR_RADIUS_;//2000*mm;//3000*mm;
      Box mirrorFvSolid(mifvWidth/2, mifvWidth/2, mifvLength/2);
      Sphere mirrorSphere(miRadius - miThickness, miRadius, 0.0, M_PI, -M_PI, M_PI);
      IntersectionSolid mirrorSolid(mirrorFvSolid, mirrorSphere, Position(0, 0, mifvLength/2 - miRadius));

      // FIXME: hardcoded to a semi-transparent pfRICH attribute;
      auto mirrorVis  = desc.visAttributes("PFRICH_mirror_vis");
  
      // mirror volume, attributes, and placement
      Volume mirrorVol(detName + "_mirror_", mirrorSolid, vesselMaterial);
      mirrorVol.setVisAttributes(mirrorVis);
      auto mirrorPV = gasvolVol.placeVolume(mirrorVol, Position(0,0,gvLength/2 - mifvLength/2));
      
      // properties
      DetElement mirrorDE(det, "mirror_de_", 0);//isec);
      mirrorDE.setPlacement(mirrorPV);
      SkinSurface mirrorSkin(desc, mirrorDE, "mirror_optical_surface_", mirrorSurf, mirrorVol);
      mirrorSkin.isValid();

#ifdef _WITH_IRT_OPTICS_
      // NB: default is concave, which is fine;
      msurface = new SphericalSurface(sign*(1/mm)*TVector3(0,0,fvOffset + gvOffset + gvLength/2 - miRadius), (1/mm)*miRadius);
      mboundary = new OpticalBoundary(cdet->GetRadiator("GasVolume"), msurface, false);
      // Need to store it in a separate call, see a comment in CherenkovDetector.h;
      cdet->StoreOpticalBoundary(mboundary);

      // Assign gas volume rear surface (this mirror) by hand;
      //printf("@B@ createDetector() %p\n", msurface);
      cdet->GetContainerVolume()->m_Borders[0].second = msurface;//mboundary;
      //printf("@B@ createDetector() (%p) (%p)\n", cdet->GetContainerVolume()->GetRearSide(0),
      //     cdet->GetContainerVolume()->m_Borders[0].second.GetObject());//, name);
#endif
    }
  }

  // A fake HRPPD container volume, window, photocathode;
  {
#ifdef _WITH_IRT_OPTICS_
    // [0,0]: have neither access to G4VSolid nor to G4Material; IRT code does not care; fine;
    auto pd = new CherenkovPhotonDetector(0, 0);
    
    // FIXME: '0' stands for the unknown (and irrelevant) G4LogicalVolume;
    geometry->AddPhotonDetector(cdet, 0, pd);
    
    // Cannot access GEANT shapes in the reconstruction code -> store this value;
    pd->SetActiveAreaSize(_HRPPD_ACTIVE_AREA_SIZE_/mm);
    
    pd->SetGeometricEfficiency(_SENSOR_PLANE_GEOMETRIC_EFFICIENCY_);// * _SAFETY_FACTOR_);
#endif
    
    // Mimic BuildFakePhotonDetectorMatrix() used in the standalone code;
    const double fvzOffset = fvOffset, _wzOffset = -fvLength/2 + _SENSOR_AREA_LENGTH_;
    const double pdthick = 0.01*mm, _zpdc = _wzOffset - _HRPPD_WINDOW_THICKNESS_ - pdthick/2;
    
    // For now assume it is a unique surface, same for all HRPPDs;
#ifdef _WITH_IRT_OPTICS_
    {	
      TVector3 nx(-1*sign,0,0), ny(0,-1,0);
      
      auto surface = 
	new FlatSurface(sign*(1/mm)*TVector3(0,0,fvzOffset + _wzOffset - _HRPPD_WINDOW_THICKNESS_/2), nx, ny);
      
      auto radiator = geometry->AddFlatRadiator(cdet, "QuartzWindow", CherenkovDetector::Downstream, 
						0, (G4LogicalVolume*)(0x3), 0, surface,
						_HRPPD_WINDOW_THICKNESS_/mm);
      radiator->SetAlternativeMaterialName(windowMaterialName.c_str());
      // FIXME: what is it good for in ePIC IRT 2.0 implementation?;
      geometry->AddRadiatorLogicalVolume(radiator, (G4LogicalVolume*)(0x3));
    }	
#endif
        
    for(unsigned ix=0; ix<_HRPPD_MATRIX_DIM_; ix++) {
      double xOffset = _HRPPD_INSTALLATION_PITCH_*(ix - (_HRPPD_MATRIX_DIM_-1)/2.);
      
      for(unsigned iy=0; iy<_HRPPD_MATRIX_DIM_; iy++) {
	unsigned id = iy*_HRPPD_MATRIX_DIM_ + ix;
	double yOffset = _HRPPD_INSTALLATION_PITCH_*(iy - (_HRPPD_MATRIX_DIM_-1)/2.);

	uint64_t sensorID = 0x0;
	
	auto hrppdBox =
	  Box(_HRPPD_TILE_SIZE_/2, _HRPPD_TILE_SIZE_/2, _HRPPD_CONTAINER_VOLUME_HEIGHT_/2);
	// FIXME: think (OLD: well, let it be vacuum inside);
	TString hrppdName; hrppdName.Form("%s-hrppd-%02d%02d", detName.c_str(), ix, iy);
	Volume hrppdVol(hrppdName.Data(), hrppdBox, gasvolMaterial);//vacuumMat);
	
	// Full size quartz window; 
	auto wndBox = Box(_HRPPD_TILE_SIZE_/2, _HRPPD_TILE_SIZE_/2, _HRPPD_WINDOW_THICKNESS_/2);
	TString wndName; wndName.Form("%s-window-%02d%02d", detName.c_str(), ix, iy);
	Volume wndVol(wndName.Data(), wndBox, windowMaterial);
	
	auto pcBox  = Box(_HRPPD_ACTIVE_AREA_SIZE_/2,
			  _HRPPD_ACTIVE_AREA_SIZE_/2, pdthick/2);
	TString pcName; pcName.Form("%s-photocathode-%02d%02d", detName.c_str(), ix, iy);
	Volume pcVol(pcName.Data(), pcBox, pcMaterial);
	pcVol.setSensitiveDetector(sens);
	auto pcVis  = desc.visAttributes("FRICH_sensor_vis");
	pcVol.setVisAttributes(pcVis);
	
	{
	  double accu = -_HRPPD_CONTAINER_VOLUME_HEIGHT_/2;
	  
	  // Window layer;
	  hrppdVol.placeVolume(wndVol, Position(0.0, 0.0, accu + _HRPPD_WINDOW_THICKNESS_/2));
	  accu += _HRPPD_WINDOW_THICKNESS_;
	  
	  // Photocathode layer;
	  {
	    auto pcPV = hrppdVol.placeVolume(pcVol, Position(0.0, 0.0, accu + pdthick/2));
	    // sensor readout // NOTE: follow `sensorIDfields`
	    pcPV.addPhysVolID("hrppd", id);
	    
	    // sensor DetElement
	    sensorID = encodeSensorID(pcPV.volIDs());
	    TString deName; deName.Form("%s-sensor-%02d%02d", detName.c_str(), ix, iy);
	    DetElement pcDE(det, deName.Data(), sensorID);
	    pcDE.setPlacement(pcPV);
	  }
	  accu += pdthick;
	  
	  // A fake absorber layer behind the photocathode; FIXME: make sure that reflection
	  // on the window and photocathode boundary still works correctly (no fake volume as
	  // in a standalone code);
	  {
	    TString absName; absName.Form("%s-absorber-%02d%02d", detName.c_str(), ix, iy);
	    // Recycle the same pcBox shape;
	    Volume absVol(absName.Data(), pcBox, absorberMaterial);
	    hrppdVol.placeVolume(absVol, Position(0.0, 0.0, accu + pdthick/2));
	  }
	}

	{
	  Rotation3D rot(RotationZYX(0, M_PI, 0));
	  Transform3D transform(rot,  Position(xOffset, yOffset, _wzOffset - _HRPPD_CONTAINER_VOLUME_HEIGHT_/2));
	  vesselVolume.placeVolume(hrppdVol, transform);
	}
	
#ifdef _WITH_IRT_OPTICS_
	{
	  // Photocathode surface;
	  auto surface = new FlatSurface((1/mm)*TVector3(sign*xOffset, yOffset, sign*(fvzOffset + _zpdc)),
					 TVector3(-1*sign,0,0), TVector3(0,-1,0));

	  {
	    // Let it be a single sector;
	    unsigned sector = 0;
	    
	    auto irt = pd->AllocateIRT(sector, sensorID);
	      
	    // Aerogel and acrylic;
	    if (cdet->m_OpticalBoundaries[CherenkovDetector::Upstream].find(sector) != 
		cdet->m_OpticalBoundaries[CherenkovDetector::Upstream].end())
	      for(auto boundary: cdet->m_OpticalBoundaries[CherenkovDetector::Upstream][sector])
		irt->AddOpticalBoundary(boundary);

	    // Mirror;
	    irt->AddOpticalBoundary(mboundary);
	      
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
      } //for iy
    } //for ix
  }
      
  // Write out optical configuration;
#ifdef _WITH_IRT_OPTICS_
  geometry->Write();
  output_file->Close();
#endif
  
  return det;
}

DECLARE_DETELEMENT(epic_FRICH, createDetector)
