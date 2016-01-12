/**
 *  Author:
 *    
 *
 *  allpix Authors:
 *   John Idarraga <idarraga@cern.ch>
 *   Mathieu Benoit <benoit@lal.in2p3.fr>
 */

#include "AllPixMCTruthDigitizer.hh"
#include "AllPixTrackerHit.hh"
#include "G4DigiManager.hh"
#include "G4VDigitizerModule.hh"
#include "AllPixGeoDsc.hh"

#include "TMath.h"
#include "TF1.h"

AllPixMCTruthDigitizer::AllPixMCTruthDigitizer(G4String modName, G4String hitsColName, G4String digitColName) 
  : AllPixDigitizerInterface (modName) {

  // Registration of digits collection name
  collectionName.push_back(digitColName);
  m_hitsColName.push_back(hitsColName);

  // threshold
  m_digitIn.thl = 0.;

}

AllPixMCTruthDigitizer::~AllPixMCTruthDigitizer(){

}

void AllPixMCTruthDigitizer::Digitize(){

  // create the digits collection
  m_digitsCollection = new AllPixMCTruthDigitsCollection("AllPixMCTruthDigitizer", collectionName[0] );

  // get the digiManager
  G4DigiManager * digiMan = G4DigiManager::GetDMpointer();

  // BoxSD_0_HitsCollection
  G4int hcID = digiMan->GetHitsCollectionID(m_hitsColName[0]);

  AllPixTrackerHitsCollection * hitsCollection = 0;
  hitsCollection = (AllPixTrackerHitsCollection*)(digiMan->GetHitsCollection(hcID));

  // temporary data structure
  map<pair<G4int, G4int>, G4double > pixelsContent;
  pair<G4int, G4int> tempPixel;

  G4int nEntries = hitsCollection->entries();
  G4cout << "nEntries=" << nEntries << G4endl;
  
  // Example of detector description handle
  // provided by the interface
  AllPixGeoDsc * gD = GetDetectorGeoDscPtr();
  gD->GetNPixelsX();
  Double_t pitchX=gD->GetPixelX();
  Double_t pitchY=gD->GetPixelY();
  
  G4double MC_deposited_energy=0.0;

  int pixelX=-10;
  int pixelY=-10;
  
  for(G4int itr  = 0 ; itr < nEntries ; itr++) 
    {
      
      //  tempPixel.first  = (*hitsCollection)[itr]->GetPixelNbX();
      //tempPixel.second = (*hitsCollection)[itr]->GetPixelNbY();
      //pixelsContent[tempPixel] += (*hitsCollection)[itr]->GetEdep();
      Double_t tempEnergy=(*hitsCollection)[itr]->GetEdep();
      MC_deposited_energy+=(*hitsCollection)[itr]->GetEdep();
      
      G4double xpos=(*hitsCollection)[itr]->GetPosWithRespectToPixel().x();
      G4double ypos=(*hitsCollection)[itr]->GetPosWithRespectToPixel().y();
      G4double zpos=(*hitsCollection)[itr]->GetPosWithRespectToPixel().z(); // [mm]; zpos=thickness corresponds to the sensor side and zpos=0 corresponds to the pixel side     
      
      xpos=(*hitsCollection)[itr]->GetPixelNbX()*pitchX+xpos+pitchX/2.0;
      ypos=(*hitsCollection)[itr]->GetPixelNbY()*pitchY+ypos+pitchY/2.0;
      
      tempPixel.first=xpos;
      tempPixel.second=zpos;
      pixelsContent[tempPixel] += (*hitsCollection)[itr]->GetEdep();
      
      // G4cout << "xpos=" << xpos << ", zpos=" << zpos << G4endl;
      // G4cout << "tempEnergy=" << tempEnergy << ", zpos=" << zpos << G4endl;
      //G4cout << "x=" << (*hitsCollection)[itr]->GetPixelNbX() << ", y=" << (*hitsCollection)[itr]->GetPixelNbY() << G4endl;
      pixelX=(*hitsCollection)[itr]->GetPixelNbX();
      pixelY=(*hitsCollection)[itr]->GetPixelNbY();
    }
  

  G4cout << "x=" << pixelX << ", y=" << pixelY << G4endl;

  G4cout << "MC_deposited_energy=" << MC_deposited_energy/keV << G4endl;
  // Now create digits, one per pixel
  // Second entry in the map is the energy deposit in the pixel
  map<pair<G4int, G4int>, G4double >::iterator pCItr = pixelsContent.begin();
  /*for( ; pCItr != pixelsContent.end() ; pCItr++)
    {
      if(((*pCItr).second)/keV > 0) // over threshold !                                                                                                             
        {
          pair<G4int, G4int> tempPixel;
          tempPixel.first=(*pCItr).first.first;
          tempPixel.second=(*pCItr).first.second;
	  
	}
	}*/
  // NOTE that there is a nice interface which provides useful info for hits.
  // For instance, the following member gives you the position of a hit with respect
  //  to the center of the pixel.
  // G4ThreeVector vec = (*hitsCollection)[itr]->GetPosWithRespectToPixel();
  // See class AllPixTrackerHit !

  // Also, you have full access to the detector geometry from this scope
  // AllPixGeoDsc * GetDetectorGeoDscPtr()
  // provides you with a pointer to the geometry description.
  // See class AllPixGeoDsc !

  // Create one digit per pixel, I need to look at all the pixels first

  if(MC_deposited_energy>0)
    {
      AllPixMCTruthDigit * digit = new AllPixMCTruthDigit;
      digit->SetPixelEnergyDep(MC_deposited_energy/keV);
      digit->SetPixelIDX(pixelX); // An ugly hack to have the same coordinate system as the test-beam
      digit->SetPixelIDY(pixelY);	 
      digit->IncreasePixelCounts(); // Counting mode
      m_digitsCollection->insert(digit);
	
    }

  G4int dc_entries = m_digitsCollection->entries();
  if(dc_entries > 0)
    {
      G4cout << "--------> Digits Collection : " << collectionName[0]
	     << "(" << m_hitsColName[0] << ")"
	     << " contains " << dc_entries
	     << " digits" << G4endl;
    }
	
  StoreDigiCollection(m_digitsCollection);
}
