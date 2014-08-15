/**
 *  Author:
 *    nalipour@cern.ch
 *
 *  allpix Authors:
 *   John Idarraga <idarraga@cern.ch>
 *   Mathieu Benoit <benoit@lal.in2p3.fr>
 */

#include "AllPixTMPXDigitizer.hh"
#include "AllPixTrackerHit.hh"
#include "G4DigiManager.hh"
#include "G4VDigitizerModule.hh"
#include "AllPixGeoDsc.hh"

#include "TMath.h"
#include "TF1.h"

AllPixTMPXDigitizer::AllPixTMPXDigitizer(G4String modName, G4String hitsColName, G4String digitColName) 
  : AllPixDigitizerInterface (modName) {

  // Registration of digits collection name
  collectionName.push_back(digitColName);
  m_hitsColName.push_back(hitsColName);

  // threshold
  m_digitIn.thl = 0.;

}

AllPixTMPXDigitizer::~AllPixTMPXDigitizer(){

}

void AllPixTMPXDigitizer::Digitize()
{
  m_digitsCollection = new AllPixTMPXDigitsCollection("AllPixTMPXDigitizer", collectionName[0] );

  // get the digiManager
  G4DigiManager * digiMan = G4DigiManager::GetDMpointer();

  // BoxSD_0_HitsCollection
  G4int hcID = digiMan->GetHitsCollectionID(m_hitsColName[0]);

  AllPixTrackerHitsCollection * hitsCollection = 0;
  hitsCollection = (AllPixTrackerHitsCollection*)(digiMan->GetHitsCollection(hcID));

  // temporary data structure
  map<pair<G4int, G4int>, MC_content> pixelsContent_MC; //contains information with only MC
  map<pair<G4int, G4int>, G4double > pixelsContent; //contains information with charge sharing
  pair<G4int, G4int> tempPixel;
  G4int nEntries = hitsCollection->entries();

  AllPixGeoDsc * gD = GetDetectorGeoDscPtr();
  gD->GetNPixelsX();

  for(G4int itr  = 0 ; itr < nEntries ; itr++)
    {

      tempPixel.first  = (*hitsCollection)[itr]->GetPixelNbX();
      tempPixel.second = (*hitsCollection)[itr]->GetPixelNbY();
      pixelsContent_MC[tempPixel].MC_energy += (*hitsCollection)[itr]->GetEdep(); //MC
      if(pixelsContent_MC[tempPixel].posX_WithRespectToPixel == -11.0 || pixelsContent_MC[tempPixel].posY_WithRespectToPixel == -11.0)
	{
	  pixelsContent_MC[tempPixel].posX_WithRespectToPixel=(*hitsCollection)[itr]->GetPosWithRespectToPixel().x();
	  pixelsContent_MC[tempPixel].posY_WithRespectToPixel=(*hitsCollection)[itr]->GetPosWithRespectToPixel().y();
	  pixelsContent_MC[tempPixel].posZ_WithRespectToPixel=(*hitsCollection)[itr]->GetPosWithRespectToPixel().z();
	}
      pixelsContent[tempPixel] = pixelsContent_MC[tempPixel].MC_energy; //+ charge sharing + noise
    }

  //------------------ RECORD DIGITS ------------------//
  // With charge sharing
  map<pair<G4int, G4int>, G4double >::iterator pCItr = pixelsContent.begin();
  for( ; pCItr != pixelsContent.end() ; pCItr++)
    {
      if((*pCItr).second > m_digitIn.thl) // over threshold !
	{
	  pair<G4int, G4int> tempPixel;
	  tempPixel.first=(*pCItr).first.first;
	  tempPixel.second=(*pCItr).first.second;

	  AllPixTMPXDigit * digit = new AllPixTMPXDigit;

	  digit->SetPixelIDX((*pCItr).first.first);
	  digit->SetPixelIDY((*pCItr).first.second);
	  digit->SetPixelEnergyDep((*pCItr).second/keV); //Energy with charge sharing
	  digit->SetPixelCounts((*pCItr).second/eV); //TOT value
	  if (pixelsContent_MC.count(tempPixel))
	    {
	      digit->SetPixelEnergyMC(pixelsContent_MC[tempPixel].MC_energy/keV); //MC value
	      digit->Set_posX_WithRespectoToPixel(pixelsContent_MC[tempPixel].posX_WithRespectToPixel);
	      digit->Set_posY_WithRespectoToPixel(pixelsContent_MC[tempPixel].posY_WithRespectToPixel);
	      digit->Set_posZ_WithRespectoToPixel(pixelsContent_MC[tempPixel].posZ_WithRespectToPixel);
	      pixelsContent_MC.erase(tempPixel);
	    }
	  digit->IncreasePixelCounts(); // Counting mode

	  m_digitsCollection->insert(digit);
	}
    }
  //MC only
  map<pair<G4int, G4int>, MC_content >::iterator pCItr_MC = pixelsContent_MC.begin();
  for( ; pCItr_MC != pixelsContent_MC.end() ; pCItr_MC++)
    {
      pair<G4int, G4int> tempPixel;
      tempPixel.first=(*pCItr_MC).first.first;
      tempPixel.second=(*pCItr_MC).first.second;

      AllPixTMPXDigit * digit = new AllPixTMPXDigit;

      digit->SetPixelIDX((*pCItr_MC).first.first);
      digit->SetPixelIDY((*pCItr_MC).first.second);
      digit->SetPixelEnergyMC(pixelsContent_MC[tempPixel].MC_energy/keV); //MC value
      digit->Set_posX_WithRespectoToPixel(pixelsContent_MC[tempPixel].posX_WithRespectToPixel);
      digit->Set_posY_WithRespectoToPixel(pixelsContent_MC[tempPixel].posY_WithRespectToPixel);
      digit->Set_posZ_WithRespectoToPixel(pixelsContent_MC[tempPixel].posZ_WithRespectToPixel);
      //digit->IncreasePixelCounts(); // Counting mode

      m_digitsCollection->insert(digit);
    }
  //----------------------------------------------------//


  G4int dc_entries = m_digitsCollection->entries();
  if(dc_entries > 0){
    G4cout << "--------> Digits Collection : " << collectionName[0]
	   << "(" << m_hitsColName[0] << ")"
	   << " contains " << dc_entries
	   << " digits" << G4endl;
  }

  StoreDigiCollection(m_digitsCollection);
}
