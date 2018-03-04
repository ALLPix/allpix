/**
 *  Author:
 *    
 *
 *  allpix Authors:
 *   John Idarraga <idarraga@cern.ch>
 *   Mathieu Benoit <benoit@lal.in2p3.fr>
 */

#include "AllPixTestBichselDigitizer.hh"
#include "AllPixTrackerHit.hh"
#include "G4DigiManager.hh"
#include "G4VDigitizerModule.hh"
#include "AllPixGeoDsc.hh"

#include "TMath.h"
#include "TF1.h"

AllPixTestBichselDigitizer::AllPixTestBichselDigitizer(G4String modName, G4String hitsColName, G4String digitColName) 
: AllPixDigitizerInterface (modName) {

  std::cout << "in the init" << std::endl;

	// Registration of digits collection name
	collectionName.push_back(digitColName);
	m_hitsColName.push_back(hitsColName);

	// threshold
	m_digitIn.thl = 0.;

	/***********************  Bichsel  *****************************/
	doBichsel=true;
	G4int n_Col=10;

	AllPixGeoDsc * gD_init = GetDetectorGeoDscPtr();
	
	std::cout << "just about to start the bichsel init" << std::endl;

	Bichsel.Initialize(n_Col, gD_init->GetPixelX(), gD_init->GetPixelY(), gD_init->GetPixelZ());
	/***************************************************************/

	std::cout << "done with the init" << std::endl;
}

AllPixTestBichselDigitizer::~AllPixTestBichselDigitizer(){

}

void AllPixTestBichselDigitizer::Digitize(){

  std::cout << "Here " << std::endl;

	// create the digits collection
	m_digitsCollection = new AllPixTestBichselDigitsCollection("AllPixTestBichselDigitizer", collectionName[0] );

	// get the digiManager
	G4DigiManager * digiMan = G4DigiManager::GetDMpointer();

	// BoxSD_0_HitsCollection
	G4int hcID = digiMan->GetHitsCollectionID(m_hitsColName[0]);

	AllPixTrackerHitsCollection * hitsCollection = 0;
	hitsCollection = (AllPixTrackerHitsCollection*)(digiMan->GetHitsCollection(hcID));

	/***************************************************************/
	//BichselHisCollection:
	bool doBichsel=false;
	std::vector<BichselHits> BichselhitsCollection;

	if(doBichsel)
	  BichselhitsCollection=Bichsel.GetBichselhitsCollection(hitsCollection,false);
	G4int nBichselEntries = BichselhitsCollection.size();
	map<pair<G4int, G4int>, G4double > pixelsContent;
	pair<G4int, G4int> tempBichselPixel;
	G4ThreeVector pos_respect_pixel;
	// for(G4int itr  = 0 ; itr < nEntries ; itr++) {
	//   tempBichselPixel.first=BichselhitsCollection[itr].PixelNb.first;
	//   tempBichselPixel.second=BichselhitsCollection[itr].PixelNb.second;
	//   pixelsContent[tempBichselPixel]+=BichselhitsCollection[itr].energy*1e6;//eV
	//   // pos_respect_pixel=BichselhitsCollection[itr].pos;
	//   // trackID=BichselhitsCollection[itr].trackid;
	// }
	/***************************************************************/




	// temporary data structure
	//map<pair<G4int, G4int>, G4double > pixelsContent;
	pair<G4int, G4int> tempPixel;

	G4int nEntries = hitsCollection->entries();

	// Example of detector description handle
	// provided by the interface
	AllPixGeoDsc * gD = GetDetectorGeoDscPtr();
	gD->GetNPixelsX();

	for(G4int itr  = 0 ; itr < nEntries ; itr++) {

		tempPixel.first  = (*hitsCollection)[itr]->GetPixelNbX();
		tempPixel.second = (*hitsCollection)[itr]->GetPixelNbY();
		pixelsContent[tempPixel] += (*hitsCollection)[itr]->GetEdep();

	}

	// Now create digits, one per pixel
	// Second entry in the map is the energy deposit in the pixel
	map<pair<G4int, G4int>, G4double >::iterator pCItr = pixelsContent.begin();

	// NOTE that there is a nice interface which provides useful info for hits.
	// For instance, the following member gives you the position of a hit with respect
	//  to the center of the pixel.
	// G4ThreeVector vec = (*hitsCollection)[itr]->GetPosWithRespectToPixel();
	// See class AllPixTrackerHit !

	// Also, you have full access to the detector geometry from this scope
	// AllPixGeoDsc * GetDetectorGeoDscPtr()
	// provides you with a pointer to the geometry description.
	// See class AllPixGeoDsc !

	for( ; pCItr != pixelsContent.end() ; pCItr++)
	{

		if((*pCItr).second > 0) // over threshold !
		{
			// Create one digit per pixel, I need to look at all the pixels first
			AllPixTestBichselDigit * digit = new AllPixTestBichselDigit;
			digit->SetPixelIDX((*pCItr).first.first);
			digit->SetPixelIDY((*pCItr).first.second);
			digit->SetPixelCounts((*pCItr).second/keV);

			G4cout << "dEdX : " << (*pCItr).second/keV << endl;

			m_digitsCollection->insert(digit);
		}
	}

	G4int dc_entries = m_digitsCollection->entries();
	if(dc_entries > 0){
		G4cout << "--------> Digits Collection : " << collectionName[0]
		                                           << "(" << m_hitsColName[0] << ")"
		                                           << " contains " << dc_entries
		                                           << " digits" << G4endl;
	}

	StoreDigiCollection(m_digitsCollection);

}
