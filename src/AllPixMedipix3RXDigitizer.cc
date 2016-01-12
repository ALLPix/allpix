/**
 *  Author:
 *    idarraga@cern.ch
 *
 *  allpix Authors:
 *   John Idarraga <idarraga@cern.ch>
 *   Mathieu Benoit <benoit@lal.in2p3.fr>
 */

#include "AllPixMedipix3RXDigitizer.hh"
#include "AllPixTrackerHit.hh"
#include "G4DigiManager.hh"
#include "G4VDigitizerModule.hh"
#include "AllPixGeoDsc.hh"

#include "TMath.h"
#include "TF1.h"

AllPixMedipix3RXDigitizer::AllPixMedipix3RXDigitizer(G4String modName, G4String hitsColName, G4String digitColName) 
: AllPixDigitizerInterface (modName) {

	// Registration of digits collection name
	collectionName.push_back(digitColName);
	m_hitsColName.push_back(hitsColName);

	// threshold
	m_digitIn.thl = 0.;

}

AllPixMedipix3RXDigitizer::~AllPixMedipix3RXDigitizer(){

}

void AllPixMedipix3RXDigitizer::Digitize(){

	// number of threholds
	int nThresholds = 8;
	// create a threshold map
	map<int, double> thlMap;
	for ( int thlI = 0 ; thlI < nThresholds ; thlI++) {
		thlMap[thlI] = 2*keV + 4*thlI*keV; // ex. starting at 2keV in steps if 4keV
	}

	// create the digits collection
	m_digitsCollection = new AllPixMedipix3RXDigitsCollection("AllPixMedipix3RXDigitizer", collectionName[0] );

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

	// Example of detector description handle
	// provided by the interface
	AllPixGeoDsc * gD = GetDetectorGeoDscPtr();
	gD->GetNPixelsX();

	G4double kinEParent;
	for(G4int itr  = 0 ; itr < nEntries ; itr++) {

		tempPixel.first  = (*hitsCollection)[itr]->GetPixelNbX();
		tempPixel.second = (*hitsCollection)[itr]->GetPixelNbY();
		pixelsContent[tempPixel] += (*hitsCollection)[itr]->GetEdep();
		kinEParent = (*hitsCollection)[itr]->GetKinEParent();

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
		// Create one digit per pixel, I need to look at all the pixels first

		AllPixMedipix3RXDigit * digit = new AllPixMedipix3RXDigit( nThresholds );

		// Iterate over the thresholds
		for ( int thlI = 0 ; thlI < nThresholds ; thlI++) {

			// Compare to the corresponding threshold level
			if ( (*pCItr).second > thlMap[thlI] ) {
				digit->SetPixelIDX((*pCItr).first.first);
				digit->SetPixelIDY((*pCItr).first.second);
				digit->SetPixelEnergyDep( (*pCItr).second );
				digit->IncreasePixelCountsMultiTHL( thlI );
			}

		}

		// kinEnergy parent.  Once per event.
		digit->SetKinEParent( kinEParent );

		m_digitsCollection->insert(digit);

	}

	//G4int dc_entries = m_digitsCollection->entries();
	//if ( dc_entries > 0 ) {
	//	G4cout << "--------> Digits Collection : " << collectionName[0]
	//	                                                             << "(" << m_hitsColName[0] << ")"
	//	                                                             << " contains " << dc_entries
	//	                                                             << " digits" << G4endl;
	//}

	StoreDigiCollection(m_digitsCollection);

}
