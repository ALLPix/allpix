/**
 *  Author:
 *    Mathieu.Benoit@CERN.CH
 *
 *  allpix Authors:
 *   John Idarraga <idarraga@cern.ch>
 *   Mathieu Benoit <benoit@lal.in2p3.fr>
 */

#include "AllPixEdepHistogrammerDigitizer.hh"
#include "AllPixTrackerHit.hh"
#include "G4DigiManager.hh"
#include "G4VDigitizerModule.hh"
#include "AllPixGeoDsc.hh"

#include "TMath.h"
#include "TF1.h"

AllPixEdepHistogrammerDigitizer::AllPixEdepHistogrammerDigitizer(G4String modName, G4String hitsColName, G4String digitColName) 
: AllPixDigitizerInterface (modName) {

	// Registration of digits collection name
	collectionName.push_back(digitColName);
	m_hitsColName.push_back(hitsColName);

	// threshold
	m_digitIn.thl = 0.;
	Edeposition = new TH1D("","",2000,-10,10);


}

AllPixEdepHistogrammerDigitizer::~AllPixEdepHistogrammerDigitizer(){

	G4cout << "[EdepHistogrammer] writing Edep histo" << endl;
	Edeposition->SaveAs("EdepHisto.root");

}

void AllPixEdepHistogrammerDigitizer::Digitize(){

	// create the digits collection
	m_digitsCollection = new AllPixEdepHistogrammerDigitsCollection("AllPixEdepHistogrammerDigitizer", collectionName[0] );

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

	for(G4int itr  = 0 ; itr < nEntries ; itr++) {
		G4double zpos = (*hitsCollection)[itr]->GetPosWithRespectToPixel().z()/um;
		// Dose in Rad , Energy in joule / kg(V*rho)
		double doseInRad = 1.60217733e-16*((*hitsCollection)[itr]->GetEdep()/keV)/(0.001*0.001*0.01e-6*2532.59);
		Edeposition->Fill(zpos,doseInRad);
		//G4cout << "Deposition of " << (*hitsCollection)[itr]->GetEdep()/keV << "keV at x,y,z : " << xpos << " " << ypos << " " << zpos << endl;
		tempPixel.first  = (*hitsCollection)[itr]->GetPixelNbX();
		tempPixel.second = (*hitsCollection)[itr]->GetPixelNbY();
		pixelsContent[tempPixel] += (*hitsCollection)[itr]->GetEdep();

	}
	
	if (nEntries>0) { 
	//G4cout << "[EdepHistogrammer] writing Edep histo" << endl;
	//Edeposition->SaveAs("EdepHisto.root");
	};

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

		if((*pCItr).second > m_digitIn.thl) // over threshold !
		{
			double doseInRad = 1.60217733e-16*((*pCItr).second/keV)/(0.001*0.001*0.00005*2532.59);
			G4cout << "E Deposited in this pixel : " << (*pCItr).second/keV << " keV" << " " << doseInRad << " Rad" << endl;
			
			// Create one digit per pixel, I need to look at all the pixels first
			AllPixEdepHistogrammerDigit * digit = new AllPixEdepHistogrammerDigit;
			digit->SetPixelIDX((*pCItr).first.first);
			digit->SetPixelIDY((*pCItr).first.second);
			digit->IncreasePixelCounts(); // Counting mode

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
