/**
 *  Author John Idarraga <idarraga@cern.ch>
 */
//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: AllPixTrackerSD.cc,v 1.9 2006/06/29 17:48:27 gunter Exp $
// GEANT4 tag $Name: geant4-09-02 $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "AllPixTrackerSD.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4RunManager.hh"
#include "G4ios.hh"
#include "G4Track.hh"
#include "G4VProcess.hh"
#include "G4DecayTable.hh"
#include "G4LogicalVolume.hh"

#include "AllPixGeoDsc.hh"

#include "TMath.h"
#include "TString.h"

// temp
G4double g_temp_edep = 0.;
G4int g_temp_pdgId = 0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

AllPixTrackerSD::AllPixTrackerSD(G4String name,
		G4ThreeVector absPosWrapper,
		G4ThreeVector relPosSD,
		AllPixGeoDsc * gD,
		G4RotationMatrix * rot)
:G4VSensitiveDetector(name)
{

	m_thisHitsCollectionName = name + "_HitsCollection";
	collectionName.insert( m_thisHitsCollectionName ); // only one hit collection per detector for now

	// Absolute position of SD
	m_absolutePosOfWrapper = absPosWrapper;
	m_relativePosOfSD = relPosSD;
	m_rotationOfWrapper = rot;
	m_gD = gD; // Geo description
	m_thisIsAPixelDetector = true;

	m_globalTrackId_Dump = 0;
	firstStrikePrimary = false;
	_totalEdep = 0;

}
/*
 * Second constructor for sensitive devices which are
 * not actual pixel detectors
 */
AllPixTrackerSD::AllPixTrackerSD(G4String name, G4ThreeVector absPos, G4RotationMatrix * rot)
:G4VSensitiveDetector(name)
{

	m_thisHitsCollectionName = name + "_HitsCollection";
	collectionName.insert( m_thisHitsCollectionName ); // only one hit collection per detector for now

	// reduced set of parameters for a SD which is not a pixel detector
	m_absolutePosOfWrapper = absPos;
	m_relativePosOfSD = G4ThreeVector(0,0,0);
	m_rotationOfWrapper = rot;
	m_gD = 0x0; // Geo description
	m_thisIsAPixelDetector = false;

	m_globalTrackId_Dump = 0;
	_totalEdep = 0;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

AllPixTrackerSD::~AllPixTrackerSD(){ 

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/**
 * Runs once per event.  The hitsCollections pointers are retrieved
 *	I have a set of pointers to control that I get them right
 */
void AllPixTrackerSD::Initialize(G4HCofThisEvent* HCE)
{

	// Create the hit collection
	hitsCollection = new AllPixTrackerHitsCollection
			(SensitiveDetectorName, collectionName[0]);

	static G4int HCID = -1;
	HCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);

	// Add to hits collection of this event
	HCE->AddHitsCollection( HCID, hitsCollection );

	// Insert the pointer in a set to check its existence later
	// Normally there is only one instance of AllPixTrackerSD per sensitive volume
	m_hitsCollectionSet.insert(hitsCollection);

}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool AllPixTrackerSD::ProcessHits(G4Step * aStep, G4TouchableHistory *)
{

	// Check first if we are working in a valid HitCollection
	// This should never happen if the user don't replicate detector Id in the macro.
	// A check is done in ReadGeoDescription ... I probably don't need to check this here.
	try {
		if( m_hitsCollectionSet.find(hitsCollection) == m_hitsCollectionSet.end()) // not found !
			throw hitsCollection;
	} catch (AllPixTrackerHitsCollection * h) {
		G4cout << "[WARNING] The following pointer to hitsCollection is invalid : " << h << G4endl;
		G4cout << "          The available set contains " << m_hitsCollectionSet.size() << " collections." << G4endl;
		G4cout << "          Trying to recover by ignoring this hit where this hitCollection" << G4endl;
		G4cout << "          !!! You may have an error in your macro, bad detector ID !!!" << G4endl;
		G4cout << "          has been given. AllPixTrackerSD::ProcessHits returns false." << G4endl;
		return false;
	}

	// track
	G4Track * aTrack = aStep->GetTrack();
	// particle
	G4ParticleDefinition * aParticle = aTrack->GetDefinition();

	// not used for now, not pretty.
	// I need to know where is this hit and tie it to a detector ID
	// G4String detId_S = touchablepre->GetSolid(2)->GetName();
	// detId_S.remove(0,4); // remove the part "Box_", take the ID

	// Isolate the information of the incoming particle when doing the first hit in the SD
	//G4cout << aParticle->GetPDGEncoding() << G4endl;

	if ( aTrack->GetTrackID() == 1 && ! firstStrikePrimary ) {
		_kinEPrimary = aTrack->GetKineticEnergy()/keV;
		_kinEPrimary -= (aStep->GetDeltaEnergy()/keV);
		firstStrikePrimary = true;
	}

	// Work with the Hit
	G4double edep = aStep->GetTotalEnergyDeposit();
	if(edep==0.) return false;

	G4StepPoint * preStepPoint = aStep->GetPreStepPoint();
	G4StepPoint * postStepPoint = aStep->GetPostStepPoint();

	const G4TouchableHandle touchablepre = preStepPoint->GetTouchableHandle();
	const G4TouchableHandle touchablepost = postStepPoint->GetTouchableHandle();

	G4int copyIDy_pre  = -1;
	G4int copyIDx_pre  = -1;
	G4int copyIDy_post = -1;
	G4int copyIDx_post = -1;
	G4ThreeVector correctedPos(0,0,0);

	if (m_thisIsAPixelDetector) {
		// This positions are global, I will bring them to pixel-centered frame
		// I can use the physical volumes for that
		G4ThreeVector prePos = preStepPoint->GetPosition();


		// Find the inverse rotation
		//G4RotationMatrix invRot2 = CLHEP::inverseOf(*m_rotationOfWrapper);
		G4RotationMatrix invRot = m_rotationOfWrapper->inverse().inverse();

		// Absolute center of Si wafer
		G4ThreeVector absCenterOfDetector = m_absolutePosOfWrapper ;
		//G4cout << "Absolute position of Wrapper : " << absCenterOfDetector << endl;

		// Bring the detector (Si layer) to the Origin
		correctedPos = prePos;
		correctedPos -= absCenterOfDetector;
		// apply rotation !
		correctedPos = invRot * correctedPos;

		// Now let's finally provide pixel-centered coordinates for each hit
		// Build the center of the Pixel
		G4ThreeVector centerOfPixel(
				m_gD->GetPixelX()*TMath::FloorNint(correctedPos.x() / m_gD->GetPixelX()) + m_gD->GetHalfPixelX(),
				m_gD->GetPixelY()*TMath::FloorNint(correctedPos.y() / m_gD->GetPixelY()) + m_gD->GetHalfPixelY(),
				0.); // in the middle of the tower

		// The position within the pixel !!!

                // 20160316: The line below wrong and only works if the relative postion is (0,0,0)   
		// correctedPos = correctedPos - centerOfPixel - m_relativePosOfSD;
		correctedPos = correctedPos - centerOfPixel;
		
		//G4cout << "uncorrectedPos : " << prePos.x()/um << " " << prePos.y()/um
		//	   << " " << prePos.z()/um << " [um]" << G4endl;

		//G4cout << "correctedPos : " << correctedPos.x()/um << " " << correctedPos.y()/um
		//	   << " " << correctedPos.z()/um << " [um]" << G4endl;

		//G4cout << "(" << shi << ") "<<	 correctedPos.z()/um << " " ;
		//G4cout << TString::Format("(%02.0f) %02.1f ",shi,correctedPos.z()/um);

		// depth 1 --> x
		// depth 0 --> y
		copyIDy_pre  = touchablepre->GetCopyNumber();
		copyIDx_pre  = touchablepre->GetCopyNumber(1);

	}

	// Look at the touchablepost only if in the same volume, i.e. in the sensitive Si Box
	// If the hit is in a different pixel, it is still the same phys volume
	// The problem is that if I the postStep is in a different volume, GetCopyNumber(1)
	//  doesn't make sense anymore.
	G4ThreeVector postPos(0,0,0);
	if(preStepPoint->GetPhysicalVolume() == postStepPoint->GetPhysicalVolume()){
		postPos = postStepPoint->GetPosition();
		copyIDy_post = touchablepost->GetCopyNumber();
		copyIDx_post = touchablepost->GetCopyNumber(1);
	}

	// process
	const G4VProcess * aProcessPointer = aStep->GetPostStepPoint()->GetProcessDefinedStep();


	// create a hit instance
	AllPixTrackerHit * newHit = new AllPixTrackerHit();
	//newHit->SetDetId(atoi(detId_S.c_str()));
	newHit->SetTrackID(aTrack->GetTrackID());
	newHit->SetParentID(aTrack->GetParentID());
	newHit->SetPixelNbX(copyIDx_pre);
	newHit->SetPixelNbY(copyIDy_pre);
	newHit->SetPostPixelNbX(copyIDx_post);
	newHit->SetPostPixelNbY(copyIDy_post);
	newHit->SetEdep(edep);
	_totalEdep += edep;
	newHit->SetPos(aStep->GetPostStepPoint()->GetPosition());

	newHit->SetPosWithRespectToPixel( correctedPos );

	newHit->SetProcessName(aProcessPointer->GetProcessName());
	newHit->SetTrackPdgId(aParticle->GetPDGEncoding());

	newHit->SetKinEParent( _kinEPrimary );

	/////////////////////////////////////////////
	g_temp_edep = edep;
	g_temp_pdgId = aParticle->GetPDGEncoding();
	/////////////////////////////////////////////

	newHit->SetTrackVolumeName(aTrack->GetVolume()->GetName());
	newHit->SetParentVolumeName(aTrack->GetLogicalVolumeAtVertex()->GetName());

	//G4cout << "hitsCollection : " << hitsCollection << G4endl;
	//G4cout << "     entries --> " << hitsCollection->entries() << G4endl;
	hitsCollection->insert(newHit);
	//newHit->Print();
	//newHit->Draw();

	if ( _totalEdep > _kinEPrimary ) {
		cout << "[WARNING] totalEdep = " << _totalEdep << ", kinEPrimary = " << _kinEPrimary << endl;
	}

	return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void AllPixTrackerSD::EndOfEvent(G4HCofThisEvent*)
{

	//G4int NbHits = hitsCollection->entries();
	//if(NbHits > 0)
	//	G4cout << "--------> Hits Collection : " << collectionName[0] << " has " << NbHits << " hits " << G4endl;

	// clear the Set of pointers to hitCollection used for verification
	m_hitsCollectionSet.clear();
	firstStrikePrimary = false;
	_totalEdep = 0.;
	_kinEPrimary = 0.;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......



/*
G4cout << " ----------------------------- " << G4endl;
G4cout << "prePos       : " << prePos.x()/mm << " " << prePos.y()/mm << " " << prePos.z()/mm << " [mm]" << G4endl;
G4cout << "detectorPos  : " << absCenterOfDetector.x()/mm << " " << absCenterOfDetector.y()/mm << " " << absCenterOfDetector.z()/mm << " [mm]" << G4endl;
G4cout << "previous correctedPos : " << correctedPos.x()/um << " " << correctedPos.y()/um
			<< " " << correctedPos.z()/um << " [um]" << G4endl << G4endl;
 */

