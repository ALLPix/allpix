/**
 *  Author John Idarraga <idarraga@cern.ch>
 */

#include "AllPixSteppingAction.hh"
//#include "G4DigiManager.hh"
//#include "G4SDManager.hh"
//#include "G4PrimaryVertex.hh"
//#include "AllPixMimosa26Digitizer.hh"
//#include "AllPixFEI3StandardDigitizer.hh"

// geometry
//#include "ReadGeoDescription.hh"

AllPixSteppingAction::AllPixSteppingAction(){

	//m_run_action = run;
	//m_nHC = 0;
	//m_nDigitizers = 0;


	trackID.clear();
	trackFate.clear();
    eOut.clear();

}

void AllPixSteppingAction::Clear() {
	trackID.clear();
	trackFate.clear();
	eOut.clear();
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

AllPixSteppingAction::~AllPixSteppingAction()
{

}

void AllPixSteppingAction::UserSteppingAction(const G4Step* step) {

	G4String name0="";
	G4String name1="";

	if (step->GetPostStepPoint()->GetPhysicalVolume()) {
	  // G4cout << "1" << ' ' << step->GetPostStepPoint()->GetPhysicalVolume() << G4endl;
	  name1 = step->GetPostStepPoint()->GetPhysicalVolume()->GetName();
	}
	if (step->GetPreStepPoint()->GetPhysicalVolume()) {
	  // G4cout << "2" << ' ' << step->GetPreStepPoint()->GetPhysicalVolume() << G4endl;
	  name0 = step->GetPreStepPoint()->GetPhysicalVolume()->GetName();
	}

	//G4cout <<  " step: " << step->GetTrack()->GetTrackID() << ' '
	//	<< step->GetTrack()->GetVolume()->GetName() << G4endl;
	//G4cout << step->GetTrack()->GetPosition() << G4endl;

	int tmp_trackID = step->GetTrack()->GetTrackID();
	G4String name2 = step->GetTrack()->GetVolume()->GetName();

	// todo: avoid hard coded volume names

	if (name2=="Box_500") {

		bool done = false;

		for (unsigned i=0; i<trackID.size(); i++) {
		  if (trackID[i]==tmp_trackID) {
			done = true;
		  }
		}

		if (!done) {
		  //G4cout << "-->first" << G4endl;
		  trackID.push_back(tmp_trackID);
		  eOut.push_back(0);
		  trackFate.push_back(0);
		}
	}

	if (name0!=name1 && name0=="Box_500") {

		//int tmp_trackID = step->GetTrack()->GetTrackID();

		bool track_done = false;

		// only use first exit of sensitive volume

		for (unsigned i=0; i<trackID.size(); i++) {
		  if (trackFate[i]!=0) {
		   track_done = true;
		  }
		}

		if (!track_done) {

		  int pos = eOut.size()-1;

		  float tmp_E = step->GetTrack()->GetKineticEnergy();
		  eOut[pos]=tmp_E;

		  //G4cout << name0 << ' ' << name1 << ' ' << G4endl;

		  //G4cout << step->GetTrack()->GetTrackID() << G4endl;
		  //G4cout << step->GetTrack()->GetMomentumDirection() << G4endl;
		  //G4cout << step->GetTrack()->GetPosition() << G4endl;
		  //G4cout << step->GetTrack()->GetKineticEnergy() << G4endl;
		  //G4cout << step->GetTrack()->GetParticleDefinition()->GetPDGEncoding() << G4endl;
		  //G4cout << step->GetPostStepPoint()->GetKineticEnergy() << G4endl;
          if (name1=="Coverlayer_500_phys") {
            trackFate[pos]=1;
          }
          // in simulation without bumps and/or chip the wrapper or the chip is the next volume
          else if (name1=="BumpBox_500_phys" || name1=="Bump_500phys" || name1=="wrapper_500_phys" || name1=="Chip_500_phys") {
        	trackFate[pos]=2;
          }
          else if (name1=="GuardRings_500_phys") {
            trackFate[pos]=3;
          }
          else {
            G4cout << "(AllPixSteppingAction) entry into unspecified Volume: " << name1 << G4endl;
	      }

	   }
	}


}





//G4String AllPixEventAction::GetNewName(G4String oldName, G4String toErase, G4String toAppend){
//
//	string oldName_s = oldName.data();
//	string toErase_s = toErase.data();
//	string toAppend_s = toAppend.data();
//	string fixedString;
//
//	size_t pos = oldName_s.find_last_of(toErase) - toErase_s.size();
//	if(pos == string::npos) // couldn't find it
//		return oldName;
//
//
//	fixedString = oldName_s.substr(0,pos);
//	fixedString.append(toAppend_s);
//
//	return G4String(fixedString);
//}

//void AllPixEventAction::SetDetectorDigitInput(G4double thl, G4int detId){
//
//	m_digiPtrs[detId]->SetDetectorDigitInputs(thl);
//
//}

//void AllPixEventAction::BeginOfEventAction(const G4Event * /*evt*/)
//{
//
//	//G4PrimaryVertex * pv = evt->GetPrimaryVertex();
//
//}

//void AllPixEventAction::EndOfEventAction(const G4Event * evt)
//{
//
//	G4DigiManager * fDM = G4DigiManager::GetDMpointer();
//	// find digitizer module and digitize
//	AllPixMimosa26Digitizer * myDM;
//
//	G4PrimaryVertex * pv = evt->GetPrimaryVertex();
//
//	for(G4int itr = 0 ; itr < m_nDigitizers ; itr++){
//		myDM = (AllPixMimosa26Digitizer*)fDM->FindDigitizerModule( digitizerModulesNames[itr] );
//		myDM->SetPrimaryVertex(pv);
//		myDM->Digitize();
//	}
//
//	// digits will be retrieved at the end of the event in AllPixRun.
//
//}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

