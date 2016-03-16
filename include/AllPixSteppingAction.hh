#ifndef AllPixSteppingAction_h
#define AllPixSteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "G4Step.hh"

//#include "G4Event.hh"

//#include "AllPixMimosa26Digitizer.hh"
//#include "AllPixDigitizerInterface.hh"
//#include "G4VDigitizerModule.hh"
#include "AllPixRunAction.hh"

#include <vector>
#include <string>

using namespace std;

//class AllPixGeoDsc;

class AllPixSteppingAction : public G4UserSteppingAction {

public:

	AllPixSteppingAction();
	virtual ~AllPixSteppingAction();

	void UserSteppingAction(const G4Step*);

	//void BeginOfEventAction(const G4Event*);
	//void   EndOfEventAction(const G4Event*);

	//void SetDetectorDigitInput(G4double, G4int);
	//void SetupDigitizers(G4String, AllPixGeoDsc *);
	//void SetupDigitizers();
	//G4String GetNewName(G4String, G4String, G4String);
	//G4int GetNumberOfDigitizers() { return m_nDigitizers; };
	//G4int GetNumberOfHC() { return m_nHC; };
	void Clear();

	vector<int> GetTrackID() { return trackID;};
	vector<int> GetTrackFate() { return trackFate;};
	vector<float> GetEOut() { return eOut;};

private:

	//AllPixRunAction * m_run_action;
	//vector<AllPixDigitizerInterface *> m_digiPtrs;

	// digitizer modules names
	//vector<G4String> digitizerModulesNames;

	// number of hit collections
	//G4int m_nHC;

	// number of digitizers
	//G4int m_nDigitizers;

	vector<int> trackID;
	vector<int> trackFate;
	vector<float> eOut;

};

#endif
