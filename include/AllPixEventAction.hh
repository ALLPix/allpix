#ifndef AllPixEventAction_h
#define AllPixEventActino_h 1

#include "G4UserEventAction.hh"
#include "G4Event.hh"
#include "AllPixSteppingAction.hh"

//#include "AllPixMimosa26Digitizer.hh"
#include "AllPixDigitizerInterface.hh"
#include "G4VDigitizerModule.hh"
#include "AllPixRunAction.hh"

#include <vector>
#include <string>

using namespace std;

class AllPixGeoDsc;

class AllPixEventAction : public G4UserEventAction {

public:

	AllPixEventAction(AllPixRunAction *, AllPixSteppingAction *);
	virtual ~AllPixEventAction();

	void BeginOfEventAction(const G4Event*);
	void   EndOfEventAction(const G4Event*);

	void SetDetectorDigitInput(G4double, G4int);
	//void SetupDigitizers(G4String, AllPixGeoDsc *);
	void SetupDigitizers();
	G4String GetNewName(G4String, G4String, G4String);
	G4int GetNumberOfDigitizers() { return m_nDigitizers; };
	G4int GetNumberOfHC() { return m_nHC; };

private:

	AllPixRunAction * m_run_action;
	AllPixSteppingAction * m_stepping_action;
	vector<AllPixDigitizerInterface *> m_digiPtrs;

	// digitizer modules names
	vector<G4String> digitizerModulesNames;

	// number of hit collections
	G4int m_nHC;

	// number of digitizers
	G4int m_nDigitizers;

};

#endif
