/**
 *  Author John Idarraga <idarraga@cern.ch>
 */

#ifndef AllPixDigitizerInterface_h
#define AllPixDigitizerInterface_h 1

#include "G4VDigitizerModule.hh"
#include "AllPixMimosa26Digit.hh"
#include "G4PrimaryVertex.hh"
#include "ReadGeoDescription.hh"

#include <map>
#include <vector>

#include "TString.h"

using namespace std;

typedef struct {
	G4double thl;
} digitInput;

class AllPixGeoDsc;

/**
 *  An Interface for Digitizer
 */

class AllPixDigitizerInterface : public G4VDigitizerModule {

public:

	AllPixDigitizerInterface(G4String modName) : G4VDigitizerModule(modName) {

		// Pickup the right index
		TString theIndex_S = modName.data();
		theIndex_S.Remove(0,6); // remove BoxSD_
		int lastpos = theIndex_S.Index("_", 1, 0, TString::kExact);
		theIndex_S.Remove(lastpos, theIndex_S.Length()); // remove suffix
		int detId = atoi(theIndex_S.Data());

		// Get the detectors
		ReadGeoDescription * geoDsc = ReadGeoDescription::GetInstance();
		map<int, AllPixGeoDsc *> * detMap = geoDsc->GetDetectorsMap();
		m_gD = (*detMap)[detId];

		//ReadGeoDescription * geoDsc = ReadGeoDescription::GetInstance();
		//vector<AllPixGeoDsc *> * detVector = geoDsc->GetDetectorsVector();
		//m_gD = (*detVector)[0];

	};
	virtual ~AllPixDigitizerInterface(){};

	virtual void SetPrimaryVertex(G4PrimaryVertex *) = 0;
	virtual void Digitize () = 0;
	virtual void SetDetectorDigitInputs(G4double) = 0;

	void SetDetectorGeoDscPtr(AllPixGeoDsc * gD){ m_gD = gD; };

protected:
	AllPixGeoDsc * GetDetectorGeoDscPtr(){ return m_gD; }; // first detector

private:
	AllPixGeoDsc * m_gD;

};

#endif
