/**
 *  Authors:
 *   John Idarraga <idarraga@cern.ch>
 *   Mathieu Benoit <Mathieu.Benoit@cern.ch>
 */

#include "AllPixEventAction.hh"
#include "G4DigiManager.hh"
#include "G4SDManager.hh"
#include "G4PrimaryVertex.hh"
#include "AllPixMedipix2Digitizer.hh"
#include "AllPixFEI3StandardDigitizer.hh"
#include "AllPixMimosa26Digitizer.hh"
#include "AllPixTimepixDigitizer.hh"
#include "AllPixMCTruthDigitizer.hh"
#include "AllPixLETCalculatorDigitizer.hh"
#include "AllPixFEI4RadDamageDigitizer.hh"

// Included by newdigitizer.sh script --> Medipix3RX
#include "AllPixMedipix3RXDigitizer.hh"
// __endofheader__

// geometry
#include "ReadGeoDescription.hh"

void AllPixEventAction::SetupDigitizers(){

	// Digit manager
	G4DigiManager * fDM = G4DigiManager::GetDMpointer();

	// geo description
	extern ReadGeoDescription * g_GeoDsc; // already loaded ! :)
	map<int, AllPixGeoDsc *> * geoMap = g_GeoDsc->GetDetectorsMap();

	// I need as many digitizer as detectors.
	// I decide that through the hitCollections
	// There is one hit collection by detector
	//  but not all of them are pixel detectors
	G4SDManager * SDman = G4SDManager::GetSDMpointer();
	G4HCtable * HCTable = SDman->GetHCtable();
	m_nHC = HCTable->entries();
	map<int, AllPixGeoDsc *>::iterator geoItr;
	G4String digitizerName;
	G4int detectorId;

	// loop over all hit collections and identify those which need a digitizer, i.e. only detectors
	for(G4int itr = 0 ; itr < m_nHC ; itr++)
	{

		G4String hcName = HCTable->GetHCname(itr);

		// Check if this is a Box.  If it is not a Si wafer it doesn't need to be digitized
		G4String searchS = "BoxSD";

		G4String digitSuffix = "_";
		digitSuffix += digitizerName;
		digitSuffix += "Digitizer";
		G4String digitizerModName = GetNewName(hcName, "HitsCollection", digitSuffix);

		bool notADetector = false;
		if(!digitizerModName.contains("Box")){
			G4cout << "SD [" << hcName
					<< "] --> This is not a sensitive detector needing a Digitizer ! ... skipping"
					<< G4endl;
			notADetector = true;
		}
		if(notADetector) continue;

		// find the detector id for this digitizer in the detector db
		bool found = false;
		for( geoItr = geoMap->begin() ; geoItr != geoMap->end() ; geoItr++ ){
			if( hcName == (*geoItr).second->GetHitsCollectionName() ){
				digitizerName = (*geoItr).second->GetSensorDigitizer();
				detectorId = (*geoItr).first;
				G4cout << "id : " << detectorId << G4endl;
				found = true;
				break;
			}
		}
		if(!found){
			cout << "Couldn't find the detector in the list !!!??? "<< endl;
			exit(1);
		}

		// Load the name of the collection
		(*geoMap)[detectorId]->SetDigitCollectionName(digitizerModName);

		// first check if this digitizer module is already there
		vector<G4String>::iterator it;
		it = find (digitizerModulesNames.begin(), digitizerModulesNames.end(), digitizerModName);

		// If the digitizer module is in the list,
		//  it means the digitizer was already created for this detector, jump !
		G4cout << "SD [" << hcName << "] Attemping to build digitizer \""
				<< digitizerModName << "\"" << G4endl;

		if(it != digitizerModulesNames.end()) {
			G4cout << "WARNING : Digitizer module \"" << digitizerModName
					<< "\" was already created ... skipping." << G4endl;
			continue;
		}

		digitizerModulesNames.push_back(digitizerModName);

		// Now build the digit Collection name
		G4String digitColectionName = GetNewName(hcName,"HitsCollection","_DigitCollection");

		// Creating an instance of the actual digitizer, and keep pointer through the interface
		AllPixDigitizerInterface * dmPtr;
        // __beginofdigitlist__
		if (digitizerName == "FEI3Standard") {
			AllPixFEI3StandardDigitizer * dp = new AllPixFEI3StandardDigitizer(digitizerModulesNames[itr] , hcName, digitColectionName);
			dmPtr = static_cast<AllPixDigitizerInterface *> (dp);
			cout << "    Setting up a " << digitizerName << " digitizer for det : " << detectorId << endl;
		} else if (digitizerName == "Medipix2") {
			AllPixMedipix2Digitizer * dp = new AllPixMedipix2Digitizer(digitizerModulesNames[itr] , hcName, digitColectionName);
			dmPtr = static_cast<AllPixDigitizerInterface *> (dp);
			cout << "    Setting up a " << digitizerName << " digitizer for det : " << detectorId << endl;
		} else if (digitizerName == "Mimosa26") {
			AllPixMimosa26Digitizer * dp = new AllPixMimosa26Digitizer(digitizerModulesNames[itr] , hcName, digitColectionName);
			dmPtr = static_cast<AllPixDigitizerInterface *> (dp);
			cout << "    Setting up a " << digitizerName << " digitizer for det : " << detectorId << endl;
		} else if (digitizerName == "Timepix") {
			G4cout << digitizerModulesNames[itr] << "  " << hcName << "  " << digitColectionName << G4endl;
			AllPixTimepixDigitizer * dp = new AllPixTimepixDigitizer(digitizerModulesNames[itr] , hcName, digitColectionName);
			dmPtr = static_cast<AllPixDigitizerInterface *> (dp);
			cout << "    Setting up a " << digitizerName << " digitizer for det : " << detectorId << endl;
		}else if (digitizerName == "MCTruth") {
			AllPixMCTruthDigitizer * dp = new AllPixMCTruthDigitizer(digitizerModulesNames[itr] , hcName, digitColectionName);
			dmPtr = static_cast<AllPixDigitizerInterface *> (dp);
			cout << "    Setting up a " << digitizerName << " digitizer for det : " << detectorId << endl;
		}// Included by newdigitizer.sh script --> LETCalculator
		else if (digitizerName == "LETCalculator") {
					AllPixLETCalculatorDigitizer * dp = new AllPixLETCalculatorDigitizer(digitizerModulesNames[itr] , hcName, digitColectionName);
					dmPtr = static_cast<AllPixDigitizerInterface *> (dp);
					cout << "    Setting up a " << digitizerName << " digitizer for det : " << detectorId << endl;
				}
		else if (digitizerName == "FEI4RadDamage") {
			AllPixFEI4RadDamageDigitizer * dp = new AllPixFEI4RadDamageDigitizer(digitizerModulesNames[itr] , hcName, digitColectionName);
			dmPtr = static_cast<AllPixDigitizerInterface *> (dp);
			cout << "    Setting up a " << digitizerName << " digitizer for det : " << detectorId << endl;
		}
 
// Included by newdigitizer.sh script --> Medipix3RX
else if (digitizerName == "Medipix3RX") {
			AllPixMedipix3RXDigitizer * dp = new AllPixMedipix3RXDigitizer(digitizerModulesNames[itr] , hcName, digitColectionName);
			dmPtr = static_cast<AllPixDigitizerInterface *> (dp);
			cout << "    Setting up a " << digitizerName << " digitizer for det : " << detectorId << endl;
		}
 
        // __endofdigitlist__
   	    else {
			G4cout << "    can't find digitizer with name : " << digitizerName << G4endl;
			exit(1);
		}

		///////////////////////////////////////////////////////////
		// Common task to all digitizers provided in the interface
		// pass here the AllPixGeoDsc ptr
		dmPtr->SetDetectorGeoDscPtr((*geoMap)[detectorId]);

		// push back the digitizer
		m_digiPtrs.push_back( dmPtr );
		fDM->AddNewModule(m_digiPtrs[itr]);
		m_nDigitizers++;

	}

}
