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
//
// $Id: ExN06SteppingVerbose.cc,v 1.4 2006/06/29 17:54:37 gunter Exp $
// GEANT4 tag $Name: geant4-08-01-patch-01 $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "AllPixSteppingVerbose.hh"
#include "AllPixPostDetConstruction.hh"

#include "G4SteppingManager.hh"
#include "G4UnitsTable.hh"
#include "CLHEP/Units/SystemOfUnits.h"

using namespace CLHEP;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

AllPixSteppingVerbose::AllPixSteppingVerbose() {

	// Special plain output used for connectivity purposes
	// e.x. connection with the theoryicists ;)
	m_simplePlainOutput = new ofstream("HitsListChargeDiffusionModel.txt");

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

AllPixSteppingVerbose::~AllPixSteppingVerbose() {

	// Special plain output used for connectivity purposes
	// e.x. connection with the theoryicists ;)
	m_simplePlainOutput->close();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void AllPixSteppingVerbose::StepInfo() {

	CopyState();

	AllPixPostDetConstruction::GetInstance()->AddPoint(fTrack->GetPosition().x(),
			fTrack->GetPosition().y(),
			fTrack->GetPosition().z());

	// particle
	//G4cout << "--> " << fTrack->GetDefinition()->GetParticleName() << G4endl;

	G4int prec = G4cout.precision(3);

	if( verboseLevel >= 2 ){
		if( verboseLevel >= 4 ) VerboseTrack();
		if( verboseLevel >= 3 ){
			G4cout << G4endl;
			G4cout << std::setw( 5) << "#Step#"     << " "
					<< std::setw( 6) << "X"          << "    "
					<< std::setw( 6) << "Y"          << "    "
					<< std::setw( 6) << "Z"          << "    "
					<< std::setw( 9) << "KineE"      << " "
					<< std::setw( 9) << "dEStep"     << " "
					<< std::setw(10) << "StepLeng"
					<< std::setw(10) << "TrakLeng"
					<< std::setw(10) << "Volume"    << "  "
					<< std::setw(10) << "Process"   << G4endl;
		}

		G4cout << std::setw(5) << fTrack->GetCurrentStepNumber() << " "
				<< std::setw(6) << G4BestUnit(fTrack->GetPosition().x(),"Length")
				<< std::setw(6) << G4BestUnit(fTrack->GetPosition().y(),"Length")
				<< std::setw(6) << G4BestUnit(fTrack->GetPosition().z(),"Length")
				<< std::setw(6) << G4BestUnit(fTrack->GetKineticEnergy(),"Energy")
				<< std::setw(6) << G4BestUnit(fStep->GetTotalEnergyDeposit(),"Energy")
				<< std::setw(6) << G4BestUnit(fStep->GetStepLength(),"Length")
				<< std::setw(6) << G4BestUnit(fTrack->GetTrackLength(),"Length")
				<< "  ";

		// Special plain output used for connectivity purposes
		// e.x. connection with the theoryicists ;)
		// This line prints out information for all hits, primaries and secondaries
		// Only if the volume is a Box (Timepix Sensor)
		if( fTrack->GetNextVolume() != 0 ) {
			if( fTrack->GetVolume()->GetName().find("Box") != string::npos ) {
				*m_simplePlainOutput
				<< fTrack->GetPosition().x()/um << "\t"            // X
				<< fTrack->GetPosition().y()/um << "\t"            // Y
				<< fTrack->GetPosition().z()/um << "\t"            // Z
				<< fStep->GetTotalEnergyDeposit()/keV << "\t"      // dEStep
				<< fTrack->GetStepLength()/um << "\t"             // StepLeng
				<< G4endl;
			}
		}

		// if( fStepStatus != fWorldBoundary){
		if( fTrack->GetNextVolume() != 0 ) {
			G4cout << std::setw(10) << fTrack->GetVolume()->GetName();
		} else {
			G4cout << std::setw(10) << "OutOfWorld";
		}

		if(fStep->GetPostStepPoint()->GetProcessDefinedStep() != 0){
			G4cout << "  "
					<< std::setw(10)
			<< fStep->GetPostStepPoint()->GetProcessDefinedStep()
			->GetProcessName();
		} else {
			G4cout << "   UserLimit";
		}

		G4cout << G4endl;

		if( verboseLevel == 2 ){
			G4int tN2ndariesTot = fN2ndariesAtRestDoIt +
					fN2ndariesAlongStepDoIt +
					fN2ndariesPostStepDoIt;
			if(tN2ndariesTot>0){
				G4cout << "    :----- List of 2ndaries - "
						<< "#SpawnInStep=" << std::setw(3) << tN2ndariesTot
						<< "(Rest="  << std::setw(2) << fN2ndariesAtRestDoIt
						<< ",Along=" << std::setw(2) << fN2ndariesAlongStepDoIt
						<< ",Post="  << std::setw(2) << fN2ndariesPostStepDoIt
						<< "), "
						<< "#SpawnTotal=" << std::setw(3) << (*fSecondary).size()
						<< " ---------------"
						<< G4endl;

				for(size_t lp1=(*fSecondary).size()-tN2ndariesTot;
						lp1<(*fSecondary).size(); lp1++){
					G4cout << "    : "
							<< std::setw(6)
					<< G4BestUnit((*fSecondary)[lp1]->GetPosition().x(),"Length")
					<< std::setw(6)
					<< G4BestUnit((*fSecondary)[lp1]->GetPosition().y(),"Length")
					<< std::setw(6)
					<< G4BestUnit((*fSecondary)[lp1]->GetPosition().z(),"Length")
					<< std::setw(6)
					<< G4BestUnit((*fSecondary)[lp1]->GetKineticEnergy(),"Energy")
					<< std::setw(10)
					<< (*fSecondary)[lp1]->GetDefinition()->GetParticleName();
					G4cout << G4endl;
				}

				G4cout << "    :-----------------------------"
						<< "----------------------------------"
						<< "-- EndOf2ndaries Info ---------------"
						<< G4endl;
			}
		}

	}
	G4cout.precision(prec);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void AllPixSteppingVerbose::TrackingStarted()
{
	CopyState();
	AllPixPostDetConstruction::GetInstance()->NewTrack(fTrack->GetPosition().x(),
			fTrack->GetPosition().y(),
			fTrack->GetPosition().z());
	G4int prec = G4cout.precision(3);
	if( verboseLevel > 1 ){

		G4cout << std::setw( 5) << "Step#"      << " "
				<< std::setw( 6) << "X"          << "    "
				<< std::setw( 6) << "Y"          << "    "
				<< std::setw( 6) << "Z"          << "    "
				<< std::setw( 9) << "KineE"      << " "
				<< std::setw( 9) << "dEStep"     << " "
				<< std::setw(10) << "StepLeng"
				<< std::setw(10) << "TrakLeng"
				<< std::setw(10) << "Volume"     << "  "
				<< std::setw(10) << "Process"    << G4endl;

		G4cout << std::setw( 5) << fTrack->GetCurrentStepNumber() << " "
				<< std::setw( 6) << G4BestUnit(fTrack->GetPosition().x(),"Length")
				<< std::setw( 6) << G4BestUnit(fTrack->GetPosition().y(),"Length")
				<< std::setw( 6) << G4BestUnit(fTrack->GetPosition().z(),"Length")
				<< std::setw( 6) << G4BestUnit(fTrack->GetKineticEnergy(),"Energy")
				<< std::setw( 6) << G4BestUnit(fStep->GetTotalEnergyDeposit(),"Energy")
				<< std::setw( 6) << G4BestUnit(fStep->GetStepLength(),"Length")
				<< std::setw( 6) << G4BestUnit(fTrack->GetTrackLength(),"Length")
				<< "  ";

		if(fTrack->GetNextVolume()){
			G4cout << std::setw(10) << fTrack->GetVolume()->GetName();
		} else {
			G4cout << "OutOfWorld";
		}
		G4cout << "    initStep" << G4endl;
	}
	G4cout.precision(prec);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
