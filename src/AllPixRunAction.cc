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
// $Id: ExAllPixRunAction.cc,v 1.10 2006/06/29 17:54:31 gunter Exp $
// GEANT4 tag $Name: geant4-09-02 $
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Make this appear first!
#include "G4Timer.hh"
#include "G4Run.hh"

#include "AllPixRunAction.hh"
#include "AllPixRun.hh"

#include "AllPixDetectorConstruction.hh"
#include "AllPix_Frames_WriteToEntuple.h"
#include "allpix_dm.h"
#include "AllPixPrimaryGeneratorMessenger.hh"

#include <vector>
#include <string>
using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// this constructor is called only once in the whole program
AllPixRunAction::AllPixRunAction(AllPixDetectorConstruction * det, TString ds, TString td, TString lciofn, TString lciofn_dut)
{

	m_detectorPtr = det;
	m_dataset = ds;
	m_tempdir = td;

	timer = new G4Timer;

	// file for lcio format conversion
	m_lciobridge_f = new ofstream(lciofn);
	m_lciobridge_dut_f = new ofstream(lciofn_dut);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

AllPixRunAction::~AllPixRunAction()
{

	m_lciobridge_f->close();
	m_lciobridge_dut_f->close();

	delete timer;
}

/**
 *  Concrete implementation
 */
G4Run * AllPixRunAction::GenerateRun(){

	
	m_writeTPixTelescopeFilesFlag = AllPixMessenger->GetTimepixTelescopeWriteFlag(); //pass it to veto RecordTelescopeDigits
	m_AllPixRun = new AllPixRun(m_detectorPtr, m_detectorPtr->GetOutputFilePrefix(),
	m_dataset, m_tempdir, m_writeTPixTelescopeFilesFlag); // keep this pointer
	m_AllPixRun->SetLCIOBridgeFileDsc(m_lciobridge_f, m_lciobridge_dut_f);

	return m_AllPixRun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void AllPixRunAction::BeginOfRunAction(const G4Run* aRun)
{
	G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;
	timer->Start();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void AllPixRunAction::EndOfRunAction(const G4Run* aRun)
{   

	// at the end of the run
	G4cout << "Filling frames ntuple" << G4endl;
	m_AllPixRun->FillFramesNtuple(aRun);

	/*
	 * write telescope files
	 * these 4 variables are set via the /allpix/timepixtelescope messenger
	 * in AllPixPrimaryGeneratorMessenger.cc
	 * and passed here through Get methods
	 */
	G4bool writeFlag    = AllPixMessenger->GetTimepixTelescopeWriteFlag();
	G4String folderName = AllPixMessenger->GetTimepixTelescopeFolderName();
	G4bool eventIDflag  = AllPixMessenger->GetTimepixTelescopeDoEventFlag();
	G4bool sumTOTflag   = AllPixMessenger->GetTimepixTelescopeSumTOTFlag();
	if (writeFlag) {
		G4cout << "Filling telescope files" << G4endl;
		m_AllPixRun->FillTelescopeFiles(aRun,folderName,eventIDflag,sumTOTflag);
	}

	timer->Stop();
	G4cout << "event Id = " << aRun->GetNumberOfEvent()
        				 << " " << *timer << G4endl;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
