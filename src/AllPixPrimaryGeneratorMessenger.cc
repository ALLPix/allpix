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
// $Id: AllPixPrimaryGeneratorMessenger.cc,v 1.3 2006/06/29 17:54:29 gunter Exp $
// GEANT4 tag $Name: geant4-09-02 $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "AllPixPrimaryGeneratorMessenger.hh"
#include "AllPixDetectorMessenger.hh"

#include "AllPixPrimaryGeneratorAction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcommand.hh"

#include "G4RunManager.hh"

#include "CLHEP/Random/RandGauss.h"
#include "CLHEP/Random/RandPoisson.h"
#include "CLHEP/Random/RandExponential.h"

#include "TFile.h"
#include "TH1.h"

#include "sys/types.h"
#include "sys/stat.h"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"


using namespace CLHEP;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

AllPixPrimaryGeneratorMessenger::AllPixPrimaryGeneratorMessenger(
								 AllPixPrimaryGeneratorAction* AllPixGun)
 :AllPixAction(AllPixGun)
{
  m_hits            = 1;
  m_frames          = 1;
  m_beamTypeHitFunc = "gauss";
  m_beamTypePar1    = 1.;
  m_beamTypePar2    = 1.;
  m_TimepixTelescopeWriteFlag = false;
  gunDir = new G4UIdirectory("/N06/gun/");
  gunDir->SetGuidance("PrimaryGenerator control");
  m_Write_MC_FilesFlag=false; //nalipour: MC hits

  polarCmd = new G4UIcmdWithADoubleAndUnit("/N06/gun/optPhotonPolar",this);
  polarCmd->SetGuidance("Set linear polarization");
  polarCmd->SetGuidance("  angle w.r.t. (k,n) plane");
  polarCmd->SetParameterName("angle",true);
  polarCmd->SetUnitCategory("Angle");
  polarCmd->SetDefaultValue(-360.0);
  polarCmd->SetDefaultUnit("deg");
  polarCmd->AvailableForStates(G4State_Idle);

  m_MCPrefix = new G4UIcmdWithAString("/allpix/config/HEPEvtFile", this);
  m_MCPrefix->SetGuidance("Set HEPEvt file path (path can be included).");
  m_MCPrefix->SetGuidance("If no path is specified, the file will be written to ./");
  m_MCPrefix->SetParameterName("MCFile", false);
  m_MCPrefix->SetDefaultValue("/afs/cern.ch/user/m/mbenoit/scratch0/Full_Tracker_Model/HEPEVT_files/AllPairs3TeV8MeV6Deg.HEPEvt00");
  m_MCPrefix->AvailableForStates(G4State_PreInit, G4State_Idle);

  m_userBeamNumberOfFramesCmd = new G4UIcmdWithAnInteger("/allpix/beam/frames",this);
  m_userBeamNumberOfFramesCmd->SetGuidance("Set number of frames");
  m_userBeamNumberOfFramesCmd->SetParameterName("Frames",false,false);
  m_userBeamNumberOfFramesCmd->SetDefaultValue(1);
  m_userBeamNumberOfFramesCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  m_userBeamOnCmd = new G4UIcmdWithoutParameter("/allpix/beam/on",this);
  m_userBeamOnCmd->SetGuidance("Set beam ON");
  m_userBeamOnCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  m_userBeamTypeCmd = new G4UIcommand("/allpix/beam/type",this);
  m_userBeamTypeCmd->SetGuidance("Select hits distribution function.");
  m_userBeamTypeCmd->SetGuidance("Current possible functions:");
  m_userBeamTypeCmd->SetGuidance("const <hits>");
  m_userBeamTypeCmd->SetGuidance("gauss <mean> <sigma>");
  m_userBeamTypeCmd->SetGuidance("poisson <mean>");
  m_userBeamTypeCmd->SetGuidance("expo <tau>");
  m_userBeamTypeCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  G4UIparameter * func = new G4UIparameter("hitFunction",'s',true);
  func->SetDefaultValue("const");
  m_userBeamTypeCmd->SetParameter(func);
  G4UIparameter * p1 = new G4UIparameter("par1",'d',true);
  p1->SetDefaultValue(1.);
  p1->SetParameterRange("par1 >= 0");
  m_userBeamTypeCmd->SetParameter(p1);
  G4UIparameter * p2 = new G4UIparameter("par2",'d',true);
  p2->SetDefaultValue(1.);
  m_userBeamTypeCmd->SetParameter(p2);
  p2->SetParameterRange("par2 > 0");

  m_TimepixTelescopeWriteCmd = new G4UIcmdWithABool("/allpix/timepixtelescope/write",this);
  m_TimepixTelescopeWriteCmd->SetGuidance("Switch on/off writing Timepix Telescope files. Default OFF.");
  m_TimepixTelescopeWriteCmd->SetDefaultValue(false);
  m_TimepixTelescopeWriteCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  m_TimepixTelescopeFolderNameCmd = new G4UIcmdWithAString("/allpix/timepixtelescope/setFolderPath",this);;
  m_TimepixTelescopeFolderNameCmd->SetGuidance("Set Timepix Telescope files folder path.");
  m_TimepixTelescopeFolderNameCmd->SetGuidance("Directory structure will be created if it does not exist.");
  m_TimepixTelescopeFolderNameCmd->SetGuidance("Default is ./");
  m_TimepixTelescopeFolderNameCmd->SetDefaultValue("./");
  m_TimepixTelescopeFolderNameCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  m_TimepixTelescopeDoEventCmd = new G4UIcmdWithABool("/allpix/timepixtelescope/setEventIDcolumn",this);;
  m_TimepixTelescopeDoEventCmd->SetGuidance("Dump EventID column in Timepix Telescope files. Default OFF.");
  m_TimepixTelescopeDoEventCmd->SetDefaultValue(false);
  m_TimepixTelescopeDoEventCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  m_TimepixTelescopeSumTOTCmd = new G4UIcmdWithABool("/allpix/timepixtelescope/setSumTOT",this);;
  m_TimepixTelescopeSumTOTCmd->SetGuidance("Sum TOT of pixels with multiple hits/frame. Default ON.");
  m_TimepixTelescopeSumTOTCmd->SetDefaultValue(true);
  m_TimepixTelescopeSumTOTCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  //nalipour: MC hits
  m_Write_MC_FilesCmd = new G4UIcmdWithABool("/allpix/WriteROOTFiles/write",this);
  m_Write_MC_FilesCmd->SetGuidance("Switch on/off writing ROOT files containing the MC and the AllPix information. Default OFF.");
  m_Write_MC_FilesCmd->SetDefaultValue(false);
  m_Write_MC_FilesCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  m_Write_MC_FolderNameCmd = new G4UIcmdWithAString("/allpix/WriteROOTFiles/setFolderPath",this);;
  m_Write_MC_FolderNameCmd->SetGuidance("Set folder path for the ROOT Files containing MC information.");
  m_Write_MC_FolderNameCmd->SetGuidance("Directory structure will be created if it does not exist.");
  m_Write_MC_FolderNameCmd->SetGuidance("Default is ./");
  m_Write_MC_FolderNameCmd->SetDefaultValue("./");
  m_Write_MC_FolderNameCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

AllPixPrimaryGeneratorMessenger::~AllPixPrimaryGeneratorMessenger()
{
  delete polarCmd;
  delete gunDir;
  delete m_MCPrefix;
  delete m_userBeamNumberOfFramesCmd;
  delete m_userBeamOnCmd;
  delete m_userBeamTypeCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void AllPixPrimaryGeneratorMessenger::SetNewValue(
						  G4UIcommand* command, G4String newValue)
{ 
  if( command == m_MCPrefix )
    {
      AllPixAction->SetHEPevtFile(newValue);
    }

  if( command == polarCmd )
    {
      G4double angle = polarCmd->GetNewDoubleValue(newValue);
      if ( angle == -360.0*deg )
	{
	  AllPixAction->SetOptPhotonPolar();
	}
      else
	{
	  AllPixAction->SetOptPhotonPolar(angle);
	}
    }

  if ( command == m_userBeamNumberOfFramesCmd )
    {
      m_frames = m_userBeamNumberOfFramesCmd->GetNewIntValue(newValue);
    }

  if ( command == m_userBeamTypeCmd )
    {
      const char* nv = (const char*)newValue;
      std::istringstream is(nv);
      is >> m_beamTypeHitFunc >> m_beamTypePar1 >> m_beamTypePar2;
    }

  if (command == m_userBeamOnCmd)
    {
      // Get run manager
      G4RunManager * runManager = G4RunManager::GetRunManager();

      // Write number of hits sent to check the distribution
      TFile * f = new TFile("hitFunction.root","recreate");
      TH1I * h = new TH1I("h","h",1000,0,1000);

      for (G4int ii = 0; ii < m_frames; ii++)
	{
	  if ( m_beamTypeHitFunc == "gauss" )
	    {
	      m_hits = CLHEP::RandGauss::shoot(m_beamTypePar1,m_beamTypePar2);
	    }
	  else if ( m_beamTypeHitFunc == "poisson" )
	    {
	      m_hits = CLHEP::RandPoisson::shoot(m_beamTypePar1);
	    }
	  else if ( m_beamTypeHitFunc == "expo" )
	    {
	      m_hits = CLHEP::RandExponential::shoot(m_beamTypePar1);
	    }
	  else if ( m_beamTypeHitFunc == "const")
	    {
	      m_hits = m_beamTypePar1;
	    }
	  else
	    {
	      G4cout << "============> Unknown parameter: " <<  m_beamTypeHitFunc  << "... sending "<< m_hits << " hit(s)..." << G4endl;
	    }

	  // same call as in /run/beamOn
	  runManager->BeamOn(m_hits);
	  h->Fill(m_hits);
	}

      f->Write();
      f->Close();
    }

  if ( command == m_TimepixTelescopeWriteCmd )
    {
      m_TimepixTelescopeWriteFlag = m_TimepixTelescopeWriteCmd->GetNewBoolValue(newValue);		
    }

  if ( command == m_TimepixTelescopeFolderNameCmd )
    {
      m_TimepixTelescopeFolderName = newValue;

      // check if folder exists, otherwise create it
      struct stat st;
      if ( stat(m_TimepixTelescopeFolderName,&st) != 0 )
	{
	  G4cout << "folder " << m_TimepixTelescopeFolderName << " does not exist, creating..." << G4endl;
	  system(TString::Format("mkdir -p %s",m_TimepixTelescopeFolderName.data()));
	}
    }

  if ( command == m_TimepixTelescopeDoEventCmd )
    {
      m_TimepixTelescopeDoEventFlag = m_TimepixTelescopeDoEventCmd->GetNewBoolValue(newValue);
    }

  if ( command == m_TimepixTelescopeSumTOTCmd )
    {
      m_TimepixTelescopeSumTOTFlag = m_TimepixTelescopeSumTOTCmd->GetNewBoolValue(newValue);
    }



  //nalipour: MC hits
  if (command == m_Write_MC_FilesCmd)
    {
      m_Write_MC_FilesFlag = m_Write_MC_FilesCmd->GetNewBoolValue(newValue);		
    }
  if (command == m_Write_MC_FolderNameCmd)
    {
      m_Write_MC_FolderName = newValue;
      
      // check if folder exists, otherwise create it
      struct stat st;
      if ( stat(m_Write_MC_FolderName,&st) != 0 )
	{
	  G4cout << "folder " <<m_Write_MC_FolderName  << " does not exist, creating..." << G4endl;
	  system(TString::Format("mkdir -p %s",m_Write_MC_FolderName.data()));
	}
    }
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

