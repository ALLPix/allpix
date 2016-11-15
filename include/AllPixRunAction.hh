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
// $Id: AllPixRunAction.hh,v 1.9 2006/06/29 17:54:10 gunter Exp $
// GEANT4 tag $Name: geant4-09-02 $
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef AllPixRunAction_h
#define AllPixRunAction_h 1

#include "globals.hh"
#include "G4UserRunAction.hh"

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>

#include "AllPixPrimaryGeneratorAction.hh"
#include "AllPixWriteROOTFile.hh" //nalipour
#include "AllPixWriter.hh"

#ifdef HAVE_LCIO
#include "AllPixLCIOwriter.hh"
#define lcio_h 1
#else
#define lcio_h 0
#endif

#include <iostream>
#include <fstream>
using namespace std;



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4Timer;
class G4Run;
class AllPixRun;
class AllPixDetectorConstruction;
class AllPixPrimaryGeneratorMessenger;
class AllPixWriteROOTFile; //nalipour
//class FramesHandler;
//class WriteToNtuple;

// simplified hits;
typedef struct {
  float edep;
  float stepL;
  float x;
  float y;
  float z;
} simplifiedHitsInfo;

class AllPixRunAction : public G4UserRunAction
{

public:
  AllPixRunAction(AllPixDetectorConstruction *, TString, TString, TString, TString);
  ~AllPixRunAction();
  
public:
  void BeginOfRunAction(const G4Run* aRun);
  void EndOfRunAction(const G4Run* aRun);
  G4Run * GenerateRun();
  void GetPrimaryGeneratorMessenger(AllPixPrimaryGeneratorAction * AllPixAction) 
  {
    this->AllPixMessenger = AllPixAction->GetPrimaryGeneratorMessenger();
  }
  AllPixRun* ReturnAllPixRun(); //nalipour

  AllPixWriteROOTFile** writeROOTFile; //nalipour: To write MC in a ROOT file
  
private:

  AllPixDetectorConstruction * m_detectorPtr;
  TString m_dataset;
  TString m_tempdir;

  G4Timer* timer;
  AllPixRun * m_AllPixRun;

  // file for lcio
  ofstream * m_lciobridge_f;
  ofstream * m_lciobridge_dut_f;

  AllPixPrimaryGeneratorMessenger * AllPixMessenger;


  G4bool m_writeTPixTelescopeFilesFlag;
  G4bool m_writeEUTelescopeFilesFlag;
  G4bool m_writeMCROOTFilesFlag; //nalipour: Flag to write MC hits in a ROOT file

  string m_writeEUTelescopeFolder;
  G4int m_EUTelescopeRunNumber;
  
  G4int euRunNr;

  vector <AllPixWriter*> writerList;

  void initWriters();
  void writeEventToWriters(const G4Run* aRun);
  void closeWriters();

#ifdef HAVE_LCIO
  AllPixLCIOwriter * lcio;
#endif

  
};


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif /*AllPixRunAction_h*/
