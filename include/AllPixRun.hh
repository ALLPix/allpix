/**
 *  Author John Idarraga <idarraga@cern.ch>
 */
#ifndef AllPixRun_h
#define AllPixRun_h 1

#include "G4Run.hh"

#include "AllPixTrackerHit.hh"
#include "AllPixRunAction.hh"
#include "AllPixMimosa26Digit.hh"
#include "AllPixDigitInterface.hh"
#include "AllPixWriteROOTFile.hh" //nalipour: Write MC hits in a ROOT file
#include "ROOTDataFormat.hh" //nalipour: Write MC hits in a ROOT file

#include <TROOT.h>
#include <TFile.h>
#include <TBranch.h>
#include <TTree.h>
#include <TChain.h>
#include <TString.h>

#include <vector>
#include <string>
#include <map>

#include <iostream>
#include <fstream>

using namespace std;

class FramesHandler;
class WriteToNtuple;
class SimpleHits;
class AllPixDetectorConstruction;

#ifdef _EUTELESCOPE
#define __magic_trigger_cntr_ack 4 // 4 scintillators with 4 primary particle hits
#endif

class AllPixRun : public G4Run {

public:

  AllPixRun(AllPixDetectorConstruction *, TString, TString, TString, G4bool, G4bool MCWriteFlag);
  virtual ~AllPixRun();

  virtual void RecordEvent(const G4Event*);
  void RecordHits(const G4Event*);
  void RecordDigits(const G4Event*);
  void RecordTelescopeDigits(const G4Event*);
  void SetLCIOBridgeFileDsc(ofstream * f, ofstream * fdut)
  { m_lciobridge_f = f; m_lciobridge_dut_f = fdut;};

  // filling frames
  void FillFramesNtuple(const G4Run *);
  // fill timepix telescope files
  void FillTelescopeFiles(const G4Run *, G4String, G4bool, G4bool);


  //nalipour: MC hits
  void FillROOTFiles(AllPixWriteROOTFile** rootFiles); //fills the ROOT files
  G4int return_detIdToIndex(G4int detId)
  {
    return m_detIdToIndex[detId];
  }
  //void RecordHitsForROOTFiles(const G4Event* evt); //Record the MC hits to save in the ROOT file
  //void RecordHitsForROOTFiles_withChargeSharing(const G4Event* evt); // Hits with charge sharing values
  void RecordDigits_all(const G4Event* evt);
  
  void GetEUTelescopeEvent(G4int runnr, const G4Run* aRun);


private:

  AllPixDetectorConstruction * m_detectorPtr;
  AllPixTrackerHitsCollection * m_hitsCollection;

  // Hits
  SimpleHits ** m_storableHits;
  vector<ROOTDataFormat*> MC_ROOT_data; //nalipour: MC hits (for each detector) 

  // map index in frames handler to det Id
  map<int, int> m_detIdToIndex;
  // map pixel X, Y, TOT per run
  map<int,vector<vector<vector<int> > > > m_data;

  // Frames ntuple  --> not storing whole Digits
  //  building frames from digits
  FramesHandler ** m_frames;
  TString m_datasetDigits;
  TString m_datasetHits;
  TString m_tempdir;

  // Information about ntuples

  // Same as number of digitizers
  G4int m_nOfDetectors;

  // Number of sensitive detectors.  Can be more than the number of
  // digitizers if the user decides to have other volumes as SDs
  G4int m_nOfSD;

  // coming from AllPixRunAction
  ofstream * m_lciobridge_f;
  ofstream * m_lciobridge_dut_f;

  G4String m_outputFilePrefix;

  time_t m_runTime;

  G4bool m_writeTPixTelescopeFilesFlag;
  G4bool m_writeMCFilesFlag; //nalipour: MC hits

};

#endif
