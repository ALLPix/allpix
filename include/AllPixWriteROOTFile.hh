/**
 *  Author Niloufar Alipour Tehrani <nalipour@cern.ch>
 */
#ifndef AllPixWriteROOTFile_h
#define AllPixWriteROOTFile_h 1

#include "G4Run.hh"

#include "AllPixTrackerHit.hh"
#include "AllPixRunAction.hh"
#include "AllPixMimosa26Digit.hh"
#include "AllPixDigitInterface.hh"
#include "ROOTDataFormat.hh"

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

class AllPixWriteROOTFile
{

public:

  AllPixWriteROOTFile(Int_t detID, TString path);
  void AllPixWriteROOTFillTree();
  void AllPixCloseROOTFile();
  void SetVectors(ROOTDataFormat* d);
  virtual ~AllPixWriteROOTFile();


  
  TFile* file;
  TTree* tree;
  TBranch* branch;

  Int_t detectorID;

  std::vector<Int_t> posX;
  std::vector<Int_t> posY;
  std::vector<Double_t> energyTotal;
  std::vector<Int_t> TOT;
  std::vector<Double_t> energyMC;
  std::vector<Double_t> posX_WithRespectToPixel;
  std::vector<Double_t> posY_WithRespectToPixel;
  std::vector<Double_t> posZ_WithRespectToPixel;
  /*
  Int_t nHits_MC;
  Int_t eventNB_MC;
  std::vector<Int_t> posX_MC;
  std::vector<Int_t> posY_MC;
  std::vector<Double_t> energy_MC;

  Int_t nHits;
  Int_t eventNB;
  std::vector<Int_t> posX;
  std::vector<Int_t> posY;
  std::vector<Double_t> energy;
  std::vector<Int_t> TOT;
  */

private:
};

#endif
