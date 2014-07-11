/**
 * Author: Niloufar Alipour Tehrani <nalipour@cern.ch> , 2014
 *
 */

#ifndef AllPix_Frames_WriteToEntuple_RD53_h
#define AllPix_Frames_WriteToEntuple_RD53_h 1

#include <map>
#include <vector>
#include <TROOT.h>
#include <TChain.h>
#include <TString.h>
#include <TFile.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TBranch.h>
#include <TObjArray.h>
#include <TMath.h>

#include "allpix_dm_consts.h"

/** Implements what is needed to write the 
 *  frames info and the MetaData 
 *  to an Ntuple ROOT file.
 */

class FramesHandler;
class FrameStruct;

class WriteToNtupleRD53 {

public:

  WriteToNtupleRD53(TString, TString, TString, Int_t, TString om="RECREATE");
  ~WriteToNtupleRD53();
  void fillVarsRD53(FramesHandler *);
  void closeNtuple();
  TString GetNtupleFileName(){return m_ntupleFileName;};
  static WriteToNtupleRD53 * GetInstanceRD53(TString, TString, TString, Int_t, Int_t, TString openmode = "RECREATE");
  TTree * GetTree(){return t2;};
  //void SetBranchAddress(FrameStruct * fs){b2->SetAddress(fs);};
  Int_t GetDetectorId(){return m_detID;};

private:

  FrameStruct * m_frame;

  TH2 * h1;
  TFile * nt;
  TTree * t2;
  TBranch * b2;
  TString m_MPXDataSetNumber;
  TString m_ntupleFileName;
  Int_t m_detID;

  std::vector<Int_t> posX;
  std::vector<Int_t> posY;
  std::vector<Double_t> energy;
  Int_t nHits;
  std::vector<int> MC_posX;
  std::vector<int> MC_posY;
  std::vector<double> MC_energy;
  Int_t MC_nHits;


  ClassDef(WriteToNtupleRD53, 1)
};

#endif
