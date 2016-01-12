/**
 * Author: John Idarraga <idarraga@cern.ch> , 2009
 *
 */

#ifndef AllPix_Frames_WriteToEntuple_h
#define AllPix_Frames_WriteToEntuple_h 1

#include <map>

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

class WriteToNtuple {

public:

  WriteToNtuple(TString, TString, TString, Int_t, TString om="RECREATE");
  virtual ~WriteToNtuple();
  void fillVars(FramesHandler *);
  void closeNtuple();
  TString GetNtupleFileName(){return m_ntupleFileName;};
  static WriteToNtuple * GetInstance(TString, TString, TString, Int_t, Int_t, TString openmode = "RECREATE");
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


  ClassDef(WriteToNtuple,1)
};

#endif
