/**
 * Author: John Idarraga <idarraga@cern.ch> , 2009
 *
 */

#ifndef AllPix_Hits_WriteToEntuple_h
#define AllPix_Hits_WriteToEntuple_h 1

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
#include <TVector3.h>

#include <vector>
#include <string>

using namespace std;

/**
 *  Simple hits holder
 *   The idea is not to store the HitCollection
 *   but something simpler to read
 */

class SimpleHits {

 public:
  SimpleHits();
  virtual ~SimpleHits(){};
  void Rewind();

  // total energy (all hits in event)
  double edepTotal;
  //
  double kinEParent;
  // interactions in hit (ioni, msc, phot, etc ..)
  vector<string> interactions;
  // position of hit
  vector<TVector3> pos;
  // pdgId of the particle producing the hit
  vector<int> pdgId;
  // energy deposition per hit
  vector<float> edep;
  // Id in the shower
  vector<int> trackId;
  // Parent in the shower
  vector<int> parentId;
  // Name of volume for this hit
  vector<string> trackVolumeName;
  // Name of volume where particle was created
  vector<string> parentVolumeName;
  // event
  int event;
  // run
  int run;

  ClassDef(SimpleHits, 1)
};

/** Implements what is needed to write 
 *  hits to an Ntuple ROOT file.
 */

class Hits_WriteToNtuple {

public:

  Hits_WriteToNtuple(TString, TString, TString, Int_t, TString om="RECREATE");
  virtual ~Hits_WriteToNtuple(){};
  void fillVars(SimpleHits *);
  void closeNtuple();
  TString GetNtupleFileName(){return m_ntupleFileName;};
  static Hits_WriteToNtuple * GetInstance(TString, TString, TString, Int_t, Int_t, TString openmode = "RECREATE");
  TTree * GetTree(){return t2;};
  //void SetBranchAddress(FrameStruct * fs){b2->SetAddress(fs);};
  //Int_t GetDetectorId(){return m_detID;};

private:
  
  TFile * nt;
  TTree * t2;

  //Int_t m_detID;
  TString m_MPXDataSetNumber;
  TString m_ntupleFileName;

  // to ntuple
  SimpleHits * m_storableHits;

  //ClassDef(Hits_WriteToNtuple,1)
};

#endif
