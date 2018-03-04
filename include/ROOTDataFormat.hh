/**
 *  Author Niloufar Alipour Tehrani <nalipour@cern.ch>
 */
#ifndef ROOTDataFormat_h
#define ROOTDataFormat_h 1


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

class ROOTDataFormat
{

public:

  ROOTDataFormat(Int_t detId);
  virtual ~ROOTDataFormat();

  void add_posX(Int_t pos) {posX.push_back(pos);};
  void add_posY(Int_t pos) {posY.push_back(pos);};
  void add_energyTotal(Double_t energy) {energyTotal.push_back(energy);};
  void add_TOT(Int_t tot) {TOT.push_back(tot);};
  void add_energyMC(Double_t energy) {energyMC.push_back(energy);};
  void add_posX_WithRespectToPixel(Double_t pos) {posX_WithRespectToPixel.push_back(pos);};
  void add_posY_WithRespectToPixel(Double_t pos) {posY_WithRespectToPixel.push_back(pos);};
  void add_posZ_WithRespectToPixel(Double_t pos) {posZ_WithRespectToPixel.push_back(pos);};
  void add_POLangle(Double_t ang) {POLang.push_back(ang);};
  void add_AZMangle(Double_t ang) {AZMang.push_back(ang);};
  void add_TruthEntryLocalY(Double_t pos) {TruthEntryLocalY.push_back(pos);};
  void add_TruthEntryLocalX(Double_t pos) {TruthEntryLocalX.push_back(pos);};
  void add_TruthExitLocalX(Double_t pos) {TruthExitLocalX.push_back(pos);};
  void add_TruthExitLocalY(Double_t pos) {TruthExitLocalY.push_back(pos);};
  void add_deltaEfrac(Double_t energy) {deltaEfrac.push_back(energy);};
  void add_path_length_first_pixel(Double_t pos) {path_length_first_pixel.push_back(pos);};
  void add_InitE(Double_t pos) {InitE.push_back(pos);};
  void add_InitID(Int_t pos) {InitID.push_back(pos);};

  vector<Double_t> get_POLangle() {return POLang;};
  vector<Double_t> get_AZMangle() {return AZMang;};
  vector<Double_t> get_TruthEntryLocalY() {return TruthEntryLocalY;};
  vector<Double_t> get_TruthEntryLocalX() {return TruthEntryLocalX;};
  vector<Double_t> get_TruthExitLocalY() {return TruthExitLocalY;};
  vector<Double_t> get_TruthExitLocalX() {return TruthExitLocalX;};
  vector<Int_t> get_posX() {return posX;};
  vector<Int_t> get_posY() {return posY;};
  vector<Double_t> get_energyTotal() {return energyTotal;};
  vector<Int_t> get_TOT() {return TOT;};
  vector<Double_t> get_energyMC() {return energyMC;};
  vector<Double_t> get_posX_WithRespectToPixel() {return posX_WithRespectToPixel;};
  vector<Double_t> get_posY_WithRespectToPixel() {return posY_WithRespectToPixel;};
  vector<Double_t> get_posZ_WithRespectToPixel() {return posZ_WithRespectToPixel;};
  vector<Double_t> get_deltaEfrac(){ return deltaEfrac;};
  vector<Double_t> get_path_length_first_pixel(){ return path_length_first_pixel;};
  vector<Double_t> get_InitE(){ return InitE;};
  vector<Int_t> get_InitID() {return InitID;};
/*
  void set_posX_MC(vector<Int_t> vec);
  void set_posY_MC(vector<Int_t> vec);
  void set_energy_MC(vector<Double_t> vec);

  vector<Int_t> get_posX_MC();
  vector<Int_t> get_posY_MC();
  vector<Double_t> get_energy_MC();

  void add_posX_MC(Int_t pos);
  void add_posY_MC(Int_t pos);
  void add_energy_MC(Double_t energy);


  void add_posX(Int_t pos);
  void add_posY(Int_t pos);
  void add_energy(Double_t energy);
  void add_TOT(Int_t TOT);

  vector<Int_t> get_posX();
  vector<Int_t> get_posY();
  vector<Double_t> get_energy();
  vector<Int_t> get_TOT();
  */

private:  
  //MC hits
  Int_t detectorID;
  std::vector<Double_t> deltaEfrac;
  std::vector<Double_t> POLang;
  std::vector<Double_t> AZMang;
  std::vector<Int_t> posX;
  std::vector<Int_t> posY;
  std::vector<Double_t> energyTotal;
  std::vector<Int_t> TOT;
  std::vector<Double_t> energyMC;
  std::vector<Double_t> posX_WithRespectToPixel;
  std::vector<Double_t> posY_WithRespectToPixel;
  std::vector<Double_t> posZ_WithRespectToPixel;
  std::vector<Double_t> TruthEntryLocalY;
  std::vector<Double_t> TruthEntryLocalX; 
  std::vector<Double_t> TruthExitLocalY;
  std::vector<Double_t> TruthExitLocalX;
  std::vector<Double_t> path_length_first_pixel;
  std::vector<Double_t> InitE;
  std::vector<Int_t> InitID;
/*  Int_t trackID;
  Int_t pdgIdTrack;
  Int_t parentID;
  Double_t posWithRespectToPixel_X;
  Double_t posWithRespectToPixel_Y;
  Double_t posWithRespectToPixel_Z;*/
  /*
  Int_t nHits_MC;
  Int_t eventNB_MC;
  std::vector<Int_t> posX_MC;
  std::vector<Int_t> posY_MC;
  std::vector<Double_t> energy_MC;
  std::vector<Double_t> TruthEntryLocalY;
  std::vector<Double_t> TruthEntryLocalX;

  //Hits considering charge sharing
  Int_t nHits;
  Int_t eventNB;
  std::vector<Int_t> posX;
  std::vector<Int_t> posY;
  std::vector<Double_t> energy;
  std::vector<Int_t> TOT;
  */
};

#endif
