/**
 *  Author Niloufar Alipour Tehrani <nalipour@cern.ch>
 */

#include "AllPixWriteROOTFile.hh"

AllPixWriteROOTFile::AllPixWriteROOTFile(Int_t detID, TString path)
{
  detectorID=detID;
  //G4cout << "nalipour AllPixWriteROOTFile" << G4endl;
  TString fileName=path+"/RD53_"+(TString)Form("%d", detID)+".root";
  file = new TFile(fileName, "RECREATE");
  tree = new TTree("tree","tree data");

  tree->Branch("posX", &posX);
  tree->Branch("posY", &posY);
  tree->Branch("energyTotal", &energyTotal);
  tree->Branch("TOT", &TOT);
  tree->Branch("energyMC", &energyMC);
  tree->Branch("posX_WithRespectToPixel", &posX_WithRespectToPixel);
  tree->Branch("posY_WithRespectToPixel", &posY_WithRespectToPixel);
  tree->Branch("posZ_WithRespectToPixel", &posZ_WithRespectToPixel);
/*
  tree->Branch("nHits_MC", &nHits_MC);
  tree->Branch("posX_MC", &posX_MC);
  tree->Branch("posY_MC", &posY_MC);
  tree->Branch("energy_MC", &energy_MC);

  tree->Branch("nHits", &nHits);
  tree->Branch("posX", &posX);
  tree->Branch("posY", &posY);
  tree->Branch("energy", &energy);
  tree->Branch("TOT", &TOT);
  */
}

void AllPixWriteROOTFile::AllPixWriteROOTFillTree()
{
  file->cd();
  tree->Fill();

  posX.clear();
  posY.clear();
  energyTotal.clear();
  TOT.clear();
  energyMC.clear();
  posX_WithRespectToPixel.clear();
  posY_WithRespectToPixel.clear();
  posZ_WithRespectToPixel.clear();
  /*
  posX_MC.clear();
  posY_MC.clear();
  energy_MC.clear();
  

  posX.clear();
  posY.clear();
  energy.clear();
  TOT.clear();
  */
}

void AllPixWriteROOTFile::SetVectors(ROOTDataFormat* d)
{
	//G4cout << "nalipour: SetVectors " << G4endl;
	posX=d->get_posX();
	posY=d->get_posY();
	energyTotal=d->get_energyTotal();
	energyMC=d->get_energyMC();
	TOT=d->get_TOT();
	posX_WithRespectToPixel=d->get_posX_WithRespectToPixel();
	posY_WithRespectToPixel=d->get_posY_WithRespectToPixel();
	posZ_WithRespectToPixel=d->get_posZ_WithRespectToPixel();
}

void AllPixWriteROOTFile::AllPixCloseROOTFile()
{
  G4cout << "nalipour AllPixCloseROOTFile: detectorID=" << detectorID << G4endl;
  file->cd();
  tree->Write();
  file->Close(); 
}
AllPixWriteROOTFile::~AllPixWriteROOTFile()
{
}
