/**
 * Author: Niloufar Alipour Tehrani <nalipour@cern.ch> , 2014
 */


#include <iostream>

#include "AllPix_Frames_WriteToEntuple_RD53.h"
#include "allpix_dm.h"

// geometry
#include "ReadGeoDescription.hh"

static WriteToNtupleRD53 ** instanceRD53 = 0;
static Int_t * indexToDetectorIdMapRD53 = 0;

WriteToNtupleRD53::WriteToNtupleRD53(TString prefix, TString dataSet, TString tempScratchDir, Int_t detID, TString openmode /* default "RECREATE" */){

	m_MPXDataSetNumber = dataSet;
	m_ntupleFileName = "";
	m_detID = detID;

	// change dir if requested
	if(tempScratchDir.Length() > 0){
		m_ntupleFileName += tempScratchDir;
		m_ntupleFileName += "/";
	}


	if(prefix.Length() > 0){
	  m_ntupleFileName += prefix;
	  m_ntupleFileName += "_";
	  m_ntupleFileName += m_MPXDataSetNumber;
	} else {
	  m_ntupleFileName += "RD53_"+m_MPXDataSetNumber;	  
	}
	m_ntupleFileName += ".root";

	nt = new TFile(m_ntupleFileName, openmode);
	t2 = new TTree("tree","tree data");

	// This instance of FrameStruct won't be stored
	// it'll be overriden when calling WriteToNtuple::fillVars
	m_frame = new FrameStruct(m_MPXDataSetNumber);
	//t2->Branch("FramesData", "FrameStruct", &m_frame, 128000, 2); //nalipour
	// t2->Branch("nHits", &nHits, "nHits/I");
	// t2->Branch("posX", &posX, "posX[nHits]/I");
	// t2->Branch("posY", &posY, "posY[nHits]/I");
	// t2->Branch("energy", &energy, "energy[nHits]/I");


	// t2->Branch("MC_nHits", &MC_nHits, "MC_nHits/I");
	// t2->Branch("MC_posX", &MC_posX, "MC_posX[nHits]/I");
	// t2->Branch("MC_posY", &MC_posY, "MC_posY[nHits]/I");
	// t2->Branch("MC_energy", &MC_energy, "MC_energy[nHits]/I");

	t2->Branch("nHits", &nHits);
	t2->Branch("posX", &posX);
	t2->Branch("posY", &posY);
	t2->Branch("energy", &energy);


	t2->Branch("MC_nHits", &MC_nHits);
	t2->Branch("MC_posX", &MC_posX);
	t2->Branch("MC_posY", &MC_posY);
	t2->Branch("MC_energy", &MC_energy);



}

WriteToNtupleRD53::~WriteToNtupleRD53(){

	// deleting FrameStruct instance
	delete m_frame;

}

/**
 *  Delivers one instance per detector
 *
 */
WriteToNtupleRD53 * WriteToNtupleRD53::GetInstanceRD53(TString prefix, TString dataset, TString tempdir, Int_t nOfDetectors, Int_t detID, TString openmode /* RECREATE */) {

	TString tempDataset = dataset;

	if (instanceRD53 == 0){

		instanceRD53 = new WriteToNtupleRD53 * [nOfDetectors];
		indexToDetectorIdMapRD53 = new Int_t [nOfDetectors];

		// Geo description
		extern ReadGeoDescription * g_GeoDsc; // already loaded ! :)
		map<int, AllPixGeoDsc *> * geoMap = g_GeoDsc->GetDetectorsMap();
		map<int, AllPixGeoDsc *>::iterator detItr;

		int cntr = 0;
		for( detItr = geoMap->begin() ; detItr != geoMap->end() ; detItr++) {

			tempDataset = dataset;
			tempDataset += (*detItr).first; // append detector id

			instanceRD53[cntr] = new WriteToNtupleRD53(prefix, tempDataset, tempdir, (*detItr).first, openmode);
			indexToDetectorIdMapRD53[cntr] = (*detItr).first;

			cntr++;
		}

		if(nOfDetectors != cntr){
			std::cout << "[OUCH] ! number of detectors don't match WriteToNtupleRD53::GetInstance " << std::endl;
			exit(1);
		}

	}

	// Search for the right index with detID (can't use a map here).
	for(Int_t i = 0 ; i < nOfDetectors ; i++){
		if(indexToDetectorIdMapRD53[i] == detID)
			return instanceRD53[i];
	}

	// If I get here is because I couldn't find the instance associated to the detID
	std::cout << "[OUCH] ! det " << detID << " couldn't be found ... giving up." << std::endl;
	exit(1);

	return 0x0;
}

void WriteToNtupleRD53::fillVarsRD53(FramesHandler * frameHandlerObj) {

	// Variables(class) in the Tree
	m_frame = frameHandlerObj->getFrameStructObject();


	std::map<std::pair<Int_t, Int_t>, Double_t > Energy_MC_map=m_frame->return_RD53_E_MC();
	std::map<std::pair<Int_t, Int_t>, Double_t > Energy_map=m_frame->return_RD53_E();

	nHits=Energy_map.size();
	G4cout << "Nilou: nHits=" << nHits << G4endl;
	map<pair<Int_t, Int_t>, Double_t >::iterator pCItr = Energy_map.begin();
	for( ; pCItr != Energy_map.end() ; pCItr++)
	  {
	    G4cout << "Nilou: posX=" << (*pCItr).first.first << ", posY=" << (*pCItr).first.second << ", energy=" << (*pCItr).second << G4endl;
	    posX.push_back((*pCItr).first.first);
	    posY.push_back((*pCItr).first.second);
	    energy.push_back((*pCItr).second);
	  }


	MC_nHits=Energy_MC_map.size();
	map<pair<Int_t, Int_t>, Double_t >::iterator pCItrMC = Energy_MC_map.begin();
	for( ; pCItrMC != Energy_MC_map.end() ; pCItrMC++)
	  {
	    MC_posX.push_back((*pCItrMC).first.first);
	    MC_posY.push_back((*pCItrMC).first.second);
	    MC_energy.push_back((*pCItrMC).second);
	  }


	// fill the Tree
	nt->cd();
	t2->Fill();
	// clean up
	frameHandlerObj->RewindAll();

	posX.clear();
	posY.clear();
	energy.clear();
	
	MC_posX.clear();
	MC_posY.clear();
	MC_energy.clear();
}

void WriteToNtupleRD53::closeNtuple()
{

	nt->cd();
	t2->Write();
	nt->Close();

}
