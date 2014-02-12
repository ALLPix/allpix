/**
 * Author: John Idarraga <idarraga@cern.ch> , 2009
 */

#include <vector>
#include <iostream>

#include "AllPix_Frames_WriteToEntuple.h"
#include "allpix_dm.h"

// geometry
#include "ReadGeoDescription.hh"

static WriteToNtuple ** instance = 0;
static Int_t * indexToDetectorIdMap = 0;

WriteToNtuple::WriteToNtuple(TString prefix, TString dataSet, TString tempScratchDir, Int_t detID, TString openmode /* default "RECREATE" */){

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
	  m_ntupleFileName += "MPXNtuple_"+m_MPXDataSetNumber;	  
	}
	m_ntupleFileName += ".root";

	nt = new TFile(m_ntupleFileName, openmode);
	t2 = new TTree("MPXTree","Medi/TimePix data");

	// This instance of FrameStruct won't be stored
	// it'll be overriden when calling WriteToNtuple::fillVars
	m_frame = new FrameStruct(m_MPXDataSetNumber);
	t2->Branch("FramesData", "FrameStruct", &m_frame, 128000, 2);

}

WriteToNtuple::~WriteToNtuple(){

	// deleting FrameStruct instance
	delete m_frame;

}

/**
 *  Delivers one instance per detector
 *
 */
WriteToNtuple * WriteToNtuple::GetInstance(TString prefix, TString dataset, TString tempdir, Int_t nOfDetectors, Int_t detID, TString openmode /* RECREATE */) {

	TString tempDataset = dataset;

	if (instance == 0){

		instance = new WriteToNtuple * [nOfDetectors];
		indexToDetectorIdMap = new Int_t [nOfDetectors];

		// Geo description
		extern ReadGeoDescription * g_GeoDsc; // already loaded ! :)
		map<int, AllPixGeoDsc *> * geoMap = g_GeoDsc->GetDetectorsMap();
		map<int, AllPixGeoDsc *>::iterator detItr;

		int cntr = 0;
		for( detItr = geoMap->begin() ; detItr != geoMap->end() ; detItr++) {

			tempDataset = dataset;
			tempDataset += (*detItr).first; // append detector id

			instance[cntr] = new WriteToNtuple(prefix, tempDataset, tempdir, (*detItr).first, openmode);
			indexToDetectorIdMap[cntr] = (*detItr).first;

			cntr++;
		}

		if(nOfDetectors != cntr){
			std::cout << "[OUCH] ! number of detectors don't match WriteToNtuple::GetInstance " << std::endl;
			exit(1);
		}

	}

	// Search for the right index with detID (can't use a map here).
	for(Int_t i = 0 ; i < nOfDetectors ; i++){
		if(indexToDetectorIdMap[i] == detID)
			return instance[i];
	}

	// If I get here is because I couldn't find the instance associated to the detID
	std::cout << "[OUCH] ! det " << detID << " couldn't be found ... giving up." << std::endl;
	exit(1);

	return 0x0;
}

void WriteToNtuple::fillVars(FramesHandler * frameHandlerObj) {

	// Variables(class) in the Tree
	m_frame = frameHandlerObj->getFrameStructObject();
	// fill the Tree
	nt->cd();
	t2->Fill();
	// clean up
	frameHandlerObj->RewindAll();

}

void WriteToNtuple::closeNtuple()
{

	nt->cd();
	t2->Write();
	nt->Close();

}
