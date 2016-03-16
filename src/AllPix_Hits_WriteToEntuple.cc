/**
 * Author: John Idarraga <idarraga@cern.ch> , 2009
 *
 */

#include <vector>
#include <iostream>

#include "AllPix_Hits_WriteToEntuple.h"

static Hits_WriteToNtuple ** instance_hit = 0;
int g_instance_hit_Cntr = 0;

Hits_WriteToNtuple::Hits_WriteToNtuple(TString prefix, TString dataSet, TString tempScratchDir,
		Int_t /*detID*/, TString openmode /* default "RECREATE" */){

	m_MPXDataSetNumber = dataSet;
	m_ntupleFileName = "";

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
		m_ntupleFileName += "Hits_"+m_MPXDataSetNumber;
	}
	m_ntupleFileName += ".root";

	nt = new TFile(m_ntupleFileName, openmode);
	t2 = new TTree("AllPixHits", dataSet);

	m_storableHits = new SimpleHits;
	m_storableHits->Rewind();
	t2->Branch("SimpleHits", "SimpleHits", &m_storableHits);

}

/**
 *  Delivers one instance per detector
 *
 */
Hits_WriteToNtuple * Hits_WriteToNtuple::GetInstance(TString prefix, TString dataset, TString tempdir, Int_t nOfInstances, Int_t detID, TString openmode /* RECREATE */){

	if(detID >= nOfInstances){
		std::cout << "[ERROR] (in Hits_WriteToNtuple::GetInstance) you are asking for the " << detID << "th hits ntuple"
				<< " out of " << nOfInstances << " aborting ..." << std::endl;
		std::cout << "        (note: Indexing starts at 0)" << std::endl;
		exit(1);
	}

	TString tempDataset = dataset;
	if (instance_hit == 0){
		// create pointers for all instances
		instance_hit = new Hits_WriteToNtuple * [nOfInstances];
	}

	// but instanciate only once at a time
	if(g_instance_hit_Cntr < nOfInstances)
	{
		std::cout << "Creating hits file : "
				<< "\"" << tempDataset
				<< "\"" << std::endl;
		instance_hit[detID] = new Hits_WriteToNtuple(prefix, tempDataset, tempdir, detID, openmode);
		g_instance_hit_Cntr++;
	}

	return instance_hit[detID];
}

void Hits_WriteToNtuple::fillVars(SimpleHits * hits_i){

	// Variables(class) in the Tree
	m_storableHits = hits_i;

	// fill the Tree
	nt->cd();
	t2->Fill();

	// clean up
	m_storableHits->Rewind();

}

void Hits_WriteToNtuple::closeNtuple()
{

	nt->cd();
	t2->Write();
	nt->Close();

}


ClassImp(SimpleHits)

SimpleHits::SimpleHits(){
	Rewind();
}

void SimpleHits::Rewind(){

	edepTotal = 0.;

	interactions.clear();
	pos.clear();
	pdgId.clear();
	edep.clear();
	trackId.clear();
	parentId.clear();
	trackVolumeName.clear();
	parentVolumeName.clear();

	event = 0;
	run = 0;
	kinEParent = 0.;

	trackID.clear();
	trackFate.clear();
	eOut.clear();
}
