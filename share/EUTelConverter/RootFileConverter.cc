#include <fstream>
#include <iostream>
#include <sstream>
#include <cmath>
#include <vector>
#include <string>

//#include <TVector3.h>
#include <TROOT.h>
//#include <TChain.h>
#include <TFile.h>
//#include <TH2.h>
//#include <TStyle.h>
//#include <TCanvas.h>
#include <TTree.h>
#include <TBranch.h>
//#include <TObjArray.h>
//#include <TMath.h>

#include <lcio.h>
#include <IO/LCWriter.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/LCEventImpl.h>
#include <IMPL/LCRunHeaderImpl.h>
#include <IMPL/LCGenericObjectImpl.h>
#include <IMPL/TrackerDataImpl.h>
#include <UTIL/CellIDEncoder.h>

#include "AllPix_Hits_WriteToEntuple.h"
//#include "AllPix_Hits_WriteToEntuple.cc"
#include "allpix_dm.h"
//#include "AllPix_Frames_WriteToEntuple.h"
//#include "AllPix_Frames_WriteToEntuple.cc"

using namespace lcio;

void print_info();

int main(int argc, char* argv[]) {

	std::string output_file, file_path;
	int n_sensors, n_DUTs;
	std::vector<int> DUT_IDs;

	if (argc < 4) print_info(); //if not enough arguments are given, print user info
	else {

		sscanf(argv[1], "%d", &n_sensors);
		sscanf(argv[2], "%d", &n_DUTs);
		sscanf(argv[3], "%s", &output_file);

		if ((argc > 4) && (n_DUTs <= 0)) sscanf(argv[4], "%s", &file_path);
		else {

			DUT_IDs = std::vector<int>(n_DUTs);

			for (int i = 0; i < n_DUTs; i++) {

				sscanf(argv[4+i], "%d", &DUT_IDs[i]);
			}

			if (argc > 4+n_DUTs) sscanf(argv[4+n_DUTs], "%s", &file_path);
		}
	}

	TFile** true_hit_files = new TFile*[n_sensors+n_DUTs];
	TTree** true_hit_trees = new TTree*[n_sensors+n_DUTs];
	TFile** frame_files = new TFile*[n_sensors+n_DUTs];
	TTree** frame_trees = new TTree*[n_sensors+n_DUTs];

	std::ostringstream oss;

	for (int i = 0; i < n_sensors; i++) {

		oss.str(std::string());//resets oss to an empty string
		oss << file_path << "Hits_BoxSD_" << i+300 << "_HitsCollection.root";

		true_hit_files[i] = new TFile(oss.str().c_str());
		true_hit_trees[i] = (TTree*) true_hit_files[i]->Get("AllPixHits");

		oss.str(std::string());
		oss << file_path << "MPXNtuple_allPix_det_" << i+300 << ".root";

		frame_files[i] = new TFile(oss.str().c_str());
		frame_trees[i] = (TTree*) frame_files[i]->Get("MPXTree");
	}

	for (int i = 0; i < n_DUTs; i++) {

		oss.str(std::string());
		oss << file_path << "Hits_BoxSD_" << DUT_IDs[i] << "_HitsCollection.root";

		true_hit_files[i+n_sensors] = new TFile(oss.str().c_str());
		true_hit_trees[i+n_sensors] = (TTree*)true_hit_files[i+n_sensors]->Get("AllPixHits");

		oss.str(std::string());
		oss << file_path << "MPXNtuple_allPix_det_" << DUT_IDs[i] << ".root";

		frame_files[i+n_sensors] = new TFile(oss.str().c_str());
		frame_trees[i+n_sensors] = (TTree*) frame_files[i+n_sensors]->Get("MPXTree");
	}

	int n_entries = true_hit_trees[0]->GetEntries();

	SimpleHits root_hit;
	true_hit_trees[0]->SetBranchAddress("SimpleHits", &root_hit);		
	true_hit_trees[0]->GetEntry(0);
	int run_number = root_hit.run;

	std::string detector_name = "EUTelescope";

	LCWriter* writer = LCFactory::getInstance()->createLCWriter();
	writer->open(output_file);

	LCRunHeaderImpl* run_header = new LCRunHeaderImpl; 
	run_header->setRunNumber(run_number);
	run_header->setDetectorName(detector_name);
	run_header->parameters().setValue  ("GeoID"           , 0);
	run_header->parameters().setValues ("MaxX"            , std::vector<int>(6,1151));
	run_header->parameters().setValues ("MaxY"            , std::vector<int>(6,575));
	run_header->parameters().setValues ("MinX"            , std::vector<int>(6,0));
	run_header->parameters().setValues ("MinY"            , std::vector<int>(6,0));
	run_header->parameters().setValue  ("NoOfDetector"    , 6);
	run_header->parameters().setValues ("AppliedProcessor", std::vector<std::string>(1,""));
	run_header->parameters().setValue  ("DAQHWName"       , "EUDRB");
	run_header->parameters().setValue  ("DAQSWName"       , "EUDAQ");
	run_header->parameters().setValue  ("DataType"        , "SimData");
	run_header->parameters().setValue  ("DateTime"        , "24.12.2000  23:59:59.000000000");
	run_header->parameters().setValue  ("EUDRBDet"        , "MIMOSA26");
	run_header->parameters().setValue  ("EUDRBMode"       , "ZS2");
   	writer->writeRunHeader(run_header);
	delete run_header;

	for (int i = 0; i < n_entries; i++) {

		LCEventImpl* event = new LCEventImpl();
		event->setRunNumber(run_number);
		event->setEventNumber(i);
		event->setDetectorName(detector_name);
		event->setTimeStamp( (int) time(NULL)*1000000000.0);
		event->parameters().setValue("EventType", 2);

		if (i == 0) {

			LCCollectionVec* eudrb_setup = new LCCollectionVec(LCIO::LCGENERICOBJECT);

			//collection parameters
			eudrb_setup->parameters().setValue("DataDescription", "type:i,mode:i,spare1:i,spare2:i,spare3:i");
			eudrb_setup->parameters().setValue("TypeName", "Setup Description");

			//create one setup object per Telescope plane
			for (int n = 0; n < n_sensors; n++) {

				LCGenericObjectImpl* setup_obj = new LCGenericObjectImpl(5,0,0);
				setup_obj->setIntVal(0, 102);
				setup_obj->setIntVal(1, 101);
				eudrb_setup->addElement(setup_obj);
			}

			event->addCollection(eudrb_setup, "eudrbSetup");
		}

		//ID encoder info
		std::string encoding_string = "sensorID:7,sparsePixelType:5";

		//Telescope data collection
		LCCollectionVec* tracker_data_coll = new LCCollectionVec(LCIO::TRACKERDATA);
		CellIDEncoder<TrackerDataImpl> id_encoder_telescope(encoding_string, tracker_data_coll);// = new CellIDEncoder<TrackerDataImpl>(encoding_string, tracker_data_coll);

		//fil telescope collection
		for (int j = 0; j < n_sensors; j++) {

			//to-do: initiate the planeData and idEncoder_Telesope

			true_hit_trees[j]->SetBranchAddress("SimpleHits", &root_hit);
			true_hit_trees[j]->GetEntry(i);

			TrackerDataImpl* plane_data = new TrackerDataImpl();

			id_encoder_telescope.reset();
			id_encoder_telescope["sensorID"] = j + 300;
			id_encoder_telescope["sparsePixelType"] = 2;
			id_encoder_telescope.setCellID(plane_data);

			//loop over hits
			std::vector<float> charge_vec;
			for (int k = 0; k < root_hit.pos.size(); k++) {

				charge_vec.push_back(root_hit.pos[k].X());
				charge_vec.push_back(root_hit.pos[k].Y());
				charge_vec.push_back(root_hit.edep[k]);
				charge_vec.push_back(0);
			}

			plane_data->setChargeValues(charge_vec);
			tracker_data_coll->addElement(plane_data);

			charge_vec.clear();
		}

		event->addCollection(tracker_data_coll, "true_hits_m26");

		//write DUT true hit data collection if there are DUTs present
		if (n_DUTs > 0) {

			LCCollectionVec* DUT_data_coll = new LCCollectionVec(LCIO::TRACKERDATA);
			CellIDEncoder<TrackerDataImpl> id_encoder_DUT(encoding_string, DUT_data_coll);

			for (int j = 0; j < n_DUTs; j++) {

				true_hit_trees[j+n_sensors]->SetBranchAddress("SimpleHits", &root_hit);
				true_hit_trees[j+n_sensors]->GetEntry(i);

				TrackerDataImpl* plane_data = new TrackerDataImpl();

				id_encoder_DUT.reset();
				id_encoder_DUT["sensorID"] = DUT_IDs[j];
				id_encoder_DUT["sparsePixelType"] = 2;
				id_encoder_DUT.setCellID(plane_data);

				std::vector<float> charge_vec;
				for (int k = 0; k < root_hit.pos.size(); k++) {

					charge_vec.push_back(root_hit.pos[k].X());
					charge_vec.push_back(root_hit.pos[k].Y());
					charge_vec.push_back(root_hit.edep[k]);
					charge_vec.push_back(0);
				}

				plane_data->setChargeValues(charge_vec);
				DUT_data_coll->addElement(plane_data);

				charge_vec.clear();
			}

			event->addCollection(DUT_data_coll, "true_hits_DUT");
		}


		writer->writeEvent(event);

		delete event;
	}

	writer->flush();
	writer->close();

	for (int i = 0; i < n_sensors; i++) {

		true_hit_files[i]->Close();
		delete true_hit_trees[i];
		delete true_hit_files[i];
	}

	delete [] true_hit_trees;
	delete [] true_hit_files;

	return 0;
}

void print_info() {

	std::cout << "Converts allpix generated root files into LCIO format" << std::endl;
	std::cout << "Usage: compile and run with the arguments in the following order:" << std::endl;
	std::cout << "<number of sensors> <number of DUTs> <output file name> <list of DUT sensor IDs> <path to folder containing input root files>" << std::endl;
	std::cout << "if no DUTs are present, inputing 0 for the number of DUTs means that the list of DUTs is no longer required. The fourth argument will be interpreted as the path to folder containing input root files." << std::endl;
}
