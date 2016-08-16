#include <fstream>
#include <iostream>
#include <sstream>
#include <cmath>
#include <vector>
#include <string>
#include <map>

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>

#include <lcio.h>
#include <IO/LCWriter.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/LCEventImpl.h>
#include <IMPL/LCRunHeaderImpl.h>
#include <IMPL/LCGenericObjectImpl.h>
#include <IMPL/TrackerDataImpl.h>
#include <IMPL/TrackerHitImpl.h>
#include <UTIL/CellIDEncoder.h>

#include "AllPix_Hits_WriteToEntuple.h"
#include "AllPix_Frames_WriteToEntuple.h"
#include "allpix_dm.h"

using namespace lcio;

void print_info();

int main(int argc, char* argv[]) {

	std::string output_file, file_path;
	int n_sensors, n_DUTs;
	std::vector<int> DUT_IDs;

	if (argc < 4) {

		print_info();//if not enough arguments are given, print user info
		return 0;
	}
	else {

		sscanf(argv[1], "%d", &n_sensors);
		sscanf(argv[2], "%d", &n_DUTs);
		output_file = argv[3];

		if ((argc > 4) && (n_DUTs <= 0)) file_path = argv[4];
		else {

			DUT_IDs = std::vector<int>(n_DUTs);

			for (int i = 0; i < n_DUTs; i++) {

				sscanf(argv[4+i], "%d", &DUT_IDs[i]);
			}

			if (argc > 4+n_DUTs) file_path = argv[4+n_DUTs];
		}
	}

	cout << "succesfully initiated variables:\noutput_file = " << output_file << "\nfile_path = " << file_path << "\nn_sensors = " << n_sensors << "\nn_DUTs = " << n_DUTs << endl;

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

		std::cout << "succesfully opened file: " << oss.str() << std::endl;

		oss.str(std::string());
		oss << file_path << "MPXNtuple_allPix_det_" << i+300 << ".root";

		frame_files[i] = new TFile(oss.str().c_str());
		frame_trees[i] = (TTree*) frame_files[i]->Get("MPXTree");

		std::cout << "succesfully opened file: " << oss.str() << std::endl;
	}

	for (int i = 0; i < n_DUTs; i++) {

		oss.str(std::string());
		oss << file_path << "Hits_BoxSD_" << DUT_IDs[i] << "_HitsCollection.root";

		true_hit_files[i+n_sensors] = new TFile(oss.str().c_str());
		true_hit_trees[i+n_sensors] = (TTree*)true_hit_files[i+n_sensors]->Get("AllPixHits");
		std::cout << "succesfully opened file: " << oss.str() << std::endl;

		oss.str(std::string());
		oss << file_path << "MPXNtuple_allPix_det_" << DUT_IDs[i] << ".root";

		frame_files[i+n_sensors] = new TFile(oss.str().c_str());
		frame_trees[i+n_sensors] = (TTree*) frame_files[i+n_sensors]->Get("MPXTree");
		std::cout << "succesfully opened file: " << oss.str() << std::endl;
	}

	int n_entries = true_hit_trees[0]->GetEntries();

	cout << "succesfully got number of entries: " << n_entries << endl;

	FrameStruct* root_frame = NULL;
	SimpleHits* root_hit = NULL;
	true_hit_trees[0]->SetBranchAddress("SimpleHits", &root_hit);		
	true_hit_trees[0]->GetEntry(0);
	int run_number = root_hit->run;

	cout << "succesfully got run number: " << run_number << endl;

	std::string detector_name = "EUTelescope";

	LCWriter* writer = LCFactory::getInstance()->createLCWriter();
	writer->open(output_file, LCIO::WRITE_NEW);

	cout << "succesfully opened instance of LCWriter\n";

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
		LCCollectionVec* sensor_hit_coll = new LCCollectionVec(LCIO::TRACKERHIT);
		//CellIDEncoder<TrackerHitImpl> true_hit_encoder_sensor(encoding_string, sensor_hit_coll);
		LCCollectionVec* sensor_frame_coll = new LCCollectionVec(LCIO::TRACKERDATA);
		CellIDEncoder<TrackerDataImpl> frame_encoder_sensor(encoding_string, sensor_frame_coll);

		//fil telescope collection
		for (int j = 0; j < n_sensors; j++) {

			TrackerDataImpl* frame_data = new TrackerDataImpl();

			true_hit_trees[j]->SetBranchAddress("SimpleHits", &root_hit);
			true_hit_trees[j]->GetEntry(i);
			frame_trees[j]->SetBranchAddress("FramesData", &root_frame);
			frame_trees[j]->GetEntry(i);

			frame_encoder_sensor.reset();
			frame_encoder_sensor["sensorID"] = j;
			frame_encoder_sensor["sparsePixelType"] = 2;
			frame_encoder_sensor.setCellID(frame_data);

			//loop over true hits
			double* hit_pos = new double[3];
			for (unsigned int k = 0; k < root_hit->pos.size(); k++) {

				TrackerHitImpl* true_hit_data = new TrackerHitImpl();

				hit_pos[0] = root_hit->pos[k].X();
				hit_pos[1] = root_hit->pos[k].Y();
				hit_pos[2] = root_hit->pos[k].Z();

				true_hit_data->setPosition(hit_pos);
				true_hit_data->setEDep(root_hit->edep[k]);

				sensor_hit_coll->addElement(true_hit_data);
			}

			//loop over hits from frame
			std::vector<float> charge_vec;
			for (auto it = root_frame->GetHitMap().begin(); it != root_frame->GetHitMap().end(); ++it) {

				//the key contains position of hit in the format key = y*width + x
				float x = (it->first)%root_frame->GetWidth();
				float y = (it->first)/root_frame->GetWidth();
				float TOT = it->second;

				charge_vec.push_back(x);
				charge_vec.push_back(y);
				charge_vec.push_back(TOT);
				charge_vec.push_back(0);
			}

			frame_data->setChargeValues(charge_vec);
			sensor_frame_coll->addElement(frame_data);

			charge_vec.clear();
			delete [] hit_pos;
		}

		event->addCollection(sensor_hit_coll, "true_hits_m26");
		event->addCollection(sensor_frame_coll, "zsdata_m26");

		//write DUT true hit data collection if there are DUTs present
		if (n_DUTs > 0) {

			LCCollectionVec* DUT_hit_coll = new LCCollectionVec(LCIO::TRACKERHIT);
			CellIDEncoder<TrackerHitImpl> true_hit_encoder_DUT(encoding_string, DUT_hit_coll);
			LCCollectionVec* DUT_frame_coll = new LCCollectionVec(LCIO::TRACKERDATA);
			CellIDEncoder<TrackerDataImpl> frame_encoder_DUT(encoding_string, DUT_frame_coll);

			for (int j = 0; j < n_DUTs; j++) {

				TrackerDataImpl* frame_data = new TrackerDataImpl();

				true_hit_trees[j+n_sensors]->SetBranchAddress("SimpleHits", &root_hit);
				true_hit_trees[j+n_sensors]->GetEntry(i);
				frame_trees[j+n_sensors]->SetBranchAddress("FramesData", &root_frame);
				frame_trees[j+n_sensors]->GetEntry(i);

				frame_encoder_DUT.reset();
				frame_encoder_DUT["sensorID"] = j + 6;
				frame_encoder_DUT["sparsePixelType"] = 2;
				frame_encoder_DUT.setCellID(frame_data);

				double* hit_pos = new double[3];
				for (unsigned int k = 0; k < root_hit->pos.size(); k++) {

					TrackerHitImpl* true_hit_data = new TrackerHitImpl();

					hit_pos[0] = root_hit->pos[k].X();
					hit_pos[1] = root_hit->pos[k].Y();
					hit_pos[2] = root_hit->pos[k].Z();

					true_hit_data->setPosition(hit_pos);
					true_hit_data->setEDep(root_hit->edep[k]);

					DUT_hit_coll->addElement(true_hit_data);
				}

				std::vector<float> charge_vec;
				for (auto it = root_frame->GetHitMap().begin(); it != root_frame->GetHitMap().end(); ++it) {

					float x = (it->first)%root_frame->GetWidth();
					float y = (it->first)/root_frame->GetWidth();
					float TOT = it->second;

					charge_vec.push_back(x);
					charge_vec.push_back(y);
					charge_vec.push_back(TOT);
					charge_vec.push_back(0);
				}

				frame_data->setChargeValues(charge_vec);
				DUT_frame_coll->addElement(frame_data);

				charge_vec.clear();
				delete [] hit_pos;
			}

			event->addCollection(DUT_hit_coll, "true_hits_DUT");
			event->addCollection(DUT_frame_coll, "zsdata_DUT");
		}


		writer->writeEvent(event);
		delete event;

		if (((i+1)%1000 == 0) || (i == 0)) cout << "succesfully recorded event " << i+1 << endl;
	}

	cout << "succesfully read in all data\n";

	writer->flush();
	writer->close();

	cout << "succesfully closed the writer\n";

	for (int i = 0; i < n_sensors+n_DUTs; i++) {

		true_hit_files[i]->Close();
		frame_files[i]->Close();

		cout << "succesfully closed files " << i << endl;

		delete true_hit_files[i];
		delete frame_files[i];

		cout << "succesfully deleted trees and files " << i << endl;
	}

	delete [] true_hit_trees;
	delete [] true_hit_files;

	cout << "succesfully deleted trees and files arrays\n";

	delete root_hit;
	delete root_frame;

	cout << "succesfully deleted root_hit and root_frame\n";

	return 0;
}

void print_info() {

	std::cout << "Converts allpix generated root files into LCIO format" << std::endl;
	std::cout << "Usage: compile and run with the arguments in the following order:" << std::endl;
	std::cout << "<number of sensors> <number of DUTs> <output file name> <list of DUT sensor IDs> <path to folder containing input root files>" << std::endl;
	std::cout << "if no DUTs are present, inputing 0 for the number of DUTs means that the list of DUTs is no longer required. The fourth argument will be interpreted as the path to folder containing input root files." << std::endl;

	return;
}
