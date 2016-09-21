#include <fstream>
#include <iostream>
#include <sstream>
#include <cmath>
#include <vector>
#include <string>
#include <array>
#include <map>
#include <memory>

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
#include "allpix_dm.h"

using namespace lcio;

enum HitProperties {
	kHitInGlobalCoord  = 1L << 0,
	kFittedHit         = 1L << 1,
	kSimulatedHit      = 1L << 2
};

void print_info();

int main(int argc, char* argv[]) {

	std::string input_file_path, output_file;
	size_t n_sensors, n_DUTs;
	std::vector<int> DUT_IDs;

	if (argc < 4) {

		print_info();//if not enough arguments are given, print user info
		return 0;
	}
	else {

		input_file_path = argv[1];
		output_file = argv[2];
		n_sensors = static_cast<size_t>(atoi(argv[3]));//atoi converts the input string to an int

		if (argc > 4) {

			for (int i = 0; i < argc-4; i++) DUT_IDs.push_back(atoi(argv[4+i]));
		}
	}

	n_DUTs = DUT_IDs.size();

	auto totSensorSize = n_DUTs+n_sensors;
	
	auto true_hit_files = std::vector<TFile*>(totSensorSize, nullptr);
	auto true_hit_trees = std::vector<TTree*>(totSensorSize, nullptr);
	auto frame_files = std::vector<TFile*>(totSensorSize, nullptr);
	auto frame_trees = std::vector<TTree*>(totSensorSize, nullptr);

	std::ostringstream oss;

	for (size_t i = 0; i < n_sensors; i++) {

		oss.str(std::string());//resets oss to an empty string
		oss << input_file_path;
		if (input_file_path.back() == '/') oss << "Hits";
		oss << "_BoxSD_" << i+300 << "_HitsCollection.root";

		true_hit_files[i] = new TFile(oss.str().c_str());
		true_hit_trees[i] = dynamic_cast<TTree*>(true_hit_files[i]->Get("AllPixHits"));

		oss.str(std::string());
		oss << input_file_path;
		if (input_file_path.back() == '/') oss << "MPXNtuple";
		oss << "_allPix_det_" << i+300 << ".root";

		frame_files[i] = new TFile(oss.str().c_str());
		frame_trees[i] = dynamic_cast<TTree*>(frame_files[i]->Get("MPXTree"));
	}

	for (size_t i = 0; i < n_DUTs; i++) {

		oss.str(std::string());
		oss << input_file_path;
		if (input_file_path.back() == '/') oss << "Hits";
		oss << "_BoxSD_" << DUT_IDs[i] << "_HitsCollection.root";

		true_hit_files[i+n_sensors] = new TFile(oss.str().c_str());
		true_hit_trees[i+n_sensors] = dynamic_cast<TTree*>(true_hit_files[i+n_sensors]->Get("AllPixHits"));

		oss.str(std::string());
		oss << input_file_path;
		if (input_file_path.back() == '/') oss << "MPXNtuple";
		oss << "_allPix_det_" << DUT_IDs[i] << ".root";

		frame_files[i+n_sensors] = new TFile(oss.str().c_str());
		frame_trees[i+n_sensors] = dynamic_cast<TTree*>(frame_files[i+n_sensors]->Get("MPXTree"));
	}

	auto root_frames = std::vector<FrameStruct*>(totSensorSize, nullptr);
	auto root_hits = std::vector<SimpleHits*>(totSensorSize, nullptr);
	for (size_t i = 0; i < n_sensors+n_DUTs; i++) {
		frame_trees[i]->SetBranchAddress("FramesData", &root_frames[i]);
		true_hit_trees[i]->SetBranchAddress("SimpleHits", &root_hits[i]);
	}

	size_t n_entries = true_hit_trees[0]->GetEntries();
	true_hit_trees[0]->GetEntry(n_entries-1);
	size_t n_runs = root_hits[0]->run;

	std::vector<int> hit_incrs(totSensorSize, 0);

	size_t run_number = 0;
	std::string detector_name = "EUTelescope";

	LCWriter* writer = LCFactory::getInstance()->createLCWriter();
	writer->open(output_file, LCIO::WRITE_NEW);

	//set up run header as in TelescopeConverter.py
	auto run_header = std::unique_ptr<LCRunHeaderImpl>(new LCRunHeaderImpl); 
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
	writer->writeRunHeader(run_header.get());

	for (size_t i = 0; i < n_runs+1; i++) {

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
			for (size_t n = 0; n < n_sensors; n++) {

				LCGenericObjectImpl* setup_obj = new LCGenericObjectImpl(5,0,0);
				setup_obj->setIntVal(0, 102);
				setup_obj->setIntVal(1, 101);
				eudrb_setup->addElement(setup_obj);
			}

			event->addCollection(eudrb_setup, "eudrbSetup");

			//check if there is a DUT
			if (n_DUTs > 0) {

				LCCollectionVec* DUT_setup = new LCCollectionVec(LCIO::LCGENERICOBJECT);
				event->addCollection(DUT_setup, "DUTSetup");
			}
		}

		//ID encoder info
		std::string encoding_string_frame = "sensorID:7,sparsePixelType:5";
		std::string encoding_string_hit = "sensorID:7,properties:7";

		//Telescope data collection
		LCCollectionVec* hit_coll = new LCCollectionVec(LCIO::TRACKERHIT);
		CellIDEncoder<TrackerHitImpl> true_hit_encoder(encoding_string_hit, hit_coll);
		LCCollectionVec* sensor_frame_coll = new LCCollectionVec(LCIO::TRACKERDATA);
		CellIDEncoder<TrackerDataImpl> frame_encoder_sensor(encoding_string_frame, sensor_frame_coll);

		//fil telescope collection
		for (size_t j = 0; j < n_sensors; j++) {

			TrackerDataImpl* frame_data = new TrackerDataImpl();

			frame_trees[j]->GetEntry(i);

			frame_encoder_sensor.reset();
			frame_encoder_sensor["sensorID"] = j;
			frame_encoder_sensor["sparsePixelType"] = 2;
			frame_encoder_sensor.setCellID(frame_data);

			true_hit_encoder.reset();
			true_hit_encoder["sensorID"] = j;
			true_hit_encoder["properties"] = kHitInGlobalCoord + kSimulatedHit;

			//loop over true hits
			//there are multiple entries with the same run (aka event) number,
			//hence must loop over all the entries with the same run number as
			//the first for loop (looping over run numbers)
			while (true_hit_trees[j]->GetEntry(hit_incrs[j]) and root_hits[j]->run == (int)i) {

				auto hit_pos = std::array<double,3>();
				if (root_hits[j]->trackId.back() == 1) {//only a single track in this event

					for (size_t k = 0; k < root_hits[j]->pos.size(); k++) {

						if (k == root_hits[j]->pos.size()/2 - 1) {

							TrackerHitImpl* true_hit_data = new TrackerHitImpl();

							hit_pos[0] = root_hits[j]->pos[k].X();
							hit_pos[1] = root_hits[j]->pos[k].Y();
							hit_pos[2] = root_hits[j]->pos[k].Z();

							true_hit_data->setTime(root_hits[j]->event);
							true_hit_data->setType(root_hits[j]->parentId[k]);
							true_hit_data->setQuality(root_hits[j]->trackId[k]);
							true_hit_data->setPosition(hit_pos.data());
							true_hit_data->setEDep(root_hits[j]->edep[k]);
							true_hit_data->setEDepError(root_hits[j]->edepTotal);

							true_hit_encoder.setCellID(true_hit_data);
							hit_coll->addElement(true_hit_data);
						}
					}
				}
				else {//more than one track in this event

					vector<SimpleHits> hits;
					int hits_index = -1;
					for (size_t k = 0; k < root_hits[j]->pos.size(); k++) {

						if ((k == 0) || (root_hits[j]->trackId[k] != root_hits[j]->trackId[k-1]) || (root_hits[j]->interactions[k-1] == "Transportation")) {

							hits.push_back(SimpleHits());
							hits_index++;
						}

						hits[hits_index].edepTotal = root_hits[j]->edepTotal;
						hits[hits_index].pos.push_back(root_hits[j]->pos[k]);
						hits[hits_index].edep.push_back(root_hits[j]->edep[k]);
						hits[hits_index].trackId.push_back(root_hits[j]->trackId[k]);
						hits[hits_index].parentId.push_back(root_hits[j]->parentId[k]);
					}

					for (size_t k = 0; k < hits.size(); k++) {

						for (size_t l = 0; l < hits[k].pos.size(); l++) {

							if ((l == hits[k].pos.size()/2 - 1) || (hits[k].pos.size() == 1)) {

								TrackerHitImpl* true_hit_data = new TrackerHitImpl();

								hit_pos[0] = hits[k].pos[l].X();
								hit_pos[1] = hits[k].pos[l].Y();
								hit_pos[2] = hits[k].pos[l].Z();

								true_hit_data->setTime(root_hits[j]->event);
								true_hit_data->setType(hits[k].parentId[l]);
								true_hit_data->setQuality(hits[k].trackId[l]);
								true_hit_data->setPosition(hit_pos.data());
								true_hit_data->setEDep(hits[k].edep[l]);
								true_hit_data->setEDepError(hits[k].edepTotal);

								true_hit_encoder.setCellID(true_hit_data);
								hit_coll->addElement(true_hit_data);
							}
						}
					}
				}
						
				hit_incrs[j]++;
			}

			//loop over hits from frame
			std::vector<float> charge_vec;
			for (auto& hit: root_frames[j]->GetHitMap()) {

				//the key contains position of hit in the format key = y*width + x
				float x = (hit.first)%root_frames[j]->GetWidth();
				float y = (hit.first)/root_frames[j]->GetWidth();
				float TOT = hit.second;

				charge_vec.push_back(x);
				charge_vec.push_back(y);
				charge_vec.push_back(TOT);
				charge_vec.push_back(0);
			}

			frame_data->setChargeValues(charge_vec);
			sensor_frame_coll->addElement(frame_data);
		}

		//write DUT true hit data collection if there are DUTs present
		if (n_DUTs > 0) {

			LCCollectionVec* DUT_frame_coll = new LCCollectionVec(LCIO::TRACKERDATA);
			CellIDEncoder<TrackerDataImpl> frame_encoder_DUT(encoding_string_frame, DUT_frame_coll);

			for (size_t j = 0; j < n_DUTs; j++) {

				TrackerDataImpl* frame_data = new TrackerDataImpl();

				frame_trees[j+n_sensors]->GetEntry(i);

				std::string DUT_ID = to_string(DUT_IDs[j]);
				std::string DUT_ID_shortened = DUT_ID.erase(DUT_ID.find_last_of('0'), 1);

				frame_encoder_DUT.reset();
				frame_encoder_DUT["sensorID"] = atoi(DUT_ID_shortened.data());
				frame_encoder_DUT["sparsePixelType"] = 2;
				frame_encoder_DUT.setCellID(frame_data);

				true_hit_encoder.reset();
				true_hit_encoder["sensorID"] = atoi(DUT_ID_shortened.data());
				true_hit_encoder["properties"] = kHitInGlobalCoord + kSimulatedHit;

				while (true_hit_trees[j+n_sensors]->GetEntry(hit_incrs[j+n_sensors]) and root_hits[j+n_sensors]->run == (int)i) {

					auto hit_pos = std::array<double,3>();
					if (root_hits[j+n_sensors]->trackId.back() == 1) {

						for (size_t k = 0; k < root_hits[j+n_sensors]->pos.size(); k++) {

							if (k == root_hits[j+n_sensors]->pos.size()/2 - 1) {

								TrackerHitImpl* true_hit_data = new TrackerHitImpl();

								hit_pos[0] = root_hits[j+n_sensors]->pos[k].X();
								hit_pos[1] = root_hits[j+n_sensors]->pos[k].Y();
								hit_pos[2] = root_hits[j+n_sensors]->pos[k].Z();

								true_hit_data->setTime(root_hits[j+n_sensors]->event);
								true_hit_data->setType(root_hits[j+n_sensors]->parentId[k]);
								true_hit_data->setQuality(root_hits[j+n_sensors]->trackId[k]);
								true_hit_data->setPosition(hit_pos.data());
								true_hit_data->setEDep(root_hits[j+n_sensors]->edep[k]);
								true_hit_data->setEDepError(root_hits[j+n_sensors]->edepTotal);

								true_hit_encoder.setCellID(true_hit_data);
								hit_coll->addElement(true_hit_data);
							}
						}
					}
					else {

						vector<SimpleHits> hits;
						int hits_index = -1;
						for (size_t k = 0; k < root_hits[j+n_sensors]->pos.size(); k++) {

							if ((k == 0) || (root_hits[j+n_sensors]->trackId[k] != root_hits[j+n_sensors]->trackId[k-1]) || (root_hits[j+n_sensors]->interactions[k-1] == "Transportation")) {

								hits.push_back(SimpleHits());
								hits_index++;
							}

							hits[hits_index].edepTotal = root_hits[j+n_sensors]->edepTotal;
							hits[hits_index].pos.push_back(root_hits[j+n_sensors]->pos[k]);
							hits[hits_index].edep.push_back(root_hits[j+n_sensors]->edep[k]);
							hits[hits_index].trackId.push_back(root_hits[j+n_sensors]->trackId[k]);
							hits[hits_index].parentId.push_back(root_hits[j+n_sensors]->parentId[k]);
						}

						for (size_t k = 0; k < hits.size(); k++) {

							for (size_t l = 0; l < hits[k].pos.size(); l++) {

								if ((l == hits[k].pos.size()/2 - 1) || (hits[k].pos.size() == 1)) {

									TrackerHitImpl* true_hit_data = new TrackerHitImpl();

									hit_pos[0] = hits[k].pos[l].X();
									hit_pos[1] = hits[k].pos[l].Y();
									hit_pos[2] = hits[k].pos[l].Z();

									true_hit_data->setTime(root_hits[j+n_sensors]->event);
									true_hit_data->setType(hits[k].parentId[l]);
									true_hit_data->setQuality(hits[k].trackId[l]);
									true_hit_data->setPosition(hit_pos.data());
									true_hit_data->setEDep(hits[k].edep[l]);
									true_hit_data->setEDepError(hits[k].edepTotal);

									true_hit_encoder.setCellID(true_hit_data);
									hit_coll->addElement(true_hit_data);
								}
							}
						}
					}

					hit_incrs[j+n_sensors]++;
				}

				std::vector<float> charge_vec;
				for (auto& hit: root_frames[j+n_sensors]->GetHitMap()) {

					float x = (hit.first)%root_frames[j+n_sensors]->GetWidth();
					float y = (hit.first)/root_frames[j+n_sensors]->GetWidth();
					float TOT = hit.second;

					charge_vec.push_back(x);
					charge_vec.push_back(y);
					charge_vec.push_back(TOT);
					charge_vec.push_back(0);
				}

				frame_data->setChargeValues(charge_vec);
				DUT_frame_coll->addElement(frame_data);
			}

			event->addCollection(DUT_frame_coll, "zsdata_DUT");
		}

		event->addCollection(hit_coll, "true_hits");
		event->addCollection(sensor_frame_coll, "zsdata_m26");

		writer->writeEvent(event);
		delete event;

		if (((i+1)%1000 == 0) || (i == 0)) std::cout << "succesfully recorded event " << i+1 << endl;
	}

	std::cout << "succesfully recorded all events\n";

	writer->flush();
	writer->close();

	for (size_t i = 0; i < n_sensors+n_DUTs; i++) {

		true_hit_files[i]->Close();
		frame_files[i]->Close();

		delete true_hit_files[i];
		delete frame_files[i];
		delete root_hits[i];
		delete root_frames[i];
	}

	return 0;
}

void print_info() {

	std::cout << "Converts allpix generated root files into LCIO format" << std::endl;
	std::cout << "Usage: compile and run with the arguments in the following order:" << std::endl;
	std::cout << "<input file path+prefix> <output file name> <number of telescope planes> <list of DUT sensor IDs>" << std::endl << std::endl;
	std::cout << "Here, input file path+prefix refers to the path following the command /allpix/config/setOutputPrefixWithPath in the allpix macro. If this is a directory, it must end with a '/' or it will be misinterpreted." << std::endl << std::endl;
	std::cout << "The DUT sensor IDs that will be input into the lcio file in the CellID fields via the CellIDEncoder are the given DUT sensor IDs with the last occurence of the character '0' removed, so as to fit in the limited bit space of the CellID" << std::endl << std::endl;

	return;
}
