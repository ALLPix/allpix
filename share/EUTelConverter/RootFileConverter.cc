#include <fstream>
#include <iostream>
#include <sstream>
#include <cmath>
#include <vector>
#include <string>
#include <array>
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

/*struct TrueHit {

	double edep_total;
	double sensor_ID;

	std::vector<double*> pos;
	std::vector<int> pdgId;
	std::vector<float> edep;

	~TrueHit() {

		for (unsigned int i = 0; i < pos.size(); i++) delete [] pos[i];
	}
};

struct Event {

	std::vector<TrueHit*> hits;

	~Event() {

		for (unsigned int i = 0; i < hits.size(); i++) delete hits[i];
	}
};*/

void print_info();

int main(int argc, char* argv[]) {

	std::string input_file_path, output_file;
	int n_sensors, n_DUTs;
	std::vector<int> DUT_IDs;

	if (argc < 4) {

		print_info();//if not enough arguments are given, print user info
		return 0;
	}
	else {

		input_file_path = argv[1];
		output_file = argv[2];
		sscanf(argv[3], "%d", &n_sensors);

		if (argc > 4) {

			for (int i = 0; i < 4-argc; i++) {

				int temp = 0;
				sscanf(argv[4+i], "%d", &temp);
				DUT_IDs.push_back(temp);
			}
		}
	}

	n_DUTs = DUT_IDs.size();

	TFile** true_hit_files = new TFile*[n_sensors+n_DUTs];
	TTree** true_hit_trees = new TTree*[n_sensors+n_DUTs];
	TFile** frame_files = new TFile*[n_sensors+n_DUTs];
	TTree** frame_trees = new TTree*[n_sensors+n_DUTs];

	std::ostringstream oss;

	for (int i = 0; i < n_sensors; i++) {

		oss.str(std::string());//resets oss to an empty string
		oss << input_file_path;
		if (input_file_path.back() == '/') oss << "Hits";
		oss << "_BoxSD_" << i+300 << "_HitsCollection.root";

		true_hit_files[i] = new TFile(oss.str().c_str());
		true_hit_trees[i] = (TTree*) true_hit_files[i]->Get("AllPixHits");

		oss.str(std::string());
		oss << input_file_path;
		if (input_file_path.back() == '/') oss << "MPXNtuple";
		oss << "_allPix_det_" << i+300 << ".root";

		frame_files[i] = new TFile(oss.str().c_str());
		frame_trees[i] = (TTree*) frame_files[i]->Get("MPXTree");
	}

	for (int i = 0; i < n_DUTs; i++) {

		oss.str(std::string());
		oss << input_file_path;
		if (input_file_path.back() == '/') oss << "Hits";
		oss << "_BoxSD_" << DUT_IDs[i] << "_HitsCollection.root";

		true_hit_files[i+n_sensors] = new TFile(oss.str().c_str());
		true_hit_trees[i+n_sensors] = (TTree*)true_hit_files[i+n_sensors]->Get("AllPixHits");

		oss.str(std::string());
		oss << input_file_path;
		if (input_file_path.back() == '/') oss << "MPXNtuple";
		oss << "_allPix_det_" << DUT_IDs[i] << ".root";

		frame_files[i+n_sensors] = new TFile(oss.str().c_str());
		frame_trees[i+n_sensors] = (TTree*) frame_files[i+n_sensors]->Get("MPXTree");
	}

	FrameStruct** root_frames = new FrameStruct*[n_sensors+n_DUTs];
	SimpleHits** root_hits = new SimpleHits*[n_sensors+n_DUTs];
	for (int i = 0; i < n_sensors+n_DUTs; i++) {

		root_frames[i] = NULL;
		root_hits[i] = NULL;

		frame_trees[i]->SetBranchAddress("FramesData", &root_frames[i]);
		true_hit_trees[i]->SetBranchAddress("SimpleHits", &root_hits[i]);

	}

	int n_entries = true_hit_trees[0]->GetEntries();
	true_hit_trees[0]->GetEntry(n_entries-1);
	int n_runs = root_hits[0]->run;

	/*Event** hit_events = new Event*[n_runs];
	for (int i = 0; i < n_runs; i++) hit_events[i] = new Event();*/

	std::vector<int> hit_incrs(n_sensors+n_DUTs);
	for (int i = 0; i < n_sensors+n_DUTs; i++) hit_incrs[i] = 0;

	int run_number = 0;
	std::string detector_name = "EUTelescope";

	LCWriter* writer = LCFactory::getInstance()->createLCWriter();
	writer->open(output_file, LCIO::WRITE_NEW);

	//set up run header as in TelescopeConverter.py
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

	for (int i = 0; i < n_runs; i++) {

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
		LCCollectionVec* sensor_frame_coll = new LCCollectionVec(LCIO::TRACKERDATA);
		//LCCollectionVec* sensor_vertex_coll = new LCCollectionVec(LCIO::TRACKERHIT);
		CellIDEncoder<TrackerDataImpl> frame_encoder_sensor(encoding_string, sensor_frame_coll);

		//fil telescope collection
		for (int j = 0; j < n_sensors; j++) {

			TrackerDataImpl* frame_data = new TrackerDataImpl();

			frame_trees[j]->GetEntry(i);

			frame_encoder_sensor.reset();
			frame_encoder_sensor["sensorID"] = j;
			frame_encoder_sensor["sparsePixelType"] = 2;
			frame_encoder_sensor.setCellID(frame_data);

			//loop over true hits
			//there are multiple entries with the same run (aka event) number,
			//hence must loop over all the entries with the same run number as
			//the first for loop (looping over run numbers)
			while (true_hit_trees[j]->GetEntry(hit_incrs[j]) and root_hits[j]->run == i) {

				/*hit_events[i]->hits.push_back(new TrueHit());
				hit_events[i]->hits.back()->sensor_ID = 300+j;
				hit_events[i]->hits.back()->edep_total = root_hits[j]->edepTotal;
				hit_events[i]->hits.back()->pdgId = root_hits[j]->pdgId;
				hit_events[i]->hits.back()->edep = root_hits[j]->edep;*/

				double* hit_pos = new double[3];
				for (unsigned int k = 0; k < root_hits[j]->pos.size(); k++) {

					//hit_events[i]->hits.back()->pos.push_back(hit_pos);

					if (k == root_hits[j]->pos.size()/2 - 1) {

						TrackerHitImpl* true_hit_data = new TrackerHitImpl();

						hit_pos[0] = root_hits[j]->pos[k].X();
						hit_pos[1] = root_hits[j]->pos[k].Y();
						hit_pos[2] = root_hits[j]->pos[k].Z();

						true_hit_data->setCellID0(300+j);
						true_hit_data->setQuality(root_hits[j]->event);
						true_hit_data->setType(root_hits[j]->pdgId[k]);
						true_hit_data->setPosition(hit_pos);
						true_hit_data->setEDep(root_hits[j]->edepTotal);

						sensor_hit_coll->addElement(true_hit_data);

					}
				}
				//if ((i%10 == 0) || (i == 0)) cout << root_hits[j]->event << endl;

				hit_incrs[j]++;

				delete [] hit_pos;
			}

			//loop over hits from frame
			std::vector<float> charge_vec;
			for (auto it = root_frames[j]->GetHitMap().begin(); it != root_frames[j]->GetHitMap().end(); ++it) {

				//the key contains position of hit in the format key = y*width + x
				float x = (it->first)%root_frames[j]->GetWidth();
				float y = (it->first)/root_frames[j]->GetWidth();
				float TOT = it->second;

				charge_vec.push_back(x);
				charge_vec.push_back(y);
				charge_vec.push_back(TOT);
				charge_vec.push_back(0);
			}

			frame_data->setChargeValues(charge_vec);
			sensor_frame_coll->addElement(frame_data);
			charge_vec.clear();

			//loop over the primary vertices stored in the frame file
			/*double* vertex_pos = new double[3];
			for (unsigned int k = 0; k < root_frames[j]->GetVertex_x().size(); k++) {

				TrackerHitImpl* frame_vertex = new TrackerHitImpl();

				vertex_pos[0] = root_frames[j]->GetVertex_x()[k];
				vertex_pos[1] = root_frames[j]->GetVertex_y()[k];
				vertex_pos[2] = root_frames[j]->GetVertex_z()[k];

				frame_vertex->setCellID0(300+j);
				frame_vertex->setPosition(vertex_pos);

				sensor_vertex_coll->addElement(frame_vertex);
			}

			delete [] vertex_pos;*/
		}

		/*TrackerHitImpl* test = new TrackerHitImpl();

		test->setCellID0(hit_events[i]->hits.back()->sensor_ID);
		test->setType(hit_events[i]->hits.back()->pdgId[hit_events[i]->hits.back()->pdgId.size()/2-1]);
		test->setPosition(hit_events[i]->hits.back()->pos[hit_events[i]->hits.back()->pos.size()/2-1]);
		test->setEDep(hit_events[i]->hits.back()->edep_total);

		sensor_hit_coll->addElement(test);*/

		event->addCollection(sensor_hit_coll, "true_hits_m26");
		event->addCollection(sensor_frame_coll, "zsdata_m26");
		//event->addCollection(sensor_vertex_coll, "MC_vertices_m26");

		//write DUT true hit data collection if there are DUTs present
		if (n_DUTs > 0) {

			LCCollectionVec* DUT_hit_coll = new LCCollectionVec(LCIO::TRACKERHIT);
			LCCollectionVec* DUT_frame_coll = new LCCollectionVec(LCIO::TRACKERDATA);
			//LCCollectionVec* DUT_vertex_coll = new LCCollectionVec(LCIO::TRACKERHIT);
			CellIDEncoder<TrackerDataImpl> frame_encoder_DUT(encoding_string, DUT_frame_coll);

			for (int j = 0; j < n_DUTs; j++) {

				TrackerDataImpl* frame_data = new TrackerDataImpl();

				frame_trees[j+n_sensors]->GetEntry(i);

				frame_encoder_DUT.reset();
				frame_encoder_DUT["sensorID"] = j + 6;
				frame_encoder_DUT["sparsePixelType"] = 2;
				frame_encoder_DUT.setCellID(frame_data);

				while (true_hit_trees[j+n_sensors]->GetEntry(hit_incrs[j+n_sensors]) and root_hits[j+n_sensors]->run == i) {

					/*hit_events[i]->hits.push_back(new TrueHit());
					hit_events[i]->hits.back()->sensor_ID = DUT_IDs[j];
					hit_events[i]->hits.back()->edep_total = root_hits[j]->edepTotal;
					hit_events[i]->hits.back()->pdgId = root_hits[j]->pdgId;
					hit_events[i]->hits.back()->edep = root_hits[j]->edep;*/

					double* hit_pos = new double[3];
					for (unsigned int k = 0; k < root_hits[j+n_sensors]->pos.size(); k++) {

						//hit_events[i]->hits.back()->pos.push_back(hit_pos);

						if (k == root_hits[j+n_sensors]->pos.size()/2 - 1) {

							TrackerHitImpl* true_hit_data = new TrackerHitImpl();

							hit_pos[0] = root_hits[j+n_sensors]->pos[k].X();
							hit_pos[1] = root_hits[j+n_sensors]->pos[k].Y();
							hit_pos[2] = root_hits[j+n_sensors]->pos[k].Z();

							true_hit_data->setCellID0(DUT_IDs[j]);
							true_hit_data->setQuality(root_hits[j+n_sensors]->event);
							true_hit_data->setType(root_hits[j+n_sensors]->pdgId[k]);
							true_hit_data->setPosition(hit_pos);
							true_hit_data->setEDep(root_hits[j+n_sensors]->edepTotal);

							DUT_hit_coll->addElement(true_hit_data);
						}
					}

					hit_incrs[j+n_sensors]++;

					delete [] hit_pos;
				}

				std::vector<float> charge_vec;
				for (auto it = root_frames[j]->GetHitMap().begin(); it != root_frames[j]->GetHitMap().end(); ++it) {

					float x = (it->first)%root_frames[j]->GetWidth();
					float y = (it->first)/root_frames[j]->GetWidth();
					float TOT = it->second;

					charge_vec.push_back(x);
					charge_vec.push_back(y);
					charge_vec.push_back(TOT);
					charge_vec.push_back(0);
				}

				frame_data->setChargeValues(charge_vec);
				DUT_frame_coll->addElement(frame_data);
				charge_vec.clear();


				/*double* vertex_pos = new double[3];
				for (unsigned int k = 0; k < root_frames[j]->GetVertex_x().size(); k++) {

					TrackerHitImpl* frame_vertex = new TrackerHitImpl();

					vertex_pos[0] = root_frames[j]->GetVertex_x()[k];
					vertex_pos[1] = root_frames[j]->GetVertex_y()[k];
					vertex_pos[2] = root_frames[j]->GetVertex_z()[k];

					frame_vertex->setCellID0(DUT_IDs[j]);
					frame_vertex->setPosition(vertex_pos);

					DUT_vertex_coll->addElement(frame_vertex);
				}

				delete [] vertex_pos;*/
			}

			event->addCollection(DUT_hit_coll, "true_hits_DUT");
			event->addCollection(DUT_frame_coll, "zsdata_DUT");
			//event->addCollection(DUT_vertex_coll, "MC_vertices_DUT");
		}


		writer->writeEvent(event);
		delete event;

		if (((i+1)%1000 == 0) || (i == 0)) std::cout << "succesfully recorded event " << i+1 << endl;
	}

	std::cout << "succesfully recorded all events\n";

	writer->flush();
	writer->close();

	for (int i = 0; i < n_sensors+n_DUTs; i++) {

		true_hit_files[i]->Close();
		frame_files[i]->Close();

		delete true_hit_files[i];
		delete frame_files[i];
		delete root_hits[i];
		delete root_frames[i];
	}

	/*for (int i = 0; i < n_runs; i++) {

		delete hit_events[i];
	}

	delete [] hit_events;*/

	delete [] true_hit_trees;
	delete [] true_hit_files;
	delete [] root_hits;
	delete [] root_frames;

	return 0;
}

void print_info() {

	std::cout << "Converts allpix generated root files into LCIO format" << std::endl;
	std::cout << "Usage: compile and run with the arguments in the following order:" << std::endl;
	std::cout << "<input file path+prefix> <output file name> <number of telescope planes> <list of DUT sensor IDs>" << std::endl;
	std::cout << "Here, input file path+prefix refers to the path following the command /allpix/config/setOutputPrefixWithPath in the allpix macro. If this is a directory, it must end with a '/' or it will be misinterpreted." << std::endl;

	return;
}
