/**
 * Allpix
 * Author: John Idarraga <idarraga@cern.ch> , 2010
 */

#include "ReadGeoDescription.hh"
#include "ReadGeoDescription_defs.hh"

#include <Riostream.h>
#include <TDOMParser.h>
#include <TXMLNode.h>
#include <TXMLAttr.h>
#include <TList.h>
#include <TString.h>

#include <set>
using namespace std;


// to be used as extern later
ReadGeoDescription * g_GeoDsc;

ReadGeoDescription::ReadGeoDescription(string xmlFile){

	if(g_GeoDsc){
		std::cout << "Can't instantiate this object twice ... ReadGeoDescription::ReadGeoDescription " << std::endl;
		exit(1);
	}

	// Initialize my vector
	//m_detsGeo = new map<int, AllPixGeoDsc *>;

	// Units map
	// Initializing units from strings to
	// map units from the xml file.
	m_unitsMap["nm"] = CLHEP::nm;
	m_unitsMap["um"] = CLHEP::um;
	m_unitsMap["mm"] = CLHEP::mm;
	m_unitsMap["cm"] = CLHEP::cm;
	m_unitsMap["m"] = CLHEP::m;

	// list of expected tags
	TDOMParser *domParser = new TDOMParser();
	m_firstIndx = -1;

	domParser->SetValidate(false); // do not validate with DTD for now
	domParser->ParseFile(xmlFile.c_str());

	TXMLNode *node = domParser->GetXMLDocument()->GetRootNode();

	// parse
	ParseContext(node);
	// apply requested replicas
	ReplicateDetectors();

	cout << "Summary: read "<< (int)m_detsGeo.size() << " detectors from xml database " << endl;
	map<int, AllPixGeoDsc *>::iterator itr = m_detsGeo.begin();
	for( ; itr != m_detsGeo.end() ; itr++){
		(*itr).second->Dump();
	}

	// keep this pointer
	g_GeoDsc = this;

}

/**
 * Keep only the detectors that will be used
 */
G4int ReadGeoDescription::UseTheseDetectorsOnly(vector<G4int> useDetectors){

	// Check for replicated sensors
	set<G4int> checkIds;

	// First verify if all detectors are in the db
	vector<G4int>::iterator itr = useDetectors.begin();
	for( ; itr != useDetectors.end() ; itr++ ){

		if(m_detsGeo.find(*itr) == m_detsGeo.end()){ // not found
			G4cout << "[OOPS] detector with Id " << *itr << " requested in the macro but not found in the db." << G4endl;
			G4cout << "       Allpix can't recover ... so long, and thanks for all the fish ;)" << G4endl;
			exit(1);
		}

		if( checkIds.find(*itr) == checkIds.end()) { // not in the list all good
			checkIds.insert( *itr );
		} else {
			G4cout << "[OOPS] detector with Id " << *itr << " duplicated in the macro ! please check your macro." << G4endl;
			G4cout << "       Allpix can't recover ... so long, and thanks for all the fish ;)" << G4endl;
			exit(1);
		}

	}

	// now erase what's not needed
	G4int nErased = 0;
	map<int, AllPixGeoDsc *>::iterator detItr = m_detsGeo.begin();
	bool found = false;
	vector<G4int> scheduledErase;

	for( ; detItr != m_detsGeo.end() ; detItr++){

		G4int aDet = (*detItr).first;
		found = false;

		for( itr = useDetectors.begin() ; itr != useDetectors.end() ; itr++ ){


			if(aDet == *itr){
				found = true;
				break;
			}
		}


		if(!found){
			G4cout << " ----> schedule for erasing " << aDet << G4endl;
			// don't erase right away, schedule and delete later, otherwise the
			// vector over which I am looping shrinks.
			scheduledErase.push_back(aDet);
			//m_detsGeo.erase(aDet);
			nErased++;
		}

	}

	// and finally erase
	vector<G4int>::iterator eraseItr = scheduledErase.begin();
	for( ; eraseItr != scheduledErase.end() ; eraseItr++ ) {
		m_detsGeo.erase(*eraseItr);
	}

	return nErased;
}

void ReadGeoDescription::BuildListOfExpectedTags(){


}

ReadGeoDescription * ReadGeoDescription::GetInstance(){

	if(!g_GeoDsc){
		std::cout << "This object has to be intantiated by DetectorConstruction first ... ReadGeoDescription::GetInstance()" << std::endl;
		exit(1);
	}

	return g_GeoDsc;
}

void ReadGeoDescription::ParseContext(TXMLNode *node)
{

	string tempContent;
	string tempAtt1;

	for ( ; node ; node = node->GetNextNode()) {

		if (node->GetNodeType() == TXMLNode::kXMLElementNode) { // Element Node

			m_currentNodeName = string(node->GetNodeName());

			//cout << m_currentNodeName << endl;

			/*
			if(m_currentNodeName == __pixeldet_node_S) {
				cout << "Creating configuration for device with id : ";
				//m_detsGeoIndx = 0;
			}
			 */

			// Catch properties first if any
			// This should be the right attribute each time
			// guaranteed by the dtd file
			if (node->HasAttributes()) {

				TList * attrList = node->GetAttributes();
				TIter next(attrList);
				TXMLAttr *attr;

				while ((attr =(TXMLAttr*)next())) {

					// verifying attributes names
					tempContent = string(attr->GetName());       // att name

					if(tempContent == __pixeldet_node_ATT_id_S && m_currentNodeName == __pixeldet_node_S) { // check if this is the right attribute "id"

						tempAtt1 = string(attr->GetValue());      // fetch the value

						if(StringIsRelevant(tempAtt1)){

							// Analyse the id string first
							// list of detector to create with this info
							vector<int> indexes = ProcessIdString(tempAtt1.c_str());
							// save first index
							m_firstIndx = indexes[__FIRST_DET_INDX];
							// get rid of it
							indexes.erase(indexes.begin());

							// put the rest in the map
							m_detsGeoIndx[m_firstIndx] = indexes;

							// create the first
							m_detsGeo[m_firstIndx] = new AllPixGeoDsc;
							m_detsGeo[m_firstIndx]->SetID(atoi(tempAtt1.c_str()));

							cout << "Creating configuration for device with id : " << m_firstIndx << endl; // display the id
							//cout << m_detsGeoIndx[m_firstIndx].size() << " copies requested" << endl;

						}
					}else if(tempContent == __pixeldet_global_ATT_units_S){
						tempAtt1 = string(attr->GetValue());      // fetch the value
						if(StringIsRelevant(tempAtt1)){
							// use this units when processing the contents
							m_currentAtt = tempAtt1;
						}
					}

				}

			}


		}
		if (node->GetNodeType() == TXMLNode::kXMLTextNode) { // Text node

			// get the other nodes into the GeoDsc for the current detector
			tempContent = string(node->GetContent());

			//if(m_detsGeoIndx[__FIRST_DET_INDX] > -1 && StringIsRelevant(tempContent)){
			if(StringIsRelevant(tempContent)){

				if(m_currentNodeName == __npix_x_S) {

					int val = atoi(tempContent.c_str());
					m_detsGeo[m_firstIndx]->SetNPixelsX(val);

				}else if(m_currentNodeName == __npix_y_S){

					int val = atoi(tempContent.c_str());
					m_detsGeo[m_firstIndx]->SetNPixelsY(val);

				}else if(m_currentNodeName == __npix_z_S){

					int val = atoi(tempContent.c_str());
					m_detsGeo[m_firstIndx]->SetNPixelsZ(val);

				}else if(m_currentNodeName == __chip_hx_S){

					float val = atof(tempContent.c_str());

					m_detsGeo[m_firstIndx]->SetChipHX(val*m_unitsMap[m_currentAtt]);

				}else if(m_currentNodeName == __chip_hy_S){

					float val = atof(tempContent.c_str());
					m_detsGeo[m_firstIndx]->SetChipHY(val*m_unitsMap[m_currentAtt]);

				}else if(m_currentNodeName == __chip_hz_S){

					float val = atof(tempContent.c_str());
					m_detsGeo[m_firstIndx]->SetChipHZ(val*m_unitsMap[m_currentAtt]);

				}else if(m_currentNodeName == __chip_posx_S){

					float val = atof(tempContent.c_str());
					m_detsGeo[m_firstIndx]->SetChipPosX(val*m_unitsMap[m_currentAtt]);

				}else if(m_currentNodeName == __chip_posy_S){

					float val = atof(tempContent.c_str());
					m_detsGeo[m_firstIndx]->SetChipPosY(val*m_unitsMap[m_currentAtt]);

				}else if(m_currentNodeName == __chip_posz_S){

					float val = atof(tempContent.c_str());
					m_detsGeo[m_firstIndx]->SetChipPosZ(val*m_unitsMap[m_currentAtt]);

				}else if(m_currentNodeName == __chip_offsetx_S){

					float val = atof(tempContent.c_str());
					m_detsGeo[m_firstIndx]->SetChipOffsetX(val*m_unitsMap[m_currentAtt]);

				}else if(m_currentNodeName == __chip_offsety_S){

					float val = atof(tempContent.c_str());
					m_detsGeo[m_firstIndx]->SetChipOffsetY(val*m_unitsMap[m_currentAtt]);

				}else if(m_currentNodeName == __chip_offsetz_S){

					float val = atof(tempContent.c_str());
					m_detsGeo[m_firstIndx]->SetChipOffsetZ(val*m_unitsMap[m_currentAtt]);

				}else if(m_currentNodeName == __pixsize_x_S){

					float val = atof(tempContent.c_str());
					m_detsGeo[m_firstIndx]->SetPixSizeX(val*m_unitsMap[m_currentAtt]);

				}else if(m_currentNodeName == __pixsize_y_S){

					float val = atof(tempContent.c_str());
					m_detsGeo[m_firstIndx]->SetPixSizeY(val*m_unitsMap[m_currentAtt]);

				}else if(m_currentNodeName == __pixsize_z_S){

					float val = atof(tempContent.c_str());
					m_detsGeo[m_firstIndx]->SetPixSizeZ(val*m_unitsMap[m_currentAtt]);

				}else if(m_currentNodeName == __sensor_hx_S){

					float val = atof(tempContent.c_str());
					m_detsGeo[m_firstIndx]->SetSensorHX(val*m_unitsMap[m_currentAtt]);

				}else if(m_currentNodeName == __sensor_hy_S){

					float val = atof(tempContent.c_str());
					m_detsGeo[m_firstIndx]->SetSensorHY(val*m_unitsMap[m_currentAtt]);

				}else if(m_currentNodeName == __sensor_hz_S){

					float val = atof(tempContent.c_str());
					m_detsGeo[m_firstIndx]->SetSensorHZ(val*m_unitsMap[m_currentAtt]);

                }else if(m_currentNodeName == __coverlayer_hz_S){

                    float val = atof(tempContent.c_str());
                    m_detsGeo[m_firstIndx]->SetCoverlayerHZ(val*m_unitsMap[m_currentAtt]);

                }else if(m_currentNodeName == __coverlayer_mat_S){

                    G4String valS(tempContent);
                    m_detsGeo[m_firstIndx]->SetCoverlayerMat(valS);

                }else if(m_currentNodeName == __sensor_posx_S){

					float val = atof(tempContent.c_str());
					m_detsGeo[m_firstIndx]->SetSensorPosX(val*m_unitsMap[m_currentAtt]);

				}else if(m_currentNodeName == __sensor_posy_S){

					float val = atof(tempContent.c_str());
					m_detsGeo[m_firstIndx]->SetSensorPosY(val*m_unitsMap[m_currentAtt]);

				}else if(m_currentNodeName == __sensor_posz_S){

					float val = atof(tempContent.c_str());
					m_detsGeo[m_firstIndx]->SetSensorPosZ(val*m_unitsMap[m_currentAtt]);

				}else if(m_currentNodeName == __pcb_hx_S){

					float val = atof(tempContent.c_str());
					m_detsGeo[m_firstIndx]->SetPCBHX(val*m_unitsMap[m_currentAtt]);

				}else if(m_currentNodeName == __pcb_hy_S){

					float val = atof(tempContent.c_str());
					m_detsGeo[m_firstIndx]->SetPCBHY(val*m_unitsMap[m_currentAtt]);

				}else if(m_currentNodeName == __pcb_hz_S){

					float val = atof(tempContent.c_str());
					m_detsGeo[m_firstIndx]->SetPCBHZ(val*m_unitsMap[m_currentAtt]);

				}else if(m_currentNodeName == __sensor_gr_excess_htop_S){

					float val = atof(tempContent.c_str());
					m_detsGeo[m_firstIndx]->SetSensorExcessHTop(val*m_unitsMap[m_currentAtt]);

				}else if(m_currentNodeName == __sensor_gr_excess_hbottom_S){

					float val = atof(tempContent.c_str());
					m_detsGeo[m_firstIndx]->SetSensorExcessHBottom(val*m_unitsMap[m_currentAtt]);

				}else if(m_currentNodeName == __sensor_gr_excess_hright_S){

					float val = atof(tempContent.c_str());
					m_detsGeo[m_firstIndx]->SetSensorExcessHRight(val*m_unitsMap[m_currentAtt]);

				}else if(m_currentNodeName == __sensor_gr_excess_hleft_S){

					float val = atof(tempContent.c_str());
					m_detsGeo[m_firstIndx]->SetSensorExcessHLeft(val*m_unitsMap[m_currentAtt]);

				}else if(m_currentNodeName == __digitizer_S){

					G4String valS(tempContent.c_str());
					m_detsGeo[m_firstIndx]->SetSensorDigitizer(valS);

				}
	
				else if(m_currentNodeName == __sensor_Resistivity){

					float val = atof(tempContent.c_str());
					m_detsGeo[m_firstIndx]->SetResistivity(val);

				}				
				
				
				else if(m_currentNodeName == __MIP_Tot_S){

					float val = atoi(tempContent.c_str());
					m_detsGeo[m_firstIndx]->SetMIPTot(val);

				}
				
				else if(m_currentNodeName == __MIP_Charge_S){

					float val = atof(tempContent.c_str());
					m_detsGeo[m_firstIndx]->SetMIPCharge(val);

				}
								
				else if(m_currentNodeName == __Counter_Depth_S){

					float val = atoi(tempContent.c_str());
					m_detsGeo[m_firstIndx]->SetCounterDepth(val);

				}				

				else if(m_currentNodeName == __Clock_Unit_S){

					float val = atof(tempContent.c_str());
					m_detsGeo[m_firstIndx]->SetClockUnit(val);

				}					
				
				else if(m_currentNodeName == __Chip_Noise_S){

					float val = atof(tempContent.c_str());
					m_detsGeo[m_firstIndx]->SetChipNoise(val);

				}

				else if(m_currentNodeName == __Chip_Threshold_S){

					float val = atof(tempContent.c_str());
					m_detsGeo[m_firstIndx]->SetThreshold(val);

				}
				
				else if(m_currentNodeName == __Cross_Talk_S){

					float val = atof(tempContent.c_str());
					m_detsGeo[m_firstIndx]->SetCrossTalk(val);

				}
								

				else if(m_currentNodeName == __Saturation_Energy_S){

					float val = atof(tempContent.c_str());
					m_detsGeo[m_firstIndx]->SetSaturationEnergy(val);

				}
								
				else if(m_currentNodeName == __Bump_Radius_S){

					float val = atof(tempContent.c_str());
					m_detsGeo[m_firstIndx]->SetBumpRadius(val*m_unitsMap[m_currentAtt]);

				}
				
				else if(m_currentNodeName == __Bump_Height_S){

					float val = atof(tempContent.c_str());
					m_detsGeo[m_firstIndx]->SetBumpHeight(val*m_unitsMap[m_currentAtt]);

				}
				
				else if(m_currentNodeName == __Bump_OffsetX_S){

					float val = atof(tempContent.c_str());
					m_detsGeo[m_firstIndx]->SetBumpOffsetX(val*m_unitsMap[m_currentAtt]);

				}

				else if(m_currentNodeName == __Bump_OffsetY_S){

					float val = atof(tempContent.c_str());
					m_detsGeo[m_firstIndx]->SetBumpOffsetY(val*m_unitsMap[m_currentAtt]);

				}

				else if(m_currentNodeName == __Bump_Dr_S){

					float val = atof(tempContent.c_str());
					m_detsGeo[m_firstIndx]->SetBumpDr(val*m_unitsMap[m_currentAtt]);

				}
				
			}
			
			/*
			if(StringIsRelevant(tempContent))
				cout << "+" << tempContent << "-" << endl;
			 */
		}
		if (node->GetNodeType() == TXMLNode::kXMLCommentNode) { //Comment node
			//cout << "Comment: " << node->GetContent();
		}

		ParseContext(node->GetChildren());
	}


}

vector<int> ReadGeoDescription::ProcessIdString(const char * st){

	TString idString(st);
	vector<int> res;

	if(!idString.Contains(',',TString::kExact)){

		res.push_back(atoi(st));
		return res;
	}

	TString tempS;
	for(int i = 0 ; i < idString.Sizeof() ; i++){

		if(idString[i] == ',' && tempS.Sizeof() != 0){
			res.push_back( atoi(tempS.Data()) );
			tempS.Clear();
		}else{
			tempS.Append(idString[i]);
		}

	}
	// last number
	if(tempS.Sizeof() != 0) res.push_back( atoi(tempS.Data()) );

	/*
	for(int i = 0 ; i < (int)res.size() ; i++){
		cout << res[i] << endl;
	}
	 */
	return res;
}

/**
 * Replicate detectors if needed
 * i.e. create new instances an copy from the first object in the vector
 */

void ReadGeoDescription::ReplicateDetectors(){

	map<int, vector<int> >::iterator itr = m_detsGeoIndx.begin(); // step in the second instance
	vector<int> replicas;
	vector<int>::iterator repItr;

	for( ; itr != m_detsGeoIndx.end() ; itr++) {

		replicas = (*itr).second;

		for(repItr = replicas.begin() ; repItr != replicas.end() ; repItr++){
			m_detsGeo[*repItr] = new AllPixGeoDsc;
			*(m_detsGeo[*repItr]) = *(m_detsGeo[(*itr).first]); // ! a copy of object contents !
			cout << "  --> Creating configuration replica with Id : " << *repItr << " (copy of " << (*itr).first << ")" << endl;
			m_detsGeo[*repItr]->SetID(*repItr);
		}

	}

}

bool ReadGeoDescription::StringIsRelevant(string s){

	TString t(s);
	if(t.Contains('\n') || t.Contains('\t'))
		return false;

	return true;
}
