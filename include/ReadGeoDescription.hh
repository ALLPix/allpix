/**
 * Allpix
 * Author: John Idarraga <idarraga@cern.ch> , 2010
 */

#ifndef ReadGeoDescriptions_hh
#define ReadGeoDescriptions_hh

#include <string>
#include <vector>
#include <map>

#include "AllPixGeoDsc.hh"
#include "G4Types.hh"

using namespace std;

class TXMLNode;



#define __FIRST_DET_INDX 0

class ReadGeoDescription {

public:
	ReadGeoDescription(){};
	ReadGeoDescription(string);
	~ReadGeoDescription(){};

	static ReadGeoDescription * GetInstance();

	map<int, AllPixGeoDsc *> * GetDetectorsMap(){return &m_detsGeo;};
	void BuildListOfExpectedTags();
	void ParseContext(TXMLNode *);
	bool StringIsRelevant(string);
	vector<int> ProcessIdString(const char *);
	void ReplicateDetectors();
	G4int UseTheseDetectorsOnly(vector<G4int>);

private:
	string m_xmlfile;
	string m_currentNodeName;
	string m_currentAtt;

	map<int, AllPixGeoDsc *> m_detsGeo;
	//vector<int> m_detsGeoIndx;
	map<int, vector<int> > m_detsGeoIndx;
	int m_firstIndx;

	// map for replicated detectors
	//vector<vector<int> > m_replicationVector;

	map<string, double> m_unitsMap;

};



#endif

