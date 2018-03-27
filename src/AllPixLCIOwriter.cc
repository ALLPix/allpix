/*
 * LCIO output for AllPix
 *
 * 
 * Author: Paul Schuetze (paul.schuetze@desy.de)
 */


#include "AllPixLCIOwriter.hh"

using namespace std;
using namespace lcio;

AllPixLCIOwriter::AllPixLCIOwriter(){
}


AllPixLCIOwriter::~AllPixLCIOwriter(){
}

void AllPixLCIOwriter::Initialize(G4String folder, int runnr){

  LCFactory * lcf = LCFactory::getInstance();
  m_lcWriter = lcf->createLCWriter();

  char filename[50];
  sprintf(filename, "run%06d-converter.slcio", runnr);
  
  folder.append(G4String(filename));

  m_lcWriter->open(folder.data(), LCIO::WRITE_NEW);

  WriteRunHeader(runnr);
  
}

void AllPixLCIOwriter::WriteRunHeader(int runnr){

  LCRunHeaderImpl * runHdr = new LCRunHeaderImpl;
  runHdr->setRunNumber(runnr);
  runHdr->setDetectorName("EUTelescope");

  m_lcWriter->writeRunHeader(runHdr);
  
  delete runHdr;
  
}

void AllPixLCIOwriter::SetWriteEventID(G4bool flag) { 

   this->m_writeEventID = flag;

}

void AllPixLCIOwriter::WriteEvent(int runnr, int eventID, map<int,vector<vector<vector<int>>>> data){

  LCEventImpl * event = new LCEventImpl();
  event->setRunNumber(runnr);
  event->setEventNumber(eventID);
  event->setDetectorName("AllPix");
  event->parameters().setValue("EventType",2);

  
  // Let's have a look at the data

  // Prerequisites
  
  int nDetectors = data.size();
  int telescopePlanes = 0;
  int DUTPlanes = 0;

  map<int,vector<vector<vector<int> > > >::iterator itr;
  for (itr = data.begin(); itr != data.end(); itr++){ // Loop over the planes
    
    int detId = itr->first;
    if(detId >= 300 && detId < 350) telescopePlanes++;
    else DUTPlanes++;

  }

  // Setup vectors for LCIO files
  
  if(eventID==0){

    cout << "Got " << DUTPlanes << " DUTs and " << telescopePlanes << " telescope planes." << endl;
  
    LCCollectionVec* setupVec = new LCCollectionVec( LCIO::LCGENERICOBJECT );

    setupVec->parameters().setValue("DataDescription","type:i,mode:i,spare1:i,spare2:i,spare3:i");
    setupVec->parameters().setValue("TypeName","Setup Description");

    for(int i=0; i<telescopePlanes; i++){

      LCGenericObjectImpl* setupObject = new LCGenericObjectImpl(5, 0, 0);
      setupObject->setIntVal(0, 102);
      setupObject->setIntVal(1, 101);

      setupVec->addElement(setupObject);
    }
    event->addCollection(setupVec, "TelescopeSetup");

    if(DUTPlanes){
      LCCollectionVec* setupVecDUT = new LCCollectionVec( LCIO::LCGENERICOBJECT );
      event->addCollection(setupVecDUT, "DUTSetup");
    }

  } // Setup


  // Now the actual data
  
  char encodingString[50];
  sprintf(encodingString, "sensorID:7,sparsePixelType:5");

  int planeNr = 0;
  int dutNr = telescopePlanes;

  LCCollectionVec* telescopeDataCollection = new LCCollectionVec( LCIO::TRACKERDATA );
  CellIDEncoder<TrackerDataImpl> telescopeEncoder(encodingString, telescopeDataCollection);

  bool inTelescope=false;
  
  for (itr = data.begin(); itr != data.end(); itr++){ // Loop over the planes
    
    int detId = itr->first;

    if(detId >= 300 && detId < 350) inTelescope=true;
    else inTelescope=false;
    
    LCCollectionVec* dataCollection = new LCCollectionVec( LCIO::TRACKERDATA );

    CellIDEncoder<TrackerDataImpl> encoder(encodingString, dataCollection);

    TrackerDataImpl* planeData = new TrackerDataImpl();


    if(inTelescope){
      telescopeEncoder["sensorID"] = planeNr;
      telescopeEncoder["sparsePixelType"] = 2;  //type 1 works for EUTelescope reconstruction
      
      telescopeEncoder.setCellID(planeData);
    }else{
      encoder["sensorID"] = dutNr;
      encoder["sparsePixelType"] = 2;      //type 1 works for EUTelescope reconstruction

      encoder.setCellID(planeData);
    }
    
    // maps to store multiple hits per pixel
    map<int,int> map_xyTOT;
    map<int,vector<int> > map_xyEvent;

    // loop on events
    vector<vector<vector<int> > > v_events = itr->second;
    for ( unsigned int i=0; i<v_events.size(); i++ )
    {
      // loop on digits for this event
      vector<vector<int> > v_digits = v_events[i];
      for ( unsigned int j=0; j<v_digits.size(); j++ )
      {
        vector<int> v_data = v_digits[j];
        int x    = v_data[0];
        int y    = v_data[1];
        int tot  = v_data[2];
        int evID = v_data[3];
        
        map_xyTOT[x*1000+y] += tot;
        map_xyEvent[x*1000+y].push_back(evID);    
      }
    }

    vector<float> pixelData;
    
    // loop on maps in order to sum up digits in same pixels
    map<int,int>::iterator map_itr;
    for (map_itr = map_xyTOT.begin(); map_itr != map_xyTOT.end(); map_itr++)
    {
      int index = map_itr->first;
      div_t division = div(index,1000);
      int x = division.quot; // pixel x coordinate
      int y = division.rem;  // pixel y coordinate
      int tot = map_itr->second;

      
      if(tot>=0){

	pixelData.push_back(x);
	pixelData.push_back(y);
	pixelData.push_back(tot);
	
	if(m_writeEventID){
	pixelData.push_back(0);
	}
        
      }
    }

    planeData->setChargeValues(pixelData);

    char collectionName[50];

    if(detId == 900){
      sprintf(collectionName, "CMSPixelDUT");
    }else if(detId == 901){
      sprintf(collectionName, "CMSPixelREF");
    }else{
      //sprintf(collectionName, "Det%d",detId);
      sprintf(collectionName, "zsdata_apix"); //EUTelescope reco
    }

    if(inTelescope){
      telescopeDataCollection->addElement(planeData);
      planeNr++;
    }else{
      dataCollection->addElement(planeData);
      event->addCollection(dataCollection, collectionName);
      dutNr++;
    }

    
  }

  if(telescopePlanes){
    event->addCollection(telescopeDataCollection, "zsdata_m26");
  }


  // Write out event to LCIO writer
  
  m_lcWriter->writeEvent(event);
  
  m_lcWriter->flush();
  cout << "Wrote Event to LCIO" << endl;  
    
}

void AllPixLCIOwriter::Close(){

  m_lcWriter->close();

}
