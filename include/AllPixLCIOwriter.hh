/*
 * LCIO output for AllPix
 *
 * 
 * Author: Paul Schuetze (paul.schuetze@desy.de)
 */


#ifndef AllPixLCIOwriter_h
#define AllPixLCIOwriter_h 1

#include "AllPixWriter.hh"

#include "lcio.h"
#include "IMPL/LCTOOLS.h"
#include "IO/LCWriter.h"
#include "IMPL/LCGenericObjectImpl.h"
#include "IMPL/LCEventImpl.h"
#include "IMPL/LCRunHeaderImpl.h"
#include "IMPL/LCCollectionVec.h"
#include "IMPL/TrackerDataImpl.h"
#include "EVENT/TrackerData.h"
#include "EVENT/LCCollection.h"
#include "EVENT/LCEvent.h"
#include "EVENT/LCIO.h"
#include "IOIMPL/LCFactory.h"
#include "UTIL/CellIDEncoder.h"

#include "G4String.hh"

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>

using namespace std;
using namespace lcio;


class AllPixLCIOwriter : public AllPixWriter{
public:
  
  AllPixLCIOwriter();
  ~AllPixLCIOwriter();

  void Initialize(G4String folder, int runnr);

  void WriteEvent(int runnr, int eventID, map<int,vector<vector<vector<int>>>> data);
  
  void SetWriteEventID(G4bool flag);

  void Close();

private:

  void WriteRunHeader(int runnr);

  LCWriter* m_lcWriter;
  
  G4bool m_writeEventID=false;


};

#endif
