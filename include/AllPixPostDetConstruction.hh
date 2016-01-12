/**
 *  Author John Idarraga <idarraga@cern.ch>
 */

#ifndef AllPixPostDetConstruction_h
#define AllPixPostDetConstruction_h 1

//Example of post detector construction user class
//#include "TG4RootDetectorConstruction.h"

#include <TROOT.h>

#include "CLHEP/Units/SystemOfUnits.h"
using namespace CLHEP;

class TPolyLine3D;
class TObjArray;

class AllPixPostDetConstruction { //: public TVirtualUserPostDetConstruction {
  
private:
  TObjArray            *fTracks;  // Array of tracks
  TPolyLine3D          *fCurrent; // Current track
  
  AllPixPostDetConstruction();
  static AllPixPostDetConstruction *fgInstance; // Self pointer
public:
  virtual ~AllPixPostDetConstruction();
  
  static AllPixPostDetConstruction *GetInstance();
  
  void                  NewTrack(Double_t x, Double_t y, Double_t z);
  void                  AddPoint(Double_t x, Double_t y, Double_t z);
  void                  WriteTracks(const char *filename);
  
  virtual void          Initialize();
};

#endif
