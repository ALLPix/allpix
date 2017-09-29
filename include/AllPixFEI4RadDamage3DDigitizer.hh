/**
 * Author:
 *    Veronica Wallangen veronica.wallangen@cern.ch
 *
 *  allpix Authors:
 *   John Idarraga <idarraga@cern.ch>
 *   Mathieu Benoit <Mathieu.Benoit@cern.ch>
 */

#ifndef AllPixFEI4RadDamage3DDigitizer_h
#define AllPixFEI4RadDamage3DDigitizer_h 1

// allpix Interface
#include "AllPixDigitizerInterface.hh"
// digits for this digitizer
#include "AllPixFEI4RadDamage3DDigit.hh"
#include "G4PrimaryVertex.hh"

#include <map>
#include <vector>

#include "TH2F.h"
#include "TH3F.h"

using namespace std;

/**
 *  Digitizer AllPixFEI4RadDamage3D implementation
 */
class AllPixFEI4RadDamage3DDigitizer : public  AllPixDigitizerInterface {

public:

  AllPixFEI4RadDamage3DDigitizer(G4String, G4String, G4String);
  virtual ~AllPixFEI4RadDamage3DDigitizer();

  void SetPrimaryVertex(G4PrimaryVertex * pv) {m_primaryVertex = pv;};
  void Digitize ();
  void SetDetectorDigitInputs(G4double){};

private:

  // digitInput typedef is defined in AllPixDigitizerInterface.hh
  digitInput m_digitIn;

  AllPixFEI4RadDamage3DDigitsCollection * m_digitsCollection;
  vector<G4String> m_hitsColName;
  G4PrimaryVertex * m_primaryVertex; // information from EventAction

  TH2F *eFieldMap;
  TH2F *timeMap_e;
  TH2F *timeMap_h;
  TH3F *xPositionMap_e;
  TH3F *xPositionMap_h;
  TH3F *yPositionMap_e;
  TH3F *yPositionMap_h;
  TH2F *ramoPotentialMap;
  TH2F *avgChargeMap_e;
  TH2F *avgChargeMap_h;

  G4double elec;
  G4double fluence;
  G4double temperature;
  G4double trappingTimeElectrons;
  G4double trappingTimeHoles;
  G4double betaElectrons;
  G4double tuning;
  G4double threshold;
  G4int precision;

  // Physics process switches
  G4bool doDiff;
  G4bool doRamo;
  G4bool doRamoNeighbors;
  G4bool doChunkCorrection;
  G4bool isHole;

  // Geometry-related constants
  G4double pitchX;
  G4double pitchY;
  G4int nPixX;
  G4int nPixY;
  G4double detectorThickness;

  G4double GetElectricField(G4double x, G4double y);
  G4double GetMobility(G4double electricFiled, G4bool isHoleBit);
  G4double GetDriftTime(G4bool isHoleBit);
  G4double GetTimeToElectrode(G4double x, G4double y, G4bool isHoleBit);
  G4double GetTrappingPositionX(G4double initX, G4double initY, G4double driftTime, G4bool isHoleBit);
  G4double GetTrappingPositionY(G4double initX, G4double initY, G4double driftTime, G4bool isHoleBit);
  G4double GetRamoPotential(G4double x, G4double y);

};

#endif
