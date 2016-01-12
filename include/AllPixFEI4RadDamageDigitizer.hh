/**
 * Author:
 *    Callie Bertsche <c.bertsche@cern.ch>
 *
 *  allpix Authors:
 *   John Idarraga <idarraga@cern.ch>
 *   Mathieu Benoit <Mathieu.Benoit@cern.ch>
 */

#ifndef AllPixFEI4RadDamageDigitizer_h
#define AllPixFEI4RadDamageDigitizer_h 1

// allpix Interface
#include "AllPixDigitizerInterface.hh"
// digits for this digitizer
#include "AllPixFEI4RadDamageDigit.hh"
#include "G4PrimaryVertex.hh"

#include <map>
#include <vector>

// added for radiation damage
#include "AllPixTrackerHit.hh"
#include "AllPixGeoDsc.hh"
#include "TString.h"
#include "TH2D.h"
#include "TH3F.h"
#include "TH1F.h"

using namespace std;

class AllPixFEI4RadDamageDigitizer : public  AllPixDigitizerInterface {

public:
  AllPixFEI4RadDamageDigitizer(G4String, G4String, G4String);
  virtual ~AllPixFEI4RadDamageDigitizer();

  void SetPrimaryVertex(G4PrimaryVertex * pv) {m_primaryVertex = pv;};
  void Digitize ();
//  void SetDetectorDigitInputs(G4double){};
  void SetDetectorDigitInputs(G4double);


private:

  // digitInput typedef is defined in AllPixDigitizerInterface.hh
  digitInput m_digitIn;

  AllPixFEI4RadDamageDigitsCollection * m_digitsCollection;
  vector<G4String> m_hitsColName;
  G4PrimaryVertex * m_primaryVertex; // information from EventAction

//***Code added for radiation damage***//

  TH3F *ramoPotentialMap;
  TH1F *eFieldMap;
  TH1F *timeMap_e;
  TH1F *timeMap_h;

  G4double GetElectricField(G4double z);
  G4double GetMobility(G4double electricField, G4double Temperature, G4bool isHoleBit);
  G4double GetDriftVelocity(G4double electricField, G4double mobility, G4bool isHoleBit);
  G4double GetMeanFreePath(G4double driftVelocity, G4bool isHoleBit);
  G4double GetTrappingProbability(G4double z, G4double meanFreePath,G4bool isHoleBit);
  G4double GetDriftTime(G4bool isHoleBit);
  G4double GetTimeToElectrode(G4double z, G4bool isHoleBit);

  G4int EnergyToTOT(G4double Energy, G4double threshold);
  G4double SlimEdgeEffect(G4int nX,G4double xpos,G4double eHit);
  G4bool isSlimEdge(G4int nX, G4int nY);

  G4double elec;

  // Variables for the charge sharing computation
  G4double mobility;
  G4double resistivity;
  G4bool   bulkType;

  G4double detectorThickness;
  G4double biasVoltage;
  G4double temperature;
  G4double fluence;
  G4double trappingTimeElectrons;
  G4double trappingTimeHoles;
  G4double betaElectrons;
  G4double betaHoles;
  G4double bField;
  G4double chipNoise;

  G4double epsilon;
  G4double echarge;
  G4int precision;

  // Physics process switches
  G4bool doTrapping;
  G4bool doRamo;
  G4bool doSlimEdge;
  G4bool isHole;

  // Geometry-related constants
  G4double pitchX;
  G4double pitchY;
  G4int nPixX;
  G4int nPixY;
  G4double chargeSharingConstant;
  G4double GRShift;
  G4int FEIX;
  G4int Sensor;

  //Tuning of the chip and counters characteristics
  G4int MipTOT;
  G4int CounterDepth;
  G4int MipCharge;
  G4double Lv1Unit;

  G4bool doDrift;
  G4double diffusion_length;
  G4double threshold;
  G4double tuning;

};

#endif
