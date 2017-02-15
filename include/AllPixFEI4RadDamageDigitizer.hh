/*
 *  Authors:
 *    Callie Bertsche <c.bertsche@cern.ch>
 *    Benjamin Nachman <bnachman@cern.ch>
 *
 *  allpix Authors:
 *   John Idarraga <idarraga@cern.ch>
 *   Mathieu Benoit <benoit@lal.in2p3.fr>
 */

#ifndef AllPixFEI4RadDamageDigitizer_h
#define AllPixFEI4RadDamageDigitizer_h 1

#include "AllPixDigitizerInterface.hh"
#include "AllPixFEI4RadDamageDigit.hh"
#include "G4PrimaryVertex.hh"
#include <map>
#include <vector>
#include "AllPixTrackerHit.hh"
#include "AllPixGeoDsc.hh"
#include "TString.h"
#include "TH2D.h"
#include "TH3F.h"
#include "TTree.h"
#include "TFile.h"
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

  TH3F *ramoPotentialMap;
  TH3F *eFieldMap;
  TH1F *m_eFieldMap1D;
  TH1F *timeMap_e;
  TH1F *timeMap_h;
  TH2F *distancemap_e;
  TH2F *distancemap_h;
  TH2F *lorentz_map_e;
  TH2F *lorentz_map_h;
  TH3F *charge_chunk_map_e;
  TH3F *charge_chunk_map_h;

  G4double GetElectricField(G4double z);
  G4double GetElectricField(G4double x, G4double y, G4double z);
  G4double GetMobility(G4double electricField, G4bool isHole);
  G4double GetDriftVelocity(G4double electricField, G4double mobility, G4bool isHole);
  G4double GetMeanFreePath(G4double driftVelocity, G4bool isHole);
  G4double GetTrappingProbability(G4double z, G4double meanFreePath);
  G4double GetDriftTime(G4bool isHole);
  G4double GetTimeToElectrode(G4double z, G4bool isHole);
  G4double GetTanLorentz(G4double electricField, G4bool isHole);
  G4double GetTanLorentz(G4double z1, G4double z2, G4bool isHole);
  G4double Phi(G4double x, G4double z, G4double Lx, G4double Lz);
  G4double N(G4double z, G4double L);
  G4double elec;
  G4int count;

  // Variables for the charge sharing computation
  G4double mobility;
  G4double resistivity;
  G4bool   bulkType;

  G4double depVoltage;
  G4double deplationLenght;

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
  G4double tuning;
  G4double threshold;
  G4double diffusion_length;
  
  G4double epsilon;
  G4double echarge;
  G4int precision;

  // Physics process switches
  G4bool doTrapping;
  G4bool doRamo;
  G4bool doDiff;
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
};

#endif
