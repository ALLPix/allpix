/**
 * Author:
 *    Ben Nachman <bnachman@cern.ch>
 *
 *  allpix Authors:
 *   John Idarraga <idarraga@cern.ch>
 *   Mathieu Benoit <Mathieu.Benoit@cern.ch>
 */

#ifndef AllPixTruthWithLorentzDigitizer_h
#define AllPixTruthWithLorentzDigitizer_h 1

// allpix Interface
#include "AllPixDigitizerInterface.hh"
// digits for this digitizer
#include "AllPixTruthWithLorentzDigit.hh"
#include "G4PrimaryVertex.hh"

#include <map>
#include <vector>
#include "TH1F.h"

using namespace std;

/**
 *  Digitizer AllPixTruthWithLorentz implementation
 */
class AllPixTruthWithLorentzDigitizer : public  AllPixDigitizerInterface {

public:

  AllPixTruthWithLorentzDigitizer(G4String, G4String, G4String);
  virtual ~AllPixTruthWithLorentzDigitizer();

  void SetPrimaryVertex(G4PrimaryVertex * pv) {m_primaryVertex = pv;};
  G4double GetMobility(G4double electricField, G4bool isHole, G4double temperature);
  void Digitize ();
  void SetDetectorDigitInputs(G4double){};

  TH1F* eFieldMap;
  TH1F* timeMap_e;
  G4double biasVoltage;
  G4double Temperature;
  G4double bfield;
  G4double tanLorentz_e;
  G4double mobility_e;
  bool doPropagateFundamental;
  bool doConstantEfield;
  bool doDiffusion;
  bool doDigitize;
  int nToTbits;
  int chargeThreshold;
  int chargeTuning;
  int method;

private:

  // digitInput typedef is defined in AllPixDigitizerInterface.hh
  digitInput m_digitIn;

  AllPixTruthWithLorentzDigitsCollection * m_digitsCollection;
  vector<G4String> m_hitsColName;
  G4PrimaryVertex * m_primaryVertex; // information from EventAction

};

#endif
