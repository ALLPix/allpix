/**
 * Author:
 *    bnachman@cern.ch
 *
 *  allpix Authors:
 *   John Idarraga <idarraga@cern.ch>
 *   Mathieu Benoit <Mathieu.Benoit@cern.ch>
 */

#ifndef AllPixTruthWithLorentzAndTOTDigitizer_h
#define AllPixTruthWithLorentzAndTOTDigitizer_h 1

// allpix Interface
#include "AllPixDigitizerInterface.hh"
// digits for this digitizer
#include "AllPixTruthWithLorentzAndTOTDigit.hh"
#include "G4PrimaryVertex.hh"

#include <map>
#include <vector>

using namespace std;

/**
 *  Digitizer AllPixTruthWithLorentzAndTOT implementation
 */
class AllPixTruthWithLorentzAndTOTDigitizer : public  AllPixDigitizerInterface {

public:

  AllPixTruthWithLorentzAndTOTDigitizer(G4String, G4String, G4String);
  virtual ~AllPixTruthWithLorentzAndTOTDigitizer();

  void SetPrimaryVertex(G4PrimaryVertex * pv) {m_primaryVertex = pv;};
  void Digitize ();
  void SetDetectorDigitInputs(G4double){};

private:

  // digitInput typedef is defined in AllPixDigitizerInterface.hh
  digitInput m_digitIn;

  AllPixTruthWithLorentzAndTOTDigitsCollection * m_digitsCollection;
  vector<G4String> m_hitsColName;
  G4PrimaryVertex * m_primaryVertex; // information from EventAction

};

#endif
