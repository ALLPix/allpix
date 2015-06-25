/**
 * Author:
 *    idarraga@cern.ch
 *
 *  allpix Authors:
 *   John Idarraga <idarraga@cern.ch>
 *   Mathieu Benoit <Mathieu.Benoit@cern.ch>
 */

#ifndef AllPixMedipix3RXDigitizer_h
#define AllPixMedipix3RXDigitizer_h 1

// allpix Interface
#include "AllPixDigitizerInterface.hh"
// digits for this digitizer
#include "AllPixMedipix3RXDigit.hh"
#include "G4PrimaryVertex.hh"

#include <map>
#include <vector>

using namespace std;

/**
 *  Digitizer AllPixMedipix3RX implementation
 */
class AllPixMedipix3RXDigitizer : public  AllPixDigitizerInterface {

public:

  AllPixMedipix3RXDigitizer(G4String, G4String, G4String);
  virtual ~AllPixMedipix3RXDigitizer();

  void SetPrimaryVertex(G4PrimaryVertex * pv) {m_primaryVertex = pv;};
  void Digitize ();
  void SetDetectorDigitInputs(G4double){};

private:

  // digitInput typedef is defined in AllPixDigitizerInterface.hh
  digitInput m_digitIn;

  AllPixMedipix3RXDigitsCollection * m_digitsCollection;
  vector<G4String> m_hitsColName;
  G4PrimaryVertex * m_primaryVertex; // information from EventAction

};

#endif
