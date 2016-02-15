/**
 * Author:
 *    Paul Schuetze <paul.schuetze@desy.de>
 *
 *  allpix Authors:
 *   John Idarraga <idarraga@cern.ch>
 *   Mathieu Benoit <Mathieu.Benoit@cern.ch>
 */

#ifndef AllPixCMSp1Digitizer_h
#define AllPixCMSp1Digitizer_h 1

// allpix Interface
#include "AllPixDigitizerInterface.hh"
// digits for this digitizer
#include "AllPixCMSp1Digit.hh"
#include "G4PrimaryVertex.hh"

#include <map>
#include <vector>

using namespace std;

/**
 *  Digitizer AllPixCMSp1 implementation
 */
class AllPixCMSp1Digitizer : public  AllPixDigitizerInterface {

public:

  AllPixCMSp1Digitizer(G4String, G4String, G4String);
  virtual ~AllPixCMSp1Digitizer();

  void SetPrimaryVertex(G4PrimaryVertex * pv) {m_primaryVertex = pv;};
  void Digitize ();
  void SetDetectorDigitInputs(G4double){};

private:

  // digitInput typedef is defined in AllPixDigitizerInterface.hh
  digitInput m_digitIn;

  AllPixCMSp1DigitsCollection * m_digitsCollection;
  vector<G4String> m_hitsColName;
  G4PrimaryVertex * m_primaryVertex; // information from EventAction

};

#endif
