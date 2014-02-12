/**
 * Author:
 *    AAAuthor
 *
 *  allpix Authors:
 *   John Idarraga <idarraga@cern.ch>
 *   Mathieu Benoit <Mathieu.Benoit@cern.ch>
 */

#ifndef __Digitizer_h
#define __Digitizer_h 1

// allpix Interface
#include "AllPixDigitizerInterface.hh"
// digits for this digitizer
#include "__Digit.hh"
#include "G4PrimaryVertex.hh"

#include <map>
#include <vector>

using namespace std;

/**
 *  Digitizer __ implementation
 */
class __Digitizer : public  AllPixDigitizerInterface {

public:

  __Digitizer(G4String, G4String, G4String);
  virtual ~__Digitizer();

  void SetPrimaryVertex(G4PrimaryVertex * pv) {m_primaryVertex = pv;};
  void Digitize ();
  void SetDetectorDigitInputs(G4double){};

private:

  // digitInput typedef is defined in AllPixDigitizerInterface.hh
  digitInput m_digitIn;

  __DigitsCollection * m_digitsCollection;
  vector<G4String> m_hitsColName;
  G4PrimaryVertex * m_primaryVertex; // information from EventAction

};

#endif
