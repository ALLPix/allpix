/**
 * Author:
 *    
 *
 *  allpix Authors:
 *   John Idarraga <idarraga@cern.ch>
 *   Mathieu Benoit <benoit@lal.in2p3.fr>
 */

#ifndef AllPixTestBichselDigitizer_h
#define AllPixTestBichselDigitizer_h 1

// allpix Interface
#include "AllPixDigitizerInterface.hh"
// digits for this digitizer
#include "AllPixTestBichselDigit.hh"
#include "G4PrimaryVertex.hh"
#include "AllPixBichselTool.hh"

#include <map>
#include <vector>

using namespace std;

/**
 *  Digitizer AllPixTestBichsel implementation
 */
class AllPixTestBichselDigitizer : public  AllPixDigitizerInterface {

public:

  AllPixTestBichselDigitizer(G4String, G4String, G4String);
  virtual ~AllPixTestBichselDigitizer();

  void SetPrimaryVertex(G4PrimaryVertex * pv) {m_primaryVertex = pv;};
  void Digitize ();
  void SetDetectorDigitInputs(G4double){};

private:

  // digitInput typedef is defined in AllPixDigitizerInterface.hh
  digitInput m_digitIn;

  AllPixTestBichselDigitsCollection * m_digitsCollection;
  vector<G4String> m_hitsColName;
  G4PrimaryVertex * m_primaryVertex; // information from EventAction
  bool doBichsel;
  AllPixBichselTool Bichsel = AllPixBichselTool();
  //MyBichselTool Bichsel;
};

#endif
