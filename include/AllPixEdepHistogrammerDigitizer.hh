/**
 * Author:
 *    Mathieu.Benoit@CERN.CH
 *
 *  allpix Authors:
 *   John Idarraga <idarraga@cern.ch>
 *   Mathieu Benoit <Mathieu.Benoit@cern.ch>
 */

#ifndef AllPixEdepHistogrammerDigitizer_h
#define AllPixEdepHistogrammerDigitizer_h 1

// allpix Interface
#include "AllPixDigitizerInterface.hh"
// digits for this digitizer
#include "AllPixEdepHistogrammerDigit.hh"
#include "G4PrimaryVertex.hh"
#include "TH1D.h"
#include <map>
#include <vector>

using namespace std;

/**
 *  Digitizer AllPixEdepHistogrammer implementation
 */
class AllPixEdepHistogrammerDigitizer : public  AllPixDigitizerInterface {

public:

  AllPixEdepHistogrammerDigitizer(G4String, G4String, G4String);
  ~AllPixEdepHistogrammerDigitizer();

  void SetPrimaryVertex(G4PrimaryVertex * pv) {m_primaryVertex = pv;};
  void Digitize ();
  void SetDetectorDigitInputs(G4double){};

private:

  // digitInput typedef is defined in AllPixDigitizerInterface.hh
  digitInput m_digitIn;

  AllPixEdepHistogrammerDigitsCollection * m_digitsCollection;
  vector<G4String> m_hitsColName;
  G4PrimaryVertex * m_primaryVertex; // information from EventAction
  TH1D * Edeposition;

};

#endif
