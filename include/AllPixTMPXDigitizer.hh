/**
 * Author:
 *    nalipour@cern.ch
 *
 *  allpix Authors:
 *   John Idarraga <idarraga@cern.ch>
 *   Mathieu Benoit <Mathieu.Benoit@cern.ch>
 */

#ifndef AllPixTMPXDigitizer_h
#define AllPixTMPXDigitizer_h 1

// allpix Interface
#include "AllPixDigitizerInterface.hh"
// digits for this digitizer
#include "AllPixTMPXDigit.hh"
#include "G4PrimaryVertex.hh"

#include <map>
#include <vector>

using namespace std;

/**
 *  Digitizer AllPixTMPX implementation
 */
class AllPixTMPXDigitizer : public  AllPixDigitizerInterface {

public:

  AllPixTMPXDigitizer(G4String, G4String, G4String);
  virtual ~AllPixTMPXDigitizer();

  void SetPrimaryVertex(G4PrimaryVertex * pv) {m_primaryVertex = pv;};
  void Digitize ();
  void SetDetectorDigitInputs(G4double){};

private:

  class MC_content
  {
  public:
	  MC_content()
	  {
		  posX_WithRespectToPixel=-11.0;
		  posY_WithRespectToPixel=-11.0;
		  posZ_WithRespectToPixel=-11.0;
		  MC_energy=0.0;
	  };
	  Double_t posX_WithRespectToPixel;
	  Double_t posY_WithRespectToPixel;
	  Double_t posZ_WithRespectToPixel;
	  Double_t MC_energy;
  };
  // digitInput typedef is defined in AllPixDigitizerInterface.hh
  digitInput m_digitIn;

  AllPixTMPXDigitsCollection * m_digitsCollection;
  vector<G4String> m_hitsColName;
  G4PrimaryVertex * m_primaryVertex; // information from EventAction

};

#endif
