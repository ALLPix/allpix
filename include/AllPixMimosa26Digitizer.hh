/**
 *  Author John Idarraga <idarraga@cern.ch>
 */

#ifndef AllPixMimosa26Digitizer_h
#define AllPixMimosa26Digitizer_h 1

// interface
#include "AllPixDigitizerInterface.hh"
// digits for this digitizer
#include "AllPixMimosa26Digit.hh"

#include "G4PrimaryVertex.hh"

#include <map>
#include <vector>

using namespace std;

#define __SIG_SIZE	40

/**
 *  A simple Digitizer implementation
 */

class AllPixMimosa26Digitizer : public  AllPixDigitizerInterface {

public:
  AllPixMimosa26Digitizer(G4String, G4String, G4String);
  virtual ~AllPixMimosa26Digitizer();

  void SetPrimaryVertex(G4PrimaryVertex * pv) {m_primaryVertex = pv;};
  void Digitize ();
  void SetDetectorDigitInputs(G4double);
  int  indexofSmallestElement(double array[], int size);
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
  digitInput m_digitIn;
  AllPixMimosa26DigitsCollection * m_digitsCollection;
  vector<G4String> m_hitsColName;
  G4PrimaryVertex * m_primaryVertex; // information from EventAction
  G4int m_randomNoise;

  //////////////////////////////////////////////////////
  // Geometry Related constants
  G4double pitchX ;
  G4double pitchY ;
  G4double hpitchX ;
  G4double hpitchY ;
  G4int nPixX;
  G4int nPixY;
  //TRandom2 * m_trand;
  G4double depletionDepth;
  G4double thickness;
  G4double crosstalk_prob;

  G4double * xsig; //G4double xsig[5];
  G4double * ysig; //G4double ysig[5];


};

#endif
