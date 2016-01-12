/**
 * Author:
 *    nalipour
 *
 *  allpix Authors:
 *   John Idarraga <idarraga@cern.ch>
 *   Mathieu Benoit <Mathieu.Benoit@cern.ch>
 */

#ifndef AllPixNoThreshDigitizer_h
#define AllPixNoThreshDigitizer_h 1

// allpix Interface
#include "AllPixDigitizerInterface.hh"
// digits for this digitizer
#include "AllPixNoThreshDigit.hh"
#include "G4PrimaryVertex.hh"

#include <map>
#include <vector>


#include "TFile.h" //calibration
#include "TTree.h" //calibration

using namespace std;

/**
 *  Digitizer AllPixNoThresh implementation
 */
class AllPixNoThreshDigitizer : public  AllPixDigitizerInterface {

public:

  AllPixNoThreshDigitizer(G4String, G4String, G4String);
  virtual ~AllPixNoThreshDigitizer();

  void SetPrimaryVertex(G4PrimaryVertex * pv) {m_primaryVertex = pv;};
  void Digitize ();
  void SetDetectorDigitInputs(G4double){};
  G4double IntegrateGaussian(G4double xhit, G4double yhit, G4double Sigma, G4double x1, G4double x2, G4double y1, G4double y2, G4double Energy);

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

  AllPixNoThreshDigitsCollection * m_digitsCollection;
  vector<G4String> m_hitsColName;
  G4PrimaryVertex * m_primaryVertex; // information from EventAction

  G4double thickness; //detector thickness
  G4double pitchX;
  G4double pitchY;
  G4int nPixX;
  G4int nPixY;

  //nalipour
  Double_t resistivity;
  Double_t V_B;
  G4String SensorType; // p-in-n or n-in-p
  G4String CalibrationFile;

  G4double **ThresholdMatrix; 
  G4double **SurrogateA;
  G4double **SurrogateB;
  G4double **SurrogateC;
  G4double **SurrogateD;
  G4int size_calibX;
  G4int size_calibY;

  void initialise_calibrationMatrices(int size_calibX, int size_calibY);


};

#endif
