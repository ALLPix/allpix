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

#include "TFile.h" //calibration
#include "TTree.h" //calibration

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

  AllPixTMPXDigitsCollection * m_digitsCollection;
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
  Double_t Threshold;

  G4double **ThresholdMatrix; 
  G4double **SurrogateA;
  G4double **SurrogateB;
  G4double **SurrogateC;
  G4double **SurrogateD;
  G4int size_calibX;
  G4int size_calibY;


  //Charge sharing parameters
  // G4double epsilon = 11.8*8.854187817e-14; // [F/cm] -> F=As/V (silicon)
  // G4double echarge=1.60217646e-19; //[C=As]
  // G4double Default_Hole_Mobility=480.0; //[cm2/Vs] Hole mobility
  // G4double Default_Hole_D=12; //;// Hole diffusion [cm2/s]

  // G4double Default_Electron_Mobility=1415.0; //[cm2/Vs] Electron mobility
  // G4double Default_Electron_D=36; //;// Electron diffusion [cm2/s]
  // G4double temperature_T=300.0;

  G4double epsilon; // [F/cm] -> F=As/V (silicon)
  G4double echarge; //[C=As]
  G4double Default_Hole_Mobility; //[cm2/Vs] Hole mobility
  G4double Default_Hole_D; //;// Hole diffusion [cm2/s]
  G4double Default_Electron_Mobility; //[cm2/Vs] Electron mobility
  G4double Default_Electron_D; //;// Electron diffusion [cm2/s]
  G4double temperature_T;
  G4double elec;//=3.6*eV;

  G4double vm;
  G4double Ec;
  G4double beta;

  G4double mobility_const;
  G4double diffusion_const;
  G4double V_D;
  G4double depletionWidth;

  void initialise_calibrationMatrices(int size_calibX, int size_calibY);
  
};

#endif
