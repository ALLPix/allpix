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
#include <fstream>

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

  AllPixGeoDsc * gD;

  // digitInput typedef is defined in AllPixDigitizerInterface.hh
  digitInput m_digitIn;

  AllPixCMSp1DigitsCollection * m_digitsCollection;
  vector<G4String> m_hitsColName;
  G4PrimaryVertex * m_primaryVertex; // information from EventAction

  G4double threshold;
  G4double elec;
  G4double Temperature;
  G4ThreeVector bfield;
  G4double detectorThickness;
  
  G4double Target_Spatial_Precision;
  G4double Timestep_max;
  G4double Timestep_min;


  // Silicon electron and hole transport constants
  G4double Electron_Mobility;
  G4double Electron_HallFactor;
  G4double Electron_Beta;
  G4double Boltzmann_kT;
  G4double Electron_ec;
  
  G4double flux;
  G4double Electron_Trap_beta0;
	G4double Electron_Trap_kappa;
	G4double Electron_Trap_T0;
	G4double Electron_Trap_TauNoFluence;
  G4double Electron_Trap_TauEff;
  
  void InitVariables();
  
  vector<G4double> RKF5Integration(G4ThreeVector position, G4double dt);
  G4double MobilityElectron(const G4ThreeVector efield);
  G4ThreeVector ElectronSpeed(const G4ThreeVector efield);
  G4ThreeVector DiffusionStep(const G4double timestep, const G4ThreeVector position);
  void SetDt(G4double& dt, const G4double uncertainty, const G4double z, const G4double dz);
  G4double GetTrappingTime();
  
  G4double Propagation(G4ThreeVector& pos, G4double& drifttime, G4bool& trapped);



};

#endif
