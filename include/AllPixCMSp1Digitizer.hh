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

  // Silicon electron and hole transport constants
  G4double Electron_Mobility;
  G4double Electron_Diffusion;
  G4double Boltzmann_kT;
  G4double Electron_ec;
  
  
  vector<G4double> RKF5Integration(G4ThreeVector position, G4double dt);
  G4double MobilityElectron(const G4ThreeVector efield);
  G4ThreeVector ElectronSpeed(const G4ThreeVector efield);
  G4double DiffusionWidth(const G4double timestep);
  
  G4double Propagation(G4ThreeVector& pos);



};

#endif
