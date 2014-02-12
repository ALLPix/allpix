//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef AllPixDetectorMessenger_h
#define AllPixDetectorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class AllPixDetectorConstruction;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithoutParameter;
class G4UIcmdWith3VectorAndUnit;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithAString;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class AllPixDetectorMessenger: public G4UImessenger
{
public:
  AllPixDetectorMessenger(AllPixDetectorConstruction* );
  ~AllPixDetectorMessenger();
  
  void SetNewValue(G4UIcommand*, G4String);
  
private:
  AllPixDetectorConstruction  * m_AllPixDetector;
  
  G4UIdirectory *             m_allpixDir;
  G4UIdirectory *             m_detDir;
  G4UIdirectory *             m_configDir;
  G4UIdirectory *             m_extrasDir;
  G4UIdirectory *             m_beamDir;

  G4UIcmdWithAnInteger *      m_detIdCmd;
  G4UIcmdWith3VectorAndUnit * m_detPosCmd;
  G4UIcmdWith3VectorAndUnit * m_detRotCmd;
  G4UIcmdWith3VectorAndUnit * m_testStructPosCmd;
  G4UIcmdWith3VectorAndUnit * m_testStructRotCmd;
  G4UIcmdWithAnInteger *      m_testStructDetLinkCmd;
  G4UIcmdWithAnInteger *      m_testStructType;
  G4UIcmdWithAnInteger *      m_AppliancesType;
  G4UIcmdWith3VectorAndUnit * m_detAppliancePosCmd;
  G4UIcmdWith3VectorAndUnit * m_wrapperEnhancementCmd;
  G4UIcmdWithADoubleAndUnit * m_magFieldCmd;

  G4UIcmdWithAString * m_worldMaterial;

  G4UIcmdWithAString * m_outputPrefix;

  G4UIcmdWithADoubleAndUnit * m_HighTHLCmd;
  G4UIcmdWithADoubleAndUnit * m_LowTHLCmd;
  G4UIcmdWithADoubleAndUnit * m_AcqTimeCmd;
  G4UIcmdWithADoubleAndUnit * m_HVCmd;
  G4UIcmdWithADoubleAndUnit * m_ClockCmd;
  G4UIcmdWithADoubleAndUnit * m_StepLengthSensor;

  G4UIcmdWithoutParameter   * m_UpdateCmd;

#ifdef _EUTELESCOPE
	// Specific EUTelescope
  G4UIcmdWith3VectorAndUnit * m_scint1PosCmd;
  G4UIcmdWith3VectorAndUnit * m_scint2PosCmd;
  G4UIcmdWith3VectorAndUnit * m_scint3PosCmd;
  G4UIcmdWith3VectorAndUnit * m_scint4PosCmd;
#endif


};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

