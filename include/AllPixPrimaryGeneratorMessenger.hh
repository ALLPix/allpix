/**
 *  Author John Idarraga <idarraga@cern.ch>
 */
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
//
// $Id: AllPixPrimaryGeneratorMessenger.hh,v 1.2 2006/06/29 17:54:07 gunter Exp $
// GEANT4 tag $Name: geant4-09-02 $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef AllPixPrimaryGeneratorMessenger_h
#define AllPixPrimaryGeneratorMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

class AllPixPrimaryGeneratorAction;
class G4UIdirectory;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;
class G4UIcmdWithABool;
class G4UIcmdWithoutParameter;
class G4UIcommand;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class AllPixPrimaryGeneratorMessenger: public G4UImessenger
{
  public:
	AllPixPrimaryGeneratorMessenger(AllPixPrimaryGeneratorAction*);
	~AllPixPrimaryGeneratorMessenger();

	G4bool   GetTimepixTelescopeWriteFlag()   {return this->m_TimepixTelescopeWriteFlag;}
	G4String GetTimepixTelescopeFolderName()  {return this->m_TimepixTelescopeFolderName;}
	G4bool   GetTimepixTelescopeDoEventFlag() {return this->m_TimepixTelescopeDoEventFlag;}
	G4bool   GetTimepixTelescopeSumTOTFlag()  {return this->m_TimepixTelescopeSumTOTFlag;}

	void SetNewValue(G4UIcommand*, G4String);

  private:
    AllPixPrimaryGeneratorAction * AllPixAction;
    G4UIdirectory                * gunDir;
    G4UIcmdWithADoubleAndUnit    * polarCmd;
    G4UIcmdWithAString           * m_MCPrefix;

    G4UIcmdWithAnInteger         * m_userBeamNumberOfFramesCmd;
    G4UIcmdWithoutParameter      * m_userBeamOnCmd;
    G4UIcommand                  * m_userBeamTypeCmd;
    G4int m_hits;
    G4int m_frames;
		G4String m_beamTypeHitFunc;
		G4double m_beamTypePar1;
		G4double m_beamTypePar2;

	G4UIcmdWithABool   * m_TimepixTelescopeWriteCmd;
	G4UIcmdWithAString * m_TimepixTelescopeFolderNameCmd;
	G4UIcmdWithABool   * m_TimepixTelescopeDoEventCmd;
	G4UIcmdWithABool   * m_TimepixTelescopeSumTOTCmd;
	G4bool   m_TimepixTelescopeWriteFlag;
	G4String m_TimepixTelescopeFolderName;
	G4bool   m_TimepixTelescopeDoEventFlag;
	G4bool   m_TimepixTelescopeSumTOTFlag;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

