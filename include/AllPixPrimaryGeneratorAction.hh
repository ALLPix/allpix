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
// $Id: AllPixPrimaryGeneratorAction.hh,v 1.6 2006/06/29 17:54:04 gunter Exp $
// GEANT4 tag $Name: geant4-09-02 $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef AllPixPrimaryGeneratorAction_h
#define AllPixPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4VPrimaryGenerator.hh"
#include "G4GeneralParticleSource.hh"
#include "G4ThreeVector.hh"

#include "globals.hh"

#include <vector>

using namespace std;

class G4ParticleGun;
class G4Event;
class AllPixPrimaryGeneratorMessenger;
class G4VPrimaryGenerator;

enum SourceType {
  _ParticleGun = 0,
  _GeneralParticleSource,
  _HEPEvtInterface
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class AllPixPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
  AllPixPrimaryGeneratorAction(SourceType);
  ~AllPixPrimaryGeneratorAction();
  
public:
  void GeneratePrimaries(G4Event*);
  void SetHEPevtFile(G4String filename);
  void SetOptPhotonPolar();
  void SetOptPhotonPolar(G4double);
  AllPixPrimaryGeneratorMessenger * GetPrimaryGeneratorMessenger() {return this->m_gunMessenger;}
  
private:

  // simple particule gun case
  G4ParticleGun* m_particleGun;
  // using gps
  G4GeneralParticleSource * m_particleSource;
  // using external MC data
  G4VPrimaryGenerator * m_HEPEvt;

  AllPixPrimaryGeneratorMessenger* m_gunMessenger;

  // store temporarily particle positions
  vector<G4ThreeVector> m_primaryParticlePos;

  SourceType m_sType;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif /*AllPixPrimaryGeneratorAction_h*/
