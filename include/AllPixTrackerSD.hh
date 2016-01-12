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
// $Id: ExN02TrackerSD.hh,v 1.7 2006/06/29 17:47:56 gunter Exp $
// GEANT4 tag $Name: geant4-09-02 $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef AllPixTrackerSD_h
#define AllPixTrackerSD_h 1

#include "G4VSensitiveDetector.hh"
#include "AllPixTrackerHit.hh"
#include "G4WrapperProcess.hh"

#include <set>

using namespace std;

class G4Step;
class G4HCofThisEvent;
class G4Event;
class AllPixGeoDsc;

#define MAX_CHAMBERS_EPIX 20

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class AllPixTrackerSD : public G4VSensitiveDetector
{
  
public:
  AllPixTrackerSD(G4String, G4ThreeVector, G4ThreeVector, AllPixGeoDsc *, G4RotationMatrix *);
  AllPixTrackerSD(G4String, G4ThreeVector, G4RotationMatrix *);
  ~AllPixTrackerSD();
  
  void Initialize(G4HCofThisEvent*);
  G4bool ProcessHits(G4Step*, G4TouchableHistory*);
  void EndOfEvent(G4HCofThisEvent*);
  
  G4String GetHitsCollectionName(){ return m_thisHitsCollectionName; };

private:

  AllPixTrackerHitsCollection* hitsCollection;
  G4ThreeVector m_absolutePosOfWrapper; // Absolute position of Wrapper
  G4ThreeVector m_relativePosOfSD;      // Relative (to Wrapper) position of SD
  G4RotationMatrix * m_rotationOfWrapper;  // rotation of Wrapper
  AllPixGeoDsc * m_gD; // Geo description !
  G4String m_thisHitsCollectionName;
  bool m_thisIsAPixelDetector;
  bool firstStrikePrimary;
  G4double _kinEPrimary;
  G4double _totalEdep;

  // used to dump tracking info in special cases
  long m_globalTrackId_Dump;

  set<AllPixTrackerHitsCollection *> m_hitsCollectionSet;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

