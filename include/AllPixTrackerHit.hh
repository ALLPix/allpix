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
// $Id: ExN02TrackerHit.hh,v 1.8 2006/06/29 17:47:53 gunter Exp $
// GEANT4 tag $Name: geant4-09-02 $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef AllPixTrackerHit_h
#define AllPixTrackerHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4Track.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class AllPixTrackerHit : public G4VHit
{
public:
  
  AllPixTrackerHit();
  ~AllPixTrackerHit();
  AllPixTrackerHit(const AllPixTrackerHit&);
  const AllPixTrackerHit& operator=(const AllPixTrackerHit&);
  G4int operator==(const AllPixTrackerHit&) const;
  
  inline void* operator new(size_t);
  inline void  operator delete(void*);
  
  void Draw();
  void Print();
  
public:
  
  void SetTrackID  (G4int track)      { trackID = track; };
  void SetParentID  (G4int track)      { parentID = track; };
  void SetPixelNbX(G4int pix)      { pixelNbX = pix; };
  void SetPixelNbY(G4int pix)      { pixelNbY = pix; };
  void SetPostPixelNbX(G4int postPix) { postPixelNbX = postPix; };
  void SetPostPixelNbY(G4int postPix) { postPixelNbY = postPix; };
  void SetEdep     (G4double de)      { edep = de; };
  void SetPos      (G4ThreeVector xyz){ pos = xyz; };
  void SetPosWithRespectToPixel (G4ThreeVector pxzy) { m_posWithRespectToPixel = pxzy; };
  void SetProcessName(G4String process) { processName = process; };
  void SetTrackPdgId(G4int pdgId) { pdgIdTrack = pdgId; };
  void SetTrackVolumeName(G4String vn) { trackVolumeName = vn; };
  void SetParentVolumeName(G4String vn) { parentVolumeName = vn; } ;
  void SetDetId(G4int i){ detID = i; };
  void SetKinEParent(G4double kinE){ kinEParent = kinE; };
  //void SetKineticEnergy     (G4double kinE) {kinE = kinE; };
  void SetIncidentPOLAngle(G4double inc) {incidentPOLangle = inc;};
  void SetIncidentAZMAngle(G4double inc) {incidentAZMangle = inc;};
  void SetTruthEntryLocalX(G4double loc) {TruthEntryLocalX = loc;};
  void SetTruthEntryLocalY(G4double loc) {TruthEntryLocalY = loc;};

  G4int GetTrackID()    { return trackID; };
  G4int GetParentID()    { return parentID; };
  G4int GetPixelNbX()  { return pixelNbX; };
  G4int GetPixelNbY()  { return pixelNbY; };
  G4int GetPostPixelNbX()  { return postPixelNbX; };
  G4int GetPostPixelNbY()  { return postPixelNbY; };
  G4int GetDetId() { return detID; };
  G4double GetEdep()    { return edep; };
  G4ThreeVector GetPos(){ return pos; };
  G4ThreeVector GetPosWithRespectToPixel() { return m_posWithRespectToPixel; };
  G4String GetProcessName() { return processName; };
  G4int GetTrackPdgId() { return pdgIdTrack; };
  G4String GetTrackVolumeName() {return trackVolumeName;};
  G4String GetParentVolumeName() {return parentVolumeName;};
  G4double GetKinEParent() { return kinEParent; };
  G4double GetIncidentPOLAngle() {return incidentPOLangle;};
  G4double GetIncidentAZMAngle() {return incidentAZMangle;};
  G4double GetTruthEntryLocalX() {return TruthEntryLocalX;};
  G4double GetTruthEntryLocalY() {return TruthEntryLocalY;};
  //G4double GetKineticEnergy(){ return kinE; };
  
private:
  
  G4double TruthEntryLocalX;
  G4double TruthEntryLocalY;
  G4int         trackID;
  G4int         detID;
  G4int         parentID;
  G4int         pixelNbX;
  G4int         pixelNbY;
  G4int         postPixelNbX;
  G4int         postPixelNbY;
  G4double      edep;
  G4ThreeVector pos;
  G4ThreeVector m_posWithRespectToPixel;
  G4String      processName;
  G4int         pdgIdTrack;
  G4String      trackVolumeName;
  G4String      parentVolumeName;
  G4double		kinEParent;
  G4double incidentPOLangle;
  G4double incidentAZMangle;
  //G4double      kinE;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


typedef G4THitsCollection<AllPixTrackerHit> AllPixTrackerHitsCollection;

extern G4Allocator<AllPixTrackerHit> AllPixTrackerHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void* AllPixTrackerHit::operator new(size_t)
{
  //void *aHit;
  //aHit = (void *) AllPixTrackerHitAllocator.MallocSingle();

	AllPixTrackerHit * aHit;
	aHit = static_cast<AllPixTrackerHit*>(AllPixTrackerHitAllocator.MallocSingle());
	return aHit;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void AllPixTrackerHit::operator delete(void *aHit)
{
  AllPixTrackerHitAllocator.FreeSingle((AllPixTrackerHit*) aHit);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
