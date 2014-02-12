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
// $Id: AllPixTrackerHit.cc,v 1.10 2006/06/29 17:48:24 gunter Exp $
// GEANT4 tag $Name: geant4-09-02 $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "AllPixTrackerHit.hh"
#include "G4UnitsTable.hh"
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4Track.hh"


G4Allocator<AllPixTrackerHit> AllPixTrackerHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

AllPixTrackerHit::AllPixTrackerHit() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

AllPixTrackerHit::~AllPixTrackerHit() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

AllPixTrackerHit::AllPixTrackerHit(const AllPixTrackerHit& right)
: G4VHit()
{
	trackID   = right.trackID;
	detID = right.detID;
	parentID = right.parentID;
	pixelNbX = right.pixelNbX;
	pixelNbY = right.pixelNbY;
	postPixelNbX = right.postPixelNbX;
	postPixelNbY = right.postPixelNbY;
	edep      = right.edep;
	pos       = right.pos;
	m_posWithRespectToPixel = right.m_posWithRespectToPixel;
	processName = right.processName;
	pdgIdTrack = right.pdgIdTrack;
	trackVolumeName = right.trackVolumeName;
	parentVolumeName = right.parentVolumeName;
	//kinE      = right.kinE;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const AllPixTrackerHit& AllPixTrackerHit::operator=(const AllPixTrackerHit& right)
{
	trackID   = right.trackID;
	detID = right.detID;
	parentID = right.parentID;
	pixelNbX = right.pixelNbX;
	pixelNbY = right.pixelNbY;
	postPixelNbX = right.postPixelNbX;
	postPixelNbY = right.postPixelNbY;
	edep      = right.edep;
	pos       = right.pos;
	m_posWithRespectToPixel = right.m_posWithRespectToPixel;
	processName = right.processName;
	pdgIdTrack = right.pdgIdTrack;
	trackVolumeName = right.trackVolumeName;
	parentVolumeName = right.parentVolumeName;
	//kinE      = right.kinE;
	return *this;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int AllPixTrackerHit::operator==(const AllPixTrackerHit& right) const
		{
	return (this==&right) ? 1 : 0;
		}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void AllPixTrackerHit::Draw()
{
	G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
	if(pVVisManager)
	{
		G4Circle circle(pos);
		circle.SetScreenSize(2.);
		circle.SetFillStyle(G4Circle::filled);
		//G4Colour colour(1.,0.,0.);
		G4VisAttributes attribs(G4Colour::Green());
		circle.SetVisAttributes(attribs);
		pVVisManager->Draw(circle);
	}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void AllPixTrackerHit::Print()
{
	G4cout << "  trackID: " << trackID << "  pixelNbX: " << pixelNbX
			<< "  pixelNbY: " << pixelNbY
			<< "  energy deposit: " << G4BestUnit(edep,"Energy")
			<< "  position: " << G4BestUnit(pos,"Length") << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

