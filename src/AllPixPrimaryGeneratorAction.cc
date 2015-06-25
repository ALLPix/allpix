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
// $Id: AllPixPrimaryGeneratorAction.cc,v 1.6 2006/06/29 17:54:27 gunter Exp $
// GEANT4 tag $Name: geant4-09-02 $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "AllPixPrimaryGeneratorAction.hh"
#include "AllPixPrimaryGeneratorMessenger.hh"

#include "Randomize.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4GeneralParticleSource.hh"
#include "G4HEPEvtInterface.hh"

#include "G4RunManager.hh"
#include "CLHEP/Units/SystemOfUnits.h"

using namespace CLHEP;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

AllPixPrimaryGeneratorAction::AllPixPrimaryGeneratorAction(SourceType st)
{
	m_sType = st;
	m_particleGun = 0x0;
	m_particleSource = 0x0;
	m_gunMessenger = 0x0;
	m_HEPEvt = 0x0;

	if(m_sType == _ParticleGun)
	{
		G4int n_particle = 1;
		m_particleGun = new G4ParticleGun(n_particle);

		//create a messenger for this class
		m_gunMessenger = new AllPixPrimaryGeneratorMessenger(this);

		//default kinematic
		//
		G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
		G4ParticleDefinition* particle = particleTable->FindParticle("e+");

		m_particleGun->SetParticleDefinition(particle);
		m_particleGun->SetParticleTime(0.0*ns);
		m_particleGun->SetParticlePosition(G4ThreeVector(-55.*um, 22.5*um, 22.5*um));
		//m_particleGun->SetParticlePosition(G4ThreeVector(0.150*mm, 0.150*mm, 0.*mm));
		//m_particleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,-1.));
		m_particleGun->SetParticleMomentumDirection(G4ThreeVector(1.,0.,0.));
		m_particleGun->SetParticleEnergy(500.0*keV);
	}
	else if(m_sType == _GeneralParticleSource)
	{
		/* Using G4GeneralParticleSource */
		m_gunMessenger = new AllPixPrimaryGeneratorMessenger(this);
		m_particleSource = new G4GeneralParticleSource();
	}
	else if(m_sType == _HEPEvtInterface)
	{
		/* using external MC data */
		m_gunMessenger = new AllPixPrimaryGeneratorMessenger(this);
		m_HEPEvt = new G4HEPEvtInterface(G4String("/afs/cern.ch/user/m/mbenoit/scratch0/Full_Tracker_Model/HEPEVT_files/AllPairs3TeV8MeV6Deg.HEPEvt"));
	}

	// store temporarily the position of incoming particles
	m_primaryParticlePos.clear();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

AllPixPrimaryGeneratorAction::~AllPixPrimaryGeneratorAction()
{
	if(m_particleGun) delete m_particleGun;
	if(m_particleSource) delete m_particleSource;
	if(m_gunMessenger) delete m_gunMessenger;
	if(m_HEPEvt) delete m_HEPEvt;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void AllPixPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{

	if(m_particleGun) m_particleGun->GeneratePrimaryVertex(anEvent);

	if(m_particleSource) {
		m_particleSource->GeneratePrimaryVertex(anEvent);
		if(anEvent->GetEventID() == 0){ m_primaryParticlePos.clear(); } // clear if at first event
		G4ThreeVector pos = m_particleSource->GetParticlePosition();
		m_primaryParticlePos.push_back(pos);
	}

	if(m_HEPEvt){
		m_HEPEvt->SetParticlePosition(G4ThreeVector(0.*mm,0.*mm,0.*mm));
		m_HEPEvt->GeneratePrimaryVertex(anEvent);
	}

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void AllPixPrimaryGeneratorAction::SetOptPhotonPolar()
{
	G4double angle = G4UniformRand() * 360.0*deg;
	SetOptPhotonPolar(angle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void AllPixPrimaryGeneratorAction::SetHEPevtFile(G4String filename)
{
	delete m_HEPEvt; 
	m_HEPEvt= new G4HEPEvtInterface(filename);
	G4cout << "Charging Event file " << filename << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void AllPixPrimaryGeneratorAction::SetOptPhotonPolar(G4double angle)
{
	if (m_particleGun->GetParticleDefinition()->GetParticleName() != "opticalphoton")
	{
		G4cout << "--> warning from PrimaryGeneratorAction::SetOptPhotonPolar() :"
				"the m_particleGun is not an opticalphoton" << G4endl;
		return;
	}

	G4ThreeVector normal (1., 0., 0.);
	G4ThreeVector kphoton = m_particleGun->GetParticleMomentumDirection();
	G4ThreeVector product = normal.cross(kphoton);
	G4double modul2       = product*product;

	G4ThreeVector e_perpend (0., 0., 1.);
	if (modul2 > 0.) e_perpend = (1./std::sqrt(modul2))*product;
	G4ThreeVector e_paralle    = e_perpend.cross(kphoton);

	G4ThreeVector polar = std::cos(angle)*e_paralle + std::sin(angle)*e_perpend;
	m_particleGun->SetParticlePolarization(polar);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
