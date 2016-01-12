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
// $Id: G4EmStandardPhysics.cc 84662 2014-10-17 14:32:32Z gcosmo $
//
//---------------------------------------------------------------------------
//
// ClassName:   G4EmCompPhotoEPhysics
//
// Author:      V.Ivanchenko 09.11.2005
//
// Modified:
// 05.12.2005 V.Ivanchenko add controlled verbosity
// 13.11.2006 V.Ivanchenko use G4hMultipleScattering
// 23.11.2006 V.Ivanchenko remove mscStepLimit option and improve cout
// 13.02.2007 V.Ivanchenko use G4hMultipleScattering for muons
// 13.02.2007 V.Ivanchenko set skin=0.0
// 21.04.2008 V.Ivanchenko add long-lived D and B mesons
//
//----------------------------------------------------------------------------
//

#include "G4EmCompPhotoEPhysics.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleDefinition.hh"
#include "G4EmParameters.hh"
#include "G4LossTableManager.hh"

#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"
#include "G4PhotoElectricEffect.hh"

#include "G4eMultipleScattering.hh"
//#include "G4MuMultipleScattering.hh"
#include "G4hMultipleScattering.hh"
#include "G4CoulombScattering.hh"
#include "G4eCoulombScatteringModel.hh"
#include "G4WentzelVIModel.hh"
#include "G4UrbanMscModel.hh"

//#include "G4MuBremsstrahlungModel.hh"
//#include "G4MuPairProductionModel.hh"
//#include "G4hBremsstrahlungModel.hh"
//#include "G4hPairProductionModel.hh"

#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4eplusAnnihilation.hh"
#include "G4UAtomicDeexcitation.hh"

#include "G4hIonisation.hh"
//#include "G4ionIonisation.hh"
//#include "G4alphaIonisation.hh"

#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
//#include "G4MuonPlus.hh"
//#include "G4MuonMinus.hh"
//#include "G4PionPlus.hh"
//#include "G4PionMinus.hh"
//#include "G4KaonPlus.hh"
//#include "G4KaonMinus.hh"
//#include "G4Proton.hh"
//#include "G4AntiProton.hh"
//#include "G4Deuteron.hh"
//#include "G4Triton.hh"
//#include "G4He3.hh"
//#include "G4Alpha.hh"
//#include "G4GenericIon.hh"

#include "G4PhysicsListHelper.hh"
#include "G4BuilderType.hh"

// factory
#include "G4PhysicsConstructorFactory.hh"
//
//G4_DECLARE_PHYSCONSTR_FACTORY(G4EmStandardPhysics);

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4EmCompPhotoEPhysics::G4EmCompPhotoEPhysics(G4int ver)
  : G4VPhysicsConstructor("G4EmStandard"), verbose(ver)
{
  G4EmParameters::Instance()->SetVerbose(verbose);
  SetPhysicsType(bElectromagnetic);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4EmCompPhotoEPhysics::G4EmCompPhotoEPhysics(G4int ver, const G4String&)
  : G4VPhysicsConstructor("G4EmStandard"), verbose(ver)
{
  G4EmParameters::Instance()->SetVerbose(verbose);
  SetPhysicsType(bElectromagnetic);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4EmCompPhotoEPhysics::~G4EmCompPhotoEPhysics()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmCompPhotoEPhysics::ConstructParticle()
{
  // gamma
  G4Gamma::Gamma();

  // leptons
  G4Electron::Electron();
  G4Positron::Positron();
  //G4MuonPlus::MuonPlus();
  //G4MuonMinus::MuonMinus();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmCompPhotoEPhysics::ConstructProcess()
{
  if(verbose > 1) {
    G4cout << "### " << GetPhysicsName() << " Construct Processes " << G4endl;
  }
  G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();

  // high energy limit for e+- scattering models
  G4double highEnergyLimit = 100*MeV;

  // Add standard EM Processes
  aParticleIterator->reset();
  while( (*aParticleIterator)() ){
    G4ParticleDefinition* particle = aParticleIterator->value();
    G4String particleName = particle->GetParticleName();

    ///////////
    // gamma //
    ///////////

    if (particleName == "gamma") {

      ph->RegisterProcess(new G4PhotoElectricEffect(), particle);
      ph->RegisterProcess(new G4ComptonScattering(), particle);
      ph->RegisterProcess(new G4GammaConversion(), particle);

    //////////////
    // electron //
    //////////////

    } else if (particleName == "e-") {

      G4eMultipleScattering* msc = new G4eMultipleScattering;
      G4UrbanMscModel* msc1 = new G4UrbanMscModel();
      G4WentzelVIModel* msc2 = new G4WentzelVIModel();
      msc1->SetHighEnergyLimit(highEnergyLimit);
      msc2->SetLowEnergyLimit(highEnergyLimit);
      msc->AddEmModel(0, msc1);
      msc->AddEmModel(0, msc2);

      G4eCoulombScatteringModel* ssm = new G4eCoulombScatteringModel();
      G4CoulombScattering* ss = new G4CoulombScattering();
      ss->SetEmModel(ssm, 1);
      ss->SetMinKinEnergy(highEnergyLimit);
      ssm->SetLowEnergyLimit(highEnergyLimit);
      ssm->SetActivationLowEnergyLimit(highEnergyLimit);

      ph->RegisterProcess(msc, particle);
      ph->RegisterProcess(new G4eIonisation(), particle);
      ph->RegisterProcess(new G4eBremsstrahlung(), particle);
      ph->RegisterProcess(ss, particle);

    //////////////
    // positron //
    //////////////

    } else if (particleName == "e+") {

      G4eMultipleScattering* msc = new G4eMultipleScattering;
      G4UrbanMscModel* msc1 = new G4UrbanMscModel();
      G4WentzelVIModel* msc2 = new G4WentzelVIModel();
      msc1->SetHighEnergyLimit(highEnergyLimit);
      msc2->SetLowEnergyLimit(highEnergyLimit);
      msc->AddEmModel(0, msc1);
      msc->AddEmModel(0, msc2);

      G4eCoulombScatteringModel* ssm = new G4eCoulombScatteringModel();
      G4CoulombScattering* ss = new G4CoulombScattering();
      ss->SetEmModel(ssm, 1);
      ss->SetMinKinEnergy(highEnergyLimit);
      ssm->SetLowEnergyLimit(highEnergyLimit);
      ssm->SetActivationLowEnergyLimit(highEnergyLimit);

      ph->RegisterProcess(msc, particle);
      ph->RegisterProcess(new G4eIonisation(), particle);
      ph->RegisterProcess(new G4eBremsstrahlung(), particle);
      ph->RegisterProcess(new G4eplusAnnihilation(), particle);
      ph->RegisterProcess(ss, particle);

    }
  }

  // Deexcitation
  //
  G4VAtomDeexcitation* de = new G4UAtomicDeexcitation();
  G4LossTableManager::Instance()->SetAtomDeexcitation(de);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
