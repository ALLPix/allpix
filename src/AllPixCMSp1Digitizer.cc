/**
 *  Author:
 *    Paul Schuetze <paul.schuetze@desy.de>
 *
 *  allpix Authors:
 *   John Idarraga <idarraga@cern.ch>
 *   Mathieu Benoit <benoit@lal.in2p3.fr>
 */

#include "AllPixCMSp1Digitizer.hh"
#include "AllPixTrackerHit.hh"
#include "G4DigiManager.hh"
#include "G4VDigitizerModule.hh"
#include "AllPixGeoDsc.hh"

#include "CLHEP/Random/RandGauss.h"
#include "CLHEP/Random/RandFlat.h"

#include "TMath.h"
#include <ctime>



AllPixCMSp1Digitizer::AllPixCMSp1Digitizer(G4String modName, G4String hitsColName, G4String digitColName)
: AllPixDigitizerInterface (modName) {
	
	// Registration of digits collection name
	collectionName.push_back(digitColName);
	m_hitsColName.push_back(hitsColName);
	
	// threshold
	m_digitIn.thl = 0.;
	
	

	InitVariables();

	
	
	
	// Test some Propagations through the sensor.
	
	ofstream ofile;
	ofile.open("drifttimeAllpix.txt", std::ofstream::out | std::ofstream::app);
	
	G4ThreeVector pos = G4ThreeVector(130.*um, 130.*um, 0.5*um);
	G4ThreeVector startpos;
	G4double drifttime;
	G4bool trapped;
	
	for (size_t i = 0; i < 284; i++) {
		pos = G4ThreeVector(3130.*um, 3130.*um, (0.5+(double)i)*um);
		startpos = pos;
		Propagation(pos, drifttime, trapped);
		G4cout << "Endposition: " << pos << G4endl;
		G4cout << "Total drift time in ns: " << drifttime*1e9 << G4endl;
		
		ofile << startpos[2] << "\t" << drifttime*1e9 << std::endl;
	
	}
	
	ofile.close();

}

AllPixCMSp1Digitizer::~AllPixCMSp1Digitizer(){
	
}

void AllPixCMSp1Digitizer::InitVariables(){
	
	// gD contains info on sensor and surroundings (bfield, temperature, ...)
	gD = GetDetectorGeoDscPtr();
	
	bfield = gD->GetMagField();
	detectorThickness = gD->GetSensorZ();
	Temperature = gD->GetTemperature();
	flux = gD->GetFlux(); // 1/cm^2
	
	///////////////////////////////////////////////////
	// Silicon electron and hole transport constants //
	///////////////////////////////////////////////////
	
	//// Unit for charge in FEIX average e/h pair creation energy in Silicon
	elec = 3.64*eV;
	
	// Variables for charge drift
	Electron_Mobility = 1.53e9*TMath::Power(Temperature, -0.87)/(1.01*TMath::Power(Temperature, 1.55))*1e-4; // mu0 from pixelav, m2/volt/s
	Electron_Beta = 0.0257*TMath::Power(Temperature, 0.66); // beta from pixelav
	
	Electron_HallFactor = 1.12;
	Electron_ec = 100*1.01 * TMath::Power(Temperature, 1.55); // ec from pixelav
	
	Boltzmann_kT = 8.6173e-5*Temperature; // eV
	// Boltzmann_kT = 1.38e-23*Temperature; // J
	
	Target_Spatial_Precision = 1e-10;
	Timestep_max = 0.1e-9;
	Timestep_min = 0.005e-9;
	
	// Variables for Trapping
	Electron_Trap_beta0 = 5.65e-7; // cm2/s
	Electron_Trap_kappa = -0.86;
	Electron_Trap_T0 = 263; // K
	Electron_Trap_TauNoFluence = 1.;
	
	if(flux > 0.){
		Electron_Trap_TauEff = 1./(flux * Electron_Trap_beta0 * TMath::Power((Temperature/Electron_Trap_T0),Electron_Trap_kappa));
	}else{
		Electron_Trap_TauEff = Electron_Trap_TauNoFluence;
	}
	
}

void AllPixCMSp1Digitizer::Digitize(){
	
	// create the digits collection
	m_digitsCollection = new AllPixCMSp1DigitsCollection("AllPixCMSp1Digitizer", collectionName[0] );
	
	// get the digiManager
	G4DigiManager * digiMan = G4DigiManager::GetDMpointer();
	
	// Get the hit collection ID
	G4int hcID = digiMan->GetHitsCollectionID(m_hitsColName[0]);
	
	// And fetch the Hits Collection
	AllPixTrackerHitsCollection * hitsCollection = 0;
	hitsCollection = (AllPixTrackerHitsCollection*)(digiMan->GetHitsCollection(hcID));

	// temporary data structure
	map<pair<G4int, G4int>, G4double > pixelsContent;
	pair<G4int, G4int> tempPixel;

	G4int nEntries = hitsCollection->entries();
	
	
	for(G4int itr  = 0 ; itr < nEntries ; itr++) {

		tempPixel.first  = (*hitsCollection)[itr]->GetPixelNbX();
		tempPixel.second = (*hitsCollection)[itr]->GetPixelNbY();
		pixelsContent[tempPixel] += (*hitsCollection)[itr]->GetEdep();

	}

	// Now create digits, one per pixel
	// Second entry in the map is the energy deposit in the pixel
	map<pair<G4int, G4int>, G4double >::iterator pCItr = pixelsContent.begin();

	// NOTE that there is a nice interface which provides useful info for hits.
	// For instance, the following member gives you the position of a hit with respect
	//  to the center of the pixel.
	// G4ThreeVector vec = (*hitsCollection)[itr]->GetPosWithRespectToPixel();
	// See class AllPixTrackerHit !

	// Also, you have full access to the detector geometry from this scope
	// AllPixGeoDsc * GetDetectorGeoDscPtr()
	// provides you with a pointer to the geometry description.
	// See class AllPixGeoDsc !

	for( ; pCItr != pixelsContent.end() ; pCItr++)
	{

		if((*pCItr).second > m_digitIn.thl) // over threshold !
		{
			// Create one digit per pixel, I need to look at all the pixels first
			AllPixCMSp1Digit * digit = new AllPixCMSp1Digit;
			digit->SetPixelIDX((*pCItr).first.first);
			digit->SetPixelIDY((*pCItr).first.second);
			digit->IncreasePixelCounts(); // Counting mode

			m_digitsCollection->insert(digit);
		}
	}

	G4int dc_entries = m_digitsCollection->entries();
	if(dc_entries > 0){
		G4cout << "--------> Digits Collection : " << collectionName[0]
		<< "(" << m_hitsColName[0] << ")"
		<< " contains " << dc_entries
		<< " digits" << G4endl;
	}

	StoreDigiCollection(m_digitsCollection);
	
}



vector<G4double>  AllPixCMSp1Digitizer::RKF5Integration(G4ThreeVector position, G4double dt)
{
	// This function transport using Euler integration, for field (Ex,Ey,Ez),
	// considered constant over time dt. The movement equation are those
	// of charges in semi-conductors, sx= mu*E*dt;;
	G4ThreeVector deltaposition;
	
	G4ThreeVector electricField;
	G4ThreeVector temppos;
	
	G4ThreeVector k1, k2, k3, k4, k5, k6;
	
	electricField = 100.*gD->GetEFieldFromMap(position);
	k1 = (ElectronSpeed(electricField)*dt);
	
	// G4cout << "z: " << position[2] << "\t\t" << ElectronSpeed(electricField) << G4endl;
	// G4cout << "z: " << position[2] << "\t\t" << electricField << G4endl;

	// G4cout << "Field: " << electricField << G4endl;
	// G4cout << "Speed: " << ElectronSpeed(electricField) << G4endl;
	
	temppos = position + (1./4.)*k1;
	electricField = 100.*gD->GetEFieldFromMap(temppos);
	k2 = (ElectronSpeed(electricField)*dt);
	
	temppos = position + (3./32.)*k1 + (9./32.)*k2;
	electricField = 100.*gD->GetEFieldFromMap(temppos);
	k3 = (ElectronSpeed(electricField)*dt);
	
	temppos = position + (1932./2197)*k1 + (-7200./2197)*k2 + (7296./2197)*k3;
	electricField = 100.*gD->GetEFieldFromMap(temppos);
	k4 = (ElectronSpeed(electricField)*dt);
	
	temppos = position + (439./216)*k1 + (-8)*k2 + (3680./513)*k3 + (-845./4104)*k4;
	electricField = 100.*gD->GetEFieldFromMap(temppos);
	k5 = (ElectronSpeed(electricField)*dt);
	
	temppos = position + (-8./27)*k1 + (2)*k2 + (-3544./2565)*k3 + (1859./4104)*k4 + (-11./40)*k5;
	electricField = 100.*gD->GetEFieldFromMap(temppos);
	k6 = (ElectronSpeed(electricField)*dt);
	
	deltaposition = (16./135)*k1 + (0.)*k2 + (6656./12825)*k3 + (28561./56430)*k4 + (-9./50)*k5 + (2./55)*k6;
	
	G4ThreeVector error;
	error = (1./360)*k1 + (0.)*k2 + (-128./4275)*k3 + (-2197./75240)*k4 + (1./50)*k5 + (2./55)*k6;
	
	G4double Erreur;
	Erreur=error.mag();
	
	vector<G4double> deltapoint(4);
	deltapoint[0]=deltaposition[0];
	deltapoint[1]=deltaposition[1];
	deltapoint[2]=deltaposition[2];
	deltapoint[3]=(Erreur);
	
	return deltapoint;
	
}

G4double AllPixCMSp1Digitizer::MobilityElectron(const G4ThreeVector efield){
	
	// calculate mobility in m2/V/s
	
	G4double mobility = Electron_Mobility * TMath::Power((1.+ TMath::Power(efield.mag()/Electron_ec,Electron_Beta)),-1.0/Electron_Beta);
	// G4cout << "Mobility: " << mobility << G4endl;
	return mobility;

}

G4ThreeVector AllPixCMSp1Digitizer::ElectronSpeed(const G4ThreeVector efield){
	
	G4ThreeVector term[3];
	G4double mobility = MobilityElectron(efield);
	
	G4double rnorm = 1. + mobility*mobility*Electron_HallFactor*Electron_HallFactor*bfield.dot(bfield);
	
	term[0] = -efield;
	term[1] = mobility*Electron_HallFactor*efield.cross(bfield);
	term[2] = -mobility*mobility*Electron_HallFactor*Electron_HallFactor*efield.dot(bfield)*bfield;
	
	G4ThreeVector speed = mobility*(term[0]+term[1]+term[2])/rnorm;
	
	return speed;
	
}

G4ThreeVector AllPixCMSp1Digitizer::DiffusionStep(const G4double timestep, const G4ThreeVector position){
	
	G4ThreeVector diffusionVector;
	
	G4ThreeVector electricField = 100.*gD->GetEFieldFromMap(position);
	G4double D = Boltzmann_kT*MobilityElectron(electricField);
	G4double Dwidth = TMath::Sqrt(2.*D*timestep);
	
	for (size_t i = 0; i < 3; i++) {
		diffusionVector[i] = CLHEP::RandGauss::shoot(0.,Dwidth)*um;
	}
	
	// TMath::Sqrt(2.*D*timestep);
	return diffusionVector;
	
}

void AllPixCMSp1Digitizer::SetDt(G4double& dt, const G4double uncertainty, const G4double z, const G4double dz){
	
	G4double dt_init = dt;
	
	if(uncertainty > Target_Spatial_Precision){dt *= 0.7;}
	if(uncertainty < 0.5*Target_Spatial_Precision){dt *= 2;}
	
	if(dt > Timestep_max){dt = Timestep_max;}
	if(dt < Timestep_min){dt = Timestep_min;}
	
	if(detectorThickness-z < dz*1.2){dt = dt_init*0.7;}
	
}

G4double AllPixCMSp1Digitizer::GetTrappingTime(){
	
	return(-Electron_Trap_TauEff*log(CLHEP::RandFlat::shoot()));
	
}

/*
	This function propagates an electron through the sensor and updates the position vector.
*/

G4double AllPixCMSp1Digitizer::Propagation(G4ThreeVector& pos, G4double& drifttime, G4bool& trapped){
	
	vector<G4double> deltapoint(4);
	
	drifttime = 0.;
	G4double dt = 0.01*1e-9;
	
	G4double trappingTime = GetTrappingTime();
	
	trapped = false;
	
	G4int nsteps=0;

	while(pos[2] > 0 && pos[2] < detectorThickness)
	{
		
		if(drifttime > trappingTime){
			trapped = true;
			break;
		}
		
		deltapoint = RKF5Integration(pos,dt);
		for (size_t i = 0; i < 3; i++) {pos[i]+=deltapoint[i]/nm*um;}
		drifttime += dt;


		pos += DiffusionStep(dt, pos);


		// Adapt step size 
		SetDt(dt, deltapoint[3], pos[2], deltapoint[2]/nm*um);

		nsteps++;
	}
	if(trapped) G4cout << "Charge was trapped." << G4endl;
	// G4cout << "Endposition: " << pos << G4endl;
	// G4cout << "Total drift time in ns: " << drifttime*1e9 << G4endl;
	// G4cout << "Step count: " << nsteps << G4endl;
	
	return drifttime;
	
}
