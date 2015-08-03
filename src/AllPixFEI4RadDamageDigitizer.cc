/**
 *  Authors:
 *    Callie Bertsche <c.bertsche@cern.ch>
 *    Benjamin Nachman <bnachman@cern.ch>
 *
 *  allpix Authors:
 *   John Idarraga <idarraga@cern.ch>
 *   Mathieu Benoit <benoit@lal.in2p3.fr>
 */

#include "AllPixFEI4RadDamageDigitizer.hh"
#include "AllPixTrackerHit.hh"
#include "G4DigiManager.hh"
#include "G4VDigitizerModule.hh"
#include "AllPixGeoDsc.hh"

#include "TMath.h"
#include "TF1.h"

// Included for radiation damage
#include "math.h"
#include "TFile.h"
#include "CLHEP/Random/RandGauss.h"
#include "CLHEP/Random/RandFlat.h"

using namespace TMath;

AllPixFEI4RadDamageDigitizer::AllPixFEI4RadDamageDigitizer(G4String modName, G4String hitsColName, G4String digitColName) 
: AllPixDigitizerInterface (modName) {

	// Registration of digits collection name
	collectionName.push_back(digitColName);
	m_hitsColName.push_back(hitsColName);

	// Threshold
	m_digitIn.thl = 0.;

	// Get ramo potential mapping
	ramoPotentialMap=0;
	TFile* file=new TFile("./share/absRamo3D-map-200um-output.root");
        ramoPotentialMap=(TH3F*)file->Get("hramomap1");
	if (ramoPotentialMap == 0){
		G4cout << "Unsuccessful picking up histogram: ramoPotentialMap" << G4endl;
	}

	// Get electric field mapping
	eFieldMap=0;
	TFile* efile=new TFile("./share/eField-fei4-200um-fl0-100V-output.root");
        eFieldMap=(TH1F*)efile->Get("hefieldz");
	if (eFieldMap == 0){
		G4cout << "Unsuccessful picking up histogram: eFieldMap" << G4endl;
	}

	// Get time to electrode mapping (derived from electric field)
	timeMap_e=0;
	timeMap_h=0;
	TFile* tfile=new TFile("./share/timeZ-fei4-200um-fl0-100V-output.root");
	timeMap_e=(TH1F*)tfile->Get("etimes");
	timeMap_h=(TH1F*)tfile->Get("htimes");
	if (timeMap_e == 0 || timeMap_h == 0){
		G4cout << "Unsuccessful picking up histogram: timeMap" << G4endl;
	}

	// Unit for charge in FEIX average e/h pair creation energy in Silicon
	elec = 3.64*eV;

	// Fetching info from pixeldetector.xml
	AllPixGeoDsc * gD = GetDetectorGeoDscPtr();

	//////////////////////////
	// Bias and Temperature //
	//////////////////////////

	// Customize for different voltages used
	biasVoltage = 100.0; 			// V
	temperature = 263.2;			// K
	detectorThickness = gD->GetSensorZ();
	resistivity=gD->GetResistivity();

	bulkType =false;
	// true = p-type
	// false = n-type


	////////////////////////////////
	// Geometry Related constants //
	////////////////////////////////
	pitchX = gD->GetPixelX(); 	// total length of pixel in x
	pitchY = gD->GetPixelY();	// total length of pixel in y
	nPixX= gD->GetNPixelsX(); 	// total number of pixels in x
	nPixY= gD->GetNPixelsY();	// total number of pixels in y

	//G4cout << TString::Format("Lx=%f Ly=%f npixX=%d npixY=%d",pitchX,pitchY,nPixX,nPixY) << endl;

	//////////////////////
	// Radiation damage //
	//////////////////////

	// Customize for different fluence/voltage combinations

//	fluence = 3.8e15;				// neq/cm^2
	fluence = 0;					// neq/cm^2
	betaElectrons = 5.1E-16*cm2/ns;			// cm^2/ns
	betaHoles = 7.1E-16*cm2/ns;			// cm^2/ns
	bField = 1.6; 				// Tesla = V*s/m^2

	if(fluence!=0.0)
	{
		trappingTimeElectrons = 1.0/(betaElectrons*fluence);
		trappingTimeHoles = 1.0/(betaHoles*fluence);
	}
	else {
		trappingTimeElectrons = 1000*s;
		trappingTimeHoles = 1000*s;
	}

	echarge = 1.60217646e-19; // Coulombs
	epsilon = 11.8*8.854187817e-12/m;


	//////////////////////
 	// physics switches //
	//////////////////////

 	doTrapping = true;
	doRamo = true;
  	
 	precision = 250;

 	/////////////////////////////
 	// Detector Chip Selection //
 	/////////////////////////////

 	FEIX = 4;

 	switch (FEIX) {

 	case 3 : // FEI3 (default)
 		MipTOT=60;
 		MipCharge=19200;
 		CounterDepth=255;
 		Lv1Unit = 25*ns;
 		chipNoise = 300;
 		m_digitIn.thl = 3500;
 		chargeSharingConstant = 0.0;
 		break;

 	case 4 : // FEI4
 		MipTOT=10;
 		MipCharge=40000;
 		//MipCharge=20000;
 		CounterDepth=15;
 		Lv1Unit = 25*ns;
 		chipNoise = 125*elec;
 		m_digitIn.thl = 3200*elec;
 		chargeSharingConstant = 0.0;
 		break;

 	case 666 : // Omegapix2
 		MipTOT=5;
 		MipCharge=12000;
 		CounterDepth=32;
 		Lv1Unit = 25*ns;
 		chipNoise = 20*elec;
 		m_digitIn.thl = 800*elec;
 		chargeSharingConstant =0.01;
 		break;

 	default :
 		MipTOT=60;
 		MipCharge=22000;
 		CounterDepth=255;
 		Lv1Unit = 25*ns;
 		chipNoise = 300*elec;
 		m_digitIn.thl = 3500*elec;
 		chargeSharingConstant =0.01;

 	}


 	//////////////////////
 	// Sensor Selection //
 	//////////////////////
 	/*
 	 1=FEI3 (DEFAULT)
 	 2,3= FEI4, SlimEdge or Conservative
 	 4=Omegapix2
 	 5= FEI3 Slim edge
 	 */
 	Sensor = 2;

 	switch (Sensor) {

 	case 1 :
 	 	doSlimEdge=false;
 		GRShift =250*um;
 		break;
 	case 2 :
 	 	doSlimEdge=true;
 		GRShift =250*um;
 		break;
 	case 3 :
 	 	doSlimEdge=true;
 		GRShift =500*um;
 		break;
 	case 4 :
 	 	doSlimEdge=true;
 		GRShift =500*um;
 		break;
 	case 5 :
 	 	doSlimEdge=true;
 		GRShift =250*um;
 		break;
 	default :
 	 	doSlimEdge=false;
 		GRShift =0*um;
 	}

}// End of AllPixFEI4RadDamageDigitizer::AllPixFEI4RadDamageDigitizer definition

AllPixFEI4RadDamageDigitizer::~AllPixFEI4RadDamageDigitizer(){
}

void AllPixFEI4RadDamageDigitizer::SetDetectorDigitInputs(G4double thl){

	// set digitization input values
	// thl
	m_digitIn.thl = thl; // <-- input !
}

///////////////////////////////////
// Radiation Damage Calculations //
///////////////////////////////////

// Look up E field from TCAD mapping
G4double AllPixFEI4RadDamageDigitizer::GetElectricField(G4double z){
	// The z position is in mm, but plot uses um to return electric field in V/cm. Convert to return electric field in MV/mm
	// Is currently using z as distance from position to readout side
	int n_binz = eFieldMap->GetXaxis()->FindBin(z*1000);						
	G4double electricField = eFieldMap->GetBinContent(n_binz);	
	return electricField*1.0E-7;
}

G4double AllPixFEI4RadDamageDigitizer::GetMobility(G4double electricField, G4double temperature, G4bool isHole){ 
	
	// Initialize variables so they have the right scope
	G4double vsat = 0;
	G4double ecrit = 0;
	G4double beta = 0;	

	if(!isHole){
		vsat = 15.3*pow(temperature,-0.87);	// mm/ns
		ecrit = 1.01E-7*pow(temperature,1.55);	// MV/ns
		beta = 2.57E-2*pow(temperature,0.66);
	}
	if(isHole){
		vsat = 1.62*pow(temperature,-0.52);	// mm/ns
		ecrit = 1.24E-7*pow(temperature,1.68);	// MV/ns
		beta = 0.46*pow(temperature,0.17);
	}

	G4double mobility = (vsat/ecrit)/pow(1+pow((electricField/ecrit),beta),(1/beta));	
	return mobility;  					// mm^2/(MV*ns)
}

G4double AllPixFEI4RadDamageDigitizer::GetDriftVelocity(G4double electricField, G4double mobility, G4bool isHole){
	G4double driftVelocity = mobility*electricField;// mm/ns
	if(isHole) driftVelocity = -1*driftVelocity; // Drifts in opposite direction
	return driftVelocity;
}

G4double AllPixFEI4RadDamageDigitizer::GetMeanFreePath(G4double driftVelocity, G4bool isHole){

	G4double meanFreePath = 0;	
	if(!isHole) meanFreePath = driftVelocity*trappingTimeElectrons; // mm
	if(isHole) meanFreePath = driftVelocity*trappingTimeHoles; // mm
	return meanFreePath;
}

G4double AllPixFEI4RadDamageDigitizer::GetTrappingProbability(G4double z, G4double meanFreePath){
	if(isHole) z = 0.2 - z; // For holes, should consider opposite direction
	G4double trappingProbability = 1.0 - TMath::Exp(-TMath::Abs(z/(meanFreePath))); 
	return trappingProbability;
}

G4double AllPixFEI4RadDamageDigitizer::GetDriftTime(G4bool isHole){
	G4double u = CLHEP::RandFlat::shoot(0.,1.); // 
	G4double driftTime = 0;

	if(!isHole) driftTime = (-1.)*trappingTimeElectrons*TMath::Log(u); // ns
	if(isHole) driftTime = (-1.)*trappingTimeHoles*TMath::Log(u); // ns
	return driftTime;
}

G4double AllPixFEI4RadDamageDigitizer::GetTimeToElectrode(G4double z, G4bool isHole){
	// Uses z position in mm to return time to electrode, in ns
	// The mapping takes care of holes going in opposite direction
	G4double timeToElectrode = 0;
	if(!isHole) 
	{
		int n_binz = timeMap_e->GetZaxis()->FindBin(0.2 - z);	// using z as distance to readout side						
		timeToElectrode = timeMap_e->GetBinContent(n_binz);	
	}
	if(isHole)
	{
		int n_binz = timeMap_h->GetZaxis()->FindBin(0.2 - z);							
		timeToElectrode = timeMap_h->GetBinContent(n_binz);			
	}
	return timeToElectrode;	
}


// ***** Process for each event *****
void AllPixFEI4RadDamageDigitizer::Digitize(){
	
	// create the digits collection
	m_digitsCollection = new AllPixFEI4RadDamageDigitsCollection("AllPixFEI4RadDamageDigitizer", collectionName[0] );

	// get the digiManager
	G4DigiManager * digiMan = G4DigiManager::GetDMpointer();

	// BoxSD_0_HitsCollection
	G4int hcID = digiMan->GetHitsCollectionID(m_hitsColName[0]);

	AllPixTrackerHitsCollection * hitsCollection = 0;
	hitsCollection = (AllPixTrackerHitsCollection*)(digiMan->GetHitsCollection(hcID));

	// temporary data structure
	map<pair<G4int, G4int>, G4double > pixelsContent;		// stored energy per (countx,county) pixel	
	pair<G4int, G4int> tempPixel;					// (countx,county) which pixel


	// Loop over the whole Hits Collection - this spreads out the effect of the hit over its path through the pixel 
	G4int nEntries = hitsCollection->entries();
	for(G4int itr  = 0 ; itr < nEntries ; itr++) {

		G4double xpos = (*hitsCollection)[itr]->GetPosWithRespectToPixel().x(); // mm
		G4double ypos = (*hitsCollection)[itr]->GetPosWithRespectToPixel().y(); // mm

		// Given as distance from center of pixel, + toward electrode (where electrons are read out)
		G4double zpos = detectorThickness/2. - (*hitsCollection)[itr]->GetPosWithRespectToPixel().z();	// mm; Gives distance to readout side for electron

		G4double eHitTotal = (*hitsCollection)[itr]->GetEdep(); // Energy deposition for the hit, internal unit
		tempPixel.first  = (*hitsCollection)[itr]->GetPixelNbX();
		tempPixel.second = (*hitsCollection)[itr]->GetPixelNbY();

		//G4cout << "x : " << tempPixel.first << " ,  y : " << tempPixel.second << ", E = " << eHitTotal/keV << G4endl;

		// In case the charge moves into a neighboring pixel
		pair<G4int, G4int> extraPixel;
		extraPixel = tempPixel;

		// Split the charge into subcharges (# = precision) that are diffused separately to the electrode
		for(G4int nQ  = 0 ; nQ < precision ; nQ++) {
				
			G4double eHit = G4double(eHitTotal)/(2*precision); // eV; divide in half because we are treating holes separately	

			// Not sure if we should do slim edge, but we are using FEI4 sensors
			if(doSlimEdge){
				if (isSlimEdge(tempPixel.first,tempPixel.second)) eHit = SlimEdgeEffect(tempPixel.first,xpos,eHit);
			}

			// Loop over everything following twice, once for holes and once for electrons
			for(G4int eholes=0 ; eholes<2 ; eholes++) { // Loop over everything twice, once for electrons and once for holes

				isHole = false; // Set a condition to keep track of electron/hole-specific functions
				if (eholes == 1) isHole = true;

				// Don't use the event if electric field is zero, unless it is zero because charge is at electrode, in which case record it
				G4double electricField = GetElectricField(zpos);
				if (electricField == 0) {
					if (zpos < 1){
						pixelsContent[extraPixel] += eHit; // eV
					}
					continue;	// flushes current event, continues with next event
				}

				// Reset extraPixel coordinates each time through loop
				extraPixel = tempPixel;
	
				G4double mobility = GetMobility(electricField, temperature, isHole);
				G4double driftVelocity = GetDriftVelocity(electricField, mobility, isHole);
				G4double meanFreePath = GetMeanFreePath(driftVelocity, isHole);
				G4double trappingProbability = GetTrappingProbability(zpos, meanFreePath);
				G4double timeToElectrode = GetTimeToElectrode(zpos, isHole);
				G4double driftTime = GetDriftTime(isHole);

				// Set up variables for diffusion, per current pixel code process (I have not seen the origin of these equations)
				G4double hallEffect = 1.13 + 0.0008*(temperature - 273.0);
				G4double tanLorentz = hallEffect*mobility*bField*(1.0E-3);	// unit conversion 
				G4double coLorentz = sqrt(1+pow(tanLorentz,2));			// from pre-rad damage code: acode-browser2.usatlas.bnl.gov/lxr-rel17

				/////////////////
				// Now process the event according to whether trapping and/or ramo are included
				/////////////////

				// Trapping and ramo potential turned off
				if (!doTrapping && !doRamo){

						G4double rdif=0.007*sqrt(zpos*coLorentz/0.3); 
						// Using RandGauss for now as RandGaussZiggurat is not accessible
						G4double xposD=xpos+zpos*tanLorentz+rdif*CLHEP::RandGauss::shoot(0,1); // Is it still +zpos in the case of the holes?
						G4double yposD=ypos+rdif*CLHEP::RandGauss::shoot(0,1); 

						// Account for drifting into another pixel
						while (fabs(xposD) >= pitchX/2){
							G4double sign = xposD/(fabs(xposD)); 			// returns +1 or -1 depending on + or - x value
							extraPixel.first = extraPixel.first + 1*sign;		// increments or decrements pixel count in x
							xposD = xposD - (pitchX*sign);				// moves xpos coordinate 1 pixel over in x
						}
						while (fabs(yposD) >= pitchY/2){
							G4double sign = yposD/(fabs(yposD)); 			// returns +1 or -1 depending on + or - y value
							extraPixel.second = extraPixel.second + 1*sign;		// increments or decrements pixel count in y
							yposD = yposD - (pitchY*sign);				// moves xpos coordinate 1 pixel over in y
						}

						// Record deposit
						pixelsContent[extraPixel] += eHit; // eV
				}

				// Only trapping
				if (doTrapping && !doRamo){

					if (driftTime >= timeToElectrode){ 	// charge wasn't trapped - record energy deposit at electrode

						// Diffusion
						G4double rdif=0.007*sqrt(zpos*coLorentz/0.3); 
						G4double xposD=xpos+zpos*tanLorentz+rdif*CLHEP::RandGauss::shoot(0,1); // Is it still +zpos in the case of the holes?
						G4double yposD=ypos+rdif*CLHEP::RandGauss::shoot(0,1); 

						// Account for drifting into another pixel
						while (fabs(xposD) >= pitchX/2){
							G4double sign = xposD/(fabs(xposD)); 			// returns +1 or -1 depending on + or - x value
							extraPixel.first = extraPixel.first + 1*sign;		// increments or decrements pixel count in x
							xposD = xposD - (pitchX*sign);				// moves xpos coordinate 1 pixel over in x
						}
						while (fabs(yposD) >= pitchY/2){
							G4double sign = yposD/(fabs(yposD)); 			// returns +1 or -1 depending on + or - y value
							extraPixel.second = extraPixel.second + 1*sign;		// increments or decrements pixel count in y
							yposD = yposD - (pitchY*sign);				// moves xpos coordinate 1 pixel over in y
						}

						// Record deposit
						pixelsContent[extraPixel] += eHit; // eV

					} // end if charge not trapped
					
					// else, charge is trapped - no energy deposit, do nothing (since no ramo)

				} // end (doTrapping && !doRamo)

				// Trapping and Ramo Potential
				if (doTrapping && doRamo){

					if (driftTime >= timeToElectrode){	// charge wasn't trapped - record energy deposit at electrode

						// Diffusion
						G4double rdif=0.007*sqrt(zpos*coLorentz/0.3); 
						G4double xposD=xpos+zpos*tanLorentz+rdif*CLHEP::RandGauss::shoot(0,1); // Is it still +zpos in the case of the holes?
						G4double yposD=ypos+rdif*CLHEP::RandGauss::shoot(0,1); 

						// Account for drifting into another pixel
						while (fabs(xposD) >= pitchX/2){
							G4double sign = xposD/(fabs(xposD)); 			// returns +1 or -1 depending on + or - x value
							extraPixel.first = extraPixel.first + 1*sign;		// increments or decrements pixel count in x
							xposD = xposD - (pitchX*sign);				// moves xpos coordinate 1 pixel over in x
						}
						while (fabs(yposD) >= pitchY/2){
							G4double sign = yposD/(fabs(yposD)); 			// returns +1 or -1 depending on + or - y value
							extraPixel.second = extraPixel.second + 1*sign;		// increments or decrements pixel count in y
							yposD = yposD - (pitchY*sign);				// moves xpos coordinate 1 pixel over in y
						}

						// Record deposit
						pixelsContent[extraPixel] += eHit; // eV
						//G4cout << "Recorded hit: " << pixelsContent[extraPixel] << G4endl;

						// Also record deposit due to diff in ramo potential between (xposD, yposD, electrode) and (xpos, ypos, zpos)
						// Initial ramo potential based on (x,y,z) position in micrometers
						int nbin = 0;
						if(!isHole) nbin = ramoPotentialMap->FindBin(fabs(ypos*1000),fabs(xpos*1000),zpos*1000);
						if(isHole) nbin = ramoPotentialMap->FindBin(fabs(ypos*1000),fabs(xpos*1000),200-zpos*1000);
						G4double ramo_i = ramoPotentialMap->GetBinContent(nbin);

						// ramo potential at electrode based on (x,y,z) position in micrometers						// -- loop in the x-coordinate
						for (int i=-1; i<=1; i++){
							G4double x_neighbor = xposD + i*pitchX;
							extraPixel.first += i;
							if (i == 0){
								extraPixel.first += 1; // For the middle neighbor, still have to add one to get back to zero
							}
							extraPixel.second = extraPixel.second - 1; // to start the y-pixel count in the middle pixel each time
	  
							// -- loop in the y-coordinate							
							for (int j=-1; j<=1; j++){
								G4double y_neighbor = yposD + j*pitchY;
								extraPixel.second += j;
								if (j == 0){
									extraPixel.second += 1; // For the middle neighbor, add one to get back to zero
								}								
								// Return ramo potential based on (x,y,z) position in micrometers; in z is at electrode
								int nbin = ramoPotentialMap->FindBin(fabs(y_neighbor*1000),fabs(x_neighbor*1000),0); 
								G4double ramo = ramoPotentialMap->GetBinContent(nbin);	

								// Record deposit
								G4double eHitRamo = eHit*(ramo - ramo_i); // eV
								//G4cout << "ramodiff: " << ramo-ramo_i << G4endl;
								//pixelsContent[extraPixel] += eHitRamo; // eV
							} // end loop over y-coordinate
						} // end loop over x-coordinate
					} // end trapped charge

					else if (driftTime < timeToElectrode){	// charge is trapped - just do ramo potential calculation

						// -- Calculate trapping position taking into account diffusion in the XY plane
						// rdif/coLorentz source: acode-browser2.usatlas.bnl.gov/lxr-rel17
						
						// Review diffusion_trap and rdif_trap calculations
						// diffusion_trap varies from ~0.005 to 0.034 depending on time spent drifting in the pixel...
						G4double mu_0 = 1923.*temperature*temperature; //cm^2*(V*s)^-1
						G4double diffusion_trap = sqrt(2.*mu_0*8.62E-5*driftTime*1.0E-9);
						G4double rdif_trap = diffusion_trap*sqrt(zpos*coLorentz/0.3); 
						
						// z position is distance to readout side									
						G4double z_trap = zpos - driftTime*driftVelocity; // estimate; integral would be better...
						G4double x_trap = xpos + z_trap*tanLorentz + rdif_trap*CLHEP::RandGauss::shoot(0,1); 
						G4double y_trap = ypos + rdif_trap*CLHEP::RandGauss::shoot(0,1); 
          	  
						// Account for drifting into another pixel
						while (fabs(x_trap) >= pitchX/2){
							G4double sign = x_trap/(fabs(x_trap)); 			// returns +1 or -1 depending on + or - x value
							extraPixel.first = extraPixel.first + 1*sign;		// increments or decrements pixel count in x
							x_trap = x_trap - (pitchX*sign);			// moves xpos coordinate 1 pixel over in x
							G4cout << "Drifted to another pixel in x" << G4endl;
						}
						while (fabs(y_trap) >= pitchY/2){
							G4double sign = y_trap/(fabs(y_trap)); 			// returns +1 or -1 depending on + or - y value
							extraPixel.second = extraPixel.second + 1*sign;		// increments or decrements pixel count in y
							y_trap = y_trap - (pitchY*sign);			// moves xpos coordinate 1 pixel over in y
							G4cout << "Drifted to another pixel in y" << G4endl;
						}

						// Record deposit due to ramo potential
						// Ramo of initial position (xpos,ypos,zpos) in micrometers. Note that x and y are reversed in the TCAD ramo mappings
						int nbin = 0;						
						if(!isHole) nbin = ramoPotentialMap->FindBin(fabs(ypos*1000),fabs(xpos*1000),zpos*1000);
						if(isHole) nbin = ramoPotentialMap->FindBin(fabs(ypos*1000),fabs(xpos*1000),200-zpos*1000);
 						G4double ramo_i = ramoPotentialMap->GetBinContent(nbin);
						if(ramo_i == 0){
							G4cout << "ramo_i = 0; xpos, ypos, zpos: " << fabs(xpos*1000) << ", " << fabs(ypos*1000) << ", " << zpos*1000 << endl;
						}

						// -- Calculate signal induced in pixel where trapped as well as the 8 pixel neighbors						// -- loop in the x-coordinate
						for (int i=-1; i<=1; i++){
							G4double x_neighbor = x_trap + i*pitchX;
							extraPixel.first += i;
							if (i == 0){
								extraPixel.first += 1; // For the middle neighbor, still have to add one to get back to zero
							}
	  						extraPixel.second = extraPixel.second - 1; // to start the y-pixel count in the middle pixel each time

							// -- loop in the y-coordinate	
							for (int j=-1; j<=1; j++){
								G4double y_neighbor = y_trap + j*pitchY;
								extraPixel.second += j;
								if (j == 0){
									extraPixel.second += 1; // For the middle neighbor, add one to get back to zero
								}
								// Return ramo potential based on (x,y,z) position in micrometers
								// Note that in ramo mapping, x and y are reversed!!!
								int nbin = ramoPotentialMap->FindBin(fabs(y_neighbor*1000),fabs(x_neighbor*1000),z_trap);
 								G4double ramo = ramoPotentialMap->GetBinContent(nbin);		
								if(ramo == 0){
									G4cout << "ramo = 0; x,y,z: " << y_neighbor*1000 << ", " << x_neighbor*1000 << ", " << z_trap*1000 << G4endl;
								}

								// Calculate charge contribution
								G4double eHitRamo = eHit*(ramo - ramo_i); // eV
								//G4cout << "ramodiff: " << ramo-ramo_i << G4endl;

								// Record deposit
								//pixelsContent[extraPixel] += eHitRamo; // eV
							} // end loop over y-coordinate
						} // end loop over x-coordinate
					} // end if charge trapped

				} // end (doTrapping && doRamo)

		} // end loop over nQ charges

		} // end loop over nEntries

	} // end loop over charges/holes
				
	// Now that pixelContent is filled, create one digit per pixel
	map<pair<G4int, G4int>, G4double >::iterator pCItr = pixelsContent.begin();
				
	  for( ; pCItr != pixelsContent.end() ; pCItr++)
	    {
	      G4double deposited_energy = (*pCItr).second;
	      double deposited_charge = deposited_energy/(3.64*eV); //1e = 3.64 eV  
	      //G4double threshold = 5.824*keV; // 1600e on twiki: https://twiki.cern.ch/twiki/bin/viewauth/Atlas/IBLTBJune2011 
	      G4double threshold = 2*5.824*keV; // 2*1600e on twiki: https://twiki.cern.ch/twiki/bin/viewauth/Atlas/IBLTBJune2011
	      double TOT = 1.2*deposited_charge/1000.; //This and the next number are from parameterizations of the tables in ATHENA
	      double TOTsig = 0.5+0.02*deposited_charge/1000.;
	      TOT = CLHEP::RandGauss::shoot(TOT, TOTsig);
	      int TOT_int = TMath::FloorNint(TOT);
	      AllPixFEI4RadDamageDigit * digit = new AllPixFEI4RadDamageDigit;
	      digit->SetPixelIDX((*pCItr).first.first);
	      digit->SetPixelIDY((*pCItr).first.second);
	      digit->SetPixelCounts(TOT_int);
	      if (deposited_energy < threshold){
		digit->SetPixelCounts(0);
	      }
	      if (deposited_energy >= threshold && TOT > 0){
		m_digitsCollection->insert(digit);
	      }
	    }

	G4int dc_entries = m_digitsCollection->entries();
	if(dc_entries > 0)
	{
		G4cout << "--------> Digits Collection : " << collectionName[0]
					 << "(" << m_hitsColName[0] << ")"
					<< " contains " << dc_entries
					 << " digits" << G4endl;
	}

	StoreDigiCollection(m_digitsCollection);


} // end Digitize function



G4int AllPixFEI4RadDamageDigitizer::EnergyToTOT(G4double Energy, G4double threshold)
{

	G4double dEdt = (MipCharge- threshold)/(MipTOT*Lv1Unit);
	//G4int TOT = TMath::FloorNint((Energy-threshold)/(dEdt*Lv1Unit*ns));
	G4int TOT = TMath::FloorNint((Energy-threshold)/(dEdt*Lv1Unit));

	if(Energy<threshold) TOT=0;
	if (TOT<0) TOT=0;
	if (TOT>CounterDepth) TOT=CounterDepth;

	return TOT;
}


//G4int AllPixFEI4RadDamageDigitizer::EnergyToTOTFEI4(G4double Energy, G4double threshold)
//{
//	double ToT=(1.0/(Energy+calib_C)+ calib_A)*calib_B;

//	return TMath::Floor(ToT);
//}


G4double AllPixFEI4RadDamageDigitizer::SlimEdgeEffect(G4int nX,G4double xpos,G4double eHit)
{

	G4double eEdge;
	G4double Etemp = 0.;

	G4double xpix= xpos + pitchX/2;
	//cout << "!!SLIMEDGES!!" << MipCharge << " " << xpix << endl;


	if(nX==0)
	{
		if((xpix < GRShift)){
			eEdge = eHit*(1-(GRShift-xpix)/(300*um));
	//cout << TString::Format("!!SLIM EDGES!! left nX:%d xpos:%f xpix:%f before:%f after:%f, R:%f",nX,xpos,xpix,eHit,Etemp,eHit*(1-(GRShift-xpix)/(GRShift))) << endl;

		}
		else{eEdge=eHit;};
		Etemp=eEdge;

		}

	else if(nX==(nPixX-1)){
		if(xpix > pitchX-GRShift){
			eEdge = eHit*(1-(xpix-(pitchX-GRShift))/(300*um));
	//cout << TString::Format("!!SLIM EDGES!! right NX:%d xpos:%f xpix:%f before:%f after:%f R:%f",nX,xpos,xpix,eHit,Etemp,eHit*(1-(xpix-(pitchX-GRShift))/(GRShift))) << endl;

		}
		else{eEdge=eHit;};
		Etemp=eEdge;


	}

	return Etemp;
}

G4bool AllPixFEI4RadDamageDigitizer::isSlimEdge(G4int nX, G4int /*nY*/)
{
	if(nX==0||nX==(nPixX-1)) return true;
	else return false;

}


