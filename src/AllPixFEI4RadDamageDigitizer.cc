/*
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

	//////////////////////////////////////////////////////////////////
	/////Setup all the inputs/////////////////////////////////////////
	////////////////////////////////////////////////////////////////// 
	//Physics switches
	doTrapping = true;
        doRamo = false;
	doDrift = false;

	//Constants
	elec = 3.64*eV;//average energy required to produce a e/h pair in Si.  This has been known for a long time - see for instance http://journals.aps.org/prb/pdf/10.1103/PhysRevB.1.2945.  It depends a little on temperature, but we are ignoring that effect here.
	betaElectrons = 3.0E-16*cm2/ns;  //The charge-trapping probability is t = -tau*ln(u), for u ~ Uniform(0,1) and tau^{-1}=beta*fluence.  The value of beta might be slightly higher for holes than electrons, but it is hard to say (we ignore this effect here).  See e.g. https://cds.cern.ch/record/685542/files/indet-2003-014.pdf. 
	diffusion_length = 0.007; //Should change this to the Einstein relation - current value from the ATLAS digitizer.

	//Conditions
	fluence = 0.; // neq/cm^2
	temperature = 263.2;// K  
	bField = 0.;// Tesla = V*s/m^2 
	threshold = 1600*elec; //This is the threshold for charge collection.
	tuning = 5./(20000*elec); //for X ToT @ Y e, this is X/Y so that a deposited energy of Y gives a ToT of X.  Typical values are 5 ToT @ 20ke.

	//Phenomenological parameters
	precision = 250; //this is the number of charges to divide the G4 hit into.

	//Efield, time, ramo maps
	TFile* efile=new TFile("/afs/cern.ch/work/b/bnachman/public/TCAD_maps/efieldmapping/Converted_TCAD/eFieldfl3.8e15v400_out.root");
	TFile* tfile=new TFile("/afs/cern.ch/work/b/bnachman/public/TCAD_maps/efieldmapping/Time_Maps/timeZ_fl3.8e15v400.root");
	TFile* rfile=new TFile("./share/absRamo3D-map-200um-output.root");
	//////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////

	// Get ramo potential mapping 
	ramoPotentialMap=0;
	//Need to get the Ramo potential for 250um!  
	ramoPotentialMap=(TH3F*)rfile->Get("hramomap1");
	if (ramoPotentialMap == 0){
	  G4cout << "Unsuccessful picking up histogram: ramoPotentialMap" << G4endl;
        }
	
	// Get electric field mapping
	eFieldMap=0;
	eFieldMap=(TH1F*)efile->Get("hefieldz");
	if (eFieldMap == 0){
	  G4cout << "Unsuccessful picking up histogram: eFieldMap" << G4endl;
	}
	
	// Get time to electrode mapping (derived from electric field)
	timeMap_e=0;
	timeMap_h=0;
	timeMap_e=(TH1F*)tfile->Get("etimes");
	timeMap_h=(TH1F*)tfile->Get("htimes");
	if (timeMap_e == 0 || timeMap_h == 0){
	  G4cout << "Unsuccessful picking up histogram: timeMap" << G4endl;
	}
	
	// Fetching info from pixeldetector.xml
	AllPixGeoDsc * gD = GetDetectorGeoDscPtr();
	resistivity=gD->GetResistivity();
	
	//Geometry constants
	detectorThickness = gD->GetSensorZ();
	pitchX = gD->GetPixelX(); 	// total length of pixel in x
	pitchY = gD->GetPixelY();	// total length of pixel in y
	nPixX= gD->GetNPixelsX(); 	// total number of pixels in x
	nPixY= gD->GetNPixelsY();	// total number of pixels in y
	
	//Compute the trapping time.
	if(fluence!=0.0)
	  {
	    trappingTimeElectrons = 1.0/(betaElectrons*fluence);
	    trappingTimeHoles = 1.0/(betaElectrons*fluence);
	  }
	else { //fluence = 0 so do not trap!
	  trappingTimeElectrons = 1000*s;
	  trappingTimeHoles = 1000*s;
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

	//These parameterizations come from C. Jacoboni et al., Solid‐State Electronics 20 (1977) 77‐89. (see also https://cds.cern.ch/record/684187/files/indet-2001-004.pdf).
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
	if(isHole) z = 0.25 - z; // For holes, should consider opposite direction.  0.25 mm is the depth of the sensor bulk.
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
	    int n_binz = timeMap_e->GetZaxis()->FindBin(0.25 - z);	// using z as distance to readout side. 0.25 mm is the depth of the sensor bulk.				
	    timeToElectrode = timeMap_e->GetBinContent(n_binz);	
	  }
	if(isHole)
	  {
	    int n_binz = timeMap_h->GetZaxis()->FindBin(0.25 - z);							
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


	bool dodebug = false;
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
		
		// In case the charge moves into a neighboring pixel
		pair<G4int, G4int> extraPixel;
		extraPixel = tempPixel;

		// Split the charge into subcharges (# = precision) that are diffused separately to the electrode
		for(G4int nQ  = 0 ; nQ < precision ; nQ++) {
				
		  G4double eHit = G4double(eHitTotal)/(2*precision); // eV; divide in half because we are treating holes separately	
		  G4double eHitbefore = eHit;
		  
		  // Loop over everything following twice, once for holes and once for electrons
		  for(G4int eholes=0 ; eholes<2 ; eholes++) { // Loop over everything twice, once for electrons and once for holes
		    
		    //Need to modify to only use holes for ramo.
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
		    G4double hallEffect = 1.13 + 0.0008*(temperature - 273.0);      //Hall Scattering Factor - taken from https://cds.cern.ch/record/684187/files/indet-2001-004.pdf
		    G4double tanLorentz = hallEffect*mobility*bField*(1.0E-3);	//unit conversion 
		    
		    G4double rdif=diffusion_length;
		    if (!doDrift) rdif = 0.; 
		    G4double xposD=xpos+zpos*tanLorentz+rdif*CLHEP::RandGauss::shoot(0,1); // Is it still +zpos in the case of the holes?                                               
		    G4double yposD=ypos+rdif*CLHEP::RandGauss::shoot(0,1);
		    
		    // Account for drifting into another pixel 
		    while (fabs(xposD) >= pitchX/2){
		      G4double sign = xposD/(fabs(xposD));                        // returns +1 or -1 depending on + or - x value                                                       
		      extraPixel.first = extraPixel.first + 1*sign;               // increments or decrements pixel count in x                                                          
		      xposD = xposD - (pitchX*sign);                              // moves xpos coordinate 1 pixel over in x                                                            
		    }
		    while (fabs(yposD) >= pitchY/2){
		      G4double sign = yposD/(fabs(yposD));                        // returns +1 or -1 depending on + or - y value                                                       
		      extraPixel.second = extraPixel.second + 1*sign;             // increments or decrements pixel count in y                                                          
		      yposD = yposD - (pitchY*sign);                              // moves xpos coordinate 1 pixel over in y                                                            
		    }
		    
		    
		    if (driftTime < timeToElectrode && doTrapping){ //charge was trapped
		      if (doRamo){
			// Also record deposit due to diff in ramo potential between (xposD, yposD, electrode) and (xpos, ypos, zpos)
			// Initial ramo potential based on (x,y,z) position in micrometers
			int nbin = 0;
			if(!isHole) nbin = ramoPotentialMap->FindBin(fabs(ypos*1000),fabs(xpos*1000),zpos*1000);
			if(isHole) nbin = ramoPotentialMap->FindBin(fabs(ypos*1000),fabs(xpos*1000),250-zpos*1000);
			G4double ramo_i = ramoPotentialMap->GetBinContent(nbin);
			
			// ramo potential at electrode based on (x,y,z) position in micrometers
			// -- loop in the x-coordinate
			for (int i=-1; i<=1; i++){
			  G4double x_neighbor = xposD + i*pitchX;
			  extraPixel.first += i;
			  if (i == 0){
			    extraPixel.first += 1;// For the middle neighbor, still have to add one to get back to zero
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
			    // Return ramo potential based on (x,y,z) position in micrometers; in z is at electrode
			    int nbin2 = ramoPotentialMap->FindBin(fabs(y_neighbor*1000),fabs(x_neighbor*1000),0);
			    G4double ramo = ramoPotentialMap->GetBinContent(nbin2);
			    
			    // Record deposit
			    G4double eHitRamo = eHit*(ramo - ramo_i);  //eV
			    pixelsContent[extraPixel] += eHitRamo; //eV
			  } //loop over y
			} //loop over x
		      } //doRamo
		    } //is trapped
		    else { //charge was not trapped or charge trapping is turned off.
		      
		      // Record deposit
		      pixelsContent[extraPixel] += eHit; // eV
		    }
		    
		  } // end loop over nQ charges
		} // end loop over nEntries
	} // end loop over charges/holes
	
	// Now that pixelContent is filled, create one digit per pixel
	map<pair<G4int, G4int>, G4double >::iterator pCItr = pixelsContent.begin();
	
	for( ; pCItr != pixelsContent.end() ; pCItr++)
	  {
	    G4double deposited_energy = (*pCItr).second;
	    int ToT = TMath::FloorNint(deposited_energy*tuning);
	    if (ToT >=15) ToT = 15; //FEI4 is 4-bit.
	    AllPixFEI4RadDamageDigit * digit = new AllPixFEI4RadDamageDigit;
	    digit->SetPixelIDX((*pCItr).first.first);
	    digit->SetPixelIDY((*pCItr).first.second);
	    digit->SetPixelCounts(ToT);
	    if (deposited_energy < threshold){
	      digit->SetPixelCounts(0);
	    }
	    if (deposited_energy >= threshold && ToT > 0){
	      m_digitsCollection->insert(digit);
	      if (dodebug){
		std::cout << "inserted a digit! " << std::endl;
	      }
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
	
	if (dodebug){
	  std::cout << "\n\n\n\n\n\n\n\n\n\n squirrel \n\n\n\n\n\n\n" << std::endl;
	}
	
} // end Digitize function





