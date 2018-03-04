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

//debugging plots
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"

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
        doRamo = true;
	doDiff = false;

	//Constants
	elec = 3.64*eV;//average energy required to produce a e/h pair in Si.  This has been known for a long time - see for instance http://journals.aps.org/prb/pdf/10.1103/PhysRevB.1.2945.  It depends a little on temperature, but we are ignoring that effect here.
	betaElectrons = 3.0E-16*cm2/ns;  //The charge-trapping probability is t = -tau*ln(u), for u ~ Uniform(0,1) and tau^{-1}=beta*fluence.  The value of beta might be slightly higher for holes than electrons, but it is hard to say (we ignore this effect here).  See e.g. https://cds.cern.ch/record/685542/files/indet-2003-014.pdf. 
	
	//Conditions
	fluence = 1e15*1/cm2; // neq/cm^2
	biasVoltage = 600; // V.  This is not used if external TCAD maps are supplied.
	temperature = 263.2;// K  
	bField = 0.;// Tesla = V*s/m^2 
	threshold = 1600*elec; //This is the threshold for charge collection.
	tuning = 8./(11000*elec); //for X ToT @ Y e, this is X/Y so that a deposited energy of Y gives a ToT of X.  Typical values are 5 ToT @ 20ke.

	//Phenomenological parameters
	precision = 100; //this is the number of charges to divide the G4 hit into.

	std::cout << "Load the input maps " << std::endl;
	
	//Efield, time, ramo maps
	//TFile* efile=new TFile("/afs/cern.ch/user/b/bnachman/work/public/RadDamage/allpix/TimeMaps/planar_fixed/Converted_TCAD/EField-map-fei4-200um-fl1e15-600V.root");
	TFile* efile=new TFile("/afs/cern.ch/user/b/bnachman/work/public/RadDamage/allpix/TimeMaps/planar_fixed/Converted_TCAD/EField-map-fei4-200um-fl5e15-1000V.root");
	//TFile* efile=new TFile("/afs/cern.ch/user/b/bnachman/work/public/RadDamage/allpix/TimeMaps/planar_fixed/Converted_TCAD/EField-map-fei4-200um-fl0-80V.root"); 
	TFile* tfile=new TFile("");
	TFile* dfile = new TFile("");
	TFile* rfile=new TFile("/afs/cern.ch/user/b/bnachman/work/public/RadDamage/allpix/share/absRamo3D-map-200um-output.root");
	//////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////

	//Fetching info from pixeldetector.xml 
	AllPixGeoDsc * gD = GetDetectorGeoDscPtr();
	double L = gD->GetSensorZ()*mm*1000;
	//For all maps, 0 is at the collecting electrode and L is the far side.

	//For debugging the inputs
	TCanvas *c1 = new TCanvas("","",500,500);
	gStyle->SetOptStat(0);

	// Get ramo potential mapping 
	ramoPotentialMap=0;
	ramoPotentialMap=(TH3F*)rfile->Get("hramomap1");
	if (ramoPotentialMap == 0){
	  G4cout << "Did not find a Ramo potential map.  Will use an approximate form." << G4endl;
	  //Requirements: goes to 0 at L and 1 at 0. Should also have the property that it decreases rapidtly away from 0 and then slowly as it approaches L. 
	  ramoPotentialMap = new TH3F("hramomap1","hramomap1",12,0,75,80,0,375,180,0,L);
	  for (int k=1; k<=ramoPotentialMap->GetNbinsZ(); k++){
	    double z = ramoPotentialMap->GetZaxis()->GetBinCenter(k);
	    double norm = exp(-10.)+exp(-1.);
	    double val = exp(-z/(0.1*L))+exp(-z/L);
	    val -= norm;
	    val /= (2.-norm); //should be 1 at 0 and 0 at L.
	    for (int i=1; i<=ramoPotentialMap->GetNbinsX(); i++){ //could add in some x and y dependence later.
	      for (int j=1; j<=ramoPotentialMap->GetNbinsY(); j++){
		ramoPotentialMap->SetBinContent(i,j,k,val);
	      }
	    }
	  }
        }
	
	// Get electric field mapping
	eFieldMap=0;
	eFieldMap=(TH1F*)efile->Get("hefieldz");
	if (eFieldMap == 0){ //Note that z=0 corresponds to the collecting electrode.
	  G4cout << "Did not find an E field map.  Will use an approximate form with no rad damage.  This will make the drift times inaccurate, but all other parts of the simulation should be valid.  This may or may not be important for your application." << G4endl;
	  eFieldMap = new TH1F("hefieldz","hefieldz",200,0,L);
	  for (int i=1; i<= eFieldMap->GetNbinsX(); i++){
	    eFieldMap->SetBinContent(i,(1e4)*biasVoltage/(L/1000.)); //in V/cm
	  }
	}	
	gPad->SetLeftMargin(0.15);
	eFieldMap->GetXaxis()->SetTitle("Depth (z) [#mum]");
	eFieldMap->GetYaxis()->SetTitle("E field [V/cm]");
	eFieldMap->GetYaxis()->SetTitleSize(0.03);
	eFieldMap->GetXaxis()->SetTitleSize(0.03);
	eFieldMap->GetYaxis()->SetTitleOffset(2.2);
	eFieldMap->GetXaxis()->SetTitleOffset(1.6);
	eFieldMap->Draw();
	c1->Print("Efield.pdf");
	
	// Get the distance at the point of trap and time to electrode mapping (derived from electric field)  
        distancemap_e=0;
        distancemap_h=0;
	timeMap_e=0;
        timeMap_h=0;
        distancemap_e=(TH2F*)dfile->Get("edistance");
        distancemap_h=(TH2F*)dfile->Get("hdistance");
	timeMap_e=(TH1F*)tfile->Get("etimes");
        timeMap_h=(TH1F*)tfile->Get("htimes");
        if (distancemap_e == 0 || distancemap_h == 0 || timeMap_e == 0 || timeMap_h == 0){
	  G4cout << "Did not find any pre-computed distance maps.  Will quickly do the integration now.  This is slow, but only needs to be done once per run." << G4endl;
	  distancemap_e = new TH2F("edistance","Electron Distance Map",100,0,L/1000.,20,0,10); //mm by ns
	  distancemap_h = new TH2F("hdistance","Holes Distance Map",100,0,L/1000.,20,0,10);
	  timeMap_e = new TH1F("etimes","Electron Time Map",100,0,L/1000.); //mm
	  timeMap_h = new TH1F("htimes","Hole Time Map",100,0,L/1000.); //mm   
	  for (int i=1; i<= distancemap_e->GetNbinsX(); i++){
	    for (int j=1; j<= distancemap_e->GetNbinsY(); j++){
	      distancemap_h->SetBinContent(i,j,L/1000.); //if you travel long enough, you will reach the electrode.
	    }
	  }
	  for (int k=1; k<= distancemap_e->GetNbinsX(); k++){
	    double z = distancemap_e->GetXaxis()->GetBinCenter(k);
	    double mysum = 0.;
	    double mysum_h = 0.;
	    for (int k2=k; k2 >= 1; k2--){
	      double z2 = distancemap_e->GetXaxis()->GetBinCenter(k2);
	      double dz = distancemap_e->GetXaxis()->GetBinWidth(k2);
	      double E = eFieldMap->GetBinContent(eFieldMap->GetXaxis()->FindBin(z2*1000))/1e7; //in MV/mm
	      if (E > 0){
		double mu = GetMobility(E, 0); //mm^2/MV*ns
		mysum+=dz/(mu*E); //mm * 1/(mm/ns) = ns
		distancemap_e->SetBinContent(k,distancemap_e->GetYaxis()->FindBin(mysum),z2);
	      }
	      timeMap_e->SetBinContent(k,mysum);
	    }
	    for (int k2=k; k2 <= distancemap_e->GetNbinsX(); k2++){ //holes go the opposite direction as electrons.
	      double z2 = distancemap_e->GetXaxis()->GetBinCenter(k2);
	      double dz = distancemap_e->GetXaxis()->GetBinWidth(k2);
	      double E = eFieldMap->GetBinContent(eFieldMap->GetXaxis()->FindBin(z2*1000))/1e7;;
	      if (E > 0){
		double mu_h = GetMobility(E, 1);
		mysum_h+=dz/(mu_h*E);
		distancemap_h->SetBinContent(k,distancemap_h->GetYaxis()->FindBin(mysum_h),z2);
	      }
	      timeMap_h->SetBinContent(k,mysum_h);
	    }
	  }
	}
	gPad->SetRightMargin(0.15);
	distancemap_e->GetXaxis()->SetNdivisions(505);
	distancemap_e->GetYaxis()->SetNdivisions(505);
	distancemap_e->GetXaxis()->SetTitle("Initial Position in Z [mm]");
	distancemap_e->GetYaxis()->SetTitle("Time Traveled [ns]");
	distancemap_e->GetZaxis()->SetTitleOffset(1.6);
	distancemap_e->SetTitle("Electron Distance Map");
	distancemap_e->GetZaxis()->SetTitle("Location in Z [mm]");
	distancemap_e->Draw("colz");
	c1->Print("distancemap_e.pdf");
	
	distancemap_h->GetXaxis()->SetNdivisions(505);
        distancemap_h->GetYaxis()->SetNdivisions(505);
        distancemap_h->GetXaxis()->SetTitle("Initial Position in Z [mm]");
        distancemap_h->GetYaxis()->SetTitle("Time Traveled [ns]");
        distancemap_h->GetZaxis()->SetTitleOffset(1.6);
        distancemap_h->SetTitle("Hole Distance Map");
        distancemap_h->GetZaxis()->SetTitle("Location in Z [mm]");
        distancemap_h->Draw("colz");
        c1->Print("distancemap_h.pdf");

	timeMap_h->SetTitle("");
	timeMap_h->GetYaxis()->SetRangeUser(0,timeMap_h->GetMaximum() > timeMap_e->GetMaximum() ? timeMap_h->GetMaximum() : timeMap_e->GetMaximum());
	timeMap_h->GetXaxis()->SetTitle("Starting Pixel Depth in Z [mm]");
	timeMap_h->GetYaxis()->SetTitleOffset(1.4);
	timeMap_h->GetYaxis()->SetTitle("Projected Time to Reach Electrode [ns]");
	timeMap_h->GetXaxis()->SetNdivisions(505);
	timeMap_h->Draw();
	timeMap_e->SetLineStyle(7);
	timeMap_e->Draw("same");

	TLegend * leg = new TLegend(0.4,0.65,0.85,0.9);
	leg->SetFillColor(0);
	leg->SetBorderSize(0);
	leg->SetFillStyle(0);
	leg->AddEntry(timeMap_e,"Electrons","l");
	leg->AddEntry(timeMap_h,"Holes","l");
	leg->Draw();
	c1->Print("timemaps.pdf");
	
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
	  trappingTimeElectrons = 1000*s; //~infinity
	  trappingTimeHoles = 1000*s;
	}
	
	//Need to pre-compute a new map that allows us to properly take into account the charge chunking. 
	double steps = 1000;
	charge_chunk_map_e = new TH1F("","",int(0.9*steps),0,detectorThickness); //just make sure there are fewer bins than steps.
	charge_chunk_map_h = new TH1F("","",int(0.9*steps),0,detectorThickness); //just make sure there are fewer bins than steps.    
        for (double z = 0; z <= detectorThickness; z+= detectorThickness/steps){
	  //Could eventually make this depend on x and y.  For now, only z.
          //for (int i=1; i<= nPixX; i++){
	  //for (int j=1; j<= nPixY; j++){
	  G4double timeToElectrode = GetTimeToElectrode(z, 0);
	  G4double timeToElectrode_h = GetTimeToElectrode(z, 1);
	  double prob = 0.;
	  double prob_h=0.;
	  double prob_trap=0.;
	  double z_avg = 0.;
	  if (timeToElectrode == 0 || timeToElectrode_h == 0) continue;
	  double ramo_intial = ramoPotentialMap->GetBinContent(ramoPotentialMap->FindBin(0.,0.,z*1000));
	  for (double t = 0; t <= timeToElectrode; t+=timeToElectrode/steps){
	    double final_pos = distancemap_e->GetBinContent(distancemap_e->FindBin(z,t)); //with respect to the electrode (at 0)
	    double ramo_final = ramoPotentialMap->GetBinContent(ramoPotentialMap->FindBin(0.,0.,final_pos*1000));
	    prob+=(ramo_final-ramo_intial)*(timeToElectrode/steps)*exp(-t/trappingTimeElectrons)/trappingTimeElectrons;
	    prob_trap+=(timeToElectrode/steps)*exp(-t/trappingTimeElectrons)/trappingTimeElectrons;
	    z_avg+=final_pos*(timeToElectrode/steps)*exp(-t/trappingTimeElectrons)/trappingTimeElectrons;
	  } 
	  for (double t = 0; t <= timeToElectrode_h; t+=timeToElectrode_h/steps){
	    double final_pos_h = distancemap_h->GetBinContent(distancemap_h->FindBin(z,t));
	    double ramo_final_h = ramoPotentialMap->GetBinContent(ramoPotentialMap->FindBin(0.,0.,final_pos_h*1000));
	    prob_h+=-(ramo_final_h-ramo_intial)*(timeToElectrode_h/steps)*exp(-t/trappingTimeHoles)/trappingTimeHoles;
          }
	  //If there is no induced charge, prob = 1 - exp(-timeToElectrode/trappingTimeElectrons).
	  //std::cout << z << " " << prob << " " << prob_h << " " << exp(-timeToElectrode/trappingTimeElectrons) << " " << ramo_intial << " " << z_avg/prob_trap << std::endl;
	  charge_chunk_map_e->SetBinContent(charge_chunk_map_e->GetXaxis()->FindBin(z),prob);
	  charge_chunk_map_h->SetBinContent(charge_chunk_map_h->GetXaxis()->FindBin(z),prob_h);
	  //}
          //}
        }
	std::cout << "here " << std::endl; 

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

G4double AllPixFEI4RadDamageDigitizer::GetMobility(G4double electricField, G4bool isHole){ 
	
	// Initialize variables so they have the right scope
	G4double vsat = 0;
	G4double ecrit = 0;
	G4double beta = 0;	

	//These parameterizations come from C. Jacoboni et al., Solid‐State Electronics 20 (1977) 77‐89. (see also https://cds.cern.ch/record/684187/files/indet-2001-004.pdf).
	if(!isHole){
		vsat = 15.3*pow(temperature,-0.87);	// mm/ns
		ecrit = 1.01E-7*pow(temperature,1.55);	// MV/mm
		beta = 2.57E-2*pow(temperature,0.66);
	}
	if(isHole){
		vsat = 1.62*pow(temperature,-0.52);	// mm/ns
		ecrit = 1.24E-7*pow(temperature,1.68);	// MV/mm
		beta = 0.46*pow(temperature,0.17);
	}

	G4double mobility = (vsat/ecrit)/pow(1+pow((electricField/ecrit),beta),(1/beta));	
	return mobility;  					// mm^2/(MV*ns)
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
	    int n_binz = timeMap_e->GetXaxis()->FindBin(z);
	    timeToElectrode = timeMap_e->GetBinContent(n_binz);	
	  }
	if(isHole)
	  {
	    int n_binz = timeMap_h->GetXaxis()->FindBin(z);							
	    timeToElectrode = timeMap_h->GetBinContent(n_binz);			
	  }
	return timeToElectrode;	
}

G4double AllPixFEI4RadDamageDigitizer::GetTanLorentz(G4double electricField, G4bool isHole){
  G4double hallEffect = 1.13 + 0.0008*(temperature - 273.0); //Hall Scattering Factor - taken from https://cds.cern.ch/record/684187/files/indet-2001-004.pdf
  G4double mobility = AllPixFEI4RadDamageDigitizer::GetMobility(electricField, isHole);
  G4double tanLorentz = hallEffect*mobility*bField*(1.0E-3);  //unit conversion
  return tanLorentz;
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

		G4double xpos = (*hitsCollection)[itr]->GetPosWithRespectToPixel().x(); // mm !!!BE WARNED: THIS IS THE ETA DIRECTION.  LORENTZ ANGLE CURRENTLY WRONG!!!
		G4double ypos = (*hitsCollection)[itr]->GetPosWithRespectToPixel().y(); // mm !!!BE WARNED: THIS IS THE ETA DIRECTION.  LORENTZ ANGLE CURRENTLY WRONG!!!

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

		  //Need to determine how many elementary charges this charge chunk represents.
		  double chunk_size = eHit/(0.5*elec);
		  double kappa = 1./sqrt(chunk_size);
		  
		  // Loop over everything following twice, once for holes and once for electrons
		  for(G4int eholes=0 ; eholes<2 ; eholes++) { // Loop over everything twice, once for electrons and once for holes
		    
		    //Need to modify to only use holes for ramo.
		    isHole = false; // Set a condition to keep track of electron/hole-specific functions
		    if (eholes == 1) isHole = true;
		    
		    G4double electricField = GetElectricField(zpos);

		    // Reset extraPixel coordinates each time through loop
		    extraPixel = tempPixel;

		    G4double timeToElectrode = GetTimeToElectrode(zpos, isHole); //ns
		    //std::cout << "    ardvark   " << timeToElectrode << std::endl;
		    G4double driftTime = GetDriftTime(isHole);
		    G4double drift_time_constant = trappingTimeElectrons;
		    if (isHole) drift_time_constant = trappingTimeHoles;
		    double time_ratio = timeToElectrode/drift_time_constant;
		    double charge_correction_exp = exp(-time_ratio);
		    double average_charge = charge_chunk_map_e->GetBinContent(charge_chunk_map_e->FindBin(detectorThickness-zpos));
		    if (isHole) average_charge = charge_chunk_map_h->GetBinContent(charge_chunk_map_h->FindBin(detectorThickness-zpos));

		    G4double tanLorentz = GetTanLorentz(electricField, isHole);
		    G4double zposD = zpos; //this is the distance between the initial position and the final position.
		    if (isHole) zposD = detectorThickness - zpos;
		    if ((driftTime < timeToElectrode) && doTrapping){
		      int nbin=0;
		      if(!isHole){
			nbin = distancemap_e->FindBin(zpos,driftTime);
			zposD = zpos-distancemap_e->GetBinContent(nbin);
		      }
		      else{
			nbin = distancemap_h->FindBin(zpos,driftTime);
			zposD = distancemap_h->GetBinContent(nbin) - zpos;
		      }
		    }
		    /*
		      D = mu * kB * T / q
		      D = (mu / mm^2/MV*ns) * (T/273 K) * 0.024 microns^2 / ns
		    */
		    G4double Dt = GetMobility(electricField, isHole)*(0.024)*timeToElectrode*temperature/273.;
		    G4double rdif=sqrt(Dt)/1000; //in mm

		    if (!doDiff) rdif = 0.; 
		    G4double xposD=xpos+zposD*tanLorentz+rdif*CLHEP::RandGauss::shoot(0,1); //Both e and h move in the same direciton under B-field (q-reversed but also the velocity direction)
		    G4double yposD=ypos+rdif*CLHEP::RandGauss::shoot(0,1);
		    
		    // Account for drifting into another pixel 
		    while (fabs(xposD) > pitchX/2){
		      G4double sign = xposD/(fabs(xposD));                        // returns +1 or -1 depending on + or - x value 
		      extraPixel.first = extraPixel.first + 1*sign;               // increments or decrements pixel count in x
		      xposD = xposD - (pitchX*sign);                              // moves xpos coordinate 1 pixel over in x
		    }
		    while (fabs(yposD) > pitchY/2){
		      G4double sign = yposD/(fabs(yposD));                        // returns +1 or -1 depending on + or - y value
		      extraPixel.second = extraPixel.second + 1*sign;             // increments or decrements pixel count in y
		      yposD = yposD - (pitchY*sign);                              // moves xpos coordinate 1 pixel over in y 
		    }
		    
		    int loc_x = tempPixel.first;
                    int loc_y = tempPixel.second;
		    
		    if ((driftTime < timeToElectrode) && doTrapping){ //charge was trapped
		      if (doRamo){
			//Record the induced charge from the difference in the ramo potential between the trapped location and all electordes.			
			// ramo potential at electrode based on (x,y,z) position in micrometers
			// -- loop in the x-coordinate
			for (int i=-1; i<=1; i++){
			  G4double x_neighbor = i*pitchX;
			  extraPixel.first = loc_x + i;
			  // -- loop in the y-coordinate
			  for (int j=-1; j<=1; j++){
			    G4double y_neighbor = j*pitchY;
			    extraPixel.second = loc_y + j;
			    int nbin = 0;
			    if(!isHole) nbin = ramoPotentialMap->FindBin(fabs((ypos-y_neighbor)*1000),fabs((xpos-x_neighbor)*1000),zpos*1000); //take the absolute value because we only have 1/8 of the ramo 
			    else nbin = ramoPotentialMap->FindBin(fabs((ypos-y_neighbor)*1000),fabs((xpos-x_neighbor)*1000),zpos*1000);
			    G4double ramo_i = ramoPotentialMap->GetBinContent(nbin);
			    int nbin2 = 0;
			    //Our distance and time maps are in 1D.  Near the collecting electrode, the electrons would bend toward x=y=0.  We mimic this by just setting the final position to x=y=0.
			    if(!isHole) nbin2 = ramoPotentialMap->FindBin(fabs(y_neighbor*1000),fabs(x_neighbor*1000),max(0.,(zpos-zposD))*1000);
			    else nbin2 = ramoPotentialMap->FindBin(fabs((yposD-y_neighbor)*1000),fabs((xposD-x_neighbor)*1000),(zpos+zposD)*1000);
			    G4double ramo = ramoPotentialMap->GetBinContent(nbin2);
			    
			    // Record deposit
			    G4double eHitRamo = (1-2*isHole)*eHit*(ramo - ramo_i);  //eV
			    
			    //X' -> \mu + kappa * (X-\mu)
			    if (isHole) pixelsContent[extraPixel] += eHit*average_charge + kappa*(eHitRamo - eHit*average_charge); //eV.
			    else pixelsContent[extraPixel] += eHit*(average_charge+charge_correction_exp) + kappa*(eHitRamo - eHit*(average_charge+charge_correction_exp));
			  } //loop over y
			} //loop over x
		      } //doRamo
		      
		      // Record deposit - take into account charge chunks that represent more than one fundamental charge!
		      //X' -> \mu + kappa * (X-\mu) 
		      if (!doRamo) pixelsContent[extraPixel] += eHit*(charge_correction_exp+kappa*(0.-charge_correction_exp)); //eV

		    } //is trapped
		    else { //charge was not trapped or charge trapping is turned off.
		      
		      // Record deposit
		      //X' -> \mu + kappa * (X-\mu) 
		      if (!doTrapping && !doRamo) pixelsContent[extraPixel] += eHit; // eV
		      else if (!doRamo) pixelsContent[extraPixel] += eHit*(charge_correction_exp+kappa*(1.-charge_correction_exp)); //eV
		      else pixelsContent[extraPixel] += eHit*(average_charge+charge_correction_exp) + kappa*(eHit - eHit*(average_charge+charge_correction_exp));
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





