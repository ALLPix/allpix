/**
 *  Author: Benjamin Nachman <bnachman@cern.ch>
 *    
 *
 *  allpix Authors:
 *   John Idarraga <idarraga@cern.ch>
 *   Mathieu Benoit <benoit@lal.in2p3.fr>
 */

#include "AllPixTruthWithLorentzDigitizer.hh"
#include "AllPixTrackerHit.hh"
#include "G4DigiManager.hh"
#include "G4VDigitizerModule.hh"
#include "AllPixGeoDsc.hh"

#include "TMath.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "CLHEP/Random/RandGauss.h"
#include "CLHEP/Random/RandFlat.h"

AllPixTruthWithLorentzDigitizer::AllPixTruthWithLorentzDigitizer(G4String modName, G4String hitsColName, G4String digitColName) 
: AllPixDigitizerInterface (modName) {
  
  // Registration of digits collection name
  collectionName.push_back(digitColName);
  m_hitsColName.push_back(hitsColName);
  
  // threshold
  m_digitIn.thl = 0.;
  
  ////Global Settings////
  //Physics Flags
  doPropagateFundamental = true   ; //if true, individual electrons are propagated to the electrode.  If false, charge chunks are propagated.
  doConstantEfield = false ; //if false, then a linear E-field is used (see inline for details).
  doDiffusion = true   ;
  doDigitize = true; //will use the tuning, threshold, and bits below to convert the charge to ToT.

  //Physical Settings
  double L = 150; //The depth of the sensor in um.
  biasVoltage = 300.0   ; //in Volts
  Temperature = 273   ; //K
  bfield = 2.0 ;// T 
  nToTbits = 8  ;
  chargeThreshold = 600; //electrons.
  chargeTuning = 128  ; //ToT @ MIP = 80 e / um
  method = 0 ; //method for converting charge to ToT: 0 is linear, 1 is kinked linear, 2 is exponential, 3 is linear with asynchronous timing.
  if (nToTbits < 0) doDigitize = false;
  ///////////////////////

  //Quickly compute the Lorentz angle; needed for time-to-electrode calculation.
  double Emax = 2*10*biasVoltage/(L/1000.); //in V/cm  
  mobility_e = GetMobility((1e-7)*Emax*0.5,0,Temperature);
  G4double r_Hall_h = 0.72 - 0.0005*(Temperature-273);
  G4double r_Hall_e = 1.13 - 0.0008*(Temperature-273);
  mobility_e*=r_Hall_e;
  //1 Tesla = V*s/m^2 = 0.001 MV*ns/mm^2                                                                                                                                           
  tanLorentz_e = bfield*mobility_e*0.001;
  
  //z = 0 is at the collecting electrode and a distance L from the bias plane.
  eFieldMap = new TH1F("hefieldz","hefieldz",200,0,L);
  for (int i=1; i<= eFieldMap->GetNbinsX(); i++){
    //The field should be linear from 2V/W (W=depletion depth) starting at the bias side down to zero at the collecting electrode.  See Pixel Detectors: From Fundamentals to Applications, Chapter 2.  We are assuming that the sensor is fully depleted, i.e. W = L.
    double efield = Emax*i/200;
    if (doConstantEfield) efield = Emax/2;
    eFieldMap->SetBinContent(i,efield); //in V/cm
  }  
  timeMap_e = new TH1F("etimes","Electron Time Map",100,0,L/1000.); //mm
  for (int k=1; k<=timeMap_e->GetNbinsX(); k++){
    double mysum = 0.;
    for (int k2=k; k2 >= 1; k2--){
      double z2 = timeMap_e->GetXaxis()->GetBinCenter(k2);
      double dz = timeMap_e->GetXaxis()->GetBinWidth(k2);
      double E = eFieldMap->GetBinContent(eFieldMap->GetXaxis()->FindBin(z2*1000))/1e7; //in MV/mm                                                                               
      if (E > 0){
	double mu = GetMobility(E, 0, Temperature); //mm^2/MV*ns
	mysum+=dz/(mu*E*cos(tanLorentz_e)); //mm * 1/(mm/ns) = ns 
      }
      timeMap_e->SetBinContent(k,mysum);
    }
  }
  
  TCanvas *c1 = new TCanvas("","",500,500);
  gStyle->SetOptStat(0);
  gPad->SetLeftMargin(0.15);
  eFieldMap->SetTitle("");
  eFieldMap->GetXaxis()->SetTitle("Depth (z) [#mum]");
  eFieldMap->GetYaxis()->SetTitle("E field [V/cm]");
  eFieldMap->GetYaxis()->SetTitleSize(0.03);
  eFieldMap->GetXaxis()->SetTitleSize(0.03);
  eFieldMap->GetYaxis()->SetTitleOffset(2.2);
  eFieldMap->GetXaxis()->SetTitleOffset(1.6);
  eFieldMap->Draw();
  c1->Print("Efield.pdf");
  
  timeMap_e->SetTitle("");
  timeMap_e->GetYaxis()->SetRangeUser(0,timeMap_e->GetMaximum());
  timeMap_e->GetXaxis()->SetTitle("Starting Pixel Depth in Z [mm]");
  timeMap_e->GetYaxis()->SetTitleOffset(1.4);
  timeMap_e->GetYaxis()->SetTitle("Projected Time to Reach Electrode [ns]");
  timeMap_e->GetXaxis()->SetNdivisions(505);
  timeMap_e->Draw();
  c1->Print("timemaps.pdf");
}

AllPixTruthWithLorentzDigitizer::~AllPixTruthWithLorentzDigitizer(){

}

G4double AllPixTruthWithLorentzDigitizer::GetMobility(G4double electricField, G4bool isHole, G4double temperature){

  // Initialize variables so they have the right scope
  G4double vsat = 0;
  G4double ecrit = 0;
  G4double beta = 0;

  //These parameterizations come from C. Jacoboni et al., Solid‐State Electronics 20 (1977) 77‐89. (see also https://cds.cern.ch/record/684187/files/indet-2001-004.pdf).
  if(!isHole){
    vsat = 15.3*pow(temperature,-0.87);     // mm/ns
    ecrit = 1.01E-7*pow(temperature,1.55);  // MV/mm
    beta = 2.57E-2*pow(temperature,0.66);
  }
  if(isHole){
    vsat = 1.62*pow(temperature,-0.52);     // mm/ns
    ecrit = 1.24E-7*pow(temperature,1.68);  // MV/mm
    beta = 0.46*pow(temperature,0.17);
  }

  G4double mobility = (vsat/ecrit)/pow(1+pow((electricField/ecrit),beta),(1/beta));
  return mobility; // mm^2/(MV*ns)                                                                                                                       
}

void AllPixTruthWithLorentzDigitizer::Digitize(){

	// create the digits collection
	m_digitsCollection = new AllPixTruthWithLorentzDigitsCollection("AllPixTruthWithLorentzDigitizer", collectionName[0] );

	// get the digiManager
	G4DigiManager * digiMan = G4DigiManager::GetDMpointer();

	// BoxSD_0_HitsCollection
	G4int hcID = digiMan->GetHitsCollectionID(m_hitsColName[0]);

	AllPixTrackerHitsCollection * hitsCollection = 0;
	hitsCollection = (AllPixTrackerHitsCollection*)(digiMan->GetHitsCollection(hcID));

	// temporary data structure
	map<pair<G4int, G4int>, G4double > pixelsContent;
	map<pair<G4int, G4int>, G4double > isdeltaContent;
	pair<G4int, G4int> tempPixel;

	G4int nEntries = hitsCollection->entries();

	// Example of detector description handle
	// provided by the interface
	AllPixGeoDsc * gD = GetDetectorGeoDscPtr();
	G4double pitchX = gD->GetPixelX();       // total length of pixel in x
        G4double pitchY = gD->GetPixelY();       // total length of pixel in y
	G4double detectorThickness = gD->GetSensorZ();
	G4double exit_x = 0.;
	G4double exit_y = 0.;

	if (nEntries == 0) return;

	G4double path_length_first_pixel = 0;
	//Can we compute the path length inside the first pixel?
	G4double xpos_i = (*hitsCollection)[0]->GetPosWithRespectToPixel().x();
	G4double ypos_i = (*hitsCollection)[0]->GetPosWithRespectToPixel().y();
	G4double zpos_i = detectorThickness/2. - (*hitsCollection)[0]->GetPosWithRespectToPixel().z();
	G4double ParentEnergy = (*hitsCollection)[0]->GetKinEParent();	
	G4int ParentType = (*hitsCollection)[0]->GetTrackPdgId();

	for(G4int itr  = 1 ; itr < nEntries ; itr++) {
	  G4double xpos = (*hitsCollection)[itr-1]->GetPosWithRespectToPixel().x();
          G4double ypos = (*hitsCollection)[itr-1]->GetPosWithRespectToPixel().y();
	  G4double zpos = detectorThickness/2. - (*hitsCollection)[itr-1]->GetPosWithRespectToPixel().z();
	  G4double xpos2 = (*hitsCollection)[itr]->GetPosWithRespectToPixel().x();
          G4double ypos2 = (*hitsCollection)[itr]->GetPosWithRespectToPixel().y();
          G4double zpos2 = detectorThickness/2. - (*hitsCollection)[itr]->GetPosWithRespectToPixel().z();
	  double trans_x = (xpos2-xpos)/pitchX;
	  if (trans_x < -0.8 || itr==nEntries-1){
	    double path_length = pow(pow(xpos-xpos_i,2)+pow(ypos-ypos_i,2)+pow(zpos-zpos_i,2),0.5);
	    //path_length_first_pixel = path_length;
	    path_length_first_pixel = fabs(xpos-xpos_i);
	    break;
	  }
	}

	//In ATLAS, we have 50 steps with 10 charges per step.
	//https://svnweb.cern.ch/trac/atlasoff/browser/InnerDetector/InDetDigitization/PixelDigitization/trunk/src/PixelECChargeTool.cxx
	for(G4int itr  = 0 ; itr < nEntries ; itr++) {
	  tempPixel.first  = (*hitsCollection)[itr]->GetPixelNbX();
	  tempPixel.second = (*hitsCollection)[itr]->GetPixelNbY();
	  G4double xpos = (*hitsCollection)[itr]->GetPosWithRespectToPixel().x(); //this is the eta direction
	  G4double ypos = (*hitsCollection)[itr]->GetPosWithRespectToPixel().y();  //this is the phi direction  
	  //We are only propagating electrons; no trapping in this (balistic) model, so we can forget about holes.
	  G4double zpos_e = detectorThickness/2. - (*hitsCollection)[itr]->GetPosWithRespectToPixel().z(); //z = 0 is the collecting electrode

	  //std::cout << "         " << (*hitsCollection)[itr]->GetTrackPdgId() << " " << itr << " " << zpos_e << " " << xpos << " " << ypos << " " << (*hitsCollection)[0]->GetTruthEntryLocalX() << std::endl; 

	  if ((*hitsCollection)[itr]->GetTrackPdgId() != 11){
	    exit_x = xpos;
	    exit_y = ypos;
	  }
	  //std::cout << (*hitsCollection)[itr]->GetTrackPdgId() << " " << itr << " " << zpos_e << " " << xpos << " " << ypos << " " << (*hitsCollection)[0]->GetTruthEntryLocalX() << std::endl;
	  G4double time_e = timeMap_e->GetBinContent(timeMap_e->GetXaxis()->FindBin(zpos_e)); //in ns
	  /*
	  //Diffusion
	  D = mu * kB * T / q;
	  D = (mu / mm^2/MV*ns) * (T/273 K) * 0.024 microns^2 / ns;
	  */
	  G4double Dt = mobility_e*(0.024)*Temperature/273.*time_e;
	  G4double rdif=sqrt(2*Dt)/1000; //in mm
          //std::cout << rdif/pitchX << std::endl;

	  if (!doDiffusion) rdif=0.;
	  int chunk_factor = int((*hitsCollection)[itr]->GetEdep()/(3.6*eV));
	  int nsubcharges = 1;
	  if (doPropagateFundamental) nsubcharges = chunk_factor;
	  G4double sub_energy = (*hitsCollection)[itr]->GetEdep() / nsubcharges;
	  for (G4int sub = 0 ; sub < nsubcharges ; sub++){
	    pair<G4int, G4int> extraPixel_e = tempPixel;
	    pair<G4int, G4int> extraPixel_h = tempPixel;
	    G4double ypos_e = ypos+CLHEP::RandGauss::shoot(0,1)*rdif+zpos_e*tanLorentz_e;
	    G4double xpos_e = xpos+CLHEP::RandGauss::shoot(0,1)*rdif;
	    // Account for drifting into another pixel
	    while (fabs(ypos_e) > pitchY/2){
	      G4double sign = ypos_e/(fabs(ypos_e));                        // returns +1 or -1 depending on + or - y value
	      extraPixel_e.second = extraPixel_e.second + 1*sign;           // increments or decrements pixel count in y
	      ypos_e = ypos_e - (pitchY*sign);                              // moves xpos coordinate 1 pixel over in y
	    }
	    while (fabs(xpos_e) > pitchX/2){
	      G4double sign = xpos_e/(fabs(xpos_e)); 
	      extraPixel_e.first = extraPixel_e.first + 1*sign; 
	      xpos_e = xpos_e - (pitchX*sign);
	    }
	    pixelsContent[extraPixel_e] += sub_energy;
	    isdeltaContent[extraPixel_e] += (*hitsCollection)[itr]->GetTrackPdgId() == 11 ? sub_energy : 0;
	  }
	}

	// Now create digits, one per pixel
	// Second entry in the map is the energy deposit in the pixel
	map<pair<G4int, G4int>, G4double >::iterator pCItr = pixelsContent.begin();
	for( ; pCItr != pixelsContent.end() ; pCItr++){
	  if((*pCItr).second > 0) // over threshold !
	    {
	      // Create one digit per pixel, I need to look at all the pixels first
	      AllPixTruthWithLorentzDigit * digit = new AllPixTruthWithLorentzDigit;
	      digit->SetPixelIDX((*pCItr).first.first);
	      digit->SetPixelIDY((*pCItr).first.second);
	      double EkeV = (*pCItr).second/keV;
	      digit->SetdeltaEfrac(isdeltaContent[(*pCItr).first]/(*pCItr).second);
	      if (doDigitize){
	      
	        double Q = (*pCItr).second/(3.6*eV) - chargeThreshold;
		double a = chargeTuning / (80*detectorThickness*1000);
		int n = nToTbits;
		if (method == 0) EkeV = min(int(a*max(Q,0.)),int(pow(2,n)-2))+1-(Q < 0);
		else if (method == 1) EkeV = Q < pow(2,n-1)/a ? min(int(a*max(Q,0.)),int(pow(2,n)-2))+1-(Q < 0) : min(int(0.5*a*(Q-pow(2,n-1)/a)+pow(2,n-1)+1),int(pow(2,n)-1));
		else if (method == 2) EkeV = Q < pow(2,n-1)/a ? min(int(a*max(Q,0.)),int(pow(2,n)-2))+1-(Q < 0) : min(int(log(1+a*(Q-pow(2,n-1)/a))/log(2)+pow(2,n-1)+1),int(pow(2,n)-1));
		else if (method == 3){
		     double t1 = 0.;
		     double tnm1 = (pow(2,n)-2)/a;
		     if (Q <= t1) EkeV = 0;
		     else if (Q >= tnm1 ) EkeV = pow(2,n)-2;
		     else{   	  
		     	  int tip1 = t1;
			  int myip1 = 1;
		     	  while (Q > tip1){
			  	tip1 += 1./a;
				myip1 += 1;
			  }
			  int ti = tip1 - 1./a;
			  double p = (Q - ti) / (tip1 - ti);
			  int coin = CLHEP::RandGauss::shoot(0,1) < p ? 0 : 1;
			  if (coin) EkeV = myip1-1;
			  else EkeV = myip1-2;
			  //flip a coin.  If heads, return i.  Otherwise, return i-1.
		     }
   		}		

		//EkeV = chargeTuning*(*pCItr).second/(3.6*eV*80*detectorThickness*1000);
		//if (int(EkeV) > pow(2,nToTbits)) EkeV = pow(2,nToTbits);
		//EkeV = int(EkeV);
	      }
	      else{
	        EkeV*=1000;
              }	      
	      
	      //Might need to add one to the first pixel.
	      if (((*pCItr).first.first == (*hitsCollection)[0]->GetPixelNbX()) && ((*pCItr).first.second == (*hitsCollection)[0]->GetPixelNbY())){
		if (int(EkeV) < 1){
		  //std::cout << "yo : " << EkeV << " " << int(EkeV) << " " << path_length_first_pixel << " " << (*pCItr).first.first << " " << std::endl;
		  //path_length_first_pixel += pitchX; //remember, X and Y are flipped.
		}
	      }
	      
	      //if (doDigitize && int(EkeV) < 1) continue;

	      if (EkeV == 1 && path_length_first_pixel/0.05 > 1){
		//std::cout << "wtf " << path_length_first_pixel/0.05 << " " << (*pCItr).first.first << std::endl;
		//exit(1);
	      }

	      if (EkeV > 0){
	      	 digit->SetPixelCounts(EkeV);//this will be recorded as an int.
	      	 digit->SetPOLangle((*hitsCollection)[0]->GetIncidentPOLAngle());
	      	 digit->SetAZMangle((*hitsCollection)[0]->GetIncidentAZMAngle());
	      	 digit->SetTruthEntryLocalX((*hitsCollection)[0]->GetTruthEntryLocalX());
	      	 digit->SetTruthEntryLocalY((*hitsCollection)[0]->GetTruthEntryLocalY());
	      	 digit->SetTruthExitLocalX(exit_x);
	      	 digit->SetTruthExitLocalY(exit_y);
	      	 digit->Setpath_length_first_pixel(path_length_first_pixel);
		 digit->SetInitE(ParentEnergy);
		 digit->SetInitID(ParentType);
              	 m_digitsCollection->insert(digit);	
	      }		
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
