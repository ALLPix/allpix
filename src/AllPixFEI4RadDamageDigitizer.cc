/*
 *  Authors:
 *    Benjamin Nachman <bnachman@cern.ch>
 *    Lingxin Meng <lingxin.meng@cern.ch>
 *    Lorenzo Rossini <lorenzo.rossini@cern.ch>
 *    Rebecca Carney <rebecca.carney@cern.ch>
 *    Mathieu Benoit <Mathieu.Benoit@cern.ch>  
 *    Callie Bertsche <c.bertsche@cern.ch>
 *
 *  Allpix Authors:
 *   John Idarraga <idarraga@cern.ch>
 *   Mathieu Benoit <Mathieu.Benoit@cern.ch>
 *
 *  Brief Description:
 *   This is a digitization model for radiation damage effects.  
 *   See "Setup all the inputs" below for available options.
 *   The user can supply TCAD maps for the electric field; however,
 *   if the maps are not given, reasonable defaults are used (though
 *   use with caution).
 *
 */

#include "AllPixFEI4RadDamageDigitizer.hh"
#include "AllPixTrackerHit.hh"
#include "G4DigiManager.hh"
#include "G4VDigitizerModule.hh"
#include "AllPixGeoDsc.hh"

#include "TMath.h"
#include "TF1.h"
#include "TTree.h"
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
	doSimplifiedModel = false;

	//Physics defaults
	defaultRamo = 10; // -1 = approximate z-dependence only; 0 = exact XZ x YZ with compensating Z-dependence; > 0 = solution to Poisson's equation with defaultRamo terms in the series.  See doc/RadDamageDefaults.pdf for details.
	defaultEfield = -1; // 0 = uniform; 1 = linear.
	debug_maps = true; //if true, prints plots of all the input maps.

	//Constants
	elec = 3.64*eV;//average energy required to produce a e/h pair in Si.  This has been known for a long time - first measurement was Phys. Rev. 91, 1079 (1953).
	betaElectrons = 3.0E-16*cm2/ns;  //The charge-trapping probability is t = -tau*ln(u), for u ~ Uniform(0,1) and tau^{-1}=beta*fluence.  The value of beta might be slightly higher for holes than electrons, but it is hard to say (we ignore this effect here).  See e.g. https://cds.cern.ch/record/685542/files/indet-2003-014.pdf. 
	
	//Conditions
	fluence = 0*1/cm2; // 1 MeV neq/cm^2
	biasVoltage = 600; // V.  This is not used if external TCAD maps are supplied.
	temperature = 263.2;// K  
	bField = 0.;// Tesla = V*s/m^2 
	threshold = 2000*elec; //This is the threshold for charge collection.
	tuning = 9./(20000*elec); //for X ToT @ Y e, this is X/Y so that a deposited energy of Y gives a ToT of X.  Typical values are 5 ToT @ 20ke.

	//Phenomenological parameters
	precision = 100; //this is the number of charges to divide the G4 hit into.

	deplationLenght=0.200;   //mm IBL 
	if(depVoltage>biasVoltage) deplationLenght=deplationLenght*pow(biasVoltage/depVoltage,0.5);
	std::cout << "Load the input maps " << std::endl;
	
	//Input from TCAD
	//TFile* efile=new TFile("/afs/cern.ch/user/b/bnachman/work/public/RadDamage/allpix/TimeMaps/planar_fixed/Converted_TCAD/EField-map-fei4-200um-fl1e15-600V.root");
	//TFile* efile=new TFile("/afs/cern.ch/user/b/bnachman/work/public/RadDamage/allpix/TimeMaps/planar_fixed/Converted_TCAD/EField-map-fei4-200um-fl5e15-1000V.root");
	TFile* efile=new TFile("");// /afs/cern.ch/user/b/bnachman/work/public/RadDamage/allpix/TimeMaps/planar_fixed/Converted_TCAD/EField-map-fei4-200um-fl0-80V.root"); 
	TFile* tfile=new TFile("");
	TFile* dfile = new TFile("");
	TFile* rfile=new TFile(""); ///afs/cern.ch/user/b/bnachman/work/public/RadDamage/allpix/share/absRamo3D-map-200um-output.root");
	//////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////

	//Sensor properties
	//Fetching info from pixeldetector.xml 
	AllPixGeoDsc * gD = GetDetectorGeoDscPtr();
	double L = gD->GetSensorZ()*mm*1000; //in microns; for all maps, 0 is at the collecting electrode and L is the far side.
	int L_int = L; //for the histograms later that have one bin per micron.
	depVoltage = 30.; //V; this is highly fluece-dependent !  Should update when you change the fluence.
	deplationLenght = L/1000.; //in mm.
	if(depVoltage>biasVoltage) deplationLenght=deplationLenght*pow(biasVoltage/depVoltage,0.5); //This is the usual formula depletion depth = \sqrt{2*\epsilon_0\epsilon_{Si}*V/(eN_D)}.  See e.g. (2.26) in Pixel Detectors by L. Rossi et al.

	//Geometry constants
	detectorThickness = gD->GetSensorZ();
	pitchX = gD->GetPixelX(); 	// total length of pixel in x
	pitchY = gD->GetPixelY();	// total length of pixel in y
	nPixX= gD->GetNPixelsX(); 	// total number of pixels in x
	nPixY= gD->GetNPixelsY();	// total number of pixels in y

	//For debugging the inputs
	TCanvas *c1 = new TCanvas("","",600,600);
	gStyle->SetOptStat(0);

	// Get the Ramo potential map.  See "Shockley-Ramo Theorem": Journal of Applied Physics 9 (1938) 635 and Proceedings of the IRE 27 (1939) 584.
	ramoPotentialMap=0;
	ramoPotentialMap=(TH3F*)rfile->Get("hramomap1");
	if (ramoPotentialMap == 0){
	  G4cout << "Did not find a Ramo potential map.  Will use an approximate form." << G4endl;
	  //There are three default Ramo potentials.  The third one should be the most accurate.  See doc/RadDamageDefaults.pdf for details.
	  
	  int binxramo=2*pitchX*1000/10; //Ramo potential is evaluated for distances that are up to 2 * pitch away from the center of the primary pixel.  One bin per 10 microns.
	  int binyramo=2*pitchY*1000/10;
	  double xminramo=0.;
	  double xmaxramo=2*pitchX*1000;
	  double yminramo=0.;
	  double ymaxramo=2*pitchY*1000;

	  ramoPotentialMap = new TH3F("hramomap1","hramomap1",binxramo,xminramo,xmaxramo,binyramo,yminramo,ymaxramo,L_int/10,0.,L);
	  G4double pi=TMath::Pi();
	  for (int k=1; k<=ramoPotentialMap->GetNbinsZ(); k++){
	    G4double z = ramoPotentialMap->GetZaxis()->GetBinCenter(k)-ramoPotentialMap->GetZaxis()->GetBinWidth(k)/2.; //use the lower bin edge.
	    for (int i=1; i<=ramoPotentialMap->GetNbinsX(); i++){
	      for (int j=1; j<=ramoPotentialMap->GetNbinsY(); j++){
		G4double x = ramoPotentialMap->GetXaxis()->GetBinCenter(i)-ramoPotentialMap->GetXaxis()->GetBinWidth(i)/2.;
		G4double y = ramoPotentialMap->GetYaxis()->GetBinCenter(j)-ramoPotentialMap->GetYaxis()->GetBinWidth(j)/2.;

		if (defaultRamo == -1){
		  if (x > pitchX*1000 / 2. || y > pitchY*1000 / 2.) ramoPotentialMap->SetBinContent(i,j,k,0.01); //outside of the primary pixel.
		  else{
		    double a = 3*L/(1000*pitchY); //Y is the short direction for FEI4.
		    double norm = exp(-a)+exp(-1.);
		    double val = exp(-a*z/L)+exp(-z/L);
		    val -= norm;
		    val /= (2.-norm); //should be 1 at 0 and 0 at L.
		    ramoPotentialMap->SetBinContent(i,j,k,val);
		    std::cout << "do I ever make it khere ??? " << z << " " << val << " " << a << std::endl;
		  }
		}//end Ramo default optoin -1.

		if (defaultRamo == 0){
		  double norm = exp(-10.)+exp(-1.);
		  double val = exp(-z/(0.1*L))+exp(-z/L);
		  val -= norm;
		  val /= (2.-norm); //should be 1 at 0 and 0 at L.
		  ramoPotentialMap->SetBinContent(i,j,k,Phi(x,z,pitchX*1000,L) * Phi(y,z,pitchY*1000,L) * val / (Phi(0,z,pitchX*1000,L) * Phi(0,z,pitchY*1000,L)));
		}//end Ramo default option 0.
		
		if (defaultRamo > 0){
		  ramoPotentialMap->SetBinContent(i,j,k,Phi3D(x/L,y/L,z/L,defaultRamo,defaultRamo,4,pitchX*1000/L,pitchY*1000/L)); //N = 4 is arbitrary; just need something bigger than ~1ish.
		}//end Ramo default option > 0.
	      }//end loop over y.
	    }//end loop over x.
	  }//end loop over z.
	}//end default Ramo potential.	

	if (debug_maps){
	  ramoPotentialMap2D = new TH2F("hramomap2D","hramomap2D",ramoPotentialMap->GetNbinsX(),ramoPotentialMap->GetXaxis()->GetBinCenter(1)-0.5*ramoPotentialMap->GetXaxis()->GetBinWidth(1),ramoPotentialMap->GetXaxis()->GetBinCenter(ramoPotentialMap->GetNbinsX())+0.5*ramoPotentialMap->GetXaxis()->GetBinWidth(ramoPotentialMap->GetNbinsX()),ramoPotentialMap->GetNbinsY(),ramoPotentialMap->GetYaxis()->GetBinCenter(1)-0.5*ramoPotentialMap->GetYaxis()->GetBinWidth(1),ramoPotentialMap->GetYaxis()->GetBinCenter(ramoPotentialMap->GetNbinsY())+0.5*ramoPotentialMap->GetYaxis()->GetBinWidth(ramoPotentialMap->GetNbinsY()));
	  ramoPotentialMap1D = new TH1F("hramomap1D","hramomap1D",ramoPotentialMap->GetNbinsZ(),ramoPotentialMap->GetZaxis()->GetBinCenter(1)-0.5*ramoPotentialMap->GetZaxis()->GetBinWidth(1),ramoPotentialMap->GetZaxis()->GetBinCenter(ramoPotentialMap->GetNbinsZ())+0.5*ramoPotentialMap->GetZaxis()->GetBinWidth(ramoPotentialMap->GetNbinsZ()));
	  for (int k=1; k<=ramoPotentialMap->GetNbinsZ(); k++){
	    for (int i=1; i<=ramoPotentialMap->GetNbinsX(); i++){
	      for (int j=1; j<=ramoPotentialMap->GetNbinsY(); j++){
		ramoPotentialMap2D->SetBinContent(i,j,ramoPotentialMap->GetBinContent(i,j,k));
		if (i==1 && j==1) ramoPotentialMap1D->SetBinContent(k,ramoPotentialMap->GetBinContent(i,j,k));
	      }
	    }
	    gPad->SetRightMargin(0.15);
	    ramoPotentialMap2D->GetXaxis()->SetNdivisions(505);
	    ramoPotentialMap2D->GetYaxis()->SetNdivisions(505);
	    ramoPotentialMap2D->GetXaxis()->SetTitle("#eta position [#mum]");
	    ramoPotentialMap2D->GetYaxis()->SetTitle("#phi position [#mum]");
	    ramoPotentialMap2D->GetZaxis()->SetTitleOffset(1.6);
	    ramoPotentialMap2D->SetTitle("Ramo Potential at z = "+TString::Format("%0.2f",ramoPotentialMap->GetZaxis()->GetBinCenter(k)-0.5*ramoPotentialMap->GetZaxis()->GetBinWidth(k)));
	    ramoPotentialMap2D->GetZaxis()->SetTitle("Ramo Potential");
	    ramoPotentialMap2D->GetZaxis()->SetRangeUser(0,1);
	    ramoPotentialMap2D->Draw("colz");
	    c1->Print("DefaultRamo2D_"+TString::Format("%i",k)+".pdf");
	  }
	  gPad->SetRightMargin(0.05);
	  ramoPotentialMap1D->SetTitle("Ramo Potential at x=y=0");
	  ramoPotentialMap1D->GetYaxis()->SetRangeUser(0,1);
	  ramoPotentialMap1D->GetXaxis()->SetTitle("Pixel Depth in z [#mum]");
	  ramoPotentialMap1D->GetYaxis()->SetTitleOffset(1.4);
	  ramoPotentialMap1D->GetYaxis()->SetTitle("Ramo Potential");
	  ramoPotentialMap1D->GetXaxis()->SetNdivisions(505);
	  ramoPotentialMap1D->Draw();
	  c1->Print("DefaultRamo1D.pdf");
	}

  	m_eFieldMap1D=0;
	m_eFieldMap1D=(TH1F*)efile->Get("hefieldz");
	if (m_eFieldMap1D == 0){ //Note that z=0 corresponds to the collecting electrode.
	  G4cout << "Did not find an E-field map.  Will use an approximate form." << G4endl;
	  //WHAT ABOUT 3D E-field maps ???  This is what Lorenzo is using at the moment.
	  m_eFieldMap1D = new TH1F("m_hefieldz","m_hefieldz",L_int,0,L_int); // 1D map of the Efield in case we want to use 1D field
	  G4double electricField=0;
	  G4double pass=deplationLenght/m_eFieldMap1D->GetNbinsX();
	  for (int i=1; i<= m_eFieldMap1D->GetNbinsX()+1; i++){
	    G4double position=pass*i;
	    if(biasVoltage<depVoltage)   electricField=(biasVoltage/deplationLenght)*(1-position/deplationLenght);
   	    if(biasVoltage>=depVoltage)  electricField=(depVoltage/deplationLenght)*(1-position/deplationLenght)+(biasVoltage-depVoltage)/(2*deplationLenght);
   	    
   	    if(position>deplationLenght) electricField=0.;
	    m_eFieldMap1D->SetBinContent(i,electricField*10); //in V/cm
	    if (defaultEfield==0) m_eFieldMap1D->SetBinContent(i,(1e4)*biasVoltage/(L/1000.)); //in V/cm
	  }
	} 

	if (debug_maps){
	  gPad->SetLeftMargin(0.15);
	  m_eFieldMap1D->GetXaxis()->SetTitle("Depth (z) [#mum]");
	  m_eFieldMap1D->GetYaxis()->SetTitle("E field [V/cm]");
	  m_eFieldMap1D->GetYaxis()->SetTitleOffset(1.6);
	  m_eFieldMap1D->GetYaxis()->SetTitleSize(0.03);
	  m_eFieldMap1D->GetXaxis()->SetTitleSize(0.03);
	  m_eFieldMap1D->Draw();
	  c1->Print("Efield.pdf");
	}

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
	      double E = m_eFieldMap1D->GetBinContent(m_eFieldMap1D->GetXaxis()->FindBin(z2*1000))/1e7; //in MV/mm
	      //double E = eFieldMap->GetBinContent(eFieldMap->FindBin(1.,1.,z2*1000))/1e7; //in MV/mm
	      //double E = GetElectricField(0.002,0.002,z2); //in MV/mm
	      //double E = GetElectricField(z2); //in MV/mm
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
	      double E = m_eFieldMap1D->GetBinContent(m_eFieldMap1D->GetXaxis()->FindBin(z2*1000))/1e7;;
	      //double E = eFieldMap->GetBinContent(eFieldMap->FindBin(1.,1.,z2*1000-1.))/1e
	      //double E = GetElectricField(0.002,0.002,z2); //in MV/mm
	      //double E = GetElectricField(z2); //in MV/mm
	      if (E > 0){
		double mu_h = GetMobility(E, 1);
		mysum_h+=dz/(mu_h*E);
		distancemap_h->SetBinContent(k,distancemap_h->GetYaxis()->FindBin(mysum_h),z2);
	      }
	      timeMap_h->SetBinContent(k,mysum_h);
	    }
	  }
	}

	if (debug_maps){
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
	}

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
	// old chunk chunking map. Still here if it is needed
	/*
	double steps = 1000;
	steps=L_int;
	charge_chunk_map_e = new TH3F("","",200,-200.,200,750,-750.,750,L_int,0.,L); //just make sure there are fewer bins than steps.
	charge_chunk_map_h = new TH3F("","",200,-200.,200,750,-750.,750,L_int,0.,L);  //just make sure there are fewer bins than steps. 
	int nstepX=charge_chunk_map_e->GetNbinsX();
	int nstepY=charge_chunk_map_e->GetNbinsY();   
        //for (double z = 0; z <= detectorThickness; z+= detectorThickness/steps)
        int nstepZ=charge_chunk_map_e->GetNbinsZ();   
        for (int k = 0; k <= nstepZ; k++){
         double z=charge_chunk_map_e->GetZaxis()->GetBinCenter(k);
	  //Could eventually make this depend on x and y.  For now, only z.
         for (int i=1; i<= nstepX; i++){
          double x = charge_chunk_map_e->GetXaxis()->GetBinCenter(i);
	  for (int j=1; j<= nstepY; j++){
	  double y = charge_chunk_map_e->GetYaxis()->GetBinCenter(j);
	  G4double timeToElectrode = GetTimeToElectrode(z, 0);
	  G4double timeToElectrode_h = GetTimeToElectrode(z, 1);
	  double prob = 0.;
	  double prob_h=0.;
	  double prob_trap=0.;
	  double z_avg = 0.;
	  if (timeToElectrode == 0 || timeToElectrode_h == 0) continue;
	  double ramo_intial = ramoPotentialMap->GetBinContent(ramoPotentialMap->FindBin(x,y,z));
	  for (double t = 0; t <= timeToElectrode; t+=timeToElectrode/steps){
	    double final_pos = distancemap_e->GetBinContent(distancemap_e->FindBin(z,t)); //with respect to the electrode (at 0)
	    double ramo_final = ramoPotentialMap->GetBinContent(ramoPotentialMap->FindBin(x,y,final_pos));
	    prob+=(ramo_final-ramo_intial)*(timeToElectrode/steps)*exp(-t/trappingTimeElectrons)/trappingTimeElectrons;
	    prob_trap+=(timeToElectrode/steps)*exp(-t/trappingTimeElectrons)/trappingTimeElectrons;
	    z_avg+=final_pos*(timeToElectrode/steps)*exp(-t/trappingTimeElectrons)/trappingTimeElectrons;
	  } 
	  for (double t = 0; t <= timeToElectrode_h; t+=timeToElectrode_h/steps){
	    double final_pos_h = distancemap_h->GetBinContent(distancemap_h->FindBin(z,t));
	    double ramo_final_h = ramoPotentialMap->GetBinContent(ramoPotentialMap->FindBin(x,y,final_pos_h));
	    prob_h+=-(ramo_final_h-ramo_intial)*(timeToElectrode_h/steps)*exp(-t/trappingTimeHoles)/trappingTimeHoles;
          }
	  //If there is no induced charge, prob = 1 - exp(-timeToElectrode/trappingTimeElectrons).
	  //std::cout << z << " " << prob << " " << prob_h << " " << exp(-timeToElectrode/trappingTimeElectrons) << " " << ramo_intial << " " << z_avg/prob_trap << std::endl;
	  charge_chunk_map_e->SetBinContent(i,j,charge_chunk_map_e->GetZaxis()->FindBin(z),prob);
	  charge_chunk_map_h->SetBinContent(i,j,charge_chunk_map_h->GetZaxis()->FindBin(z),prob_h);
	  }
         }
        }
        */
        //Need to pre-compute a new map that allows us to properly take into account the charge chunking. 
	int stepsZ = int(L);
	int stepsY = int(pitchY*500); //number of bin= half the width of the sensor
	int stepsX = int(pitchX*500);
	charge_chunk_map_e = new TH3F("","",stepsX,0.,0.25,stepsY,0.,0.125,stepsZ,0,detectorThickness); 
	charge_chunk_map_h = new TH3F("","",stepsX,0.,0.25,stepsY,0.,0.125,stepsZ,0,detectorThickness); 
        //for (double z = 0; z <= detectorThickness; z+= detectorThickness/steps)
        for (int k = 0; k <= stepsZ; k++){
	  double z =  charge_chunk_map_e->GetZaxis()->GetBinCenter(k);
         //for (int i=1; i<= int(0.9*stepsX); i++)
         for (int i=1; i<= stepsX; i++){  				//not sure of the correctness of the dependency on x and y. above there is 
          double x = charge_chunk_map_e->GetXaxis()->GetBinCenter(i);   // still the version with just 1D
	  //for (int j=1; j<= int(0.9*stepsY); j++)
	  for (int j=1; j<= stepsY; j++){
	  double y= charge_chunk_map_e->GetYaxis()->GetBinCenter(j);
	  G4double timeToElectrode = GetTimeToElectrode(z, 0);
	  G4double timeToElectrode_h = GetTimeToElectrode(z, 1);
	  double prob = 0.;
	  double prob_h=0.;
	  double prob_trap=0.;
	  double z_avg = 0.;
	  if (timeToElectrode == 0 || timeToElectrode_h == 0) continue;
	  double ramo_intial = ramoPotentialMap->GetBinContent(ramoPotentialMap->FindBin(x*1000,y*1000,z*1000));
	  for (double t = 0; t <= timeToElectrode; t+=timeToElectrode/stepsZ){
	    double final_pos = distancemap_e->GetBinContent(distancemap_e->FindBin(z,t)); //with respect to the electrode (at 0)
	    double ramo_final = ramoPotentialMap->GetBinContent(ramoPotentialMap->FindBin(x*1000.,y*1000.,final_pos*1000));
	    prob+=(ramo_final-ramo_intial)*(timeToElectrode/stepsZ)*exp(-t/trappingTimeElectrons)/trappingTimeElectrons;
	    prob_trap+=(timeToElectrode/stepsZ)*exp(-t/trappingTimeElectrons)/trappingTimeElectrons;
	    z_avg+=final_pos*(timeToElectrode/stepsZ)*exp(-t/trappingTimeElectrons)/trappingTimeElectrons;
	  } 
	  for (double t = 0; t <= timeToElectrode_h; t+=timeToElectrode_h/stepsZ){
	    double final_pos_h = distancemap_h->GetBinContent(distancemap_h->FindBin(z,t));
	    double ramo_final_h = ramoPotentialMap->GetBinContent(ramoPotentialMap->FindBin(0.,0.,final_pos_h*1000));
	    prob_h+=-(ramo_final_h-ramo_intial)*(timeToElectrode_h/stepsZ)*exp(-t/trappingTimeHoles)/trappingTimeHoles;
          }
	  //If there is no induced charge, prob = 1 - exp(-timeToElectrode/trappingTimeElectrons).
	  //std::cout << z << " " << prob << " " << prob_h << " " << exp(-timeToElectrode/trappingTimeElectrons) << " " << ramo_intial << " " << z_avg/prob_trap << std::endl;
	  charge_chunk_map_e->SetBinContent(charge_chunk_map_e->FindBin(x,y,z),prob);
	  charge_chunk_map_h->SetBinContent(charge_chunk_map_h->FindBin(x,y,z),prob_h);
	  }
         }
        }
        
        
        if (lorentz_map_e == 0 || lorentz_map_h==0 ){
	  G4cout << "Did not find any pre-computed lorentz maps.  Will quickly do the integration now.  This is slow, but only needs to be done once per run." << G4endl;
	  lorentz_map_e = new TH2F("lorentz_map_e","Lorentz Map e",100,0,L/1000.,100,0,L/1000); //mm by ns
	  lorentz_map_h = new TH2F("lorentz_map_h","Lorentz Map h",100,0,L/1000.,100,0,L/1000); //mm by ns
	  
	  
	  for (int k=1; k<= lorentz_map_e->GetNbinsX(); k++){
	    double z = lorentz_map_e->GetXaxis()->GetBinCenter(k);
	    double mysum = 0.;
	    double mysum_h = 0.;
	    double lenght_e=0;
	    double lenght_h=0;
	    for (int k2=k; k2 >= 1; k2--){
	      double z2 = lorentz_map_e->GetXaxis()->GetBinCenter(k2);
	      double dz = lorentz_map_e->GetXaxis()->GetBinWidth(k2);	            
	      //double E = GetElectricField(0.002,0.002,z2); //in MV/mm
	      double E = m_eFieldMap1D->GetBinContent(m_eFieldMap1D->GetXaxis()->FindBin(z2*1000))/1e7;;
	      if (E > 0){
	         lenght_e+=dz;
		// integral tan thetaL * dz / integral dz
		double tanthetaL = GetTanLorentz(E, 0);	 
		mysum+=(tanthetaL*dz);				  // the angle is computed as the integral over the path
		lorentz_map_e->SetBinContent(k,k2,mysum/lenght_e);
		
	      }
	      //lorentz_map_e->SetBinContent(k,mysum);
	    }
	    for (int k2=k; k2 <= lorentz_map_e->GetNbinsX(); k2++){ //holes go the opposite direction as electrons.
	      double z2 = lorentz_map_e->GetXaxis()->GetBinCenter(k2);
	      double dz = lorentz_map_e->GetXaxis()->GetBinWidth(k2);
	      //double E = GetElectricField(0.002,0.002,z2); //in MV/mm
	      double E = m_eFieldMap1D->GetBinContent(m_eFieldMap1D->GetXaxis()->FindBin(z2*1000))/1e7;;
	      //double E = GetElectricField(z2); //in MV/mm
	      if (E > 0){
		double tanthetaL = GetTanLorentz(E, 1);
		lenght_h+=dz;
		mysum_h+=(tanthetaL *dz);
		lorentz_map_h->SetBinContent(k,k2,mysum_h/lenght_h);
	      }
	     // lorentz_map_h->SetBinContent(k,mysum_h);
	    }
	  }
	}      
	
	if (debug_maps){
	  gPad->SetRightMargin(0.15);
	  lorentz_map_e->GetXaxis()->SetNdivisions(505);
	  lorentz_map_e->GetYaxis()->SetNdivisions(505);
	  lorentz_map_e->GetXaxis()->SetTitle("Initial Position in Z [mm]");
	  lorentz_map_e->GetYaxis()->SetTitle("Distance Traveled in Z [mm]");
	  lorentz_map_e->GetZaxis()->SetTitleOffset(1.6);
	  lorentz_map_e->GetYaxis()->SetTitleOffset(1.3);
	  lorentz_map_e->SetTitle("Electron Lorentz Map");
	  lorentz_map_e->GetZaxis()->SetTitle("Tangent Lorentz Angle");
	  lorentz_map_e->Draw("colz");
	  c1->Print("lorentz_map_e.pdf");
	  
	  lorentz_map_h->GetXaxis()->SetNdivisions(505);
	  lorentz_map_h->GetYaxis()->SetNdivisions(505);
	  lorentz_map_h->GetXaxis()->SetTitle("Initial Position in Z [mm]");
	  lorentz_map_h->GetYaxis()->SetTitle("Distance Traveled in Z [mm]");
	  lorentz_map_h->GetYaxis()->SetTitleOffset(1.3);
	  lorentz_map_h->GetZaxis()->SetTitleOffset(1.6);
	  lorentz_map_h->SetTitle("Hole Lorentz Map");
	  lorentz_map_h->GetZaxis()->SetTitle("Tangent Lorentz Angle");
	  lorentz_map_h->Draw("colz");
	  c1->Print("lorentz_map_h.pdf");
	}

}// End of AllPixFEI4RadDamageDigitizer::AllPixFEI4RadDamageDigitizer definition

AllPixFEI4RadDamageDigitizer::~AllPixFEI4RadDamageDigitizer(){
}


G4double AllPixFEI4RadDamageDigitizer::Phi(G4double x, G4double z, G4double Lx, G4double L){
  
  if (z==0) z=0.00001;//a pathology in the definition.
  G4double pi = 4.*TMath::ATan(1.);
  G4double val = (TMath::Sin(pi*z/L)*TMath::SinH(0.5*pi*Lx/L)/(TMath::CosH(pi*x/L)-TMath::Cos(pi*z/L)*TMath::CosH(0.5*pi*Lx/L)));
  if (val > 0) return TMath::ATan(val)/pi;
  else return TMath::ATan(val)/pi+1;
}

G4double AllPixFEI4RadDamageDigitizer::N(G4double z, G4double L){
  z=z/L;
  return (pow(z+1.,2)-1.35*pow(z+1,3)+0.6*pow(z+1,4))/(0.25);
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
	// return the value of the 1D dimensional Efield -> Linear approximation
        G4int nbinEmap=m_eFieldMap1D->GetNbinsX();
        
	
	int n_binz = m_eFieldMap1D->GetXaxis()->FindBin(z*1000);
								
	G4double electricField;

	if(n_binz>nbinEmap){ 
	 
	 G4int bindist=n_binz-nbinEmap;
	 G4double binwidth=m_eFieldMap1D->GetBinWidth(nbinEmap);
	 if(m_eFieldMap1D->GetBinContent(nbinEmap)==0) nbinEmap=nbinEmap-1;
	 electricField=m_eFieldMap1D->GetBinContent(nbinEmap) - bindist*(m_eFieldMap1D->GetBinContent(nbinEmap)-m_eFieldMap1D->GetBinContent(nbinEmap-1))/m_eFieldMap1D->GetBinContent(nbinEmap);
	} else 	electricField = m_eFieldMap1D->GetBinContent(n_binz);	
		
	//electricField=biasVoltage*10/deplationLenght;
	//G4cout<<" binz "<<n_binz<<" nbin tot "<<nbinEmap<<" field "<< electricField<<endl;
	return electricField*1.0E-7;
}

G4double AllPixFEI4RadDamageDigitizer::GetElectricField(G4double x, G4double y, G4double z){
	// The position is in mm, but plot uses um to return electric field in V/cm. Convert to return electric field in MV/mm
	// Is currently using z as distance from position to readout side
	// return the value of the 3D E field
	// TCAD maps are just 1/4 of the pixel cell
	// uses fabs() in order to overcome this
	G4double electricField=0;
	if(eFieldMap==0){
		 
	 electricField=(GetElectricField(z))*1E07;  // if there is no Efield map loaded,get the "default" 1D Efield
	}else{
	 
	 while(fabs(x)>pitchX/2) x=fabs(x)-pitchX/2;  // if it is outside the limits of the cell, get the E field of the other pixel
	 while(fabs(y)>pitchY/2) y=fabs(y)-pitchY/2;
	
	 if(fabs(x)>pitchX/2-0.001 && fabs(x)<=pitchX/2) x=fabs(x)-0.001;  // maps doesn't allow for near border value, ie:
	 if(fabs(y)>pitchY/2-0.001 && fabs(y)<=pitchY/2) y=fabs(y)-0.001;  // maps goes from 1 um to 24 um (for the X sides of the IBL, and similary
         if(fabs(z)>0.199) z=fabs(z)-0.001;				  // for the other sides) so, if it is needed to evaluate the field there
        								  // it just get the field in the nearest position available
         if(fabs(x)<0.001) x=0.0012;
	 if(fabs(y)<0.001) y=0.0012;
         if(fabs(z)<0.001) z=0.0012;
        
         int n_bin = eFieldMap->FindBin(fabs(x)*1000,fabs(y)*1000,z*1000); //x y and z should be in micrometers
        
        
         electricField = eFieldMap->GetBinContent(n_bin);
         if(z<0) electricField=0;
        }
        
        
        return electricField*1.0E-7;
        //return 4000*1.0E-7;

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
  if (isHole) hallEffect = 0.72 - 0.0005*(temperature - 273.0);
  G4double mobility = AllPixFEI4RadDamageDigitizer::GetMobility(electricField, isHole);
  //if(!isHole) mobility =1100;
  //if(isHole) mobility =500;
  G4double tanLorentz = hallEffect*mobility*bField*(1.0E-3);  //unit conversion
  return tanLorentz;
}

G4double AllPixFEI4RadDamageDigitizer::GetTanLorentz(G4double z1, G4double z2, G4bool isHole){

  G4double tanLorentz = -999;
  int nbin=0;
  if(!isHole){
        if(z2>detectorThickness)z2=detectorThickness;
        if(z2<0)z2=0;
	nbin = lorentz_map_e->FindBin(z1,z2);
	tanLorentz = lorentz_map_e->GetBinContent(nbin);
  }else{
        if(z2>=detectorThickness) z2=detectorThickness-0.001;
        if(z2<0)z2=0.001; 
       
	nbin = lorentz_map_h->FindBin(z1,z2);
	tanLorentz = lorentz_map_h->GetBinContent(nbin) ;
  }
  
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
	pair<G4int, G4int> sumPixel;					// (countx,county) which pixel
	map<pair<G4int, G4int>, G4double > sumpixelsContent;		// stored energy per (countx,county) pixel	
	pair<G4int, G4int> AppotempPixel;
	
	

	bool dodebug = false;
	// Loop over the whole Hits Collection - this spreads out the effect of the hit over its path through the pixel 
	G4int nEntries = hitsCollection->entries();
	G4double sumEnToT=0;
	for(G4int itr  = 0 ; itr < nEntries ; itr++)  sumEnToT+= (*hitsCollection)[itr]->GetEdep(); 
	for(G4int itr  = 0 ; itr < nEntries ; itr++) {

		G4double xpos = (*hitsCollection)[itr]->GetPosWithRespectToPixel().x(); // mm !!!BE WARNED: THIS IS THE ETA DIRECTION.  LORENTZ ANGLE CURRENTLY WRONG!!!
		G4double ypos = (*hitsCollection)[itr]->GetPosWithRespectToPixel().y(); // mm !!!BE WARNED: THIS IS THE ETA DIRECTION.  LORENTZ ANGLE CURRENTLY WRONG!!!

		// Given as distance from center of pixel, + toward electrode (where electrons are read out)
		G4double zpos = detectorThickness/2. - (*hitsCollection)[itr]->GetPosWithRespectToPixel().z();	// mm; Gives distance to readout side for electron

		G4double eHitTotal = (*hitsCollection)[itr]->GetEdep(); // Energy deposition for the hit, internal unit
		tempPixel.first  = (*hitsCollection)[itr]->GetPixelNbX();
		tempPixel.second = (*hitsCollection)[itr]->GetPixelNbY();


		
		if(itr>0){ // this is done in order to prevent a bug where if the charge was exactly on the border of the pixel it would change the 
			   // pixel position, eg: xpos1=23 um pixelX 121 -> xpos2=25 um pixelX 122, where it should have been xpos2=25 um pixelX 121
		 if(((pitchX/2)-fabs(xpos)) <0.000001 && fabs((*hitsCollection)[itr]->GetPixelNbX() - (*hitsCollection)[itr-1]->GetPixelNbX())>1) {
		  tempPixel.first  = tempPixel.first -1;
		  G4cout<<" ok, change "<<endl;
		 }
		}
		// In case the charge moves into a neighboring pixel
		pair<G4int, G4int> extraPixel;
		extraPixel = tempPixel;

		// Split the charge into subcharges (# = precision) that are diffused separately to the electrode
		for(G4int nQ  = 0 ; nQ < precision ; nQ++) {
				
		  //G4double eHit = G4double(eHitTotal)/(2*precision);  // eV; divide in half because we are treating holes separately	
		  G4double eHit = G4double(eHitTotal)/(precision);      // eV; but charges are always "e", regardless of them beeing 
		  							//electrons or holes, isn't it?		  
		 						 	// point should be: if using ramo ->no 2. If NOT using ramo ->yes 2
		  G4double eHitbefore = eHit;			
		  
		  //Need to determine how many elementary charges this charge chunk represents.
		  double chunk_size = eHit/(0.5*elec);			// this here should not depend on the discussion above, since
		  //double chunk_size = eHit/elec;			// it is just a numerical factor
		  double kappa = 1./sqrt(chunk_size);
		  
		  // Loop over everything following twice, once for holes and once for electrons
		  for(G4int eholes=0 ; eholes<2 ; eholes++) { // Loop over everything twice, once for electrons and once for holes
		    
		    //Need to modify to only use holes for ramo.
		    isHole = false; // Set a condition to keep track of electron/hole-specific functions
		    if (eholes == 1) isHole = true;
		    //if(detectorThickness-zpos<0.001) zpos = detectorThickness -0.001;
		    
		    //G4double electricField = GetElectricField(zpos);	    		//1D Efield
		    G4double electricField = GetElectricField(xpos,ypos,zpos);		//3D EField
		    //if(fabs(electricField)<0.000001*1.0E-6) continue;

		    // Reset extraPixel coordinates each time through loop
		    extraPixel = tempPixel;

		    G4double timeToElectrode = GetTimeToElectrode(zpos, isHole); //ns
		    //if(nQ==0)G4cout << "    ardvark   " << timeToElectrode <<" electricField "<<electricField<< G4endl;
		    G4double driftTime = GetDriftTime(isHole);
		    G4double drift_time_constant = trappingTimeElectrons;
		    if (isHole) drift_time_constant = trappingTimeHoles;
		    double time_ratio = timeToElectrode/drift_time_constant;
		    double charge_correction_exp = exp(-time_ratio);
		    //double average_charge = charge_chunk_map_e->GetBinContent(charge_chunk_map_e->FindBin(zpos*1000));
		    //if (isHole) average_charge = charge_chunk_map_h->GetBinContent(charge_chunk_map_h->FindBin((detectorThickness-zpos)*1000));
		    double average_charge = charge_chunk_map_e->GetBinContent(charge_chunk_map_e->FindBin(fabs(xpos)*1000,fabs(ypos)*1000,zpos*1000));
		    if (isHole) average_charge = charge_chunk_map_h->GetBinContent(charge_chunk_map_h->FindBin(fabs(xpos)*1000,fabs(ypos)*1000,zpos*1000));
		    
		    
		   // G4double tanLorentz = GetTanLorentz(electricField, isHole);
		    //tanLorentz=0.22;
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
		    G4double tanLorentz =-1.;
		    if(!isHole) tanLorentz =GetTanLorentz(zpos,zpos-zposD, isHole);
		    else tanLorentz =GetTanLorentz(zpos,zpos+zposD, isHole);
		    /*
		      D = mu * kB * T / q
		      D = (mu / mm^2/MV*ns) * (T/273 K) * 0.024 microns^2 / ns
		    */
		    G4double Dt = GetMobility(electricField, isHole)*(0.024)*timeToElectrode*temperature/273.;
		    G4double rdif=sqrt(Dt)/1000; //in mm
		    //double rdif=0.007*sqrt(zpos/0.3); 		//ATHENA diffusion
		    if (!doDiff) rdif = 0.; 
		    G4double termRand=CLHEP::RandGauss::shoot(0,1);
		    G4double xposD=xpos+zposD*tanLorentz+termRand*rdif; //Both e and h move in the same direciton under B-field (q-reversed but also the velocity direction)
		    G4double yposD=ypos+rdif*CLHEP::RandGauss::shoot(0,1);
		    
		    
		    int loc_x = tempPixel.first;
                    int loc_y = tempPixel.second;
		    sumPixel=tempPixel;
		   //if ((driftTime < timeToElectrode) && doTrapping) //charge was trapped
		   if ( doTrapping)    //if I want to use always the Ramo potential, even when charges are not trapped
		    {
		      //G4double TrapeHit=2*eHit; 
		      //G4double TrapeHit=eHit; 
  		      if (doRamo){
			//Record the induced charge from the difference in the ramo potential between the trapped location and all electordes.			
			// ramo potential at electrode based on (x,y,z) position in micrometers
			// -- loop in the x-coordinate
			for (int i=-2; i<=2; i++){
			  G4double x_neighbor = i*pitchX;
			  extraPixel.first = loc_x + i;
			  // -- loop in the y-coordinate
			  for (int j=-2; j<=2; j++){
			    //pixelsContent[extraPixel]=0;
			    G4double y_neighbor = j*pitchY;
			    extraPixel.second = loc_y + j;
			    int nbin = 0;
			    //G4double ramo_i=0;
			    double ramo_i=0;
			    
			     int nbinmaxX=ramoPotentialMap->GetNbinsX();
			     int nbinmaxY=ramoPotentialMap->GetNbinsY();
			     double binwidthX=ramoPotentialMap->GetXaxis()->GetBinWidth(2);
			     double binwidthY=ramoPotentialMap->GetYaxis()->GetBinWidth(2);
			     
			    if(  fabs(xpos-x_neighbor)*1000<(nbinmaxX*binwidthX) && fabs(ypos-y_neighbor)*1000<(nbinmaxY*binwidthY) ) { // check if the position we are evaluating the rampo potential is inside the map, if not->ramo=0
			     if(!isHole) nbin = ramoPotentialMap->FindBin(fabs((xpos-x_neighbor)*1000),fabs((ypos-y_neighbor)*1000),zpos*1000); //take the absolute value because we only have 1/8 of the ramo 
			     else nbin = ramoPotentialMap->FindBin(fabs((xpos-x_neighbor)*1000),fabs((ypos-y_neighbor)*1000),zpos*1000);
			     ramo_i = ramoPotentialMap->GetBinContent(nbin);
			     
			    }

			    int nbin2 = 0;
			    //Our distance and time maps are in 1D.  Near the collecting electrode, the electrons would bend toward x=y=0.  We mimic this by just setting the final position to x=y=0.
			    //G4double ramo=0;
			    double ramo=0;
			    if(  fabs(xposD-x_neighbor)*1000<(nbinmaxX*binwidthX) && fabs(yposD-y_neighbor)*1000<(nbinmaxY*binwidthY) ){// check if the position we are evaluating the rampo potential is inside the map, if not->ramo_i=0
			     
			     //if(!isHole) nbin2 = ramoPotentialMap->FindBin(fabs((xposD-x_neighbor)*1000),fabs((yposD-y_neighbor)*1000),max(0.0001,(zpos-zposD))*1000);
			     if(!isHole) nbin2 = ramoPotentialMap->FindBin(fabs((xposD-x_neighbor)*1000),fabs((yposD-y_neighbor)*1000),(zpos-zposD)*1000);
			     else nbin2 = ramoPotentialMap->FindBin(fabs((xposD-x_neighbor)*1000),fabs((yposD-y_neighbor)*1000),(zpos+zposD)*1000);
			    
			     ramo = ramoPotentialMap->GetBinContent(nbin2);		     
			    }
			     if( (zpos-zposD==0.)  &&  (fabs(xposD-x_neighbor)>pitchX/2 || fabs(yposD-y_neighbor)>pitchY/2) ) ramo=0;
			     if( (zpos-zposD==0.)  &&  fabs(xposD-x_neighbor)<=pitchX/2 && fabs(yposD-y_neighbor)<=pitchY/2 ) ramo=1;
			    // Record deposit
			    //if(nQ==0 && isHole)G4cout << isHole <<" z+zD "<< zpos+zposD << " xD " << xposD-x_neighbor << " yD " << yposD-y_neighbor << " ramo_i " << ramo_i << " ramo " << ramo << " pitchX "<<pitchX<<" pitchY "<<pitchY<<G4endl; 
			    //if(nQ==0 && !isHole)G4cout << isHole <<" z-zD "<< zpos-zposD << " xD " << xposD-x_neighbor<< " pitchX/2 "<<pitchX/2<< " yD " << yposD-y_neighbor << " pitchY/2 "<<pitchY/2<< " ramo_i " << ramo_i << " ramo " << ramo << G4endl;   
			   //  if(nQ==0 && !isHole)G4cout<<"----------------"<<endl;
			    //G4double eHitRamo = (1-2*isHole)*eHit*(ramo - ramo_i);  //eV
			    double eHitRamo =(1-2*isHole)*eHit*(ramo - ramo_i);  //eV
			    //double eHitRamo =(1-2*isHole)*TrapeHit*(ramo - ramo_i);  //eV
			    double pixelprecont=pixelsContent[extraPixel];
			    pixelsContent[extraPixel] += eHitRamo;  // if you don't want to use charge chunk just uncomment here and comment afterwards
			    if(fabs(xposD-x_neighbor)<pitchX/2 && fabs(yposD-y_neighbor)<pitchY/2 ) sumpixelsContent[extraPixel] += eHitRamo;
			    //X' -> \mu + kappa * (X-\mu)
			    
			    double appo_outside=1.;			     // charge chunking part still under revision :(
		   	    //if(ramo==0. || ramo_i==0. ) appo_outside=0.;  // there shouldn't be any charge induced if we are in a neighborhood pixel
		   	    				    		    // This term "should" account for that
		   	    				    		  
			    //if(fabs(xpos-x_neighbor)>pitchX/2 && fabs(ypos-y_neighbor)>pitchY/2) 	appo_outside=0.;	   	    				    
			    /*	  
			    if (isHole) pixelsContent[extraPixel] += appo_outside*(TrapeHit*average_charge + kappa*(eHitRamo - TrapeHit*average_charge)); //eV.
			    else pixelsContent[extraPixel] += appo_outside*(TrapeHit*(average_charge+charge_correction_exp) + kappa*(eHitRamo - TrapeHit*(average_charge+charge_correction_exp)));
			    */	  
			    //if (isHole) pixelsContent[extraPixel] += appo_outside*(eHit*average_charge + kappa*(eHitRamo - eHit*average_charge)); //eV.
			    //else pixelsContent[extraPixel] += appo_outside*(eHit*(average_charge+charge_correction_exp) + kappa*(eHitRamo - eHit*(average_charge+charge_correction_exp)));
			  } //loop over y
			} //loop over x
		      } //doRamo
		      
		      // Record deposit - take into account charge chunks that represent more than one fundamental charge!
		      //X' -> \mu + kappa * (X-\mu) 
		      //if (!doRamo) pixelsContent[extraPixel] += eHit*(charge_correction_exp+kappa*(0.-charge_correction_exp)); //eV
		     
		    } //is trapped
		    else { //charge was not trapped or charge trapping is turned off.
		        // Account for drifting into another pixel 
		        //G4cout<<"why here????"<<G4endl;
		        
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
		     
		       // Record deposit
		       //X' -> \mu + kappa * (X-\mu)
		        //pixelsContent[extraPixel] += eHit; // eV
		       
		       if(!isHole){		        
		        if (!doTrapping && !doRamo) pixelsContent[extraPixel] += eHit; // eV
		        else if (!doRamo) pixelsContent[extraPixel] += eHit*(charge_correction_exp+kappa*(1.-charge_correction_exp)); //eV
		        else pixelsContent[extraPixel] += eHit*(average_charge+charge_correction_exp) + kappa*(eHit - eHit*(average_charge+charge_correction_exp));
		      }else{
		        if (!doTrapping && !doRamo) pixelsContent[extraPixel] += eHit; // eV
		        else if (!doRamo) pixelsContent[extraPixel] += eHit*(charge_correction_exp+kappa*(1.-charge_correction_exp)); //eV
		        else pixelsContent[extraPixel] += eHit*(average_charge) + kappa*(eHit - eHit*(average_charge));
		       }
		       
		       
		       	      		    
		    }		  
		    
		      
		  }  // end loop over charges/holes
		} // end loop over nQ charges
map<pair<G4int, G4int>, G4double >::iterator pCItrSUM = sumpixelsContent.begin();
		    G4double sumdeposited_energy=0;
		    for( ; pCItrSUM != sumpixelsContent.end() ; pCItrSUM++)  sumdeposited_energy+= (*pCItrSUM).second;
	            //G4cout<<"  z " <<zpos<< " sumdeposited_energy "<<sumdeposited_energy <<" endepTOT "<<sumEnToT<<" frac "<< sumdeposited_energy/sumEnToT<<G4endl;
	            //G4cout<<"  zpos " <<zpos<< "  "<< sumdeposited_energy/eHitTotal<<G4endl;
		    
	map<pair<G4int, G4int>, G4double >::iterator pCItrDEL = sumpixelsContent.begin();
	for( ; pCItrDEL != sumpixelsContent.end() ; pCItrDEL++) sumpixelsContent[ (*pCItrDEL).first]=0;
		    
	}// end loop over nEntries
		    
	// Now that pixelContent is filled, create one digit per pixel
	map<pair<G4int, G4int>, G4double >::iterator pCItr = pixelsContent.begin();
	
	for( ; pCItr != pixelsContent.end() ; pCItr++)
	  {

	    G4double deposited_energy = (*pCItr).second;
	    int ToT = TMath::FloorNint(deposited_energy*tuning);
	    //G4cout<<"  X " << (*pCItr).first.first<< " Y "<< (*pCItr).first.second<< " EN " << (*pCItr).second <<" Threshold "<<threshold<<" tot "<< ToT<< G4endl;
	    if (ToT >=15) ToT = 15; //FEI4 is 4-bit.
	    AllPixFEI4RadDamageDigit * digit = new AllPixFEI4RadDamageDigit;
	    digit->SetPixelIDX((*pCItr).first.first);
	    digit->SetPixelIDY((*pCItr).first.second);
	    digit->SetPixelCounts(ToT);
	    digit->SetPixelEnergyDep(deposited_energy);
	    if (deposited_energy < threshold){
	      digit->SetPixelCounts(0);
	      digit->SetPixelEnergyDep(0);
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
	    count++;
	    G4cout << "--------> Digits Collection : " << collectionName[0]
		   << "(" << m_hitsColName[0] << ")"
		   << " contains " << dc_entries
		   << " digits" 
		   << " count "<<count<<  G4endl;
	  }
	
	StoreDigiCollection(m_digitsCollection);
	
	if (dodebug){
	  std::cout << "\n\n\n\n\n\n\n\n\n\n squirrel \n\n\n\n\n\n\n" << std::endl;
	}
	
} // end Digitize function





