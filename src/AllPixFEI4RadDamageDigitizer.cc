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
#include "TF2.h"
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
	doDiff = false;
	doSimplifiedModel = false;
	doChunkCorrection = false; //account for the fact that we are not propagating fundamental charge carriers.

	//Physics defaults
	defaultRamo = 0; // -1 = approximate z-dependence only; 0 = exact XZ x YZ with compensating Z-dependence; > 0 = solution to Poisson's equation with defaultRamo terms in the series.  See doc/RadDamageDefaults.pdf for details.
	defaultEfield = 1; // 0 = uniform; 1 = linear.
	defaultDiffusion = -1; //mm; if positive, will use a diffusion length of dif=defaultDiffusion*sqrt(zposD/0.3); N.B. the ATLAS digitization default is defaultDiffusion = 0.007.
	
	//Debugging information
	debug_maps = true; //if true, prints plots of all the input maps.
	dodebug = true; //prints debugging statements inside the digitizer

	//Constants
	elec = 3.64*eV;//average energy required to produce a e/h pair in Si.  This has been known for a long time - first measurement was Phys. Rev. 91, 1079 (1953).
	betaElectrons = 3.0E-16*cm2/ns;  //The charge-trapping probability is t = -tau*ln(u), for u ~ Uniform(0,1) and tau^{-1}=beta*fluence.  The value of beta might be slightly higher for holes than electrons, but it is hard to say (we ignore this effect here).  See e.g. https://cds.cern.ch/record/685542/files/indet-2003-014.pdf. 
	
	//Conditions
	fluence = 1e15*1/cm2; // 1 MeV neq/cm^2
	biasVoltage = 80; // V.  This is not used if external TCAD maps are supplied.
	temperature = 263.2;// K  
	bField = 0.;// Tesla = V*s/m^2 
	threshold = 2000*elec; //This is the threshold for charge collection.
	tuning = 9./(20000*elec); //for X ToT @ Y e, this is X/Y so that a deposited energy of Y gives a ToT of X.  Typical values are 5 ToT @ 20ke.

	//Phenomenological parameters
	precision = 10; //this is the number of charges to divide the G4 hit into.

	depletionLength=0.200;   //mm IBL 
	
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
	depVoltage = 60.; //V; this is highly fluece-dependent !  Should update when you change the fluence.
	depletionLength = L/1000.; //in mm.
	if(biasVoltage<depVoltage) depletionLength=depletionLength*pow(biasVoltage/depVoltage,0.5); //This is the usual formula depletion depth = \sqrt{2*\epsilon_0\epsilon_{Si}*V/(eN_D)}.  See e.g. (2.26) in Pixel Detectors by L. Rossi et al.
		
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
	    gPad->SetLeftMargin(0.15);
	    ramoPotentialMap2D->GetXaxis()->SetNdivisions(505);
	    ramoPotentialMap2D->GetYaxis()->SetNdivisions(505);
	    ramoPotentialMap2D->GetXaxis()->SetTitle("#eta position [#mum]");
	    ramoPotentialMap2D->GetYaxis()->SetTitle("#phi position [#mum]");
	    ramoPotentialMap2D->GetZaxis()->SetTitleOffset(1.4);
	    ramoPotentialMap2D->GetYaxis()->SetTitleOffset(1.5);
	    ramoPotentialMap2D->GetXaxis()->SetTitleOffset(1.3);
	    ramoPotentialMap2D->SetTitle("Ramo Potential at z = "+TString::Format("%0.2f",ramoPotentialMap->GetZaxis()->GetBinCenter(k)-0.5*ramoPotentialMap->GetZaxis()->GetBinWidth(k))+" mm");
	    ramoPotentialMap2D->GetZaxis()->SetTitle("Ramo Potential");
	    ramoPotentialMap2D->GetZaxis()->SetRangeUser(0,1);
	    ramoPotentialMap2D->Draw("colz");
	    c1->Print("DefaultRamo2D_"+TString::Format("%i",k)+".pdf");
	  }
	  gPad->SetRightMargin(0.05);
	  ramoPotentialMap1D->SetTitle("Ramo Potential at x=y=0");
	  ramoPotentialMap1D->GetYaxis()->SetRangeUser(0,1.2);
	  ramoPotentialMap1D->GetXaxis()->SetTitle("Pixel Depth in z [#mum]");
	  ramoPotentialMap1D->GetYaxis()->SetTitleOffset(1.4);
	  ramoPotentialMap1D->GetYaxis()->SetTitle("Ramo Potential");
	  ramoPotentialMap1D->GetXaxis()->SetNdivisions(505);
	  ramoPotentialMap1D->Draw();
	  c1->Print("DefaultRamo1D.pdf");
	}

	//First, attempt to get a 3D E field. MAP IS GIVEN IN UNITS OF MICRONS AND V/cm.
	eFieldMap=0;
	eFieldMap=(TH3F*)efile->Get("hefield3d");
	Efield3D = true;
	if (eFieldMap == 0){
	  Efield3D = false;
	  m_eFieldMap1D=0;
	  m_eFieldMap1D=(TH1F*)efile->Get("hefieldz");
	  if (m_eFieldMap1D == 0){ //Note that z=0 corresponds to the collecting electrode.
	    G4cout << "Did not find an E-field map.  Will use an approximate form." << G4endl;
	    Efield3D = false;
	    m_eFieldMap1D = new TH1F("m_hefieldz","m_hefieldz",L_int,0,L_int); // 1D map of the Efield in case we want to use 1D field
	    G4double electricField=0;
	    G4double pass=L/1000/m_eFieldMap1D->GetNbinsX();
	    //G4double pass=depletionLength/m_eFieldMap1D->GetNbinsX();
	    for (int i=1; i<= m_eFieldMap1D->GetNbinsX()+1; i++){
	      G4double position=pass*i;
	      if (depletionLength != 0) {
			  if(biasVoltage<depVoltage)	electricField= 2*(biasVoltage/depletionLength)*(1-position/depletionLength);
			  if(biasVoltage>=depVoltage)  electricField=(2*depVoltage/depletionLength)*(1-position/depletionLength)+(biasVoltage-depVoltage)/(depletionLength);
			  if(position>depletionLength) electricField=0.;
		  }
			  
	      m_eFieldMap1D->SetBinContent(m_eFieldMap1D->GetNbinsX()-i+1,electricField*10); //in V/cm; n in n prior to type inversion.
	      if (defaultEfield==0) m_eFieldMap1D->SetBinContent(i,(1e4)*biasVoltage/(L/1000.)); //in V/cm
	    }
	  } 
	}	  

	if (debug_maps){
	  gPad->SetLeftMargin(0.15);
	  TFile *eOutFile = new TFile("efield_1D.root", "RECREATE");
	  if (Efield3D){
	    m_eFieldMap1D = new TH1F("hefieldz","hefieldz",eFieldMap->GetNbinsZ(),eFieldMap->GetZaxis()->GetBinCenter(1)-0.5*eFieldMap->GetZaxis()->GetBinWidth(1),eFieldMap->GetZaxis()->GetBinCenter(eFieldMap->GetNbinsZ())+0.5*eFieldMap->GetZaxis()->GetBinWidth(eFieldMap->GetNbinsZ()));
	    for (int k=1; k<=eFieldMap->GetNbinsZ(); k++){
	      double avg_num = 0.;
	      double avg_den = 0.;
	      for (int i=1; i<=eFieldMap->GetNbinsX(); i++){
          for (int j=1; j<=eFieldMap->GetNbinsY(); j++){
            avg_num+=eFieldMap->GetBinContent(i,j,k);
            avg_den+=1.;
          }
	      }
	      m_eFieldMap1D->SetBinContent(k,avg_num/avg_den);
	    }
	  }
	  m_eFieldMap1D->GetYaxis()->SetRangeUser(0.,m_eFieldMap1D->GetMaximum()*1.5);
	  m_eFieldMap1D->SetTitle("Electric Field");
	  m_eFieldMap1D->GetXaxis()->SetTitle("Depth (z) [#mum]");
	  m_eFieldMap1D->GetYaxis()->SetTitle("E field [V/cm]");
	  m_eFieldMap1D->GetYaxis()->SetTitleOffset(2.3);
	  m_eFieldMap1D->GetXaxis()->SetTitleOffset(1.3);
	  m_eFieldMap1D->GetYaxis()->SetTitleSize(0.03);
	  m_eFieldMap1D->GetXaxis()->SetTitleSize(0.03);
	  m_eFieldMap1D->Draw();
	  c1->Print("Efield.pdf");
	  m_eFieldMap1D->Write("efield_1D.root");
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
	  
	  // filling distancemap_e,h
	  /**
	  for (int i=1; i<= distancemap_e->GetNbinsX(); i++){
	    for (int j=1; j<= distancemap_e->GetNbinsY(); j++){
	      //distancemap_e->SetBinContent(i,j,L/1000.); //if you travel long enough, you will reach the electrode.
	      distancemap_h->SetBinContent(i,j,L/1000.); //if you travel long enough, you will reach the electrode.
	    }
	  }**/
	  
	  // filling timeMap_e,h
	  for (int k=1; k<= distancemap_e->GetNbinsX(); k++){
	    double z = distancemap_e->GetXaxis()->GetBinCenter(k);
	    double mysum = 0.;
	    double mysum_h = 0.;
	    for (int k2=k; k2 >= 1; k2--){
	      double z2 = distancemap_e->GetXaxis()->GetBinCenter(k2);
	      double dz = distancemap_e->GetXaxis()->GetBinWidth(k2);
	      double E = Efield3D ? GetElectricField(0.002,0.002,z2) : m_eFieldMap1D->GetBinContent(m_eFieldMap1D->GetXaxis()->FindBin(z2*1000))/1e7; //in MV/mm
	      double mu = GetMobility(E, 0); //mm^2/MV*ns
	      if (E > 0){
          mysum+=dz/(mu*E); //mm * 1/(mm/ns) = ns
          distancemap_e->SetBinContent(k,distancemap_e->GetYaxis()->FindBin(mysum),z2);
	      } else if (E == 0) {
			  mysum = 3.40282e+38; // without efield: travel almost forever to reach the electrode since we don't do actual diffusion
			  distancemap_e->SetBinContent(k,distancemap_e->GetYaxis()->FindBin(mysum),z2);
        }
	    }
	    timeMap_e->SetBinContent(k,mysum);
	    
	    for (int k2=k; k2 <= distancemap_e->GetNbinsX(); k2++){ //holes go the opposite direction as electrons.
	      double z2 = distancemap_e->GetXaxis()->GetBinCenter(k2);
	      double dz = distancemap_e->GetXaxis()->GetBinWidth(k2);
	      double E = Efield3D ? GetElectricField(0.002,0.002,z2) : m_eFieldMap1D->GetBinContent(m_eFieldMap1D->GetXaxis()->FindBin(z2*1000))/1e7; //in MV/mm
	      double mu_h = GetMobility(E, 1);
	      if (E > 0){
          mysum_h+=dz/(mu_h*E);
          distancemap_h->SetBinContent(k,distancemap_h->GetYaxis()->FindBin(mysum_h),z2);
	      } else if (E == 0) {
				mysum_h = 3.40282e+38; // without efield: travel almost forever to reach the electrode since we don't do actual diffusion
				distancemap_h->SetBinContent(k,distancemap_h->GetYaxis()->FindBin(mysum_h),z2);
        }
	    }
	    timeMap_h->SetBinContent(k,mysum_h);
	  }
	  
	  FillHoles(distancemap_e, isHole, L/1000., "distancefit_e");
	  FillHoles(distancemap_h, isHole, L/1000., "distancefit_h");
		
	}//end generating maps

	if (debug_maps){
		c1->cd();
	  gPad->SetRightMargin(0.15);
	  distancemap_e->GetXaxis()->SetNdivisions(505);
	  distancemap_e->GetYaxis()->SetNdivisions(505);
	  distancemap_e->GetXaxis()->SetTitle("Initial Position in z [mm]");
	  distancemap_e->GetYaxis()->SetTitle("Time Traveled [ns]");
	  distancemap_e->GetZaxis()->SetTitleOffset(1.6);
	  distancemap_e->SetTitle("Electron Distance Map");
	  distancemap_e->GetZaxis()->SetTitle("Location in z [mm]");
	  distancemap_e->DrawCopy("colz");
	  c1->Print("distancemap_e.pdf");
	  
	  distancemap_h->GetXaxis()->SetNdivisions(505);
	  distancemap_h->GetYaxis()->SetNdivisions(505);
	  distancemap_h->GetXaxis()->SetTitle("Initial Position in z [mm]");
	  distancemap_h->GetYaxis()->SetTitle("Time Traveled [ns]");
	  distancemap_h->GetZaxis()->SetTitleOffset(1.6);
	  distancemap_h->SetTitle("Hole Distance Map");
	  distancemap_h->GetZaxis()->SetTitle("Location in z [mm]");
	  distancemap_h->DrawCopy("colz");
	  c1->Print("distancemap_h.pdf");
	  
	  timeMap_h->SetTitle("");
	  timeMap_h->GetYaxis()->SetRangeUser(0,timeMap_h->GetMaximum() > timeMap_e->GetMaximum() ? timeMap_h->GetMaximum() : timeMap_e->GetMaximum());
	  timeMap_h->GetXaxis()->SetTitle("Starting Pixel Depth in z [mm]");
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
	if(fluence!=0.0) {
	    trappingTimeElectrons = 1.0/(betaElectrons*fluence);
	    trappingTimeHoles = 1.0/(betaElectrons*fluence);
  }
	else { //fluence = 0 so do not trap!
	  trappingTimeElectrons = 1000*s; //~infinity
	  trappingTimeHoles = 1000*s;
	}
	
        //Need to pre-compute a new map that allows us to properly take into account the charge chunking. 
	int stepsZ = int(L/2); //one bin per 10 microns
	int stepsY = int(pitchY*1000/2); //half width of sensor to take advantage of symmetry.
	int stepsX = int(pitchX*1000/2);
	charge_chunk_map_e = new TH3F("","",stepsX,0.,pitchX/2,stepsY,0.,pitchY/2,stepsZ,0.,detectorThickness); //distances are in mm.
	charge_chunk_map_h = new TH3F("","",stepsX,0.,pitchX/2,stepsY,0.,pitchY/2,stepsZ,0.,detectorThickness); 
  for (int k = 1; k <= stepsZ; k++){
    double z =  charge_chunk_map_e->GetZaxis()->GetBinCenter(k);
      for (int i=1; i<= stepsX; i++){
        double x = charge_chunk_map_e->GetXaxis()->GetBinCenter(i);
        for (int j=1; j<= stepsY; j++){
          double y= charge_chunk_map_e->GetYaxis()->GetBinCenter(j);
          G4double timeToElectrode = GetTimeToElectrode(z, 0);
          G4double timeToElectrode_h = GetTimeToElectrode(z, 1);
          double ramo_intial = ramoPotentialMap->GetBinContent(ramoPotentialMap->FindBin(fabs(x*1000),fabs(y*1000),z*1000)); //fabs because the Ramo is only 1/4 of the pixel.
          double prob = exp(-timeToElectrode/trappingTimeElectrons)*(1.-ramo_intial);
          double prob_h = exp(-timeToElectrode_h/trappingTimeHoles)*ramo_intial;
          for (double t = 0; t <= timeToElectrode; t+=timeToElectrode/stepsZ){
            double final_pos = distancemap_e->GetBinContent(distancemap_e->FindBin(z,t)); //with respect to the electrode (at 0)
            double ramo_final = ramoPotentialMap->GetBinContent(ramoPotentialMap->FindBin(fabs(x*1000),fabs(y*1000),final_pos*1000));
            prob+=(ramo_final-ramo_intial)*(timeToElectrode/stepsZ)*exp(-t/trappingTimeElectrons)/trappingTimeElectrons;
          } 
          for (double t = 0; t <= timeToElectrode_h; t+=timeToElectrode_h/stepsZ){
            double final_pos_h = distancemap_h->GetBinContent(distancemap_h->FindBin(z,t));
            double ramo_final_h = ramoPotentialMap->GetBinContent(ramoPotentialMap->FindBin(fabs(x*1000),fabs(y*1000),final_pos_h*1000));
            prob_h+=-(ramo_final_h-ramo_intial)*(timeToElectrode_h/stepsZ)*exp(-t/trappingTimeHoles)/trappingTimeHoles;
          }
          charge_chunk_map_e->SetBinContent(charge_chunk_map_e->FindBin(x,y,z),prob);
          charge_chunk_map_h->SetBinContent(charge_chunk_map_h->FindBin(x,y,z),prob_h);
      }
    }
  }
	
	if (debug_maps){
	  //as with the Ramo potential, let's get 2D and 1D maps here.
	  TH2F* charge_chunk_map_e2D = new TH2F("charge_chunk_map_e2D","charge_chunk_map_e2D",charge_chunk_map_e->GetNbinsX(),charge_chunk_map_e->GetXaxis()->GetBinCenter(1)-0.5*charge_chunk_map_e->GetXaxis()->GetBinWidth(1),charge_chunk_map_e->GetXaxis()->GetBinCenter(charge_chunk_map_e->GetNbinsX())+0.5*charge_chunk_map_e->GetXaxis()->GetBinWidth(charge_chunk_map_e->GetNbinsX()),charge_chunk_map_e->GetNbinsY(),charge_chunk_map_e->GetYaxis()->GetBinCenter(1)-0.5*charge_chunk_map_e->GetYaxis()->GetBinWidth(1),charge_chunk_map_e->GetYaxis()->GetBinCenter(charge_chunk_map_e->GetNbinsY())+0.5*charge_chunk_map_e->GetYaxis()->GetBinWidth(charge_chunk_map_e->GetNbinsY()));
    TH1F* charge_chunk_map_e1D = new TH1F("charge_chunk_map_e1D","charge_chunk_map_e1D",charge_chunk_map_e->GetNbinsZ(),charge_chunk_map_e->GetZaxis()->GetBinCenter(1)-0.5*charge_chunk_map_e->GetZaxis()->GetBinWidth(1),charge_chunk_map_e->GetZaxis()->GetBinCenter(charge_chunk_map_e->GetNbinsZ())+0.5*charge_chunk_map_e->GetZaxis()->GetBinWidth(charge_chunk_map_e->GetNbinsZ()));
	  TH2F* charge_chunk_map_h2D = new TH2F("charge_chunk_map_h2D","charge_chunk_map_h2D",charge_chunk_map_h->GetNbinsX(),charge_chunk_map_h->GetXaxis()->GetBinCenter(1)-0.5*charge_chunk_map_h->GetXaxis()->GetBinWidth(1),charge_chunk_map_h->GetXaxis()->GetBinCenter(charge_chunk_map_h->GetNbinsX())+0.5*charge_chunk_map_h->GetXaxis()->GetBinWidth(charge_chunk_map_h->GetNbinsX()),charge_chunk_map_h->GetNbinsY(),charge_chunk_map_h->GetYaxis()->GetBinCenter(1)-0.5*charge_chunk_map_h->GetYaxis()->GetBinWidth(1),charge_chunk_map_h->GetYaxis()->GetBinCenter(charge_chunk_map_h->GetNbinsY())+0.5*charge_chunk_map_h->GetYaxis()->GetBinWidth(charge_chunk_map_h->GetNbinsY()));
	  TH1F* charge_chunk_map_h1D = new TH1F("charge_chunk_map_h1D","charge_chunk_map_h1D",charge_chunk_map_h->GetNbinsZ(),charge_chunk_map_h->GetZaxis()->GetBinCenter(1)-0.5*charge_chunk_map_h->GetZaxis()->GetBinWidth(1),charge_chunk_map_h->GetZaxis()->GetBinCenter(charge_chunk_map_h->GetNbinsZ())+0.5*charge_chunk_map_h->GetZaxis()->GetBinWidth(charge_chunk_map_h->GetNbinsZ()));
	  TH1F* charge_chunk_map_1D = (TH1F*)charge_chunk_map_h1D->Clone("sum");
    for (int k=1; k<=charge_chunk_map_e->GetNbinsZ(); k++){
      for (int i=1; i<=charge_chunk_map_e->GetNbinsX(); i++){
        for (int j=1; j<=charge_chunk_map_e->GetNbinsY(); j++){
          charge_chunk_map_e2D->SetBinContent(i,j,charge_chunk_map_e->GetBinContent(i,j,k));
          charge_chunk_map_h2D->SetBinContent(i,j,charge_chunk_map_h->GetBinContent(i,j,k));
          if (i==1 && j==1){
            charge_chunk_map_e1D->SetBinContent(k,charge_chunk_map_e->GetBinContent(i,j,k));
            charge_chunk_map_h1D->SetBinContent(k,charge_chunk_map_h->GetBinContent(i,j,k));
            charge_chunk_map_1D->SetBinContent(k,charge_chunk_map_e1D->GetBinContent(k)+charge_chunk_map_h1D->GetBinContent(k));
          }
        }
      }
      gPad->SetRightMargin(0.15);
      charge_chunk_map_e2D->GetXaxis()->SetNdivisions(505);
      charge_chunk_map_e2D->GetYaxis()->SetNdivisions(505);
      charge_chunk_map_e2D->GetXaxis()->SetTitle("#eta position [#mum]");
      charge_chunk_map_e2D->GetYaxis()->SetTitle("#phi position [#mum]");
      charge_chunk_map_e2D->GetZaxis()->SetTitleOffset(1.4);
	    charge_chunk_map_e2D->GetYaxis()->SetTitleOffset(1.8);
	    charge_chunk_map_e2D->GetXaxis()->SetTitleOffset(1.3);
      charge_chunk_map_e2D->SetTitle("Average e induced charge at z = "+TString::Format("%0.2f",charge_chunk_map_e->GetZaxis()->GetBinCenter(k)-0.5*charge_chunk_map_e->GetZaxis()->GetBinWidth(k))+" mm");
      charge_chunk_map_e2D->GetZaxis()->SetTitle("Average Fractional Induced Charge");
      charge_chunk_map_e2D->GetZaxis()->SetRangeUser(0,1);
      charge_chunk_map_e2D->Draw("colz");
      c1->Print("charge_chunk_map_e2D_"+TString::Format("%i",k)+".pdf");
	    
	    charge_chunk_map_h2D->GetXaxis()->SetNdivisions(505);
      charge_chunk_map_h2D->GetYaxis()->SetNdivisions(505);
      charge_chunk_map_h2D->GetXaxis()->SetTitle("#eta position [#mum]");
      charge_chunk_map_h2D->GetYaxis()->SetTitle("#phi position [#mum]");
	    charge_chunk_map_h2D->GetXaxis()->SetTitleOffset(1.3);
	    charge_chunk_map_h2D->GetYaxis()->SetTitleOffset(1.8);
      charge_chunk_map_h2D->GetZaxis()->SetTitleOffset(1.4);
      charge_chunk_map_h2D->SetTitle("Average h induced charge at z = "+TString::Format("%0.2f",charge_chunk_map_e->GetZaxis()->GetBinCenter(k)-0.5*charge_chunk_map_e->GetZaxis()->GetBinWidth(k))+" mm");
      charge_chunk_map_h2D->GetZaxis()->SetTitle("Average Fractional Induced Charge");
      charge_chunk_map_h2D->GetZaxis()->SetRangeUser(0,1);
      charge_chunk_map_h2D->Draw("colz");
      c1->Print("charge_chunk_map_h2D_"+TString::Format("%i",k)+".pdf");
    }
    gPad->SetRightMargin(0.05);
    charge_chunk_map_e1D->SetTitle("Average induced charge at x=y=0");
    charge_chunk_map_e1D->GetYaxis()->SetRangeUser(0,1);
    charge_chunk_map_e1D->GetXaxis()->SetTitle("Pixel Depth in z [#mum]");
    charge_chunk_map_e1D->GetYaxis()->SetTitleOffset(1.4);
    charge_chunk_map_e1D->GetYaxis()->SetTitle("Average Fractional Induced Charge");
    charge_chunk_map_e1D->GetXaxis()->SetNdivisions(505);
	  charge_chunk_map_e1D->SetLineColor(2);
	  charge_chunk_map_e1D->GetYaxis()->SetRangeUser(0,1.2);
    charge_chunk_map_e1D->Draw();
	  charge_chunk_map_h1D->SetLineColor(4);
	  charge_chunk_map_h1D->SetLineStyle(3);
	  charge_chunk_map_h1D->Draw("same");
	  charge_chunk_map_1D->SetLineColor(1);
	  charge_chunk_map_1D->Draw("same");

	  TLegend * leg = new TLegend(0.55,0.35,0.9,0.6);
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->AddEntry(charge_chunk_map_e1D,"Electrons","l");
    leg->AddEntry(charge_chunk_map_h1D,"Holes","l");
	  leg->AddEntry(charge_chunk_map_1D,"Total","l");
    leg->Draw();
	  
    c1->Print("charge_chunk_map_1D.pdf");
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
	      double E = Efield3D ? GetElectricField(0.002,0.002,z2) : m_eFieldMap1D->GetBinContent(m_eFieldMap1D->GetXaxis()->FindBin(z2*1000))/1e7; //in MV/mm
	      if (E > 0){
          lenght_e+=dz;
          // integral tan thetaL * dz / integral dz
          double tanthetaL = GetTanLorentz(E, 0);	 
          mysum+=(tanthetaL*dz);				  // the angle is computed as the integral over the path
          lorentz_map_e->SetBinContent(k,k2,mysum/lenght_e);
	      }
	    }
	    for (int k2=k; k2 <= lorentz_map_e->GetNbinsX(); k2++){ //holes go the opposite direction as electrons.
	      double z2 = lorentz_map_e->GetXaxis()->GetBinCenter(k2);
	      double dz = lorentz_map_e->GetXaxis()->GetBinWidth(k2);
	      double E = Efield3D ? GetElectricField(0.002,0.002,z2) : m_eFieldMap1D->GetBinContent(m_eFieldMap1D->GetXaxis()->FindBin(z2*1000))/1e7; //in MV/mm
	      if (E > 0){
          double tanthetaL = GetTanLorentz(E, 1);
          lenght_h+=dz;
          mysum_h+=(tanthetaL *dz);
          lorentz_map_h->SetBinContent(k,k2,mysum_h/lenght_h);
	      }
	    }
	  }
	}      
	
	if (debug_maps){
	  gPad->SetRightMargin(0.18);
	  lorentz_map_e->GetXaxis()->SetNdivisions(505);
	  lorentz_map_e->GetYaxis()->SetNdivisions(505);
	  lorentz_map_e->GetXaxis()->SetTitle("Initial Position in z [mm]");
	  lorentz_map_e->GetYaxis()->SetTitle("Distance Traveled in z [mm]");
	  lorentz_map_e->GetZaxis()->SetTitleOffset(1.7);
	  lorentz_map_e->GetYaxis()->SetTitleOffset(1.5);
	  lorentz_map_e->GetXaxis()->SetTitleOffset(1.3);
	  lorentz_map_e->SetTitle("Electron Lorentz Map");
	  lorentz_map_e->GetZaxis()->SetTitle("Tangent Lorentz Angle");
	  lorentz_map_e->Draw("colz");
	  c1->Print("lorentz_map_e.pdf");
	  
	  lorentz_map_h->GetXaxis()->SetNdivisions(505);
	  lorentz_map_h->GetYaxis()->SetNdivisions(505);
	  lorentz_map_h->GetXaxis()->SetTitle("Initial Position in z [mm]");
	  lorentz_map_h->GetYaxis()->SetTitle("Distance Traveled in z [mm]");
	  lorentz_map_h->GetYaxis()->SetTitleOffset(1.5);     
	  lorentz_map_h->GetXaxis()->SetTitleOffset(1.3);
	  lorentz_map_h->GetZaxis()->SetTitleOffset(1.7);
	  lorentz_map_h->SetTitle("Hole Lorentz Map");
	  lorentz_map_h->GetZaxis()->SetTitle("Tangent Lorentz Angle");
	  lorentz_map_h->Draw("colz");
	  c1->Print("lorentz_map_h.pdf");
	}

	if (dodebug){
	  debug_inducedcharge_z_versus_time_00_num = new TH2F("","",20,0,L,20,0,10);
	  debug_inducedcharge_z_versus_time_00_den = new TH2F("","",20,0,L,20,0,10);
	  debug_inducedcharge_z_versus_time_00 = new TH2F("","",20,0,L,20,0,10);

	  debug_inducedcharge_z_versus_time_00_num_holes = new TH2F("","",20,0,L,20,0,10);
    debug_inducedcharge_z_versus_time_00_den_holes = new TH2F("","",20,0,L,20,0,10);
    debug_inducedcharge_z_versus_time_00_holes = new TH2F("","",20,0,L,20,0,10);

	  debug_inducedcharge_z_versus_time_01_num = new TH2F("","",20,0,L,20,0,10);
	  debug_inducedcharge_z_versus_time_01_den = new TH2F("","",20,0,L,20,0,10);
    debug_inducedcharge_z_versus_time_01 = new TH2F("","",20,0,L,20,0,10);

	  debug_inducedcharge_z_versus_time_01_num_holes = new TH2F("","",20,0,L,20,0,10);
    debug_inducedcharge_z_versus_time_01_den_holes = new TH2F("","",20,0,L,20,0,10);
    debug_inducedcharge_z_versus_time_01_holes = new TH2F("","",20,0,L,20,0,10);

	  debug_inducedcharge_z_versus_time_10_num = new TH2F("","",20,0,L,20,0,10);
	  debug_inducedcharge_z_versus_time_10_den = new TH2F("","",20,0,L,20,0,10);
    debug_inducedcharge_z_versus_time_10 = new TH2F("","",20,0,L,20,0,10);

	  debug_inducedcharge_z_versus_time_10_num_holes = new TH2F("","",20,0,L,20,0,10);
    debug_inducedcharge_z_versus_time_10_den_holes = new TH2F("","",20,0,L,20,0,10);
    debug_inducedcharge_z_versus_time_10_holes = new TH2F("","",20,0,L,20,0,10);

	  debug_inducedcharge_versus_time_nocorr = new TProfile("","",25,0,0.2,0.,1.2,"s");
	  debug_inducedcharge_versus_time_corr = new TProfile("","",25,0.0,0.2,0.,1.2,"s");
	  debug_inducedcharge_versus_time_nocorr_holes = new TProfile("","",25,0,0.2,0.,1.2,"s");
    debug_inducedcharge_versus_time_corr_holes = new TProfile("","",25,0.0,0.2,0.,1.2,"s");

	  debug_chunksize = new TH1F("","",100,0,100);
	}

}// End of AllPixFEI4RadDamageDigitizer::AllPixFEI4RadDamageDigitizer definition

AllPixFEI4RadDamageDigitizer::~AllPixFEI4RadDamageDigitizer(){

  if (dodebug){
    TCanvas *c1 = new TCanvas("","",600,600);
    gStyle->SetOptStat(0);
    for (int i=1; i<=debug_inducedcharge_z_versus_time_00_num->GetNbinsX(); i++){
      for (int j=1; j<=debug_inducedcharge_z_versus_time_00_num->GetNbinsY(); j++){
        if (debug_inducedcharge_z_versus_time_00_den->GetBinContent(i,j) > 0) debug_inducedcharge_z_versus_time_00->SetBinContent(i,j,debug_inducedcharge_z_versus_time_00_num->GetBinContent(i,j)/debug_inducedcharge_z_versus_time_00_den->GetBinContent(i,j));
        else debug_inducedcharge_z_versus_time_00->SetBinContent(i,j,0.);

        if (debug_inducedcharge_z_versus_time_01_den->GetBinContent(i,j) > 0) debug_inducedcharge_z_versus_time_01->SetBinContent(i,j,debug_inducedcharge_z_versus_time_01_num->GetBinContent(i,j)/debug_inducedcharge_z_versus_time_01_den->GetBinContent(i,j));
        else debug_inducedcharge_z_versus_time_01->SetBinContent(i,j,0.);

        if (debug_inducedcharge_z_versus_time_10_den->GetBinContent(i,j) > 0) debug_inducedcharge_z_versus_time_10->SetBinContent(i,j,debug_inducedcharge_z_versus_time_10_num->GetBinContent(i,j)/debug_inducedcharge_z_versus_time_10_den->GetBinContent(i,j));
        else debug_inducedcharge_z_versus_time_10->SetBinContent(i,j,0.);

        //now for holes
        if (debug_inducedcharge_z_versus_time_00_den_holes->GetBinContent(i,j) > 0) debug_inducedcharge_z_versus_time_00_holes->SetBinContent(i,j,debug_inducedcharge_z_versus_time_00_num_holes->GetBinContent(i,j)/debug_inducedcharge_z_versus_time_00_den_holes->GetBinContent(i,j));
        else debug_inducedcharge_z_versus_time_00_holes->SetBinContent(i,j,0.);

        if (debug_inducedcharge_z_versus_time_01_den_holes->GetBinContent(i,j) > 0) debug_inducedcharge_z_versus_time_01_holes->SetBinContent(i,j,debug_inducedcharge_z_versus_time_01_num_holes->GetBinContent(i,j)/debug_inducedcharge_z_versus_time_01_den_holes->GetBinContent(i,j));
        else debug_inducedcharge_z_versus_time_01_holes->SetBinContent(i,j,0.);

        if (debug_inducedcharge_z_versus_time_10_den_holes->GetBinContent(i,j) > 0) debug_inducedcharge_z_versus_time_10_holes->SetBinContent(i,j,debug_inducedcharge_z_versus_time_10_num_holes->GetBinContent(i,j)/debug_inducedcharge_z_versus_time_10_den_holes->GetBinContent(i,j));
        else debug_inducedcharge_z_versus_time_10_holes->SetBinContent(i,j,0.);
	
      }
    }
    gPad->SetRightMargin(0.18);
    debug_inducedcharge_z_versus_time_00->GetZaxis()->SetRangeUser(0,1);
    debug_inducedcharge_z_versus_time_00->GetXaxis()->SetNdivisions(505);
    debug_inducedcharge_z_versus_time_00->GetYaxis()->SetNdivisions(505);
    debug_inducedcharge_z_versus_time_00->GetXaxis()->SetTitle("Initial Position in z [#mum]");
    debug_inducedcharge_z_versus_time_00->GetYaxis()->SetTitle("Electrons Time Traveled [ns]");
    debug_inducedcharge_z_versus_time_00->GetZaxis()->SetTitleOffset(1.7);
    debug_inducedcharge_z_versus_time_00->GetYaxis()->SetTitleOffset(1.4);
    debug_inducedcharge_z_versus_time_00->GetXaxis()->SetTitleOffset(1.3);
    debug_inducedcharge_z_versus_time_00->SetTitle("Center pixel");
    debug_inducedcharge_z_versus_time_00->GetZaxis()->SetTitle("Average Induced Charge / e");
    debug_inducedcharge_z_versus_time_00->Draw("colz");
    c1->Print("debug_inducedcharge_z_versus_time_00.pdf");

    debug_inducedcharge_z_versus_time_00_holes->GetZaxis()->SetRangeUser(0,1);
    debug_inducedcharge_z_versus_time_00_holes->GetXaxis()->SetNdivisions(505);
    debug_inducedcharge_z_versus_time_00_holes->GetYaxis()->SetNdivisions(505);
    debug_inducedcharge_z_versus_time_00_holes->GetXaxis()->SetTitle("Initial Position in z [#mum]");
    debug_inducedcharge_z_versus_time_00_holes->GetYaxis()->SetTitle("Holes Time Traveled [ns]");
    debug_inducedcharge_z_versus_time_00_holes->GetZaxis()->SetTitleOffset(1.7);
    debug_inducedcharge_z_versus_time_00_holes->GetYaxis()->SetTitleOffset(1.4);
    debug_inducedcharge_z_versus_time_00_holes->GetXaxis()->SetTitleOffset(1.3);
    debug_inducedcharge_z_versus_time_00_holes->SetTitle("Center pixel");
    debug_inducedcharge_z_versus_time_00_holes->GetZaxis()->SetTitle("Average Induced Charge / e");
    debug_inducedcharge_z_versus_time_00_holes->Draw("colz");
    c1->Print("debug_inducedcharge_z_versus_time_00_holes.pdf");

    debug_inducedcharge_z_versus_time_01->GetZaxis()->SetRangeUser(-0.2,0.2);
    debug_inducedcharge_z_versus_time_01->GetXaxis()->SetNdivisions(505);
    debug_inducedcharge_z_versus_time_01->GetYaxis()->SetNdivisions(505);
    debug_inducedcharge_z_versus_time_01->GetXaxis()->SetTitle("Initial Position in z [#mum]");
    debug_inducedcharge_z_versus_time_01->GetYaxis()->SetTitle("Electrons Time Traveled [ns]");
    debug_inducedcharge_z_versus_time_01->GetZaxis()->SetTitleOffset(1.7);
    debug_inducedcharge_z_versus_time_01->GetYaxis()->SetTitleOffset(1.4);
    debug_inducedcharge_z_versus_time_01->GetXaxis()->SetTitleOffset(1.3);
    debug_inducedcharge_z_versus_time_01->SetTitle("+1 pixel in #phi");
    debug_inducedcharge_z_versus_time_01->GetZaxis()->SetTitle("Average Induced Charge (electrons) / e");
    debug_inducedcharge_z_versus_time_01->Draw("colz");
    c1->Print("debug_inducedcharge_z_versus_time_01.pdf");

    debug_inducedcharge_z_versus_time_01_holes->GetZaxis()->SetRangeUser(-0.2,0.2);
    debug_inducedcharge_z_versus_time_01_holes->GetXaxis()->SetNdivisions(505);
    debug_inducedcharge_z_versus_time_01_holes->GetYaxis()->SetNdivisions(505);
    debug_inducedcharge_z_versus_time_01_holes->GetXaxis()->SetTitle("Initial Position in z [#mum]");
    debug_inducedcharge_z_versus_time_01_holes->GetYaxis()->SetTitle("Holes Time Traveled [ns]");
    debug_inducedcharge_z_versus_time_01_holes->GetZaxis()->SetTitleOffset(1.7);
    debug_inducedcharge_z_versus_time_01_holes->GetYaxis()->SetTitleOffset(1.4);
    debug_inducedcharge_z_versus_time_01_holes->GetXaxis()->SetTitleOffset(1.3);
    debug_inducedcharge_z_versus_time_01_holes->SetTitle("+1 pixel in #phi");
    debug_inducedcharge_z_versus_time_01_holes->GetZaxis()->SetTitle("Average Induced Charge (holes) / e");
    debug_inducedcharge_z_versus_time_01_holes->Draw("colz");
    c1->Print("debug_inducedcharge_z_versus_time_01_holes.pdf");

    debug_inducedcharge_z_versus_time_10->GetZaxis()->SetRangeUser(-0.2,0.2);
    debug_inducedcharge_z_versus_time_10->GetXaxis()->SetNdivisions(505);
    debug_inducedcharge_z_versus_time_10->GetYaxis()->SetNdivisions(505);
    debug_inducedcharge_z_versus_time_10->GetXaxis()->SetTitle("Initial Position in z [#mum]");
    debug_inducedcharge_z_versus_time_10->GetYaxis()->SetTitle("Electrons Time Traveled [ns]");
    debug_inducedcharge_z_versus_time_10->GetZaxis()->SetTitleOffset(1.7);
    debug_inducedcharge_z_versus_time_10->GetYaxis()->SetTitleOffset(1.4);
    debug_inducedcharge_z_versus_time_10->GetXaxis()->SetTitleOffset(1.3);
    debug_inducedcharge_z_versus_time_10->SetTitle("+1 pixel in #eta");
    debug_inducedcharge_z_versus_time_10->GetZaxis()->SetTitle("Average Induced Charge (electrons) / e");
    debug_inducedcharge_z_versus_time_10->Draw("colz");
    c1->Print("debug_inducedcharge_z_versus_time_10.pdf");

    debug_inducedcharge_z_versus_time_10_holes->GetZaxis()->SetRangeUser(-0.2,0.2);
    debug_inducedcharge_z_versus_time_10_holes->GetXaxis()->SetNdivisions(505);
    debug_inducedcharge_z_versus_time_10_holes->GetYaxis()->SetNdivisions(505);
    debug_inducedcharge_z_versus_time_10_holes->GetXaxis()->SetTitle("Initial Position in z [#mum]");
    debug_inducedcharge_z_versus_time_10_holes->GetYaxis()->SetTitle("Electrons Time Traveled [ns]");
    debug_inducedcharge_z_versus_time_10_holes->GetZaxis()->SetTitleOffset(1.7);
    debug_inducedcharge_z_versus_time_10_holes->GetYaxis()->SetTitleOffset(1.4);
    debug_inducedcharge_z_versus_time_10_holes->GetXaxis()->SetTitleOffset(1.3);
    debug_inducedcharge_z_versus_time_10_holes->SetTitle("+1 pixel in #eta");
    debug_inducedcharge_z_versus_time_10_holes->GetZaxis()->SetTitle("Average Induced Charge (holes) / e");
    debug_inducedcharge_z_versus_time_10_holes->Draw("colz");
    c1->Print("debug_inducedcharge_z_versus_time_10_holes.pdf");

    gPad->SetRightMargin(0.05);
    debug_inducedcharge_versus_time_nocorr->GetYaxis()->SetRangeUser(0,1.1);
    debug_inducedcharge_versus_time_nocorr->GetXaxis()->SetNdivisions(505);
    debug_inducedcharge_versus_time_nocorr->GetYaxis()->SetNdivisions(505);
    debug_inducedcharge_versus_time_nocorr->GetXaxis()->SetTitle("Starting z Position [mm]");
    debug_inducedcharge_versus_time_nocorr->GetYaxis()->SetTitle("Average Induced Charge (electrons) / e");
    debug_inducedcharge_versus_time_nocorr->GetYaxis()->SetTitleOffset(1.4);
    debug_inducedcharge_versus_time_nocorr->GetXaxis()->SetTitleOffset(1.3);
    debug_inducedcharge_versus_time_nocorr->SetTitle("Induced charge versus Starting Position");
    debug_inducedcharge_versus_time_nocorr->Draw();
    debug_inducedcharge_versus_time_corr->SetLineColor(2);
    debug_inducedcharge_versus_time_corr->Draw("same");

    TLegend * leg = new TLegend(0.2,0.7,0.6,0.9);
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->AddEntry(debug_inducedcharge_versus_time_nocorr,"uncorrected","l");
    leg->AddEntry(debug_inducedcharge_versus_time_corr,"corrected","l");
    leg->Draw();

    c1->Print("debug_inducedcharge_versus_time.pdf");

    debug_inducedcharge_versus_time_nocorr_holes->GetYaxis()->SetRangeUser(0,1.1);
    debug_inducedcharge_versus_time_nocorr_holes->GetXaxis()->SetNdivisions(505);
    debug_inducedcharge_versus_time_nocorr_holes->GetYaxis()->SetNdivisions(505);
    debug_inducedcharge_versus_time_nocorr_holes->GetXaxis()->SetTitle("Starting z Position [mm]");
    debug_inducedcharge_versus_time_nocorr_holes->GetYaxis()->SetTitle("Average Induced Charge (holes) / e");
    debug_inducedcharge_versus_time_nocorr_holes->GetYaxis()->SetTitleOffset(1.4);
    debug_inducedcharge_versus_time_nocorr_holes->GetXaxis()->SetTitleOffset(1.3);
    debug_inducedcharge_versus_time_nocorr_holes->SetTitle("Induced charge versus Starting Position");
    debug_inducedcharge_versus_time_nocorr_holes->Draw();
    debug_inducedcharge_versus_time_corr_holes->SetLineColor(2);
    debug_inducedcharge_versus_time_corr_holes->Draw("same");
    leg->Draw();
    c1->Print("debug_inducedcharge_versus_time_holes.pdf");

    debug_chunksize->GetXaxis()->SetNdivisions(505);
    debug_chunksize->GetYaxis()->SetNdivisions(505);
    debug_chunksize->GetXaxis()->SetTitle("Chunk size");
    debug_chunksize->GetYaxis()->SetTitle("arbitrary units");
    debug_chunksize->GetYaxis()->SetTitleOffset(1.4);
    debug_chunksize->GetXaxis()->SetTitleOffset(1.3);
    debug_chunksize->Draw();
    c1->Print("debug_chunk_size.pdf");
  }
}  

G4double AllPixFEI4RadDamageDigitizer::Phi(G4double x, G4double z, G4double Lx, G4double L){
  
  if (z==0) z=0.00001;//a pathology in the definition.
  G4double pi = 4.*TMath::ATan(1.);
  G4double val = (TMath::Sin(pi*z/L)*TMath::SinH(0.5*pi*Lx/L)/(TMath::CosH(pi*x/L)-TMath::Cos(pi*z/L)*TMath::CosH(0.5*pi*Lx/L)));
  if (val > 0) return TMath::ATan(val)/pi;
  else return TMath::ATan(val)/pi+1;
}

void AllPixFEI4RadDamageDigitizer::SetDetectorDigitInputs(G4double thl){

	// set digitization input values
	// thl
	m_digitIn.thl = thl; // <-- input !
}

Double_t fun_e(Double_t *x, Double_t *par) {
	Double_t result = 0.0;
	for (int i = 0; i<2; i++){
		result += par[5*i]*pow(x[0]+par[5*i+1],i) + par[5*i+2]*pow(x[0]+par[5*i+3],i)*pow(x[1]+par[5*i+4],i);
	}

	return result;
	}

Double_t fun_h(Double_t *x, Double_t *par) {
	Double_t result = 0.0;
	for (int i = 0; i<2; i++){
		result += par[4*i]*pow(x[1]+par[4*i+1],i)*pow(x[0],i) + par[4*i+2]*pow(x[1]+par[4*i+3],i);
	}

	return result;
	}

//2dfit
void AllPixFEI4RadDamageDigitizer::FillHoles(TH2F *distancemap, G4bool isHole, double sensorZ, char *name) {
	std::string filename(name);
	std::string title;
	if (!isHole) {
		title +="Electron Distance Map "+std::to_string((int)biasVoltage)+"V";
  }
	else {
		title +="Hole Distance Map "+std::to_string((int)biasVoltage)+"V";
  }
	filename += std::to_string((int)biasVoltage);
	std::string rootfile(filename+".root");
	TFile *distanceOutFile = new TFile(rootfile.c_str(),"RECREATE");
	TCanvas *c1 = new TCanvas(name,"I am a canvas",600,600);
	
	const Int_t npar = 12;
	
	Double_t f2params_e[npar] = {1,0,1,0,0, 1,0,0.1,1,1, 0,0};
	Double_t f2params_h[npar] = {1,0,1,0,0, 1,0,0.1,1,1, 0,0};
	Double_t f2params[npar] = {};

	
	for (int i = 0; i< npar; i++) {
		f2params[i] = f2params_e[i];
	}
	
	if (isHole) {
		for (int i = 0; i< npar; i++) {
			f2params[i] = f2params_h[i];
		}
	}
	
	
	TF2 *f2 = new TF2("f2",fun_e, 0.,10.,0.,detectorThickness,npar);
	//if (isHole) f2 = new TF2("f2",fun_h, 0.,10.,0.,detectorThickness,npar);
	f2->SetParameters(f2params);	
	
	distancemap->Draw("lego2z");
	
	if (!isHole) {
		distancemap->Fit("f2","WEM");
	}
	else {
		distancemap->Fit("f2","EM");
	}
	
	distancemap->SetTitle(title.c_str());
	distancemap->GetXaxis()->SetNdivisions(505);
	distancemap->GetXaxis()->SetTitle("Initial Position in z [#mum]");
  distancemap->GetYaxis()->SetTitle("Electrons Time Traveled [ns]");
	distancemap->GetZaxis()->SetTitle("Location in z [mm]");
  distancemap->GetXaxis()->SetTitleOffset(1.7);
  distancemap->GetYaxis()->SetTitleOffset(1.7);
	distancemap->GetZaxis()->SetTitleOffset(1.7);
  c1->SetLeftMargin(0.12);
  c1->SetRightMargin(0.17);
  c1->SetBottomMargin(0.12);
  //f2->Draw("cont1 same");
	
	if (debug_maps) {
		std::string pdffile(filename+".pdf");
		G4cout << pdffile << G4endl;
		c1->Print(pdffile.c_str());
		distanceOutFile->Write();
		distancemap->Write(rootfile.c_str());
	}
	else {
		distanceOutFile->Delete();
	}

	for (int i=1; i<= distancemap->GetNbinsX(); i++){
		for (int j=1; j<= distancemap->GetNbinsY(); j++) {
			Double_t bincontent = distancemap->GetBinContent(i,j);
			if (bincontent == 0) {
				Double_t coords[2] = {distancemap->GetXaxis()->GetBinCenter(i),distancemap->GetYaxis()->GetBinCenter(j)};
				Double_t params[npar] = {};
				for (int k = 0; k<npar; k++) {
					params[k] = f2->GetParameter(k);
				}
				Double_t gapvalue = fun_e(coords,params);
				//if (isHole) gapvalue = fun_h(coords,params);
				if (gapvalue>sensorZ) gapvalue = sensorZ;
				else if (gapvalue<0) gapvalue = 0;
				distancemap->SetBinContent(i,j,gapvalue);
			}
		}
	}
	return;
}



///////////////////////////////////
// Radiation Damage Calculations //
///////////////////////////////////

// Look up E field from TCAD mapping
G4double AllPixFEI4RadDamageDigitizer::GetElectricField(G4double z){
	// The z position is in mm, but plot uses um to return electric field in V/cm. Convert to return electric field in MV/mm
	// Is currently using z as distance from position to readout side
  G4int nbinEmap=m_eFieldMap1D->GetNbinsX();
	int n_binz = m_eFieldMap1D->GetXaxis()->FindBin(z*1000);
	G4double zval_bin = m_eFieldMap1D->GetXaxis()->GetBinCenter(n_binz)/1000.;
	G4double electricField=0.;
	if (n_binz < 1) electricField = 0.;
	else if (n_binz > nbinEmap) electricField = 0.;
	else if (n_binz == 1) electricField = m_eFieldMap1D->GetBinContent(1);
  else if (n_binz == nbinEmap) electricField = m_eFieldMap1D->GetBinContent(nbinEmap);
	else if (z==zval_bin) electricField = m_eFieldMap1D->GetBinContent(n_binz);
	else{//attempt to interpolate.
	  double z_other = 0.;
	  if (z > zval_bin) z_other=m_eFieldMap1D->GetXaxis()->GetBinCenter(n_binz+1)/1000.;
	  else z_other=m_eFieldMap1D->GetXaxis()->GetBinCenter(n_binz-1)/1000.;
	  double Delta = fabs(z-zval_bin);
	  double Epsilon = fabs(z-z_other);
	  double EDelta = m_eFieldMap1D->GetBinContent(n_binz);
	  double EEpsilon = m_eFieldMap1D->GetBinContent(m_eFieldMap1D->GetXaxis()->FindBin(z_other*1000));
	  electricField = (pow(Delta,-1)*EDelta+pow(Epsilon,-1)*EEpsilon)/(pow(Delta,-1)+pow(Epsilon,-1));
	}
	return electricField*1.0E-7; //V/cm -> MV/mm.
}

G4double AllPixFEI4RadDamageDigitizer::GetElectricField(G4double x, G4double y, G4double z){
	// The position is in mm, but plot uses um to return electric field in V/cm. Convert to return electric field in MV/mm
	// Is currently using z as distance from position to readout side
	// return the value of the 3D E field
	// TCAD maps are just 1/4 of the pixel cell
	// uses fabs() in order to overcome this
        // defaults to the 1D field of the 3D one is not available.
	G4double electricField=0;
	if(eFieldMap==0){
	  electricField=(GetElectricField(z))*1E07;  // if there is no Efield map loaded, get the "default" 1D Efield
	}
	else{
	 
	 while(fabs(x)>pitchX/2) x=fabs(x)-pitchX/2;  // if it is outside the limits of the cell, get the E field of the other pixel
	 while(fabs(y)>pitchY/2) y=fabs(y)-pitchY/2;
	
   int n_bin = eFieldMap->FindBin(fabs(x)*1000,fabs(y)*1000,z*1000); //x y and z should be in micrometers
   int n_bin_z = eFieldMap->GetZaxis()->FindBin(z*1000);
	 electricField = eFieldMap->GetBinContent(n_bin);
   if(n_bin_z < 1 || n_bin_z > eFieldMap->GetNbinsZ()) electricField=0;
  }
	return electricField*1.0E-7; //V/cm -> MV/mm. 
}

G4double AllPixFEI4RadDamageDigitizer::GetMobility(G4double electricField, G4bool isHole){ 
  //Expects the electric field in MV/mm and returns the mobility in mm^2/(MV*ns).  Temperature is a global variable and is in K.
  
  // Initialize variables so they have the right scope
  G4double vsat = 0;
  G4double ecrit = 0;
  G4double beta = 0;	
  
  //These parameterizations come from C. Jacoboni et al., Solid‐State Electronics 20 (1977) 77‐89. (see also https://cds.cern.ch/record/684187/files/indet-2001-004.pdf).
  if(!isHole){
    vsat = 15.3*pow(temperature,-0.87);	// mm/ns
    ecrit = 1.01E-7*pow(temperature,1.55); // MV/mm
    beta = 2.57E-2*pow(temperature,0.66);
  }
  if(isHole){
    vsat = 1.62*pow(temperature,-0.52);	// mm/ns
    ecrit = 1.24E-7*pow(temperature,1.68); // MV/mm
    beta = 0.46*pow(temperature,0.17);
  }
  
  G4double mobility = (vsat/ecrit)/pow(1+pow((electricField/ecrit),beta),(1/beta));	
  return mobility; // mm^2/(MV*ns)
}

G4double AllPixFEI4RadDamageDigitizer::GetDriftTime(G4bool isHole){
  //returns the drift time in ns.
  
  G4double u = CLHEP::RandFlat::shoot(0.,1.); 
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
  //Expects the electricField in MV/mm.

  G4double hallEffect = 1.13 + 0.0008*(temperature - 273.0); //Hall Scattering Factor - taken from https://cds.cern.ch/record/684187/files/indet-2001-004.pdf
  if (isHole) hallEffect = 0.72 - 0.0005*(temperature - 273.0);
  G4double mobility = AllPixFEI4RadDamageDigitizer::GetMobility(electricField, isHole); //mobility is in mm^2/(MV*ns)
  G4double tanLorentz = hallEffect*mobility*bField*(1.0E-3);  //unit conversion; b-Field is in T = V*s/m^2 
  return tanLorentz;
}

G4double AllPixFEI4RadDamageDigitizer::GetTanLorentz(G4double z1, G4double z2, G4bool isHole){
  //Expects the starting depth z1 and the ending depth z2 in mm.
  //The GetTanLorentz function above that takes only the E field gives the instantaneous Lorentz angle.  This returns the integrated angle.
 
  G4double tanLorentz = -999;
  int nbin=0;
  if(!isHole){
    if(z2>detectorThickness)z2=detectorThickness;
    if(z2<0)z2=0;
    nbin = lorentz_map_e->FindBin(z1,z2);
    tanLorentz = lorentz_map_e->GetBinContent(nbin);
  }
  else{
    if(z2>=detectorThickness) z2=detectorThickness-0.001;
    if(z2<0)z2=0.001;        
    nbin = lorentz_map_h->FindBin(z1,z2);
    tanLorentz = lorentz_map_h->GetBinContent(nbin) ;
  }
  
  return tanLorentz;
		    
}

// ***** Process for each event *****
void AllPixFEI4RadDamageDigitizer::Digitize(){
  //This is the main function!
  
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
  
  G4int nEntries = hitsCollection->entries();
  for(G4int itr  = 0 ; itr < nEntries ; itr++) {

    G4double xpos = (*hitsCollection)[itr]->GetPosWithRespectToPixel().x(); // mm; N.B. this is the eta direction.
    G4double ypos = (*hitsCollection)[itr]->GetPosWithRespectToPixel().y(); // mm; N.B. this is the phi direction.

    // Given as distance from center of pixel, + toward electrode (where electrons are read out)
    G4double zpos = detectorThickness/2. - (*hitsCollection)[itr]->GetPosWithRespectToPixel().z();	// mm; Gives distance to readout side for electron

    G4double eHitTotal = (*hitsCollection)[itr]->GetEdep(); // Energy deposition for the hit, internal unit
    tempPixel.first  = (*hitsCollection)[itr]->GetPixelNbX();
    tempPixel.second = (*hitsCollection)[itr]->GetPixelNbY();
		
    if(itr>0){ // this is done in order to prevent a bug where if the charge was exactly on the border of the pixel it would change the 
      // pixel position, eg: xpos1=23 um pixelX 121 -> xpos2=25 um pixelX 122, where it should have been xpos2=25 um pixelX 121
      if(((pitchX/2)-fabs(xpos)) <0.000001 && fabs((*hitsCollection)[itr]->GetPixelNbX() - (*hitsCollection)[itr-1]->GetPixelNbX())>1) {
        tempPixel.first  = tempPixel.first -1;
      }
    }
    // In case the charge moves into a neighboring pixel
    pair<G4int, G4int> extraPixel;
    extraPixel = tempPixel;

    map<pair<G4int, G4int>, G4double > total_induced_uncorrected;
    map<pair<G4int, G4int>, G4double > total_induced_corrected;
    map<pair<G4int, G4int>, G4double > total_induced_den;
    map<pair<G4int, G4int>, G4double > total_induced_uncorrected_holes;
    map<pair<G4int, G4int>, G4double > total_induced_corrected_holes;
    map<pair<G4int, G4int>, G4double > total_induced_den_holes;
    // Split the charge into subcharges (# = precision) that are diffused separately to the electrode
    //precision = int(eHitTotal/elec);
    for(G4int nQ  = 0 ; nQ < precision ; nQ++) {
      
      //Each e-h chunk pair has eHitTotal / precision energy.  When the Ramo potential is used, the total induced energy is
      //e*[phi(final h) - phi(initial)] - e*[phi(final e) - phi(initial)].  Therefore, the entire chunk needs to have the
      //proper energy (i.e. should not divide by 2).  
      G4double eHit = G4double(eHitTotal)/(precision); //eV
      G4double eHitbefore = eHit;			
      
      //Need to determine how many elementary charges this charge chunk represents.
      double chunk_size = eHit/elec; 
      double kappa = 1./sqrt(chunk_size);
      if (dodebug) debug_chunksize->Fill(chunk_size); //precision = int(eHitTotal/elec);  would have chunk_size = 1.

      // Loop over everything following twice, once for holes and once for electrons
      map<pair<G4int, G4int>, G4double > total_induced;
      map<pair<G4int, G4int>, G4double > total_induced_holes;
      double driftTime_avg = 0.;
      double travelTime_avg = 0.;
      double driftTime_avg_holes = 0.;
      double travelTime_avg_holes = 0.;
      for(G4int eholes=0 ; eholes<2 ; eholes++) {
		    
        isHole = false; // Set a condition to keep track of electron/hole-specific functions
        if (eholes == 1) isHole = true;
        
        G4double electricField = GetElectricField(xpos,ypos,zpos); //in MV/mm; N.B. defaults to 1D if 3D is not provided.  

        // Reset extraPixel coordinates each time through loop
        extraPixel = tempPixel;

        G4double timeToElectrode = GetTimeToElectrode(zpos, isHole); //ns
        G4double driftTime = GetDriftTime(isHole); //ns
        G4double drift_time_constant = trappingTimeElectrons; //ns
        //for debugging plots later.
        if (eholes==0) travelTime_avg=(driftTime < timeToElectrode ? driftTime : timeToElectrode);
        if (eholes==0) driftTime_avg=timeToElectrode;
        if (eholes==1) travelTime_avg_holes=(driftTime < timeToElectrode ? driftTime : timeToElectrode);
        if (eholes==1) driftTime_avg_holes=timeToElectrode;

        if (isHole) drift_time_constant = trappingTimeHoles;
        double average_charge = charge_chunk_map_e->GetBinContent(charge_chunk_map_e->FindBin(fabs(xpos),fabs(ypos),zpos)); //expects input in mm
        if (isHole) average_charge = charge_chunk_map_h->GetBinContent(charge_chunk_map_h->FindBin(fabs(xpos),fabs(ypos),zpos));
        
        G4double zposD = zpos; //this is the distance between the initial position and the final position.
        if (isHole) zposD = detectorThickness - zpos;
        if (driftTime < timeToElectrode){
          int nbin=0;
          if(!isHole){
            nbin = distancemap_e->FindBin(zpos,min(driftTime,timeToElectrode));
            zposD = zpos-distancemap_e->GetBinContent(nbin);
          }
          else {
            nbin = distancemap_h->FindBin(zpos,min(driftTime,timeToElectrode));
            zposD = distancemap_h->GetBinContent(nbin) - zpos;
          }
        }

        G4double tanLorentz =-1.;
        if(!isHole) tanLorentz = GetTanLorentz(zpos,zpos-zposD, isHole);
        else tanLorentz = GetTanLorentz(zpos,zpos+zposD, isHole);
        
        /*
          Diffusion via the Einstein relation
          D = mu * kB * T / q
          D = (mu / mm^2/MV*ns) * (T/273 K) * 0.024 microns^2 / ns
        */

        G4double Dt = GetMobility(electricField, isHole)*(0.024)*timeToElectrode*temperature/273.;
        G4double rdif=sqrt(Dt)/1000; //in mm
        if (defaultDiffusion > 0) rdif=defaultDiffusion*sqrt(zposD/0.3); 
        if (!doDiff) rdif = 0.;
        G4double termRand=CLHEP::RandGauss::shoot(0,1);
        G4double xposD=xpos+zposD*tanLorentz+termRand*rdif; //Both e and h move in the same direciton under B-field (q-reversed but also the velocity direction)
        G4double yposD=ypos+rdif*CLHEP::RandGauss::shoot(0,1);
                
        int loc_x = tempPixel.first;
        int loc_y = tempPixel.second;
        //Record the induced charge from the difference in the ramo potential between the trapped location and all electordes.			
        // ramo potential at electrode based on (x,y,z) position in micrometers
        // -- loop in the x-coordinate
        for (int i=-2; i<=2; i++){
          G4double x_neighbor = i*pitchX;
          extraPixel.first = loc_x + i;
          // -- loop in the y-coordinate
          for (int j=-2; j<=2; j++){
            G4double y_neighbor = j*pitchY;
            extraPixel.second = loc_y + j;
            
            if (i!=0 && j!=0) continue;

            int nbin = ramoPotentialMap->FindBin(fabs((xpos-x_neighbor)*1000),fabs((ypos-y_neighbor)*1000),zpos*1000); //take the absolute value because we only have 1/4 of the ramo
            double ramo_i=0;
            if (!ramoPotentialMap->IsBinUnderflow(nbin) && !ramoPotentialMap->IsBinOverflow(nbin)){ //check if the position is inside the map, else ramo=0
              ramo_i = ramoPotentialMap->GetBinContent(nbin);
            }
            int nbin2 = ramoPotentialMap->FindBin(fabs((xposD-x_neighbor)*1000),fabs((yposD-y_neighbor)*1000),(zpos+zposD)*1000); //for this check, zpos doesn't matter.
            if(!isHole) nbin2 = ramoPotentialMap->FindBin(fabs((xposD-x_neighbor)*1000),fabs((yposD-y_neighbor)*1000),(zpos-zposD)*1000);
            double ramo=0;
            if (!ramoPotentialMap->IsBinUnderflow(nbin2) && !ramoPotentialMap->IsBinOverflow(nbin2)){ //check if the position is inside the map, else ramo=0
              ramo = ramoPotentialMap->GetBinContent(nbin2);		     
            }
            if( (zpos-zposD==0.)  &&  (fabs(xposD-x_neighbor)>pitchX/2 || fabs(yposD-y_neighbor)>pitchY/2) ) ramo=0;
            if( (zpos-zposD==0.)  &&  fabs(xposD-x_neighbor)<=pitchX/2 && fabs(yposD-y_neighbor)<=pitchY/2 ) ramo=1;
            
            // Record deposit
            double eHitRamo =(1-2*isHole)*eHit*(ramo - ramo_i);  //eV
            
            if (doChunkCorrection){
              //X' -> \mu + kappa * (X-\mu)  
              eHitRamo = (1-2*isHole)*eHit*(average_charge + kappa*((1-2*isHole)*(ramo - ramo_i)-average_charge));
            }
            pixelsContent[extraPixel] += eHitRamo; 
            if (eholes==0){
              total_induced[extraPixel] += eHitRamo/eHit;
              total_induced_den[extraPixel] += 1.;
              total_induced_uncorrected[extraPixel] += eHitRamo;
              total_induced_corrected[extraPixel] += eHit*(average_charge + kappa*((ramo - ramo_i)-average_charge));
            }
            else{
              total_induced_holes[extraPixel] += eHitRamo/eHit;
              total_induced_den_holes[extraPixel] += 1.;
              total_induced_uncorrected_holes[extraPixel] += eHitRamo;
              total_induced_corrected_holes[extraPixel] += eHit*(average_charge + kappa*(-(ramo - ramo_i)-average_charge));
            }
          } //loop over y
        } //loop over x
      }  // end loop over charges/holes
      if (dodebug){
        //N.B. this test only works if the B field is turned off.
        map<pair<G4int, G4int>, G4double >::iterator pCItr = total_induced.begin();
        for( ; pCItr != total_induced.end() ; pCItr++){
          if ((*hitsCollection)[itr]->GetPixelNbX()==(*pCItr).first.first && (*hitsCollection)[itr]->GetPixelNbY()==(*pCItr).first.second){ //central pixel
            debug_inducedcharge_z_versus_time_00_num->Fill(zpos*1000,travelTime_avg,(*pCItr).second);
            debug_inducedcharge_z_versus_time_00_den->Fill(zpos*1000,travelTime_avg);
            debug_inducedcharge_z_versus_time_00_num_holes->Fill(zpos*1000,travelTime_avg_holes,total_induced_holes[(*pCItr).first]);
            debug_inducedcharge_z_versus_time_00_den_holes->Fill(zpos*1000,travelTime_avg_holes);
          }
          if ((*hitsCollection)[itr]->GetPixelNbX()==(*pCItr).first.first && (*hitsCollection)[itr]->GetPixelNbY()+1==(*pCItr).first.second){ //shift one in phi
            debug_inducedcharge_z_versus_time_01_num->Fill(zpos*1000,travelTime_avg,(*pCItr).second);
            debug_inducedcharge_z_versus_time_01_den->Fill(zpos*1000,travelTime_avg);
            debug_inducedcharge_z_versus_time_01_num_holes->Fill(zpos*1000,travelTime_avg_holes,total_induced_holes[(*pCItr).first]);
            debug_inducedcharge_z_versus_time_01_den_holes->Fill(zpos*1000,travelTime_avg_holes);
          }
          if ((*hitsCollection)[itr]->GetPixelNbX()+1==(*pCItr).first.first && (*hitsCollection)[itr]->GetPixelNbY()==(*pCItr).first.second){ //shift one in eta
            debug_inducedcharge_z_versus_time_10_num->Fill(zpos*1000,travelTime_avg,(*pCItr).second);
            debug_inducedcharge_z_versus_time_10_den->Fill(zpos*1000,travelTime_avg);
            debug_inducedcharge_z_versus_time_10_num_holes->Fill(zpos*1000,travelTime_avg_holes,total_induced_holes[(*pCItr).first]);
            debug_inducedcharge_z_versus_time_10_den_holes->Fill(zpos*1000,travelTime_avg_holes);
          }
        }
      }
    } // end loop over nQ charges
    if (dodebug){
      map<pair<G4int, G4int>, G4double >::iterator pCItr = total_induced_corrected.begin();
      for( ; pCItr != total_induced_corrected.end() ; pCItr++){
        if ((*hitsCollection)[itr]->GetPixelNbX()==(*pCItr).first.first && (*hitsCollection)[itr]->GetPixelNbY()==(*pCItr).first.second){ //central pixel  
          debug_inducedcharge_versus_time_corr->Fill(zpos,(*pCItr).second/eHitTotal);
          debug_inducedcharge_versus_time_nocorr->Fill(zpos,total_induced_uncorrected[(*pCItr).first]/eHitTotal);
          debug_inducedcharge_versus_time_corr_holes->Fill(zpos,total_induced_corrected_holes[(*pCItr).first]/eHitTotal);
          debug_inducedcharge_versus_time_nocorr_holes->Fill(zpos,total_induced_uncorrected_holes[(*pCItr).first]/eHitTotal);
        }
      }
    }
  }// end loop over nEntries
  
  // Now that pixelContent is filled, create one digit per pixel
  map<pair<G4int, G4int>, G4double >::iterator pCItr = pixelsContent.begin();
  
  for( ; pCItr != pixelsContent.end() ; pCItr++) {

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
  
} // end Digitize function
