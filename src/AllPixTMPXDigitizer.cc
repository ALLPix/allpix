/**
 *  Author:
 *    nalipour@cern.ch
 *
 *  allpix Authors:
 *   John Idarraga <idarraga@cern.ch>
 *   Mathieu Benoit <benoit@lal.in2p3.fr>
 */

#include "AllPixTMPXDigitizer.hh"
#include "AllPixTrackerHit.hh"
#include "G4DigiManager.hh"
#include "G4VDigitizerModule.hh"
#include "AllPixGeoDsc.hh"

#include "TMath.h"
#include "TF1.h"

#include "CLHEP/Random/RandGauss.h"

AllPixTMPXDigitizer::AllPixTMPXDigitizer(G4String modName, G4String hitsColName, G4String digitColName) 
  : AllPixDigitizerInterface (modName) {

  // Registration of digits collection name
  collectionName.push_back(digitColName);
  m_hitsColName.push_back(hitsColName);

  // threshold
  //m_digitIn.thl = 1026;// For L04-W0125 //800 electrons // TO DO from gear file
  m_digitIn.thl = 900;  

  //Geometry description
  AllPixGeoDsc * gD = GetDetectorGeoDscPtr();
  thickness=gD->GetSensorZ(); //thickness[mm] 
  pitchX=gD->GetPixelX();
  pitchY=gD->GetPixelY();
  nPixX=gD->GetNPixelsX();
  nPixY=gD->GetNPixelsY();
  

//G4cout << "nalipour detector thickness=" << thickness << G4endl;
}

AllPixTMPXDigitizer::~AllPixTMPXDigitizer()
{
}

void AllPixTMPXDigitizer::Digitize()
{
  this->elec=3.64*eV;
  m_digitsCollection = new AllPixTMPXDigitsCollection("AllPixTMPXDigitizer", collectionName[0] );

  // get the digiManager
  G4DigiManager * digiMan = G4DigiManager::GetDMpointer();

  // BoxSD_0_HitsCollection
  G4int hcID = digiMan->GetHitsCollectionID(m_hitsColName[0]);

  AllPixTrackerHitsCollection * hitsCollection = 0;
  hitsCollection = (AllPixTrackerHitsCollection*)(digiMan->GetHitsCollection(hcID));

  // temporary data structure
  map<pair<G4int, G4int>, MC_content> pixelsContent_MC; //contains information with only MC
  map<pair<G4int, G4int>, G4double > pixelsContent; //contains information with charge sharing
  pair<G4int, G4int> tempPixel;
  G4int nEntries = hitsCollection->entries();

  // AllPixGeoDsc * gD = GetDetectorGeoDscPtr();
  // nPixX=gD->GetNPixelsX();
  // nPixY=gD->GetNPixelsY();




  //Parameters for Charge sharing and TOT!!!!!!
  // BE careful
//  epsilon = 11.8*8.854187817e-14; // [F/cm] -> F=As/V (silicon)
  echarge=1.60217646e-19; //[C=As]
  Default_Hole_Mobility=480.0; //[cm2/Vs] Hole mobility
  Default_Hole_D=12; //;// Hole diffusion [cm2/s]

  Default_Electron_Mobility=1415.0; //[cm2/Vs] Electron mobility
  Default_Electron_D=36; //;// Electron diffusion [cm2/s]


  //// ============PARAMETERS TO ADJUST================
  string sensorType="n-in-p"; // or n-in-p
  //Double_t V_B=35; //[V] //Run 1189
  V_B=35; //[V] //Run 2302 Vb=-35[V]

  //-------L04-W0125-------// //100um p-in-n
  //G4double a=14.2;
  //G4double b=437.2;
  //G4double c=1830;
  //G4double t=-9.26e-7;
  //-----------------------//
  // //-------B06-W0125-------// //200um n-in-p
  //B06: dopant conc.=9.88e11                                                                                                                                                        
  TString carrier_type="n";
  TString bulk_type="p";
  //Double_t nDopants=9.88e11;
  V_D=30.31;
  Double_t temperature=300; //[K]

  G4double a=29.8;
  G4double b=534.1;
  G4double c=1817;
  G4double t=0.7;
  // //-----------------------//

  resistivity=5000; //[ohm cm]
  ///=================================================
  //mobility_const=0.0;
  //diffusion_const=0.0;
  if (sensorType=="p-in-n")
    {
      mobility_const=Default_Hole_Mobility;
      diffusion_const=Default_Hole_D;
    }
  else if (sensorType=="n-in-p")
    {
      mobility_const=Default_Electron_Mobility;
      diffusion_const=Default_Electron_D;
    }
  ///=================================================

  // =========== To Correct later: If there is only one particle per frame ======
  G4double AvgPosX=0.0;
  G4double AvgPosY=0.0;
  //=====================================================================


  for(G4int itr  = 0 ; itr < nEntries ; itr++)
    {
      G4cout << "=================itr: " << itr << G4endl;
      //G4cout << "Thickness=" << thickness << G4endl;
      tempPixel.first  = (*hitsCollection)[itr]->GetPixelNbX();
      tempPixel.second = (*hitsCollection)[itr]->GetPixelNbY();
      
      G4double xpos=(*hitsCollection)[itr]->GetPosWithRespectToPixel().x();
      G4double ypos=(*hitsCollection)[itr]->GetPosWithRespectToPixel().y();
      G4double zpos=(*hitsCollection)[itr]->GetPosWithRespectToPixel().z()+thickness/2.0; // [mm]; zpos=thickness corresponds to the sensor side and zpos=0 corresponds to the pixel side

      AvgPosX+=tempPixel.first*pitchX+xpos+pitchX/2.0;
      AvgPosY+=tempPixel.second*pitchY+ypos+pitchY/2.0;
      
      depletionWidth=(thickness/(2*V_D))*(V_D+V_B)-2e-3; //TMath::Sqrt(V_B/V_D)*thickness; //[cm]
      G4cout << "test depletionWidth=" << (thickness/(2*V_D))*(V_D+V_B) << "[mm]" << G4endl; 
      G4cout << "depletionWidth=" << depletionWidth << "[mm]" << G4endl;

      pair<G4int, G4int> extraPixel;
      extraPixel = tempPixel;
      G4double hit_energy=(*hitsCollection)[itr]->GetEdep();
      
      
      if(zpos<depletionWidth && zpos>=0) // Only charge sharing for the depletion region //depletionWidth*10 [mm]
      	{
	  // zpos=TMath::Abs(zpos); // sometimes zpos=-1.38778e-17 -> Due to the step size
	  G4cout << "zpos=" << zpos << G4endl;
	  // // ======= Linear ======= //
	  // Double_t electric_field=-((V_B-V_D)/(thickness/cm)+(1-(zpos/cm)/(thickness/cm))*2*V_D/(thickness/cm));
	  // Double_t drift_time=(zpos/cm)/(mobility_const*TMath::Abs(electric_field)); //constant drift time
	  // Double_t diffusion_RMS=TMath::Sqrt(2.0*diffusion_const*drift_time); //[cm]
	  // diffusion_RMS=diffusion_RMS*10;//[mm]
	  //  G4cout << "diffusion_RMS=" << diffusion_RMS << "[mm]" << G4endl;

	  // ======= Non-linear ======= //
	  Double_t diffusion_RMS=TMath::Sqrt((TMath::K()*temperature*(thickness/cm)*(thickness/cm)/(echarge*V_D))*TMath::Log((V_B+V_D)/(V_B+V_D-2.0*V_D*(zpos/cm)/(thickness/cm))));
	  diffusion_RMS=TMath::Sqrt(2)*diffusion_RMS*10;//[mm]
	  G4cout << "diffusion_RMS=" << diffusion_RMS << "[mm]" << G4endl;
  
	  
	  for(int i=-1; i<=1; i++)
	    {
	      for(int j=-1; j<=1; j++)
		{
		  extraPixel=tempPixel;
		  extraPixel.first +=i;
		  extraPixel.second+=j;
		  if(extraPixel.first >= 0 && extraPixel.second>=0 && extraPixel.first < nPixX && extraPixel.second < nPixY)
		    {		      
		      G4double Etemp = IntegrateGaussian(xpos/nm, ypos/nm, diffusion_RMS/nm, (-pitchX/2.0 + i*pitchX)/nm, (-pitchX/2.+(i+1)*pitchX)/nm, (-pitchY/2 + j*pitchY)/nm, (-pitchY/2+(j+1)*pitchY)/nm, hit_energy);
		      pixelsContent[extraPixel]+=Etemp;
		    }
		}
	    }
	}
    }

  // // ------------ Capacitive coupling ------------ //
  // pair<G4int, G4int> extraPixel;
  // extraPixel = tempPixel;
  // G4double coupling_percentage=0.05;//TMath::Abs(CLHEP::RandGauss::shoot(0.05, 0.05));
  // G4double couplingEnergy=coupling_percentage*(G4double)pixelsContent[tempPixel];
  // G4cout << "couplingEnergy=" << couplingEnergy << G4endl;

  // pixelsContent[tempPixel]-=couplingEnergy;
  // for(int i=-1; i<=1; i++)
  //   {
  //     for(int j=-1; j<=1; j++)
  // 	{
  // 	  if(i!=0 && j!=0 && i*j==0)
  // 	    {
  // 	      extraPixel=tempPixel;
  // 	      extraPixel.first +=i;
  // 	      extraPixel.second+=j;
  // 	      pixelsContent[extraPixel]+=couplingEnergy/4.0;
  // 	    }
  // 	}
  //   }
  // // ------------ End capacitive coupling ------------ //
 

  //------------------ RECORD DIGITS ------------------//
  // With charge sharing
  map<pair<G4int, G4int>, G4double >::iterator pCItr = pixelsContent.begin();
  for( ; pCItr != pixelsContent.end() ; pCItr++)
    {
      // Double_t threshold=CLHEP::RandGauss::shoot(m_digitIn.thl, 35); // ~35 electrons noise on the threshold
      // Double_t threshold=m_digitIn.thl;
      Double_t threshold=3.836+CLHEP::RandGauss::shoot(0, 0.057);
      // G4cout << "threshold=" << threshold << G4endl;
      // G4cout << "energy=" << ((*pCItr).second)/keV << " [keV]" << G4endl;
      // //--- Electronic noise ---//
      // Double_t electronic_noise=elec*CLHEP::RandGauss::shoot(0, 200)/MeV;  //200 electrons electronic noise
      // ((*pCItr).second)+=electronic_noise;

      // // // TOT noise
      // // ((*pCItr).second)=CLHEP::RandGauss::shoot(((*pCItr).second), 5.0/100.0*((*pCItr).second)); //?????

      G4cout << "threshold*elec/keV=" << threshold*elec/keV << G4endl;

      
      //if(((*pCItr).second)/keV > threshold*elec/keV) // over threshold !
      if(((*pCItr).second)/keV > threshold) // over threshold !
	{
	  G4cout << "yes" << G4endl;
	  tempPixel.first=(*pCItr).first.first;
	  tempPixel.second=(*pCItr).first.second;

	  AllPixTMPXDigit * digit = new AllPixTMPXDigit;

	  digit->SetPixelIDX((*pCItr).first.first);
	  digit->SetPixelIDY((*pCItr).first.second);
	  digit->SetPixelEnergyDep(((*pCItr).second)/keV); //Energy with charge sharing

	  // ====== TO be corrected later==================== //
	  digit->Set_posX_WithRespectoToPixel(AvgPosX/nEntries);
	  digit->Set_posY_WithRespectoToPixel(AvgPosY/nEntries);
	  //===================================================//


	  G4int TOT=a*((*pCItr).second)/keV+b-c/(((*pCItr).second/keV)-t);
	  digit->SetPixelCounts(TOT); //TOT value

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


 G4double AllPixTMPXDigitizer::IntegrateGaussian(G4double xhit,G4double yhit,G4double Sigma, G4double x1, G4double x2, G4double y1, G4double y2, G4double Energy)
 {
   G4double Integral=(-TMath::Erf((x1-xhit)/(TMath::Sqrt(2.)*Sigma))+TMath::Erf((x2-xhit)/(TMath::Sqrt(2.)*Sigma)))*(-TMath::Erf((y1-yhit)/(TMath::Sqrt(2.)*Sigma))+TMath::Erf((y2-yhit)/(TMath::Sqrt(2.0)*Sigma)));

   G4double energybis=Integral*Energy/4.0; //*(TMath::Pi())*(TMath::Pi());
   return energybis;
}
