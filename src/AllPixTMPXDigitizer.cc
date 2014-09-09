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

AllPixTMPXDigitizer::AllPixTMPXDigitizer(G4String modName, G4String hitsColName, G4String digitColName) 
  : AllPixDigitizerInterface (modName) {

  // Registration of digits collection name
  collectionName.push_back(digitColName);
  m_hitsColName.push_back(hitsColName);

  // threshold
  m_digitIn.thl = 1026;// For L04-W0125 //800 electrons // TO DO from gear file
  
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
  G4double elec=3.64*eV;
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
  const Double_t epsilon = 11.8*8.854187817e-14; // [F/cm] -> F=As/V (silicon)
  const Double_t echarge=1.60217646e-19; //[C=As]
  const Double_t Default_Hole_Mobility=480.0; //[cm2/Vs] Hole mobility
  const Double_t Default_Hole_D=12; //;// Hole diffusion [cm2/s]

  const Double_t Default_Electron_Mobility=1415.0; //[cm2/Vs] Electron mobility
  const Double_t Default_Electron_D=36; //;// Electron diffusion [cm2/s]


  //// ============PARAMETERS TO ADJUST================
  string sensorType="p-in-n"; // or n-in-p
  //Double_t V_B=35; //[V] //Run 1189
  Double_t V_B=35; //[V] //Run 2302 Vb=-35[V]

  //-------L04-W0125-------// //100um p-in-n
  G4double a=14.2;
  G4double b=437.2;
  G4double c=1830;
  G4double t=-9.26e-7;
  //-----------------------//
  // // //-------B06-W0125-------// //200um n-in-p
  // G4double a=29.8;
  // G4double b=534.1;
  // G4double c=1817;
  // G4double t=0.7;
  // // //-----------------------//

  Double_t resistivity=5000; //[ohm cm] 
  ///=================================================
  Double_t mobility_const=0.0;
  Double_t diffusion_const=0.0;
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




  for(G4int itr  = 0 ; itr < nEntries ; itr++)
    {
      G4cout << "=================itr: " << itr << G4endl;
      //G4cout << "Thickness=" << thickness << G4endl;
      tempPixel.first  = (*hitsCollection)[itr]->GetPixelNbX();
      tempPixel.second = (*hitsCollection)[itr]->GetPixelNbY();
      
      G4double xpos=(*hitsCollection)[itr]->GetPosWithRespectToPixel().x();
      G4double ypos=(*hitsCollection)[itr]->GetPosWithRespectToPixel().y();
      G4double zpos=(*hitsCollection)[itr]->GetPosWithRespectToPixel().z()+thickness/2.0; // [mm]; zpos=thickness corresponds to the sensor side and zpos=0 corresponds to the pixel side

      //**************MC********************//
      if(pixelsContent_MC[tempPixel].MC_energy==0.0 && pixelsContent_MC[tempPixel].posX_WithRespectToPixel == -11.0 && pixelsContent_MC[tempPixel].posY_WithRespectToPixel == -11.0)
	{
	  G4cout << "nalipour set MC" << endl;
	  pixelsContent_MC[tempPixel].posX_WithRespectToPixel=xpos;//(*hitsCollection)[itr]->GetPosWithRespectToPixel().x();
	  pixelsContent_MC[tempPixel].posY_WithRespectToPixel=ypos;//(*hitsCollection)[itr]->GetPosWithRespectToPixel().y();
	  pixelsContent_MC[tempPixel].posZ_WithRespectToPixel=zpos;//(*hitsCollection)[itr]->GetPosWithRespectToPixel().z();
	}
      pixelsContent_MC[tempPixel].MC_energy+=(*hitsCollection)[itr]->GetEdep(); //MC


      //**************Charge sharing********************//

      
 
 
      Double_t Neff=1.0/(resistivity*echarge*mobility_const); //[1/cm3]
      Double_t V_D=(echarge*Neff*(thickness/cm)*(thickness/cm))/(2*epsilon);
 
      Double_t depletionWidth=TMath::Sqrt(2*epsilon*V_B/(echarge*Neff)); //[cm]
   
      pair<G4int, G4int> extraPixel;
      extraPixel = tempPixel;
      G4double hit_energy=(*hitsCollection)[itr]->GetEdep();
      
      if(zpos<depletionWidth*10) // Only charge sharing for the depletion region //depletionWidth*10 [mm]
      	{
	  Double_t electric_field=-((V_B-V_D)/(thickness/cm)+(1-(zpos/cm)/(thickness/cm))*2*V_D/(thickness/cm));
	  
	  //Double_t drift_time=(zpos/cm)/(mobility_const*TMath::Abs(electric_field)); //constant drift time
	  Double_t drift_time=((thickness/cm)*(thickness/cm))/(2*mobility_const*V_D)*TMath::Log((V_B+V_D)/(V_B+V_D-2*V_D*(zpos/cm)/(thickness/cm))); // non-constant drift velocity
	  Double_t diffusion_RMS=TMath::Sqrt(2.0*diffusion_const*drift_time); //[cm]
	  diffusion_RMS=diffusion_RMS*10;//[mm]
	  

	  G4cout << "diffusion_RMS=" << diffusion_RMS << G4endl;
	  G4cout << "pitchX/2.-3.0*diffusion_RMS=" << pitchX/2.-3.0*diffusion_RMS << G4endl;
	  G4cout << "pitchY/2.-3.0*diffusion_RMS=" << pitchY/2.-3.0*diffusion_RMS << G4endl;
	  
	  if(fabs(xpos)>=pitchX/2.-3.0*diffusion_RMS || fabs(ypos)>=pitchY/2.-3.0*diffusion_RMS)
	    {
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
			  //G4cout << "i=" << i << ", j=" << j << ", hit_energy=" << hit_energy << ", Etemp=" << Etemp << G4endl;
			  pixelsContent[extraPixel]+=Etemp;
			}
		    }
		}
	    }
	  else
	    {
	      G4cout << "No charge sharing: hit_energy=" << hit_energy << G4endl;
	      pixelsContent[extraPixel]+=hit_energy; ///////ADD ENERGY
	    }
	}
      // else
      // 	{
      // 	  G4cout << "No charge sharing: hit_energy=" << hit_energy << G4endl;
      // 	  pixelsContent[extraPixel]+=hit_energy; ///////ADD ENERGY
      // 	}
      // 	}
      // else // No charge sharing
      // 	{
      // 	  pixelsContent[extraPixel]+=hit_energy;
      // 	}
      //pixelsContent[tempPixel] = pixelsContent_MC[tempPixel].MC_energy; //+ charge sharing + noise

      // ///===============TEST===============///
      // pixelsContent[extraPixel]+=hit_energy;
      // ///===============END TEST===============///
    }

  //------------------ RECORD DIGITS ------------------//
  // With charge sharing
  map<pair<G4int, G4int>, G4double >::iterator pCItr = pixelsContent.begin();
  for( ; pCItr != pixelsContent.end() ; pCItr++)
    {
      G4cout << "Before thresh: x=" << (*pCItr).first.first << ", y=" << (*pCItr).first.second << ", energy=" << (*pCItr).second << G4endl;
      G4cout << "threshold=" <<  m_digitIn.thl*elec/keV << G4endl;
      if((*pCItr).second > m_digitIn.thl*elec) // over threshold !
	{
	  pair<G4int, G4int> tempPixel;
	  tempPixel.first=(*pCItr).first.first;
	  tempPixel.second=(*pCItr).first.second;

	  G4cout << "x=" << tempPixel.first << ", y=" << tempPixel.second << ", energy=" << (*pCItr).second/keV << G4endl;

	  AllPixTMPXDigit * digit = new AllPixTMPXDigit;

	  digit->SetPixelIDX((*pCItr).first.first);
	  digit->SetPixelIDY((*pCItr).first.second);
	  digit->SetPixelEnergyDep(((*pCItr).second)/keV); //Energy with charge sharing




	  G4int TOT=a*((*pCItr).second)/keV+b-c/(((*pCItr).second/keV)-t);
	  digit->SetPixelCounts(TOT); //TOT value
	  if (pixelsContent_MC.count(tempPixel))
	    {
	      G4cout << "x=" << tempPixel.first << ", y=" << tempPixel.second << ", energyMC=" << (pixelsContent_MC[tempPixel].MC_energy)/keV << G4endl;
	      digit->SetPixelEnergyMC((pixelsContent_MC[tempPixel].MC_energy)/keV); //MC value
	      digit->Set_posX_WithRespectoToPixel(pixelsContent_MC[tempPixel].posX_WithRespectToPixel);
	      digit->Set_posY_WithRespectoToPixel(pixelsContent_MC[tempPixel].posY_WithRespectToPixel);
	      digit->Set_posZ_WithRespectoToPixel(pixelsContent_MC[tempPixel].posZ_WithRespectToPixel);
	      pixelsContent_MC.erase(tempPixel);
	    }
	  //digit->IncreasePixelCounts(); // Counting mode

	  m_digitsCollection->insert(digit);
	}
    }
  //MC only
  map<pair<G4int, G4int>, MC_content >::iterator pCItr_MC = pixelsContent_MC.begin();
  for( ; pCItr_MC != pixelsContent_MC.end() ; pCItr_MC++)
    {
      G4cout << "x=" << tempPixel.first << ", y=" << tempPixel.second << ", energyMC=" << (pixelsContent_MC[tempPixel].MC_energy)/keV << G4endl;
      pair<G4int, G4int> tempPixel;
      tempPixel.first=(*pCItr_MC).first.first;
      tempPixel.second=(*pCItr_MC).first.second;

      AllPixTMPXDigit * digit = new AllPixTMPXDigit;

      digit->SetPixelIDX((*pCItr_MC).first.first);
      digit->SetPixelIDY((*pCItr_MC).first.second);
      digit->SetPixelEnergyMC(pixelsContent_MC[tempPixel].MC_energy/keV); //MC value
      digit->Set_posX_WithRespectoToPixel(pixelsContent_MC[tempPixel].posX_WithRespectToPixel);
      digit->Set_posY_WithRespectoToPixel(pixelsContent_MC[tempPixel].posY_WithRespectToPixel);
      digit->Set_posZ_WithRespectoToPixel(pixelsContent_MC[tempPixel].posZ_WithRespectToPixel);
      //digit->IncreasePixelCounts(); // Counting mode

      m_digitsCollection->insert(digit);
    }
  //----------------------------------------------------//


  G4int dc_entries = m_digitsCollection->entries();
  if(dc_entries > 0){
    G4cout << "--------> Digits Collection : " << collectionName[0]
	   << "(" << m_hitsColName[0] << ")"
	   << " contains " << dc_entries
	   << " digits" << G4endl;
  }

  StoreDigiCollection(m_digitsCollection);
}


//------------------------------------//
// // Uniform electric field. p-in-n
// G4double ComputeElectricField(G4double z, G4double V_B, G4double V_D, G4double thick)
// {
//   return -(V_B/thick+2*z*V_D/(thick*thick));
// }
//------------------------------------//     
 G4double AllPixTMPXDigitizer::IntegrateGaussian(G4double xhit,G4double yhit,G4double Sigma, G4double x1, G4double x2, G4double y1, G4double y2, G4double Energy)
 {
   G4double Integral=(-TMath::Erf((x1-xhit)/(TMath::Sqrt(2.)*Sigma))+TMath::Erf((x2-xhit)/(TMath::Sqrt(2.)*Sigma)))*(-TMath::Erf((y1-yhit)/(TMath::Sqrt(2.)*Sigma))+TMath::Erf((y2-yhit)/(TMath::Sqrt(2.0)*Sigma)));
   //G4cout << "integral=" << Integral << G4endl;
   G4double energybis=Integral*Energy/4.0; //*(TMath::Pi())*(TMath::Pi());
   // G4double energybis= (Energy*(-TMath::Erf((x1-xhit)/(TMath::Sqrt(2.)*Sigma))+TMath::Erf((x2-xhit)/(TMath::Sqrt(2.)*Sigma)))
   // 			*(-TMath::Erf((y1-yhit)/(TMath::Sqrt(2.)*Sigma))+TMath::Erf((y2-yhit)/(TMath::Sqrt(2.0)*Sigma))))/4.0*(TMath::Pi())*(TMath::Pi());
   return energybis;
}
