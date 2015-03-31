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

  epsilon = 11.8*8.854187817e-14; // [F/cm] -> F=As/V (silicon)
  echarge=1.60217646e-19; //[C=As]
  Default_Hole_Mobility=480.0; //[cm2/Vs] Hole mobility
  Default_Hole_D=12; //;// Hole diffusion [cm2/s]

  Default_Electron_Mobility=1415.0; //[cm2/Vs] Electron mobility
  Default_Electron_D=36; //;// Electron diffusion [cm2/s]
  temperature_T=300.0;
  elec=3.6*eV;

  // Registration of digits collection name
  collectionName.push_back(digitColName);
  m_hitsColName.push_back(hitsColName);


  //Geometry description
  AllPixGeoDsc * gD = GetDetectorGeoDscPtr();
  thickness=gD->GetSensorZ(); //thickness[mm] 
  pitchX=gD->GetPixelX();
  pitchY=gD->GetPixelY();
  nPixX=gD->GetNPixelsX();
  nPixY=gD->GetNPixelsY();

  // Parameters read from pixeldetector.xml:
  resistivity=gD->GetResistivity();
  V_B=gD->GetBiasVoltage();
  SensorType=gD->GetSensorType();
  Threshold=gD->GetThreshold()*elec;


  //-------- Calibration --------//
  CalibrationFile=gD->GetCalibrationFile();
  if (CalibrationFile.find(".root")!= std::string::npos)
    {
      G4cout << "pixel-by-pixel calibration" << G4endl;
      //Pixel-by-pixel calibration
      size_calibX=nPixX;
      size_calibY=nPixY;
      TFile* file_calib=TFile::Open(CalibrationFile);
      if (file_calib==0 || file_calib->IsZombie())
	{
	  G4cout << "Error opening calibration file!!! " << G4endl;
	  //return -1; 
	}
      else
	{
	  initialise_calibrationMatrices(size_calibX, size_calibY);
	  TTree* tree=(TTree*)file_calib->Get("fitPara"); 
	  int pixx;
	  int pixy;
	  Float_t a;
	  Float_t b;
	  Float_t c;
	  Float_t d;

	  tree->SetBranchAddress("pixx", &pixx);
	  tree->SetBranchAddress("pixy", &pixy);
	  tree->SetBranchAddress("a", &a);
	  tree->SetBranchAddress("b", &b);
	  tree->SetBranchAddress("c", &c);
	  tree->SetBranchAddress("d", &d);
	  int nEntries = tree->GetEntries();
	  for (int i=0; i<nEntries; ++i)
	    {
	      tree->GetEntry(i);
	      SurrogateA[pixx][pixy]=a;
	      SurrogateB[pixx][pixy]=b;
	      SurrogateC[pixx][pixy]=c;
	      SurrogateD[pixx][pixy]=d;
	      G4cout << "a=" << a << G4endl;
	      ThresholdMatrix[pixx][pixy]=(d*a-b+sqrt((b+d*a)*(b+d*a)+4*a*c))/(2*a);
	    }
	}
    }
  else if (CalibrationFile.find(".txt")!= std::string::npos)
    {
      G4cout << "Global calibration" << G4endl;
      //Global calibration
      size_calibX=1;
      size_calibY=1;

      initialise_calibrationMatrices(size_calibX, size_calibY);
      ifstream in;
      TString x;
      in.open(CalibrationFile);

      while (1) 
	{
	  in >> x;
	  G4cout << x << G4endl;
	  TObjArray* str = x.Tokenize(","); 
	  if((TObjString*)str->At(0) && (TObjString*)str->At(1) && (TObjString*)str->At(2)  && (TObjString*)str->At(3))
	    {
	      double a=((TObjString*)str->At(0))->String().Atof();
	      double b=((TObjString*)str->At(1))->String().Atof();
	      double c=((TObjString*)str->At(2))->String().Atof();
	      double d=((TObjString*)str->At(3))->String().Atof();
	      SurrogateA[0][0]=a;
	      SurrogateB[0][0]=b;
	      SurrogateC[0][0]=c;
	      SurrogateD[0][0]=d;
	      ThresholdMatrix[0][0]=(d*a-b+sqrt((b+d*a)*(b+d*a)+4*a*c))/(2*a);
	      G4cout << "ThresholdMatrix[0][0]=" << ThresholdMatrix[0][0] << endl;
	    }      
	  if (!in.good()) break;
	}
    }
  //-------- End Calibration --------//
  G4cout << "nilou: resistivity= " << resistivity << G4endl;

  ///=========== Mobility calculation =============///
  if (SensorType=="p-in-n")
    {
      mobility_const=Default_Hole_Mobility;
      diffusion_const=Default_Hole_D;

      
      //non const mobility
      vm=1.62*TMath::Power(10, 8)*TMath::Power(temperature_T, -0.52);
      Ec=1.24*TMath::Power(temperature_T, 1.68);
      beta=0.46*TMath::Power(temperature_T, 0.17);
    }
  else if (SensorType=="n-in-p")
    {
      mobility_const=Default_Electron_Mobility;
      diffusion_const=Default_Electron_D;

      //non const mobility
      vm=1.53*TMath::Power(10, 9)*TMath::Power(temperature_T, -0.87);
      Ec=1.01*TMath::Power(temperature_T, 1.55);
      beta=2.57*TMath::Power(10, -2)*TMath::Power(temperature_T, 1.55);
    }
  G4cout << "mobility_const=" << mobility_const << G4endl;
  G4cout << "diffusion_const=" << diffusion_const << G4endl;
  ///=================================================///

  //**************Charge sharing********************//
  V_D=((thickness/cm)*(thickness/cm))/(2*epsilon*mobility_const*resistivity);
  depletionWidth=TMath::Sqrt(2*V_B*(thickness/cm)*(thickness/cm)/(2*V_D));      
  
  G4cout << "V_D=" << V_D << G4endl;
  G4cout << "depletionWidth=" << depletionWidth*10 << "[mm]" << G4endl;
  G4cout << "V bias=" << V_B << G4endl;
 
}

AllPixTMPXDigitizer::~AllPixTMPXDigitizer()
{
}

void AllPixTMPXDigitizer::Digitize()
{
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

  // =========== To Correct later: If there is only one particle per frame ======
  G4double AvgPosX=0.0;
  G4double AvgPosY=0.0;
  //=====================================================================


  for(G4int itr  = 0 ; itr < nEntries ; itr++)
    {
      G4cout << "=================itr: " << itr << G4endl;

      tempPixel.first  = (*hitsCollection)[itr]->GetPixelNbX();
      tempPixel.second = (*hitsCollection)[itr]->GetPixelNbY();
      
      G4double xpos=(*hitsCollection)[itr]->GetPosWithRespectToPixel().x();
      G4double ypos=(*hitsCollection)[itr]->GetPosWithRespectToPixel().y();
      G4double zpos=(*hitsCollection)[itr]->GetPosWithRespectToPixel().z()+thickness/2.0; // [mm]; zpos=thickness corresponds to the sensor side and zpos=0 corresponds to the pixel side

      AvgPosX+=tempPixel.first*pitchX+xpos+pitchX/2.0;
      AvgPosY+=tempPixel.second*pitchY+ypos+pitchY/2.0;
      G4cout << "posX=" << tempPixel.first*pitchX+xpos+pitchX/2.0 << ", posY=" << tempPixel.second*pitchY+ypos+pitchY/2.0 << ", zpos=" << zpos << G4endl;
      // //**************MC********************//
      // // if(pixelsContent_MC[tempPixel].MC_energy==0.0 && pixelsContent_MC[tempPixel].posX_WithRespectToPixel == -11.0 && pixelsContent_MC[tempPixel].posY_WithRespectToPixel == -11.0)
      // if(itr==0)
      // 	{
      // 	  pixelsContent_MC[tempPixel].posX_WithRespectToPixel=tempPixel.first*pitchX+xpos+pitchX/2.0;
      // 	  pixelsContent_MC[tempPixel].posY_WithRespectToPixel=tempPixel.second*pitchY+ypos+pitchY/2.0;
      // 	  pixelsContent_MC[tempPixel].posZ_WithRespectToPixel=zpos;

      // 	  G4cout << "TMPX: nalipour set MC" << endl;
      // 	  G4cout << "pitchX=" << pitchX << ", pitchY=" << pitchY << endl; 
      // 	  G4cout << "xpos=" << xpos << ", ypos=" << ypos << endl; 
      // 	  G4cout << "posX=" << tempPixel.first*pitchX+xpos+pitchX/2.0 << ", posY=" << tempPixel.second*pitchY+ypos+pitchY/2.0 << G4endl;
      // 	}
      // pixelsContent_MC[tempPixel].MC_energy+=(*hitsCollection)[itr]->GetEdep(); //MC
      


      pair<G4int, G4int> extraPixel;
      extraPixel = tempPixel;
      G4double hit_energy=(*hitsCollection)[itr]->GetEdep();
      // //------- noise 1 ------- // //NOT CORRECT
      // G4cout << "hit_energy before=" << hit_energy/keV << "[keV]" << G4endl;
      // hit_energy = elec*CLHEP::RandGauss::shoot(hit_energy/elec,3*TMath::Sqrt(hit_energy/elec));
      // G4cout << "hit_energy after=" << hit_energy/keV << "[keV]" << G4endl;
      // // hit_energy=CLHEP::RandGauss::shoot(hit_energy, 2.35*TMath::Sqrt(0.1*hit_energy*elec*10e-6)); //intrinsic resolution of semiconductor detectors
      // //hit_energy+=elec*CLHEP::RandGauss::shoot(hit_energy/elec, hit_energy/elec);
      // //Double_t noise=elec*CLHEP::RandGauss::shoot(0, 3*TMath::Sqrt(hit_energy/elec));
      // // G4cout << "hit_energy=" << hit_energy << ", noise=" << noise << G4endl;
      // //hit_energy+=noise;
      // //------- End noise 1 ------- //


      // G4cout << "hit_energy=" << hit_energy << G4endl;
      // hit_energy = hit_energy+CLHEP::RandGauss::shoot(0, TMath::Sqrt(0.1*hit_energy*elec));//noise
      // hit_energy = hit_energy+CLHEP::RandGauss::shoot(0, 0.00572);//noise
      // G4cout << "hit_energy 2=" << hit_energy << G4endl;
      
      if(zpos<depletionWidth*10) // Only charge sharing for the depletion region //depletionWidth*10 [mm]
      	{
	  Double_t electric_field=-((V_B-V_D)/(thickness/cm)+(1-(zpos/cm)/(thickness/cm))*2*V_D/(thickness/cm));
	  
	  Double_t non_const_mobility=(vm/Ec)/(TMath::Power((1+TMath::Power(TMath::Abs(electric_field)/Ec, beta)), 1.0/beta));
	  //G4cout << "non_const_mobility=" << non_const_mobility << G4endl;
	  Double_t drift_time=((thickness/cm)*(thickness/cm))/(2*non_const_mobility*V_D)*TMath::Log((V_B+V_D)/(V_B+V_D-2*V_D*(zpos/cm)/(thickness/cm))); // non-constant drift velocity + non-constant drift velocity
	  Double_t diffusion_RMS=TMath::Sqrt(2.0*diffusion_const*drift_time); //[cm]
	  diffusion_RMS=diffusion_RMS*10;//[mm] 
	  G4cout << "diffusion_RMS=" << diffusion_RMS << " [mm]" << G4endl;
	  	  
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

  // ------------ Add gaussian noise ------------ //
  map<pair<G4int, G4int>, G4double >::iterator pCItr1 = pixelsContent.begin();
  for( ; pCItr1 != pixelsContent.end() ; pCItr1++)
    {
      pair<G4int, G4int> tempPixel;
      tempPixel.first=(*pCItr1).first.first;
      tempPixel.second=(*pCItr1).first.second;
      G4cout << "*********before noise: E=" << pixelsContent[tempPixel]/keV << "[keV]" << G4endl;
      Double_t noise=CLHEP::RandGauss::shoot(0, 200);//electrons
      G4cout << "noise=" << noise << G4endl;
      G4cout << "noise 2=" << noise*elec/MeV << G4endl;
      pixelsContent[tempPixel]+=noise*elec/MeV;
      G4cout << "after noise=" << pixelsContent[tempPixel]/keV << "[keV]" << G4endl;

      // TOT resolution:
      // pixelsContent[tempPixel]+=CLHEP::RandGauss::shoot(0, 5./100.*pixelsContent[tempPixel]); //5% noise due to TOT resolution
    }
  // ------------ End Add gaussian noise ------------ //
  // // ------------ Capacitive coupling ------------ //
  // pair<G4int, G4int> extraPixel;
  // extraPixel = tempPixel;
  // G4double coupling_percentage=0.05;//TMath::Abs(CLHEP::RandGauss::shoot(0.05, 0.05));
  // G4cout << "coupling_percentage=" << coupling_percentage << G4endl;
  // // if (coupling_percentage<0)
  // //   {
  // //     G4cout << "coupling percentage negative" << G4endl;
  // //   }
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
  G4cout << "Threshold=" <<  Threshold/keV << " [keV]" << G4endl;
  G4cout << "AvgPosX=" << AvgPosX/nEntries << ", AvgPosY=" << AvgPosY/nEntries << G4endl;
  // With charge sharing
  map<pair<G4int, G4int>, G4double >::iterator pCItr = pixelsContent.begin();
  for( ; pCItr != pixelsContent.end() ; pCItr++)
    {
      pair<G4int, G4int> tempPixel;
      tempPixel.first=(*pCItr).first.first;
      tempPixel.second=(*pCItr).first.second;

      
      G4double a=0.0;
      G4double b=0.0;
      G4double c=0.0;
      G4double t=0.0;
      
      if (size_calibX==1 && size_calibY==1)
	{
	  // m_digitIn.thl=ThresholdMatrix[0][0]; //keV
	  a=SurrogateA[0][0];
	  b=SurrogateB[0][0];
	  c=SurrogateC[0][0];
	  t=SurrogateD[0][0];
	}
      else
	{
	  // m_digitIn.thl=ThresholdMatrix[tempPixel.first][tempPixel.second]; //keV
	  a=SurrogateA[tempPixel.first][tempPixel.second];
	  b=SurrogateB[tempPixel.first][tempPixel.second];
	  c=SurrogateC[tempPixel.first][tempPixel.second];
	  t=SurrogateD[tempPixel.first][tempPixel.second];
	}
      
      if((*pCItr).second/keV > Threshold/keV) // over threshold !
	{
	  G4int TOT=a*((*pCItr).second)/keV+b-c/(((*pCItr).second/keV)-t);
	  
	  AllPixTMPXDigit * digit = new AllPixTMPXDigit;

	  digit->SetPixelIDX((*pCItr).first.first);
	  digit->SetPixelIDY((*pCItr).first.second);
	  digit->SetPixelEnergyDep(((*pCItr).second)/keV); //Energy with charge sharing
	  digit->SetPixelCounts(TOT); //TOT value

	  // ====== TO be corrected later==================== //
	  digit->Set_posX_WithRespectoToPixel(AvgPosX/nEntries);
	  digit->Set_posY_WithRespectoToPixel(AvgPosY/nEntries);
	  //===================================================//
	  // if (pixelsContent_MC.count(tempPixel))
	  //   {
	  //     digit->SetPixelEnergyMC((pixelsContent_MC[tempPixel].MC_energy)/keV); //MC value
	  //     digit->Set_posX_WithRespectoToPixel(pixelsContent_MC[tempPixel].posX_WithRespectToPixel);
	  //     digit->Set_posY_WithRespectoToPixel(pixelsContent_MC[tempPixel].posY_WithRespectToPixel);
	  //     digit->Set_posZ_WithRespectoToPixel(pixelsContent_MC[tempPixel].posZ_WithRespectToPixel);
	  //     pixelsContent_MC.erase(tempPixel);
	  //   }
	  m_digitsCollection->insert(digit);
	}
    }
  // //MC only or under-threshold pixels
  // map<pair<G4int, G4int>, MC_content >::iterator pCItr_MC = pixelsContent_MC.begin();
  // for( ; pCItr_MC != pixelsContent_MC.end() ; pCItr_MC++)
  //   {
  //     G4cout << "x=" << tempPixel.first << ", y=" << tempPixel.second << ", energyMC=" << (pixelsContent_MC[tempPixel].MC_energy)/keV << G4endl;
  //     pair<G4int, G4int> tempPixel;
  //     tempPixel.first=(*pCItr_MC).first.first;
  //     tempPixel.second=(*pCItr_MC).first.second;

  //     AllPixTMPXDigit * digit = new AllPixTMPXDigit;

  //     digit->SetPixelIDX((*pCItr_MC).first.first);
  //     digit->SetPixelIDY((*pCItr_MC).first.second);
  //     digit->SetPixelEnergyMC(pixelsContent_MC[tempPixel].MC_energy/keV); //MC value
  //     digit->Set_posX_WithRespectoToPixel(pixelsContent_MC[tempPixel].posX_WithRespectToPixel);
  //     digit->Set_posY_WithRespectoToPixel(pixelsContent_MC[tempPixel].posY_WithRespectToPixel);
  //     digit->Set_posZ_WithRespectoToPixel(pixelsContent_MC[tempPixel].posZ_WithRespectToPixel);
  //     m_digitsCollection->insert(digit);
  //   }
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

   //G4cout << "Integral=" << Integral/4.0 << G4endl;
   
   G4double energybis=Integral*Energy/4.0; //*(TMath::Pi())*(TMath::Pi());
   // G4double energybis= (Energy*(-TMath::Erf((x1-xhit)/(TMath::Sqrt(2.)*Sigma))+TMath::Erf((x2-xhit)/(TMath::Sqrt(2.)*Sigma)))
   // 			*(-TMath::Erf((y1-yhit)/(TMath::Sqrt(2.)*Sigma))+TMath::Erf((y2-yhit)/(TMath::Sqrt(2.0)*Sigma))))/4.0*(TMath::Pi())*(TMath::Pi());
   return energybis;
}


void AllPixTMPXDigitizer::initialise_calibrationMatrices(int size_calibX, int size_calibY)
{
  ThresholdMatrix = new G4double*[size_calibX];
  for(int i=0;i<size_calibX;i++)
    {
      ThresholdMatrix[i] = new G4double[size_calibY];
    }
  SurrogateA = new G4double*[size_calibX];
  for(int i=0;i<size_calibX;i++)
    {
      SurrogateA[i] = new G4double[size_calibY];
    }
  SurrogateB = new G4double*[size_calibX];
  for(int i=0;i<size_calibX;i++)
    {
      SurrogateB[i] = new G4double[size_calibY];
    }
  SurrogateC = new G4double*[size_calibX];
  for(int i=0;i<size_calibX;i++)
    {
      SurrogateC[i] = new G4double[size_calibY];
    }
  SurrogateD = new G4double*[size_calibX];
  for(int i=0;i<size_calibX;i++)
    {
      SurrogateD[i] = new G4double[size_calibY];
    }
}
