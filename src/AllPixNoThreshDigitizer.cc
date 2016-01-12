/**
 *  Author:
 *    nalipour
 *
 *  allpix Authors:
 *   John Idarraga <idarraga@cern.ch>
 *   Mathieu Benoit <benoit@lal.in2p3.fr>
 */

#include "AllPixNoThreshDigitizer.hh"
#include "AllPixTrackerHit.hh"
#include "G4DigiManager.hh"
#include "G4VDigitizerModule.hh"
#include "AllPixGeoDsc.hh"

#include "TMath.h"
#include "TF1.h"

#include "CLHEP/Random/RandGauss.h"

AllPixNoThreshDigitizer::AllPixNoThreshDigitizer(G4String modName, G4String hitsColName, G4String digitColName) 
  : AllPixDigitizerInterface (modName) {

  // Registration of digits collection name
  collectionName.push_back(digitColName);
  m_hitsColName.push_back(hitsColName);

  // threshold
  //m_digitIn.thl = 0.;

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

}

AllPixNoThreshDigitizer::~AllPixNoThreshDigitizer(){

}

void AllPixNoThreshDigitizer::Digitize()
{
  G4double elec=3.64*eV;
  // create the digits collection
  m_digitsCollection = new AllPixNoThreshDigitsCollection("AllPixNoThreshDigitizer", collectionName[0] );

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


  //Parameters for Charge sharing and TOT!!!!!!
  // BE careful
  const Double_t epsilon = 11.8*8.854187817e-14; // [F/cm] -> F=As/V (silicon)
  const Double_t echarge=1.60217646e-19; //[C=As]
  const Double_t Default_Hole_Mobility=480.0; //[cm2/Vs] Hole mobility
  const Double_t Default_Hole_D=12; //;// Hole diffusion [cm2/s]

  const Double_t Default_Electron_Mobility=1415.0; //[cm2/Vs] Electron mobility
  const Double_t Default_Electron_D=36; //;// Electron diffusion [cm2/s]



  // ---------------- Mobility depends on the electric field ---------------- //
  Double_t vm=0.0;
  Double_t Ec=0.0;
  Double_t beta=0.0;
  Double_t temperature_T=300.0;
  // ---------------- End  ---------------- //


  //// ============PARAMETERS TO ADJUST================
  G4cout << "nilou: resistivity= " << resistivity << G4endl;
  ///=================================================
  Double_t mobility_const=0.0;
  Double_t diffusion_const=0.0;
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
  ///=================================================


  G4cout << "vm=" << vm << G4endl;
  G4cout << "Ec=" << Ec << G4endl;
  G4cout << "beta=" << beta << G4endl;

  for(G4int itr  = 0 ; itr < nEntries ; itr++)
    {
      tempPixel.first  = (*hitsCollection)[itr]->GetPixelNbX();
      tempPixel.second = (*hitsCollection)[itr]->GetPixelNbY();
      //pixelsContent[tempPixel] += (*hitsCollection)[itr]->GetEdep();
      
      G4double xpos=(*hitsCollection)[itr]->GetPosWithRespectToPixel().x();
      G4double ypos=(*hitsCollection)[itr]->GetPosWithRespectToPixel().y();
      G4double zpos=(*hitsCollection)[itr]->GetPosWithRespectToPixel().z()+thickness/2.0;

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
      G4cout << "depletionWidth=" << depletionWidth << "[cm]" << G4endl;
   
      pair<G4int, G4int> extraPixel;
      extraPixel = tempPixel;
      G4double hit_energy=(*hitsCollection)[itr]->GetEdep();
      G4cout << "V bias=" << V_B << G4endl;
      
      //------- noise 1 ------- //
      G4cout << "hit_energy before=" << hit_energy/keV << "[keV]" << G4endl;
      hit_energy = elec*CLHEP::RandGauss::shoot(hit_energy/elec,3*TMath::Sqrt(hit_energy/elec));
      G4cout << "hit_energy after=" << hit_energy/keV << "[keV]" << G4endl;
      //------- End noise 1 ------- //
      
      if(zpos<depletionWidth*10) // Only charge sharing for the depletion region //depletionWidth*10 [mm]
      	{
	  Double_t electric_field=-((V_B-V_D)/(thickness/cm)+(1-(zpos/cm)/(thickness/cm))*2*V_D/(thickness/cm));
	  
	  Double_t non_const_mobility=(vm/Ec)/(TMath::Power((1+TMath::Power(TMath::Abs(electric_field)/Ec, beta)), 1.0/beta));
	  Double_t drift_time=(zpos/cm)/(non_const_mobility*TMath::Abs(electric_field)); //non-constant mobility
	  Double_t diffusion_RMS=TMath::Sqrt(2.0*diffusion_const*drift_time); //[cm]
	  diffusion_RMS=diffusion_RMS*10;//[mm]
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
      Double_t noise=CLHEP::RandGauss::shoot(0, 100);//electrons
      pixelsContent[tempPixel]+=noise*elec/MeV;
    }
  // ------------ End Add gaussian noise ------------ //


  //------------------ RECORD DIGITS ------------------//
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
	  m_digitIn.thl=ThresholdMatrix[0][0]; //keV
	  a=SurrogateA[0][0];
	  b=SurrogateB[0][0];
	  c=SurrogateC[0][0];
	  t=SurrogateD[0][0];
	}
      else
	{
	  m_digitIn.thl=ThresholdMatrix[tempPixel.first][tempPixel.second]; //keV
	  a=SurrogateA[tempPixel.first][tempPixel.second];
	  b=SurrogateB[tempPixel.first][tempPixel.second];
	  c=SurrogateC[tempPixel.first][tempPixel.second];
	  t=SurrogateD[tempPixel.first][tempPixel.second];
	}
      
      G4cout << "*************Energy =" << (*pCItr).second/keV << " [keV]" << G4endl;
      G4cout << "Threshold=" <<  m_digitIn.thl << " [keV]" << G4endl;

      G4int TOT=a*((*pCItr).second)/keV+b-c/(((*pCItr).second/keV)-t);
      G4cout << "a=" << a << ", b=" << b << ", c=" << c << ", t=" << t << G4endl;


      
      AllPixNoThreshDigit * digit = new AllPixNoThreshDigit;
      
      digit->SetPixelIDX((*pCItr).first.first);
      digit->SetPixelIDY((*pCItr).first.second);
      digit->SetPixelEnergyDep(((*pCItr).second)/keV); //Energy with charge sharing
      digit->SetPixelCounts(TOT); //TOT value
      if (pixelsContent_MC.count(tempPixel))
	{
	      //G4cout << "x=" << tempPixel.first << ", y=" << tempPixel.second << ", energyMC=" << (pixelsContent_MC[tempPixel].MC_energy)/keV << G4endl;
	  digit->SetPixelEnergyMC((pixelsContent_MC[tempPixel].MC_energy)/keV); //MC value
	  digit->Set_posX_WithRespectoToPixel(pixelsContent_MC[tempPixel].posX_WithRespectToPixel);
	  digit->Set_posY_WithRespectoToPixel(pixelsContent_MC[tempPixel].posY_WithRespectToPixel);
	  digit->Set_posZ_WithRespectoToPixel(pixelsContent_MC[tempPixel].posZ_WithRespectToPixel);
	  pixelsContent_MC.erase(tempPixel);
	}
      //digit->IncreasePixelCounts(); // Counting mode
      
      m_digitsCollection->insert(digit);
    }

  //MC only
  map<pair<G4int, G4int>, MC_content >::iterator pCItr_MC = pixelsContent_MC.begin();
  for( ; pCItr_MC != pixelsContent_MC.end() ; pCItr_MC++)
    {
      G4cout << "x=" << tempPixel.first << ", y=" << tempPixel.second << ", energyMC=" << (pixelsContent_MC[tempPixel].MC_energy)/keV << G4endl;
      pair<G4int, G4int> tempPixel;
      tempPixel.first=(*pCItr_MC).first.first;
      tempPixel.second=(*pCItr_MC).first.second;

      AllPixNoThreshDigit * digit = new AllPixNoThreshDigit;

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


void AllPixNoThreshDigitizer::initialise_calibrationMatrices(int size_calibX, int size_calibY)
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


G4double AllPixNoThreshDigitizer::IntegrateGaussian(G4double xhit,G4double yhit,G4double Sigma, G4double x1, G4double x2, G4double y1, G4double y2, G4double Energy)
{
  G4double Integral=(-TMath::Erf((x1-xhit)/(TMath::Sqrt(2.)*Sigma))+TMath::Erf((x2-xhit)/(TMath::Sqrt(2.)*Sigma)))*(-TMath::Erf((y1-yhit)/(TMath::Sqrt(2.)*Sigma))+TMath::Erf((y2-yhit)/(TMath::Sqrt(2.0)*Sigma)));
  
  G4double energybis=Integral*Energy/4.0; //*(TMath::Pi())*(TMath::Pi());
  // G4double energybis= (Energy*(-TMath::Erf((x1-xhit)/(TMath::Sqrt(2.)*Sigma))+TMath::Erf((x2-xhit)/(TMath::Sqrt(2.)*Sigma)))
  // 			*(-TMath::Erf((y1-yhit)/(TMath::Sqrt(2.)*Sigma))+TMath::Erf((y2-yhit)/(TMath::Sqrt(2.0)*Sigma))))/4.0*(TMath::Pi())*(TMath::Pi());
  return energybis;
}
