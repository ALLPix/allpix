/**
 *  Author:
 *    Mathieu Benoit <mathieu.benoit@cern.ch>
 *
 *  allpix Authors:
 *   John Idarraga <idarraga@cern.ch>
 *   Mathieu Benoit <benoit@cern.ch>
 */

#include "AllPixTimepixDigitizer.hh"
#include "AllPixTrackerHit.hh"
#include "G4DigiManager.hh"
#include "G4VDigitizerModule.hh"
#include "AllPixGeoDsc.hh"

#include "TMath.h"
#include "TF1.h"
#include "CLHEP/Random/RandGauss.h"
#include "CLHEP/Random/RandFlat.h"

#define CALIBRATION_CLOCK_UNIT 10.416666666666e-9
#define dE 5. // 50V/cm cm supposed to be equal to 10

using namespace TMath;
AllPixTimepixDigitizer::AllPixTimepixDigitizer(G4String modName, G4String hitsColName, G4String digitColName) 
: AllPixDigitizerInterface (modName) {

	
	doFastErf= true;

	if(doFastErf){
	ifstream erfFile("share/Erf.dat");

	while(!erfFile.eof()){

		G4double XT,YT;
		erfFile >> XT >> YT;
		ErrorFunction[XT]=YT;
	}
	erfFile.close();
	}

//	for(G4double i= -4 ; i< 4 ; i+=0.01){
//		cout << "ERF Test " << i << " " << MyErf(i) << endl;
//	}


	hitindex=0;
	
	nseek=1;
	
	doAnimation = false;
	
	
	// Registration of digits collection name
	collectionName.push_back(digitColName);
	m_hitsColName.push_back(hitsColName);

	//// Unit for charge in FEIX average e/h pair creation energy in Silicon
	elec = 3.64*eV;
	
	readoutType = HOLE;


	// Example of detector description handle
	// provided by the interface
	AllPixGeoDsc * gD = GetDetectorGeoDscPtr();

	//////////////////////////
	// Bias and Temperature //
	//////////////////////////

	biasVoltage=80.0; //[V]
	Temperature = 300.0;
	detectorThickness = gD->GetSensorZ();
	resistivity=gD->GetResistivity();

	// true = p-type
	// false = n-type
	bulkType =false;


	////////////////////////////////
	// Geometry Related constants //
	////////////////////////////////
	pitchX = gD->GetPixelX();
	pitchY = gD->GetPixelY();
	nPixX= gD->GetNPixelsX();
	nPixY= gD->GetNPixelsY();

	offsetX= gD->GetSensorXOffset()*mm;
	offsetY= gD->GetSensorYOffset()*mm;

	//cout << TString::Format("[Timepix Digitizer] Lx=%f Ly=%f npixX=%d npixY=%d, offsetX=%f offsetY=%f",pitchX/um,pitchY/um,nPixX,nPixY,offsetX/um,offsetY/um) << endl;


	/////////////////////////////////
	// Magnetic Field 			   //
	////////////////////////////////
	B_Field = 5.0;
	r_H_e = 1.1;
	r_H_h=0.7;	




	///////////////////////////////////////////////////
	// Silicon electron and hole transport constants //
	///////////////////////////////////////////////////

	// Default mobilities
	Default_Electron_Mobility=1415.0*cm2/s; // Electron mobility (cm2/Vs)
	Default_Hole_Mobility=450.0*cm2/s;// Hole mobility (cm2/Vs
	Default_Electron_D=36.62*cm2/s; // Electron mobility (cm2/s)
	Default_Hole_D=18*cm2/s;// Hole mobility (cm2/s



	//mobility dependence on electric field
	Electron_AlphaField = 2.4e7*cm/s ; //[um/ns]
	Electron_ThetaField = 0.8;
	Electron_TempNominal = 600.0 ; // [K]
	Electron_Beta = 2.0;
	Electron_Saturation_Velocity = Electron_AlphaField*TMath::Power(1.0+Electron_ThetaField*TMath::Exp(Temperature/Electron_TempNominal),-1.0);
	Hole_AlphaField = 2.4e7*cm/s ; //[um/ns]
	Hole_ThetaField = 0.8;
	Hole_TempNominal = 600.0 ; // [K]
	Hole_Beta = 1.0 ;
	Hole_Saturation_Velocity = Hole_AlphaField*TMath::Power(1.0+Hole_ThetaField*TMath::Exp(Temperature/Hole_TempNominal),-1.0);


	//for(double E=0;E<=100000/cm;E+=50/cm){
	for(int i = 0 ; i<=2001;i++){
		EMobilityHoleLUT[i]=i*dE;
		MobilityHoleLUT[i]=Default_Hole_Mobility*
				TMath::Power((1.0/(1.+ TMath::Power((Default_Hole_Mobility*i*dE)/
						Hole_Saturation_Velocity,Hole_Beta))),1.0/Hole_Beta);
	}




	//////////////////////
	// Radiation damage //
	//////////////////////

	fluence = 0.0/cm2;
	//Beta_electrons = 5e-16*cm2/ns;
	Beta_electrons = 5e-16*cm2/ns;
	b=7.94e-3/cm;
	c=3.54e-13*cm2;
	if(fluence!=0.0)trappingTime = 1.0/(Beta_electrons*fluence);
	else trappingTime=1000*s;

	echarge = 1.60217646e-19;
	epsilon = 11.8*8.854187817e-12/m;

	Neff=0;
	double Neff0=0;
	if(bulkType) Neff0=-1.0/(resistivity*Default_Hole_Mobility*(s/cm2)*echarge*cm3);
	else  Neff0=1.0/(resistivity*Default_Electron_Mobility*(s/cm2)*echarge*cm3);

	if(!bulkType)Neff=Neff0*TMath::Exp(-c*fluence)-b*fluence;
	else Neff=Neff0-b*fluence;

	depletionVoltage=echarge*TMath::Abs(Neff)*detectorThickness*detectorThickness/(2*epsilon);
	depletedDepth=TMath::Sqrt(2*epsilon*biasVoltage/(echarge*TMath::Abs(Neff)));
	if(depletedDepth>detectorThickness)depletedDepth=detectorThickness;
	electricFieldZ = biasVoltage/depletedDepth; // V/um
	electricFieldX = 0; // V/um
	electricFieldY = 0; // V/um

	//G4cout << "!!!!!!!Radiation Damage Report !!!!!!!!" << endl
	//	 << TString::Format("depletionVoltage : %f depleted Depth : %f",depletionVoltage,depletedDepth/um) << endl
	//	 << TString::Format("Neff : %e Neff0 : %e ",Neff*cm3,Neff0*cm3);
	//if(bulkType==false && Neff<0) cout << "Bulk inverted" << endl;






	//////////////////////
 	// physics switches //
	//////////////////////

 	doTrapping =false;
 	doFullField = true;


 	////////////////////////////////////
 	// Numerical integration accuracy //
 	////////////////////////////////////
 	Target = 1e-4 ;
 	tlow =  0.001*ns;
 	tup=	1*ns;
 	dtIni = 0.01*ns;

 	precision = 1;


 	/////////////////////////////
 	// Detector Chip Selection //
 	/////////////////////////////
 	/*


 	 3 = FEI3
 	 4 = FEI4
 	 666 = Omegapix2
 	 Default = FEI3

 	 */

// 	double calib_A=1;
// 	double calib_B=1;
// 	double calib_C=1;


 	MipTOT = gD->GetMIPTot();
 	MipCharge = gD->GetMIPCharge();
 	
 	
 
 // Initializing Pixel Per Pixel preamp characteristics 
 
 
 	Gain = new G4double*[nPixX];
 	for(int i=0;i<nPixX;i++){
 		Gain[i] = new G4double[nPixY];
 	};
 
  	RisingSlope = new G4double*[nPixX];
 	for(int i=0;i<nPixX;i++){
 		RisingSlope[i] = new G4double[nPixY];
 	};	

  	FallingSlope = new G4double*[nPixX];
 	for(int i=0;i<nPixX;i++){
 		FallingSlope[i] = new G4double[nPixY];
 	};	
 	
  	ThresholdMatrix = new G4double*[nPixX];
 	for(int i=0;i<nPixX;i++){
 		ThresholdMatrix[i] = new G4double[nPixY];
 	};	
 	
 	CounterDepth = gD->GetCounterDepth();
 	ClockUnit = gD->GetClockUnit();
 	ChipNoise = gD->GetChipNoise()*elec;
 	m_digitIn.thl = gD->GetThreshold()*elec;
 	SaturationEnergy = gD->GetSaturationEnergy()*keV;
 	chargeSharingConstant = gD->GetCrossTalk();
 	
 	DefaultGain = 10.;
 	DefaultRisingSlope = 20.0/1e-9; 
 	DefaultFallingSlope= GetIkrum(MipCharge*elec,MipTOT);
 	
 	G4double GainDispersion = 0.005;
 	G4double RisingSlopeDispersion = 0.005;
 	G4double FallingSlopeDispersion = 0.005;
 	
 	
  	for(int i=0;i<nPixX;i++){
  	for(int j=0;j<nPixY;j++){
  			Gain[i][j]=CLHEP::RandGauss::shoot(DefaultGain,GainDispersion*DefaultGain);
 			RisingSlope[i][j]=CLHEP::RandGauss::shoot(DefaultRisingSlope,RisingSlopeDispersion*DefaultRisingSlope);
 			FallingSlope[i][j]=CLHEP::RandGauss::shoot(DefaultFallingSlope,FallingSlopeDispersion*DefaultGain);
 	  		ThresholdMatrix[i][j]=CLHEP::RandGauss::shoot(m_digitIn.thl,ChipNoise);
  			};};
 	

  	A= 0.50554e+01;
  	B= 4.75150e+02;
    C=-8.54405e+02;
   	D= 4.01736e+00;

    double Surrogate_Parameter_dispertion = 0.001;


  	SurrogateA = new G4double*[nPixX];
 	for(int i=0;i<nPixX;i++){
 		SurrogateA[i] = new G4double[nPixY];
 	};
 	SurrogateB = new G4double*[nPixX];
 	for(int i=0;i<nPixX;i++){
 		SurrogateB[i] = new G4double[nPixY];
 	};
 	SurrogateC = new G4double*[nPixX];
 	for(int i=0;i<nPixX;i++){
 		SurrogateC[i] = new G4double[nPixY];
 	};
 	SurrogateD = new G4double*[nPixX];
 	for(int i=0;i<nPixX;i++){
 		SurrogateD[i] = new G4double[nPixY];
 	};

 	doRealCalibrationFile = false;


    if (doRealCalibrationFile){

      	for(int i=0;i<nPixX;i++){
      	for(int j=0;j<nPixY;j++){
      			SurrogateA[i][j]=0;
     			SurrogateB[i][j]=0;
     			SurrogateC[i][j]=0;
     	  		SurrogateD[i][j]=1e6;
      			};};

      	TFile *calib = TFile::Open("share/FitParameters_v4.root","open");
    	TTree *t =(TTree*)calib->Get("fitPara");
    	t->SetMakeClass(1);

    	Float_t At;
    	Float_t Bt;
    	Float_t Ct;
    	Float_t Dt;
    	Int_t Xt;
    	Int_t Yt;

    	Float_t At_err;
    	Float_t Bt_err;
    	Float_t Ct_err;
    	Float_t Dt_err;
    	Float_t Chi;


    	t->SetBranchAddress("pixx",&Xt);
    	t->SetBranchAddress("pixy",&Yt);
    	t->SetBranchAddress("a",&At);
    	t->SetBranchAddress("b",&Bt);
    	t->SetBranchAddress("c",&Ct);
    	t->SetBranchAddress("d",&Dt);

    	t->SetBranchAddress("a_err",&At_err);
    	t->SetBranchAddress("b_err",&Bt_err);
    	t->SetBranchAddress("c_err",&Ct_err);
    	t->SetBranchAddress("d_err",&Dt_err);

    	t->SetBranchAddress("chi2ndf",&Chi);

    	G4int nevents = t->GetEntries();

    	for (int i =0 ; i<nevents; i++){

    		t->GetEntry(i);

    		SurrogateA[Xt][Yt]=At;
    		SurrogateB[Xt][Yt]=Bt;
    		SurrogateC[Xt][Yt]=Ct;
    		SurrogateD[Xt][Yt]=Dt;

    		};

    }
    else {

  	for(int i=0;i<nPixX;i++){
  	for(int j=0;j<nPixY;j++){
  			SurrogateA[i][j]=CLHEP::RandGauss::shoot(A,Surrogate_Parameter_dispertion*A);
 			SurrogateB[i][j]=CLHEP::RandGauss::shoot(B,Surrogate_Parameter_dispertion*B);
 			SurrogateC[i][j]=CLHEP::RandGauss::shoot(C,Surrogate_Parameter_dispertion*C);
 	  		SurrogateD[i][j]=CLHEP::RandGauss::shoot(D,Surrogate_Parameter_dispertion*D);
  			};};

    };




}

AllPixTimepixDigitizer::~AllPixTimepixDigitizer(){

}

void AllPixTimepixDigitizer::SetDetectorDigitInputs(G4double thl){

	// set digitization input values
	// thl
	m_digitIn.thl = thl; // <-- input !
}



G4double AllPixTimepixDigitizer::MyErf(G4double x){

	G4double yl;

	if(doFastErf){
		if(x<ErrorFunction.begin()->first)return -1;
		else if(x>-ErrorFunction.begin()->first)return 1;

		else{

		map<G4double,G4double>::iterator itrl = ErrorFunction.lower_bound(x);
		//map<G4double,G4double>::iterator itrh = ErrorFunction.upper_bound(x);

	   // xl = itrl->first;
	   // xh= itrh->first;
		yl= itrl->second;
	   // yh= itrh->second;

	  //  cout << TString::Format("xl:%f xh:%f yl:%f yh:%f",xl,xh,yl,yh) << endl;

		//return (x-xl)*((yh-yl)/(xh-xl)) + yl;
		return yl;
		}
	}
	else{
		yl=TMath::Erf(x);
		return yl;
	}

}


G4double AllPixTimepixDigitizer::GetIkrum(double energy,double Target){
	
	G4double Max=0;
	if(energy>SaturationEnergy)  Max= (1.0/keV)*SaturationEnergy*DefaultGain;
	else	 Max= (1.0/keV)*energy*DefaultGain;

	G4double tpeak = Max/DefaultRisingSlope;
	
	G4double a = TMath::FloorNint(tpeak/ClockUnit);
	
/* 	G4cout <<endl << "[IKRum Calculation]" << endl ; 
	G4cout << "[IKRum Calculation] DefaultGain : " << DefaultGain << endl;
	G4cout << "[IKRum Calculation] energy : " << energy/keV << endl;
	G4cout << "[IKRum Calculation] DefaultRisingSlope : " << DefaultRisingSlope << endl;
	G4cout << "[IKRum Calculation] ClockUnit : " << ClockUnit << endl;
	G4cout << "[IKRum Calculation] Max : " << Max << endl;
	G4cout << "[IKRum Calculation] SE : " << SaturationEnergy << endl; */

	
	return (Max-(1.0/keV)*m_digitIn.thl*DefaultGain)/(ClockUnit*(Target-a));
}


vector<G4double>  AllPixTimepixDigitizer::ComputeDriftTimeFullField(G4double x, G4double y, G4double z,G4double energy)
{
	G4double driftTime=0;
	G4double dt=dtIni;
	G4double xtemp=x;
	G4double ytemp=y;
	G4double ztemp=z;
	//cout << "z : " << z/um << endl;
	vector<G4double> step;
	int counter=0;
	
	vector<G4double> vx,vy,vz,vt; 
	
	
	hitindex++;
	
	int iter =0;
	while( ztemp>=-detectorThickness/2 && ztemp<=detectorThickness/2 && iter<1000){

		if(doAnimation){
		vx.push_back(xtemp/um+pixelPositionWithRegardToCorner_x/um);
		vy.push_back(ytemp/um+pixelPositionWithRegardToCorner_y/um);
		vz.push_back(detectorThickness-ztemp/um);
		vt.push_back(driftTime/s);
		};
		
		driftTime+=dt;
		
		if(B_Field!=0.0){
			
			if(readoutType==ELECTRON){
			step=RKF5IntegrationBFieldElectrons(xtemp,ytemp,ztemp,dt);
			}
			else{
			step=RKF5IntegrationBFieldHoles(xtemp,ytemp,ztemp,dt);				
			}
		}
		else {
			if(readoutType==ELECTRON){
			step=RKF5IntegrationElectrons(xtemp,ytemp,ztemp,dt);
			}
			else{
			step=RKF5IntegrationHoles(xtemp,ytemp,ztemp,dt);				
			}
		};
		
		
	
		xtemp+=step[0];
		ytemp+=step[1];
		ztemp+=step[2];
		

		if(GetElectricFieldNorm(xtemp,ytemp,ztemp)==0)ztemp=detectorThickness;
		dt=SetDt(dt,step[3]);
		
		
		counter++;
		
		iter++;
		//cout << "!!!!!! Drifting youhou! !!!!! " << TString::Format("x: %f y : %f z : %f ",xtemp/um,ytemp/um,ztemp/um) << endl;
	}
	
	
	
	if(doAnimation)anim->AddTrack(vx,vy,vz,vt,ComputeDiffusionRMS(driftTime)/um,energy);
	//cout << "[animation]" << ComputeDiffusionRMS(driftTime)/um << endl;

	vector<G4double> output(4);
	output[0]=xtemp;
	output[1]=ytemp;
	output[2]=ztemp;
	output[3]=driftTime;

	return output;
}


void AllPixTimepixDigitizer::ComputeElectricField(G4double /*x*/, G4double /*y*/, G4double z){

		Efield1D(z);
}

void AllPixTimepixDigitizer::ComputeEffectiveElectricField(G4double x, G4double y, G4double z){

		Efield1D(z);
		electricFieldX =electricFieldZ*TMath::Tan(ComputeLorentzAngle(x,y,z)) ;
}

G4double AllPixTimepixDigitizer::GetElectricFieldNorm(G4double /*x*/, G4double /*y*/, G4double /*z*/){

	return TMath::Sqrt(electricFieldX*electricFieldX +electricFieldY*electricFieldY +electricFieldZ*electricFieldZ );
}

G4double  AllPixTimepixDigitizer::MobilityElectron(G4double /*x*/, G4double /*y*/, G4double /*z*/){

	//ComputeElectricField(x,y,z);
	G4double parElectricField = electricFieldZ;//GetElectricFieldNorm(x,y,z);
	G4double mobilite = Default_Electron_Mobility*
	TMath::Power((1.0/(1.+ TMath::Power((Default_Electron_Mobility*parElectricField)/
			Electron_Saturation_Velocity,Electron_Beta))),1.0/Electron_Beta);
	return mobilite;
	//cout << "!!!!!!!!!" << mobilite << " " << parElectricField*cm << endl;

}


G4double AllPixTimepixDigitizer::MobilityHole(G4double x, G4double y, G4double z){

	ComputeElectricField(x,y,z);
	G4double parElectricField = GetElectricFieldNorm(x,y,z);

	if(parElectricField<0) return Default_Hole_Mobility;
	else if (parElectricField>EMobilityHoleLUT[2000]) parElectricField=EMobilityHoleLUT[2000];


	//map<G4double,G4double>::iterator itrl =MobilityHoleLUT.lower_bound(parElectricField);
	//G4double mobilite = itrl->second;
	G4double mobilite = MobilityHoleLUT[TMath::FloorNint(parElectricField/(dE))];
	return mobilite;


}


G4double AllPixTimepixDigitizer::ComputeDriftTimeUniformField(AllPixTrackerHit * /*a_hit*/)
{
	// drift in a uniform electric field from hit position to surface
	//G4cout << "Hit position with respect to pixel (um) "<< (a_hit->GetPosWithRespectToPixel().z() + detectorThickness/2.0)/um  <<endl;
	//return (a_hit->GetPosWithRespectToPixel().z() + detectorThickness/2.0 ) /(mobility*electricFieldY);
	double drift = ( CLHEP::RandGauss::shoot((depletedDepth/detectorThickness)*detectorThickness/2,(depletedDepth/detectorThickness)*detectorThickness/4)) /(mobility*electricFieldZ);
	//cout << "drift: " << mobility*s/cm2 << endl;
	return drift;
}

G4double AllPixTimepixDigitizer::ComputeDiffusionRMS(G4double tDrift)
{
	if(readoutType==ELECTRON){
		return TMath::Sqrt(2.*Default_Electron_D*tDrift);
	}

	return TMath::Sqrt(2.*Default_Hole_D*tDrift);


}

G4double AllPixTimepixDigitizer::ApplyTrapping(G4double tDrift, G4double Energy)
{

	//cout << TString::Format("[Trapping!!!] Energy:%f after Trapping:%f",Energy, Energy*TMath::Exp(-tDrift/trappingTime)) << endl;
	return Energy*TMath::Exp(-tDrift/trappingTime);

}



G4double AllPixTimepixDigitizer::IntegrateGaussian(G4double xhit,G4double yhit,G4double Sigma, G4double x1, G4double x2, G4double y1, G4double y2, G4double Energy )
{
	//G4cout << TString::Format("Integration borns xhit=%f yhit=%f x1=%f x2=%f y1=%f y2=%f",xhit/1e3,yhit/1e3,x1/1e3,x2/1e3,y1/1e3,y2/1e3) << endl;
	//G4cout << TString::Format("TMath::Erf test %f %f %f %f",TMath::Erf((x1 - xhit)/(Sqrt(2.)*Sigma)),TMath::Erf((x2 - xhit)/(TMath::Sqrt(2.)*Sigma)),TMath::Erf((y1 - yhit)/(TMath::Sqrt(2.)*Sigma)),TMath::Erf((y2 - yhit)/(TMath::Sqrt(2.)*Sigma)))<< endl;
	double energybis= (Energy*(-MyErf((x1 - xhit)/(Sqrt(2.)*Sigma)) + MyErf((x2 - xhit)/(TMath::Sqrt(2.)*Sigma)))
			*(-MyErf((y1 - yhit)/(TMath::Sqrt(2.)*Sigma)) + MyErf((y2 - yhit)/(TMath::Sqrt(2.0)*Sigma))))/4.;

	//G4cout << "energy after integration : " << energybis << endl;
	return energybis;

}


void AllPixTimepixDigitizer::Efield1D(G4double z){


	if(readoutType==ELECTRON)
	{
		electricFieldZ=(biasVoltage-depletionVoltage)/detectorThickness+(1-z/detectorThickness)*2*depletionVoltage/detectorThickness;
		//electricFieldZ=(detectorThickness-z)*biasVoltage/depletedDepth;
	}
	else {
		electricFieldZ=-(biasVoltage-depletionVoltage)/detectorThickness+(1-z/detectorThickness)*2*depletionVoltage/detectorThickness;
		//electricFieldZ=-(detectorThickness-z)*biasVoltage/depletedDepth;
	}
	
	if(z>depletedDepth){
		electricFieldZ=0;
		//cout << "bip" << endl;
	}

}




/**
 * Digitizer for FEI3Standard
 */
void AllPixTimepixDigitizer::Digitize(){

	// Create the digits collection
	m_digitsCollection = new AllPixTimepixDigitsCollection("AllPixTimepixDigitizer", collectionName[0] );

	// Get a pointer to the digiManager
	G4DigiManager * digiMan = G4DigiManager::GetDMpointer();

	// Get the hit collection ID
	G4int hcID = digiMan->GetHitsCollectionID(m_hitsColName[0]);

	// And fetch the Hits Collection
	AllPixTrackerHitsCollection * hitsCollection = 0;
	hitsCollection = (AllPixTrackerHitsCollection*)(digiMan->GetHitsCollection(hcID));

	// Temporary data structure to store hits
	//  collection information
	map<pair<G4int, G4int>, G4double > pixelsContent;
	pair<G4int, G4int> tempPixel;

	// Loop over the whole Hits Collection
	G4int nEntries = hitsCollection->entries();


	//TF1 * f1 = new TF1("expint","exp((x-[0])*(x-[0])/2.*[1])",-100.,100.);

	G4double hitsETotal = 0.;
	
	
	G4double sx=0;
	G4double sy=0;
	
	if(doAnimation){
	//Create a Digit Graphical representation 
	anim = new AllPixDigitAnimation(20,20,detectorThickness/um,pitchX/um,pitchY/um,nEntries,hcID);
	
	
	//////////////////////////////////////////////////////////////////////	
	//Shift for event display (to recenter the hit around trhe geometry)//
	//////////////////////////////////////////////////////////////////////
	for(G4int itr  = 0 ; itr < nEntries ; itr++) {
		
		tempPixel.first  = (*hitsCollection)[itr]->GetPixelNbX();
		tempPixel.second = (*hitsCollection)[itr]->GetPixelNbY();
		G4double xpos = (*hitsCollection)[itr]->GetPosWithRespectToPixel().x();
		G4double ypos = (*hitsCollection)[itr]->GetPosWithRespectToPixel().y();
		pixelPositionWithRegardToCorner_x=tempPixel.first*pitchX -pitchX/2.;
		pixelPositionWithRegardToCorner_y=tempPixel.second*pitchY -pitchY/2. ;
			
		sx+=xpos+pixelPositionWithRegardToCorner_x;
		sy+=ypos+pixelPositionWithRegardToCorner_y;
	};
	sx/=nEntries;
	sx = TMath::FloorNint(sx/pitchX)*pitchX;
	sy/=nEntries;
	sy = TMath::FloorNint(sy/pitchY)*pitchY;	
	
	anim->SetShift(sx/(um),sy/(um));
	///////////////////////////////////////////////////////////////////////
	};
	
	for(G4int itr  = 0 ; itr < nEntries ; itr++) {

		G4double eHitTotal;

		//G4double dispersion = 4*TMath::Sqrt((*hitsCollection)[itr]->GetEdep()/elec) + 40000.0/((*hitsCollection)[itr]->GetEdep()/elec -m_digitIn.thl/elec);

		//eHitTotal = elec*CLHEP::RandGauss::shoot((*hitsCollection)[itr]->GetEdep()/elec,dispersion);
		//eHitTotal += elec*CLHEP::RandGauss::shoot(ChipNoise/elec,50);
//		if ((*hitsCollection)[itr]->GetEdep()>30*keV) {
//		 eHitTotal = elec*CLHEP::RandGauss::shoot((*hitsCollection)[itr]->GetEdep()/elec,1.2*keV);
//
//		}
//		else
//		{
		//Fe55
		// eHitTotal = elec*CLHEP::RandGauss::shoot((*hitsCollection)[itr]->GetEdep()/elec,8.6*TMath::Sqrt((*hitsCollection)[itr]->GetEdep()/elec));

		eHitTotal = elec*CLHEP::RandGauss::shoot((*hitsCollection)[itr]->GetEdep()/elec,3*TMath::Sqrt((*hitsCollection)[itr]->GetEdep()/elec));
//		}
		//G4double eHitTotal = CLHEP::RandGauss::shoot((*hitsCollection)[itr]->GetEdep(),0.6*keV);


		hitsETotal += eHitTotal;

		//under-depletion
		eHitTotal = eHitTotal*depletedDepth/detectorThickness;

		tempPixel.first  = (*hitsCollection)[itr]->GetPixelNbX();
		tempPixel.second = (*hitsCollection)[itr]->GetPixelNbY();



		G4double xpos = (*hitsCollection)[itr]->GetPosWithRespectToPixel().x();
		G4double ypos = (*hitsCollection)[itr]->GetPosWithRespectToPixel().y();
		//FIXME geom should be switched !!!
		G4double zpos = -(*hitsCollection)[itr]->GetPosWithRespectToPixel().z();
		//G4cout << "z position of hit : " << zpos/um << endl;

		pixelPositionWithRegardToCorner_x=tempPixel.first*pitchX -pitchX/2.;
		pixelPositionWithRegardToCorner_y=tempPixel.second*pitchY -pitchY/2. ;
		

		

		//G4cout << TString::Format("[TimepixDigi] hit position x,y,z : %5.5f %5.5f %5.5f",xpos/um,ypos/um,zpos/um)<<endl;
		//G4cout << TString::Format("[TimepixDigi] hit position Nx,Ny : %d %d %f ",tempPixel.first,tempPixel.second,eHitTotal/keV)<<endl;

		for(G4int nQ  = 0 ; nQ < precision ; nQ++) {

		double eHit = double(eHitTotal)/precision;

		//G4cout << "[before digi] x : " << tempPixel.first << " ,  y : " << tempPixel.second << ", E = " << eHit/keV << G4endl;

		// Hit Info
		//Ugly Hack !!!
		//if(xpos>0)xpos=-pitchX/2+xpos;
		//else xpos=pitchX/2+xpos;

		//G4double zpos = CLHEP::RandGauss::shoot(zpos,10*um);


		G4double driftTime;
		sigma = 0;
		if(doFullField==false){
			 driftTime = ComputeDriftTimeUniformField((*hitsCollection)[itr]);
			 sigma = ComputeDiffusionRMS(driftTime);
			 //cout << TString::Format("!!!!!!!!! vd/vdep : %f drift time : %f sigma : %f",depletedDepth/detectorThickness,driftTime,sigma) << endl;
		}
		else{

			// Until we get out position
			vector<G4double> data = ComputeDriftTimeFullField(xpos,ypos,zpos,eHit);
			driftTime = data[3];
			sigma = ComputeDiffusionRMS(driftTime);
			xpos=data[0];
			ypos=data[1];
			//G4cout << TString::Format("!!!!!!!!! vd/vdep : %f drift time : %f sigma : %f",depletedDepth/detectorThickness,driftTime,sigma) << endl;

		}

		// see if need more energy in neighbor
		pair<G4int, G4int> extraPixel;
		extraPixel = tempPixel;

		//loop over all hit pixel neighbors
		//G4cout << TString::Format("[MC Truth] Hit Energy=%f, x=%f, y=%f pixel(%i,%i)",eHitTruth/elec,xpos/um,ypos/um,tempPixel.first,tempPixel.second) << endl;
		//if(fabs(xpos)>=pitchX/2.-10*sigma && fabs(ypos)>pitchY/2.0-10*sigma){

		if( (fabs(xpos) >= pitchX/2.-3*sigma || fabs(ypos) >=pitchY/2.-3*sigma) ){

//		if( (fabs(xpos) >= pitchX/2.-3*sigma || fabs(ypos) >=pitchY/2.-3*sigma) ){
		  // G4cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*****!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" <<endl ;
		  for(int i=-nseek;i<=nseek;i++){
			for(int j=-nseek;j<=nseek;j++){
				extraPixel = tempPixel;
				extraPixel.first +=i;
				extraPixel.second+=j;
				if(extraPixel.first >= 0 && extraPixel.second>=0 && extraPixel.first < nPixX && extraPixel.second < nPixY)
					{
					//We compute contribution of the hit to each pixels

					//double Etemp = IntegrateGaussian(xpos/nm+offsetX/nm,ypos/nm+offsetY/nm,sigma/nm,(-pitchX/2.0 + i*pitchX)/nm,(-pitchX/2.+(i+1)*pitchX)/nm,(-pitchY/2 + j*pitchY)/nm,(-pitchY/2 + (j+1)*pitchY)/nm, eHit );
					double Etemp = IntegrateGaussian(xpos/nm,ypos/nm,sigma/nm,(-pitchX/2.0 + i*pitchX)/nm,(-pitchX/2.+(i+1)*pitchX)/nm,(-pitchY/2 + j*pitchY)/nm,(-pitchY/2 + (j+1)*pitchY)/nm, eHit );

					if(doTrapping==true) pixelsContent[extraPixel]+=ApplyTrapping(driftTime,Etemp);
					else pixelsContent[extraPixel] +=Etemp;

					//G4cout << TString::Format("[Digitizer] Pixel %i %i Energy=%f, Energy after Trapping=%f",extraPixel.first,extraPixel.second,Etemp,ApplyTrapping(driftTime,pixelsContent[extraPixel])/elec) << endl;
					//if(Etemp/keV>0)cout << TString::Format("Pixel %i %i, Energy collected = %f sigma=%f tdrift=%f",extraPixel.first,extraPixel.second,Etemp/keV,sigma/um,driftTime/ns) << endl;

					};
			};
		};
		}

		else{
		  if(doTrapping==true)pixelsContent[extraPixel] +=ApplyTrapping(driftTime,pixelsContent[tempPixel]);
		  else pixelsContent[extraPixel] +=eHit;
		  //G4cout << TString::Format("[Digitizer] Pixel %i %i Energy=%f, Energy after Trapping=%f",extraPixel.first,extraPixel.second,pixelsContent[extraPixel]/elec,ApplyTrapping(driftTime,pixelsContent[extraPixel])/elec) << endl;

		};

		// FEI3 Chip crosstalk
		double sharedCharge = chargeSharingConstant*pixelsContent[tempPixel];
		pixelsContent[tempPixel]-= sharedCharge;
		int nPixelSharing=0;
		for(int i=-1;i<=1;i++){
					for(int j=-1;j<=1;j++){
						extraPixel = tempPixel;
						extraPixel.first +=i;
						extraPixel.second+=j;
						if(extraPixel.first >= 0 && extraPixel.second>=0 && extraPixel.first < nPixX && extraPixel.second < nPixY)	{
							nPixelSharing+=1;
							};
						};
		};
		for(int i=-1;i<=1;i++){
					for(int j=-1;j<=1;j++){
						extraPixel = tempPixel;
						extraPixel.first +=i;
						extraPixel.second+=j;
						if(extraPixel.first >= 0 && extraPixel.second>=0 && extraPixel.first < nPixX && extraPixel.second < nPixY)	{
							pixelsContent[extraPixel]+=sharedCharge/nPixelSharing;
							};
						};
				};

	}

	}

	//G4cout << "total = " << hitsETotal/keV << " keV" << G4endl;

	// Now create digits.  One per pixel
	map<pair<G4int, G4int>, G4double >::iterator pCItr = pixelsContent.begin();

	for( ; pCItr != pixelsContent.end() ; pCItr++)
	{
		// If the charge in a given pixel is over the threshold
//		G4cout << "pixel : " << (*pCItr).first.first << " , " << (*pCItr).first.second
//				<< ", E =  " << ((*pCItr).second)/keV << " keV "
//				<< G4endl;
//		//if((*pCItr).second > m_digitIn.thl)
		
		int x=(*pCItr).first.first;
		int y=(*pCItr).first.second;
		
		
		
		if(doAnimation){
			if((*pCItr).second>100*elec) anim->SubThresholdPixel(x-TMath::FloorNint(sx/pitchX),y-TMath::FloorNint(sy/pitchY));
		};
		// if preamp pulse rise over trigger
		
		//G4cout << TString::Format("Energy : %f , Threshold : %f",((*pCItr).second/keV),ThresholdMatrix[x][y]/keV) << endl;
		
		
		if(((*pCItr).second) > ThresholdMatrix[x][y])
		{
			// create one digit per pixel, I need to look at all the pixels first
			AllPixTimepixDigit * digit = new AllPixTimepixDigit;
			digit->SetPixelIDX(x);
			digit->SetPixelIDY(y);
			//digit->SetPixelCounts(EnergyToTOTSurogate((*pCItr).second,x,y));
			digit->SetPixelCounts((*pCItr).second/eV);
			//digit->IncreasePixelCounts(); // Counting mode
			
			
			if(doAnimation)anim->FirePixel(x-TMath::FloorNint(sx/pitchX),y-TMath::FloorNint(sy/pitchY));
			
			//G4cout << "TOT= "<< EnergyToTOTSurogate((*pCItr).second,x,y) << endl;
			//G4cout << "Energy= "<< 1e6*(*pCItr).second/keV << endl;

			// MC only //
			// Replicating the same information in all pixels
			// FIXME !
			digit->SetPrimaryVertex(m_primaryVertex->GetPosition());
			digit->SetPixelEnergyDep((*pCItr).second/eV);

			// Finally insert the digit in the digit collection
			//  Again, one per pixel.
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
	
	if(doAnimation)delete anim;

}





vector<G4double>  AllPixTimepixDigitizer::RKF5IntegrationHoles(G4double x, G4double y, G4double z,G4double dt)
{
// This function transport using Euler integration, for field (Ex,Ey,Ez),
// considered constant over time dt. The movement equation are those
// of charges in semi-conductors, sx= mu*E*dt;;
	double k1x,k2x,k3x,k4x,k5x,k6x;
	double k1y,k2y,k3y,k4y,k5y,k6y;
	double k1z,k2z,k3z,k4z,k5z,k6z;
	double dx,dy,dz;

      ComputeElectricField(x,y,z);

      k1x=MobilityHole(x,y,z)*electricFieldX*dt;
      k1y=MobilityHole(x,y,z)*electricFieldY*dt;
      k1z=MobilityHole(x,y,z)*electricFieldZ*dt;

      ComputeElectricField(x+k1x/4,y+k1y/4,z+k1z/4);

      k2x=MobilityHole(x+k1x/4,y+k1y/4,z+k1z/4)*electricFieldX*dt;
      k2y=MobilityHole(x+k1x/4,y+k1y/4,z+k1z/4)*electricFieldY*dt;
      k2z=MobilityHole(x+k1x/4,y+k1y/4,z+k1z/4)*electricFieldZ*dt;

     ComputeElectricField(x+(9./32)*k2x+(3./32)*k1x,y+(9./32)*k2y+(3./32)*k1y,z+(9./32)*k2z+(3./32)*k1z);

      k3x=MobilityHole(x+(9./32)*k2x+(3./32)*k1x,y+(9./32)*k2y+(3./32)*k1y,z+(9./32)*k2z+(3./32)*k1z)*electricFieldX*dt;
      k3y=MobilityHole(x+(9./32)*k2x+(3./32)*k1x,y+(9./32)*k2y+(3./32)*k1y,z+(9./32)*k2z+(3./32)*k1z)*electricFieldY*dt;
      k3z=MobilityHole(x+(9./32)*k2x+(3./32)*k1x,y+(9./32)*k2y+(3./32)*k1y,z+(9./32)*k2z+(3./32)*k1z)*electricFieldZ*dt;

      ComputeElectricField(x-(7200./2197)*k2x+(1932./2197)*k1x+(7296./2197)*k3x,y-(7200./2197)*k2y+(1932./2197)*k1y+(7296./2197)*k3y,z-(7200./2197)*k2z+(1932./2197)*k1z+(7296./2197)*k3z);

      k4x=MobilityHole(x-(7200./2197)*k2x+(1932./2197)*k1x+(7296./2197)*k3x,y-(7200./2197)*k2y+(1932./2197)*k1y+(7296./2197)*k3y,z-(7200./2197)*k2z+(1932./2197)*k1z+(7296./2197)*k3z)*electricFieldX*dt;
      k4y=MobilityHole(x-(7200./2197)*k2x+(1932./2197)*k1x+(7296./2197)*k3x,y-(7200./2197)*k2y+(1932./2197)*k1y+(7296./2197)*k3y,z-(7200./2197)*k2z+(1932./2197)*k1z+(7296./2197)*k3z)*electricFieldY*dt;
      k4z=MobilityHole(x-(7200./2197)*k2x+(1932./2197)*k1x+(7296./2197)*k3x,y-(7200./2197)*k2y+(1932./2197)*k1y+(7296./2197)*k3y,z-(7200./2197)*k2z+(1932./2197)*k1z+(7296./2197)*k3z)*electricFieldZ*dt;

      ComputeElectricField(x-(8)*k2x+(439./216)*k1x+(3680./513)*k3x-(845./4104)*k4x,y-(8)*k2y+(439./216)*k1y+(3680./513)*k3y-(845./4104)*k4y,z-(8)*k2z+(439./216)*k1z+(3680./513)*k3z-(845./4104)*k4z);

      k5x=MobilityHole(x-(8)*k2x+(439./216)*k1x+(3680./513)*k3x-(845./4104)*k4x,y-(8)*k2y+(439./216)*k1y+(3680./513)*k3y-(845./4104)*k4y,z-(8)*k2z+(439./216)*k1z+(3680./513)*k3z-(845./4104)*k4z)*electricFieldX*dt;
      k5y=MobilityHole(x-(8)*k2x+(439./216)*k1x+(3680./513)*k3x-(845./4104)*k4x,y-(8)*k2y+(439./216)*k1y+(3680./513)*k3y-(845./4104)*k4y,z-(8)*k2z+(439./216)*k1z+(3680./513)*k3z-(845./4104)*k4z)*electricFieldY*dt;
      k5z=MobilityHole(x-(8)*k2x+(439./216)*k1x+(3680./513)*k3x-(845./4104)*k4x,y-(8)*k2y+(439./216)*k1y+(3680./513)*k3y-(845./4104)*k4y,z-(8)*k2z+(439./216)*k1z+(3680./513)*k3z-(845./4104)*k4z)*electricFieldZ*dt;

      ComputeElectricField(x+(2)*k2x-(8./27)*k1x-(3544./2565)*k3x-(1859./4104)*k4x-(11./40)*k5x,
      y+(2)*k2y-(8./27)*k1y-(3544./2565)*k3y-(1859./4104)*k4y-(11./40)*k5y,
      z+(2)*k2z-(8./27)*k1z-(3544./2565)*k3z-(1859./4104)*k4z-(11./40)*k5z);

      k6x=MobilityHole(x+(2)*k2x-(8./27)*k1x-(3544./2565)*k3x-(1859./4104)*k4x-(11./40)*k5x,
    	      y+(2)*k2y-(8./27)*k1y-(3544./2565)*k3y-(1859./4104)*k4y-(11./40)*k5y,
    	      z+(2)*k2z-(8./27)*k1z-(3544./2565)*k3z-(1859./4104)*k4z-(11./40)*k5z)*electricFieldX*dt;
      k6y=MobilityHole(x+(2)*k2x-(8./27)*k1x-(3544./2565)*k3x-(1859./4104)*k4x-(11./40)*k5x,
    	      y+(2)*k2y-(8./27)*k1y-(3544./2565)*k3y-(1859./4104)*k4y-(11./40)*k5y,
    	      z+(2)*k2z-(8./27)*k1z-(3544./2565)*k3z-(1859./4104)*k4z-(11./40)*k5z)*electricFieldY*dt;
      k6z=MobilityHole(x+(2)*k2x-(8./27)*k1x-(3544./2565)*k3x-(1859./4104)*k4x-(11./40)*k5x,
    	      y+(2)*k2y-(8./27)*k1y-(3544./2565)*k3y-(1859./4104)*k4y-(11./40)*k5y,
    	      z+(2)*k2z-(8./27)*k1z-(3544./2565)*k3z-(1859./4104)*k4z-(11./40)*k5z)*electricFieldZ*dt;

      dx=((16./135)*k1x+(6656./12825)*k3x+(28561./56430)*k4x-(9./50)*k5x+(2./55)*k6x);
      dy=((16./135)*k1y+(6656./12825)*k3y+(28561./56430)*k4y-(9./50)*k5y+(2./55)*k6y);
      dz=((16./135)*k1z+(6656./12825)*k3z+(28561./56430)*k4z-(9./50)*k5z+(2./55)*k6z);

      double Ex,Ey,Ez,Erreur;
      Ex=((1./360)*k1x-(128./4275)*k3x-(2197./75240)*k4x-(1./50)*k5x+(2./55)*k6x);
      Ey=((1./360)*k1y-(128./4275)*k3y-(2197./75240)*k4y-(1./50)*k5y+(2./55)*k6y);
      Ez=((1./360)*k1z-(128./4275)*k3z-(2197./75240)*k4z-(1./50)*k5z+(2./55)*k6z);
      Erreur=sqrt(Ex*Ex+Ey*Ey+Ez*Ez);

      vector<G4double> newpoint(4);
      newpoint[0]=(dx);
      newpoint[1]=(dy);
      newpoint[2]=(dz);
      newpoint[3]=(Erreur);

     return newpoint;

}

vector<G4double>  AllPixTimepixDigitizer::RKF5IntegrationBFieldHoles(G4double x, G4double y, G4double z,G4double dt)
{
// This function transport using Euler integration, for field (Ex,Ey,Ez),
// considered constant over time dt. The movement equation are those
// of charges in semi-conductors, sx= mu*E*dt;;
	double k1x,k2x,k3x,k4x,k5x,k6x;
	double k1y,k2y,k3y,k4y,k5y,k6y;
	double k1z,k2z,k3z,k4z,k5z,k6z;
	double dx,dy,dz;

      ComputeEffectiveElectricField(x,y,z);

      k1x=r_H_h*MobilityHole(x,y,z)*electricFieldX*dt;
      k1y=r_H_h*MobilityHole(x,y,z)*electricFieldY*dt;
      k1z=r_H_h*MobilityHole(x,y,z)*electricFieldZ*dt;

      ComputeEffectiveElectricField(x+k1x/4,y+k1y/4,z+k1z/4);

      k2x=r_H_h*MobilityHole(x+k1x/4,y+k1y/4,z+k1z/4)*electricFieldX*dt;
      k2y=r_H_h*MobilityHole(x+k1x/4,y+k1y/4,z+k1z/4)*electricFieldY*dt;
      k2z=r_H_h*MobilityHole(x+k1x/4,y+k1y/4,z+k1z/4)*electricFieldZ*dt;

     ComputeEffectiveElectricField(x+(9./32)*k2x+(3./32)*k1x,y+(9./32)*k2y+(3./32)*k1y,z+(9./32)*k2z+(3./32)*k1z);

      k3x=r_H_h*MobilityHole(x+(9./32)*k2x+(3./32)*k1x,y+(9./32)*k2y+(3./32)*k1y,z+(9./32)*k2z+(3./32)*k1z)*electricFieldX*dt;
      k3y=r_H_h*MobilityHole(x+(9./32)*k2x+(3./32)*k1x,y+(9./32)*k2y+(3./32)*k1y,z+(9./32)*k2z+(3./32)*k1z)*electricFieldY*dt;
      k3z=r_H_h*MobilityHole(x+(9./32)*k2x+(3./32)*k1x,y+(9./32)*k2y+(3./32)*k1y,z+(9./32)*k2z+(3./32)*k1z)*electricFieldZ*dt;

      ComputeEffectiveElectricField(x-(7200./2197)*k2x+(1932./2197)*k1x+(7296./2197)*k3x,y-(7200./2197)*k2y+(1932./2197)*k1y+(7296./2197)*k3y,z-(7200./2197)*k2z+(1932./2197)*k1z+(7296./2197)*k3z);

      k4x=r_H_h*MobilityHole(x-(7200./2197)*k2x+(1932./2197)*k1x+(7296./2197)*k3x,y-(7200./2197)*k2y+(1932./2197)*k1y+(7296./2197)*k3y,z-(7200./2197)*k2z+(1932./2197)*k1z+(7296./2197)*k3z)*electricFieldX*dt;
      k4y=r_H_h*MobilityHole(x-(7200./2197)*k2x+(1932./2197)*k1x+(7296./2197)*k3x,y-(7200./2197)*k2y+(1932./2197)*k1y+(7296./2197)*k3y,z-(7200./2197)*k2z+(1932./2197)*k1z+(7296./2197)*k3z)*electricFieldY*dt;
      k4z=r_H_h*MobilityHole(x-(7200./2197)*k2x+(1932./2197)*k1x+(7296./2197)*k3x,y-(7200./2197)*k2y+(1932./2197)*k1y+(7296./2197)*k3y,z-(7200./2197)*k2z+(1932./2197)*k1z+(7296./2197)*k3z)*electricFieldZ*dt;

      ComputeEffectiveElectricField(x-(8)*k2x+(439./216)*k1x+(3680./513)*k3x-(845./4104)*k4x,y-(8)*k2y+(439./216)*k1y+(3680./513)*k3y-(845./4104)*k4y,z-(8)*k2z+(439./216)*k1z+(3680./513)*k3z-(845./4104)*k4z);

      k5x=r_H_h*MobilityHole(x-(8)*k2x+(439./216)*k1x+(3680./513)*k3x-(845./4104)*k4x,y-(8)*k2y+(439./216)*k1y+(3680./513)*k3y-(845./4104)*k4y,z-(8)*k2z+(439./216)*k1z+(3680./513)*k3z-(845./4104)*k4z)*electricFieldX*dt;
      k5y=r_H_h*MobilityHole(x-(8)*k2x+(439./216)*k1x+(3680./513)*k3x-(845./4104)*k4x,y-(8)*k2y+(439./216)*k1y+(3680./513)*k3y-(845./4104)*k4y,z-(8)*k2z+(439./216)*k1z+(3680./513)*k3z-(845./4104)*k4z)*electricFieldY*dt;
      k5z=r_H_h*MobilityHole(x-(8)*k2x+(439./216)*k1x+(3680./513)*k3x-(845./4104)*k4x,y-(8)*k2y+(439./216)*k1y+(3680./513)*k3y-(845./4104)*k4y,z-(8)*k2z+(439./216)*k1z+(3680./513)*k3z-(845./4104)*k4z)*electricFieldZ*dt;

      ComputeEffectiveElectricField(x+(2)*k2x-(8./27)*k1x-(3544./2565)*k3x-(1859./4104)*k4x-(11./40)*k5x,
      y+(2)*k2y-(8./27)*k1y-(3544./2565)*k3y-(1859./4104)*k4y-(11./40)*k5y,
      z+(2)*k2z-(8./27)*k1z-(3544./2565)*k3z-(1859./4104)*k4z-(11./40)*k5z);

      k6x=r_H_h*MobilityHole(x+(2)*k2x-(8./27)*k1x-(3544./2565)*k3x-(1859./4104)*k4x-(11./40)*k5x,
    	      y+(2)*k2y-(8./27)*k1y-(3544./2565)*k3y-(1859./4104)*k4y-(11./40)*k5y,
    	      z+(2)*k2z-(8./27)*k1z-(3544./2565)*k3z-(1859./4104)*k4z-(11./40)*k5z)*electricFieldX*dt;
      k6y=r_H_h*MobilityHole(x+(2)*k2x-(8./27)*k1x-(3544./2565)*k3x-(1859./4104)*k4x-(11./40)*k5x,
    	      y+(2)*k2y-(8./27)*k1y-(3544./2565)*k3y-(1859./4104)*k4y-(11./40)*k5y,
    	      z+(2)*k2z-(8./27)*k1z-(3544./2565)*k3z-(1859./4104)*k4z-(11./40)*k5z)*electricFieldY*dt;
      k6z=r_H_h*MobilityHole(x+(2)*k2x-(8./27)*k1x-(3544./2565)*k3x-(1859./4104)*k4x-(11./40)*k5x,
    	      y+(2)*k2y-(8./27)*k1y-(3544./2565)*k3y-(1859./4104)*k4y-(11./40)*k5y,
    	      z+(2)*k2z-(8./27)*k1z-(3544./2565)*k3z-(1859./4104)*k4z-(11./40)*k5z)*electricFieldZ*dt;

      dx=((16./135)*k1x+(6656./12825)*k3x+(28561./56430)*k4x-(9./50)*k5x+(2./55)*k6x);
      dy=((16./135)*k1y+(6656./12825)*k3y+(28561./56430)*k4y-(9./50)*k5y+(2./55)*k6y);
      dz=((16./135)*k1z+(6656./12825)*k3z+(28561./56430)*k4z-(9./50)*k5z+(2./55)*k6z);

      double Ex,Ey,Ez,Erreur;
      Ex=((1./360)*k1x-(128./4275)*k3x-(2197./75240)*k4x-(1./50)*k5x+(2./55)*k6x);
      Ey=((1./360)*k1y-(128./4275)*k3y-(2197./75240)*k4y-(1./50)*k5y+(2./55)*k6y);
      Ez=((1./360)*k1z-(128./4275)*k3z-(2197./75240)*k4z-(1./50)*k5z+(2./55)*k6z);
      Erreur=sqrt(Ex*Ex+Ey*Ey+Ez*Ez);

      vector<G4double> newpoint(4);
      newpoint[0]=(dx);
      newpoint[1]=(dy);
      newpoint[2]=(dz);
      newpoint[3]=(Erreur);

     return newpoint;

}

vector<G4double>  AllPixTimepixDigitizer::RKF5IntegrationElectrons(G4double x, G4double y, G4double z,G4double dt)
{
// This function transport using Euler integration, for field (Ex,Ey,Ez),
// considered constant over time dt. The movement equation are those
// of charges in semi-conductors, sx= mu*E*dt;;
	double k1x,k2x,k3x,k4x,k5x,k6x;
	double k1y,k2y,k3y,k4y,k5y,k6y;
	double k1z,k2z,k3z,k4z,k5z,k6z;
	double dx,dy,dz;

      ComputeElectricField(x,y,z);

      k1x=-MobilityElectron(x,y,z)*electricFieldX*dt;
      k1y=-MobilityElectron(x,y,z)*electricFieldY*dt;
      k1z=-MobilityElectron(x,y,z)*electricFieldZ*dt;

      ComputeElectricField(x+k1x/4,y+k1y/4,z+k1z/4);

      k2x=-MobilityElectron(x+k1x/4,y+k1y/4,z+k1z/4)*electricFieldX*dt;
      k2y=-MobilityElectron(x+k1x/4,y+k1y/4,z+k1z/4)*electricFieldY*dt;
      k2z=-MobilityElectron(x+k1x/4,y+k1y/4,z+k1z/4)*electricFieldZ*dt;

     ComputeElectricField(x+(9./32)*k2x+(3./32)*k1x,y+(9./32)*k2y+(3./32)*k1y,z+(9./32)*k2z+(3./32)*k1z);

      k3x=-MobilityElectron(x+(9./32)*k2x+(3./32)*k1x,y+(9./32)*k2y+(3./32)*k1y,z+(9./32)*k2z+(3./32)*k1z)*electricFieldX*dt;
      k3y=-MobilityElectron(x+(9./32)*k2x+(3./32)*k1x,y+(9./32)*k2y+(3./32)*k1y,z+(9./32)*k2z+(3./32)*k1z)*electricFieldY*dt;
      k3z=-MobilityElectron(x+(9./32)*k2x+(3./32)*k1x,y+(9./32)*k2y+(3./32)*k1y,z+(9./32)*k2z+(3./32)*k1z)*electricFieldZ*dt;

      ComputeElectricField(x-(7200./2197)*k2x+(1932./2197)*k1x+(7296./2197)*k3x,y-(7200./2197)*k2y+(1932./2197)*k1y+(7296./2197)*k3y,z-(7200./2197)*k2z+(1932./2197)*k1z+(7296./2197)*k3z);

      k4x=-MobilityElectron(x-(7200./2197)*k2x+(1932./2197)*k1x+(7296./2197)*k3x,y-(7200./2197)*k2y+(1932./2197)*k1y+(7296./2197)*k3y,z-(7200./2197)*k2z+(1932./2197)*k1z+(7296./2197)*k3z)*electricFieldX*dt;
      k4y=-MobilityElectron(x-(7200./2197)*k2x+(1932./2197)*k1x+(7296./2197)*k3x,y-(7200./2197)*k2y+(1932./2197)*k1y+(7296./2197)*k3y,z-(7200./2197)*k2z+(1932./2197)*k1z+(7296./2197)*k3z)*electricFieldY*dt;
      k4z=-MobilityElectron(x-(7200./2197)*k2x+(1932./2197)*k1x+(7296./2197)*k3x,y-(7200./2197)*k2y+(1932./2197)*k1y+(7296./2197)*k3y,z-(7200./2197)*k2z+(1932./2197)*k1z+(7296./2197)*k3z)*electricFieldZ*dt;

      ComputeElectricField(x-(8)*k2x+(439./216)*k1x+(3680./513)*k3x-(845./4104)*k4x,y-(8)*k2y+(439./216)*k1y+(3680./513)*k3y-(845./4104)*k4y,z-(8)*k2z+(439./216)*k1z+(3680./513)*k3z-(845./4104)*k4z);

      k5x=-MobilityElectron(x-(8)*k2x+(439./216)*k1x+(3680./513)*k3x-(845./4104)*k4x,y-(8)*k2y+(439./216)*k1y+(3680./513)*k3y-(845./4104)*k4y,z-(8)*k2z+(439./216)*k1z+(3680./513)*k3z-(845./4104)*k4z)*electricFieldX*dt;
      k5y=-MobilityElectron(x-(8)*k2x+(439./216)*k1x+(3680./513)*k3x-(845./4104)*k4x,y-(8)*k2y+(439./216)*k1y+(3680./513)*k3y-(845./4104)*k4y,z-(8)*k2z+(439./216)*k1z+(3680./513)*k3z-(845./4104)*k4z)*electricFieldY*dt;
      k5z=-MobilityElectron(x-(8)*k2x+(439./216)*k1x+(3680./513)*k3x-(845./4104)*k4x,y-(8)*k2y+(439./216)*k1y+(3680./513)*k3y-(845./4104)*k4y,z-(8)*k2z+(439./216)*k1z+(3680./513)*k3z-(845./4104)*k4z)*electricFieldZ*dt;

      ComputeElectricField(x+(2)*k2x-(8./27)*k1x-(3544./2565)*k3x-(1859./4104)*k4x-(11./40)*k5x,
      y+(2)*k2y-(8./27)*k1y-(3544./2565)*k3y-(1859./4104)*k4y-(11./40)*k5y,
      z+(2)*k2z-(8./27)*k1z-(3544./2565)*k3z-(1859./4104)*k4z-(11./40)*k5z);

      k6x=-MobilityElectron(x+(2)*k2x-(8./27)*k1x-(3544./2565)*k3x-(1859./4104)*k4x-(11./40)*k5x,
    	      y+(2)*k2y-(8./27)*k1y-(3544./2565)*k3y-(1859./4104)*k4y-(11./40)*k5y,
    	      z+(2)*k2z-(8./27)*k1z-(3544./2565)*k3z-(1859./4104)*k4z-(11./40)*k5z)*electricFieldX*dt;
      k6y=-MobilityElectron(x+(2)*k2x-(8./27)*k1x-(3544./2565)*k3x-(1859./4104)*k4x-(11./40)*k5x,
    	      y+(2)*k2y-(8./27)*k1y-(3544./2565)*k3y-(1859./4104)*k4y-(11./40)*k5y,
    	      z+(2)*k2z-(8./27)*k1z-(3544./2565)*k3z-(1859./4104)*k4z-(11./40)*k5z)*electricFieldY*dt;
      k6z=-MobilityElectron(x+(2)*k2x-(8./27)*k1x-(3544./2565)*k3x-(1859./4104)*k4x-(11./40)*k5x,
    	      y+(2)*k2y-(8./27)*k1y-(3544./2565)*k3y-(1859./4104)*k4y-(11./40)*k5y,
    	      z+(2)*k2z-(8./27)*k1z-(3544./2565)*k3z-(1859./4104)*k4z-(11./40)*k5z)*electricFieldZ*dt;

      dx=((16./135)*k1x+(6656./12825)*k3x+(28561./56430)*k4x-(9./50)*k5x+(2./55)*k6x);
      dy=((16./135)*k1y+(6656./12825)*k3y+(28561./56430)*k4y-(9./50)*k5y+(2./55)*k6y);
      dz=((16./135)*k1z+(6656./12825)*k3z+(28561./56430)*k4z-(9./50)*k5z+(2./55)*k6z);

      double Ex,Ey,Ez,Erreur;
      Ex=((1./360)*k1x-(128./4275)*k3x-(2197./75240)*k4x-(1./50)*k5x+(2./55)*k6x);
      Ey=((1./360)*k1y-(128./4275)*k3y-(2197./75240)*k4y-(1./50)*k5y+(2./55)*k6y);
      Ez=((1./360)*k1z-(128./4275)*k3z-(2197./75240)*k4z-(1./50)*k5z+(2./55)*k6z);
      Erreur=sqrt(Ex*Ex+Ey*Ey+Ez*Ez);

      vector<G4double> newpoint(4);
      newpoint[0]=(dx);
      newpoint[1]=(dy);
      newpoint[2]=(dz);
      newpoint[3]=(Erreur);

     return newpoint;

}

vector<G4double>  AllPixTimepixDigitizer::RKF5IntegrationBFieldElectrons(G4double x, G4double y, G4double z,G4double dt)
{
// This function transport using Euler integration, for field (Ex,Ey,Ez),
// considered constant over time dt. The movement equation are those
// of charges in semi-conductors, sx= mu*E*dt;;
	double k1x,k2x,k3x,k4x,k5x,k6x;
	double k1y,k2y,k3y,k4y,k5y,k6y;
	double k1z,k2z,k3z,k4z,k5z,k6z;
	double dx,dy,dz;

      ComputeEffectiveElectricField(x,y,z);

      k1x=-r_H_e*MobilityElectron(x,y,z)*electricFieldX*dt;
      k1y=-r_H_e*MobilityElectron(x,y,z)*electricFieldY*dt;
      k1z=-r_H_e*MobilityElectron(x,y,z)*electricFieldZ*dt;

      ComputeEffectiveElectricField(x+k1x/4,y+k1y/4,z+k1z/4);

      k2x=-r_H_e*MobilityElectron(x+k1x/4,y+k1y/4,z+k1z/4)*electricFieldX*dt;
      k2y=-r_H_e*MobilityElectron(x+k1x/4,y+k1y/4,z+k1z/4)*electricFieldY*dt;
      k2z=-r_H_e*MobilityElectron(x+k1x/4,y+k1y/4,z+k1z/4)*electricFieldZ*dt;

     ComputeEffectiveElectricField(x+(9./32)*k2x+(3./32)*k1x,y+(9./32)*k2y+(3./32)*k1y,z+(9./32)*k2z+(3./32)*k1z);

      k3x=-r_H_e*MobilityElectron(x+(9./32)*k2x+(3./32)*k1x,y+(9./32)*k2y+(3./32)*k1y,z+(9./32)*k2z+(3./32)*k1z)*electricFieldX*dt;
      k3y=-r_H_e*MobilityElectron(x+(9./32)*k2x+(3./32)*k1x,y+(9./32)*k2y+(3./32)*k1y,z+(9./32)*k2z+(3./32)*k1z)*electricFieldY*dt;
      k3z=-r_H_e*MobilityElectron(x+(9./32)*k2x+(3./32)*k1x,y+(9./32)*k2y+(3./32)*k1y,z+(9./32)*k2z+(3./32)*k1z)*electricFieldZ*dt;

      ComputeEffectiveElectricField(x-(7200./2197)*k2x+(1932./2197)*k1x+(7296./2197)*k3x,y-(7200./2197)*k2y+(1932./2197)*k1y+(7296./2197)*k3y,z-(7200./2197)*k2z+(1932./2197)*k1z+(7296./2197)*k3z);

      k4x=-r_H_e*MobilityElectron(x-(7200./2197)*k2x+(1932./2197)*k1x+(7296./2197)*k3x,y-(7200./2197)*k2y+(1932./2197)*k1y+(7296./2197)*k3y,z-(7200./2197)*k2z+(1932./2197)*k1z+(7296./2197)*k3z)*electricFieldX*dt;
      k4y=-r_H_e*MobilityElectron(x-(7200./2197)*k2x+(1932./2197)*k1x+(7296./2197)*k3x,y-(7200./2197)*k2y+(1932./2197)*k1y+(7296./2197)*k3y,z-(7200./2197)*k2z+(1932./2197)*k1z+(7296./2197)*k3z)*electricFieldY*dt;
      k4z=-r_H_e*MobilityElectron(x-(7200./2197)*k2x+(1932./2197)*k1x+(7296./2197)*k3x,y-(7200./2197)*k2y+(1932./2197)*k1y+(7296./2197)*k3y,z-(7200./2197)*k2z+(1932./2197)*k1z+(7296./2197)*k3z)*electricFieldZ*dt;

      ComputeEffectiveElectricField(x-(8)*k2x+(439./216)*k1x+(3680./513)*k3x-(845./4104)*k4x,y-(8)*k2y+(439./216)*k1y+(3680./513)*k3y-(845./4104)*k4y,z-(8)*k2z+(439./216)*k1z+(3680./513)*k3z-(845./4104)*k4z);

      k5x=-r_H_e*MobilityElectron(x-(8)*k2x+(439./216)*k1x+(3680./513)*k3x-(845./4104)*k4x,y-(8)*k2y+(439./216)*k1y+(3680./513)*k3y-(845./4104)*k4y,z-(8)*k2z+(439./216)*k1z+(3680./513)*k3z-(845./4104)*k4z)*electricFieldX*dt;
      k5y=-r_H_e*MobilityElectron(x-(8)*k2x+(439./216)*k1x+(3680./513)*k3x-(845./4104)*k4x,y-(8)*k2y+(439./216)*k1y+(3680./513)*k3y-(845./4104)*k4y,z-(8)*k2z+(439./216)*k1z+(3680./513)*k3z-(845./4104)*k4z)*electricFieldY*dt;
      k5z=-r_H_e*MobilityElectron(x-(8)*k2x+(439./216)*k1x+(3680./513)*k3x-(845./4104)*k4x,y-(8)*k2y+(439./216)*k1y+(3680./513)*k3y-(845./4104)*k4y,z-(8)*k2z+(439./216)*k1z+(3680./513)*k3z-(845./4104)*k4z)*electricFieldZ*dt;

      ComputeEffectiveElectricField(x+(2)*k2x-(8./27)*k1x-(3544./2565)*k3x-(1859./4104)*k4x-(11./40)*k5x,
      y+(2)*k2y-(8./27)*k1y-(3544./2565)*k3y-(1859./4104)*k4y-(11./40)*k5y,
      z+(2)*k2z-(8./27)*k1z-(3544./2565)*k3z-(1859./4104)*k4z-(11./40)*k5z);

      k6x=-r_H_e*MobilityElectron(x+(2)*k2x-(8./27)*k1x-(3544./2565)*k3x-(1859./4104)*k4x-(11./40)*k5x,
    	      y+(2)*k2y-(8./27)*k1y-(3544./2565)*k3y-(1859./4104)*k4y-(11./40)*k5y,
    	      z+(2)*k2z-(8./27)*k1z-(3544./2565)*k3z-(1859./4104)*k4z-(11./40)*k5z)*electricFieldX*dt;
      k6y=-r_H_e*MobilityElectron(x+(2)*k2x-(8./27)*k1x-(3544./2565)*k3x-(1859./4104)*k4x-(11./40)*k5x,
    	      y+(2)*k2y-(8./27)*k1y-(3544./2565)*k3y-(1859./4104)*k4y-(11./40)*k5y,
    	      z+(2)*k2z-(8./27)*k1z-(3544./2565)*k3z-(1859./4104)*k4z-(11./40)*k5z)*electricFieldY*dt;
      k6z=-r_H_e*MobilityElectron(x+(2)*k2x-(8./27)*k1x-(3544./2565)*k3x-(1859./4104)*k4x-(11./40)*k5x,
    	      y+(2)*k2y-(8./27)*k1y-(3544./2565)*k3y-(1859./4104)*k4y-(11./40)*k5y,
    	      z+(2)*k2z-(8./27)*k1z-(3544./2565)*k3z-(1859./4104)*k4z-(11./40)*k5z)*electricFieldZ*dt;

      dx=((16./135)*k1x+(6656./12825)*k3x+(28561./56430)*k4x-(9./50)*k5x+(2./55)*k6x);
      dy=((16./135)*k1y+(6656./12825)*k3y+(28561./56430)*k4y-(9./50)*k5y+(2./55)*k6y);
      dz=((16./135)*k1z+(6656./12825)*k3z+(28561./56430)*k4z-(9./50)*k5z+(2./55)*k6z);

      double Ex,Ey,Ez,Erreur;
      Ex=((1./360)*k1x-(128./4275)*k3x-(2197./75240)*k4x-(1./50)*k5x+(2./55)*k6x);
      Ey=((1./360)*k1y-(128./4275)*k3y-(2197./75240)*k4y-(1./50)*k5y+(2./55)*k6y);
      Ez=((1./360)*k1z-(128./4275)*k3z-(2197./75240)*k4z-(1./50)*k5z+(2./55)*k6z);
      Erreur=sqrt(Ex*Ex+Ey*Ey+Ez*Ez);

      vector<G4double> newpoint(4);
      newpoint[0]=(dx);
      newpoint[1]=(dy);
      newpoint[2]=(dz);
      newpoint[3]=(Erreur);

     return newpoint;

}


//______________________________________________________________________________
G4double  AllPixTimepixDigitizer::SetDt(G4double dt, G4double ErreurMoy)
{

	double Dt=dt;
	//if(isnan(ErreurMoy)){Dt=tup;}
	if( ErreurMoy != ErreurMoy ){Dt=tup;}
	else if(ErreurMoy > Target){ Dt*=0.9;}
	else if(ErreurMoy < Target){ Dt*=1.1;};


	if(Dt<tlow) Dt=tlow;
	if(Dt>tup)  Dt=tup;

	return Dt;

}

G4int AllPixTimepixDigitizer::EnergyToTOT(G4double Energy,G4int x,G4int y){
	double Max=0;
	
	//G4cout << TString::Format("Th= %e RS=%e FS=%e ",ThresholdMatrix[x][y],RisingSlope[x][y],FallingSlope[x][y]) << endl;
	
	
	if(Energy>SaturationEnergy)  Max= (1.0/keV)*SaturationEnergy*Gain[x][y];
	else	 Max= (1.0/keV)*Energy*Gain[x][y];
	
	double offset = CLHEP::RandFlat::shoot(0.,1.)*ClockUnit;
	
	double t0 = ((1.0/keV)*ThresholdMatrix[x][y]*Gain[x][y])/RisingSlope[x][y] + offset;
	
	double tpeak = Max/RisingSlope[x][y] + offset;
	
	
	
	double tf = (Max-((1.0/keV)*ThresholdMatrix[x][y]*Gain[x][y]))/FallingSlope[x][y] + tpeak;
	
	
	double ToT= TMath::FloorNint((tf-t0)/ClockUnit);
	
	//G4cout << TString::Format("to=%e tpeak=%e tf=%e ",t0,tpeak,tf) << endl;
	if(ToT==0)ToT=1;
	if(Max <= ((1.0/keV)*ThresholdMatrix[x][y]*Gain[x][y])) ToT=1;
	
	return ToT;
}

G4int AllPixTimepixDigitizer::EnergyToTOTSurogate(G4double Energy,G4int x,G4int y){

		//double ToT = (SurrogateD[x][y]*SurrogateA[x][y] +Energy/keV -SurrogateB[x][y]+TMath::Sqrt((SurrogateB[x][y]+SurrogateD[x][y]*SurrogateA[x][y]-Energy/keV)*(SurrogateB[x][y]+SurrogateD[x][y]*SurrogateA[x][y]-Energy/keV)+4*SurrogateA[x][y]*SurrogateC[x][y]))/(2*SurrogateA[x][y]);

//	double zero1 = (1.0/(2*SurrogateA[x][y]))*(SurrogateB[x][y] - SurrogateA[x][y]*SurrogateD[x][y] +
//			TMath::Sqrt(SurrogateB[x][y]*SurrogateB[x][y] - 4*SurrogateA[x][y]*SurrogateC[x][y]+2*SurrogateA[x][y]*SurrogateB[x][y]*SurrogateD[x][y] + SurrogateA[x][y]*SurrogateD[x][y]*SurrogateA[x][y]*SurrogateD[x][y]));
//	double zero2 = (1.0/(2*SurrogateA[x][y]))*(-SurrogateB[x][y] + SurrogateA[x][y]*SurrogateD[x][y] +
//			TMath::Sqrt(SurrogateB[x][y]*SurrogateB[x][y] - 4*SurrogateA[x][y]*SurrogateC[x][y]+2*SurrogateA[x][y]*SurrogateB[x][y]*SurrogateD[x][y] + SurrogateA[x][y]*SurrogateD[x][y]*SurrogateA[x][y]*SurrogateD[x][y]));

//	G4cout << zero1 << " " << zero2 << endl;

	double ToT;
	//if(Energy/keV<TMath::Max(zero1,zero2)){ToT=0;}
	if(Energy/keV<SurrogateD[x][y])ToT=0;
	else {ToT = SurrogateA[x][y]*Energy/keV+SurrogateB[x][y]-SurrogateC[x][y]/(Energy/keV - SurrogateD[x][y]);}
	//G4cout << TString::Format("A= %f B=%f C=%f D=%f",SurrogateA[x][y],SurrogateB[x][y],SurrogateC[x][y],SurrogateD[x][y]) << endl;

	if(ToT<0)ToT=0;

	ToT=TMath::Floor(ToT*CALIBRATION_CLOCK_UNIT/ClockUnit);

	return ToT;


}


G4double AllPixTimepixDigitizer::ComputeLorentzAngle(G4double x, G4double y, G4double z){
	
	ComputeElectricField(x,y,z);
	G4double angle=0;
	if (readoutType==ELECTRON){
		 angle = TMath::ATan(r_H_e*(MobilityElectron(x,y,z)/(cm2/s))*B_Field*1e-4);
	}
	else { 
		 angle = TMath::ATan(r_H_h*(MobilityHole(x,y,z)/(cm2/s))*B_Field*1e-4);	
	};
			
	//cout << "[TimepixDigi] Lorentz Angle = " << angle*TMath::RadToDeg() << endl;
	return angle;
}

//
//
//G4int AllPixTimepixDigitizer::EnergyToTOT(G4double Energy, G4double threshold)
//{
//
//	G4double dEdt = (MipCharge- threshold)/(MipTOT*ClockUnit);
//	//G4int TOT = TMath::FloorNint((Energy-threshold)/(dEdt*Lv1Unit*ns));
//	
//	G4int TOT;
//	
//	if(Energy>SaturationEnergy/elec){
//	TOT = TMath::FloorNint((SaturationEnergy/elec-threshold)/(dEdt*ClockUnit));
//	}
//	else{
//	TOT = TMath::FloorNint((Energy-threshold)/(dEdt*ClockUnit));
//	}
//	if(Energy<threshold) TOT=0;
//	if (TOT<0) TOT=0;
//	if (TOT>CounterDepth) TOT=CounterDepth;
//
//	return TOT;
//}


//G4int AllPixTimepixDigitizer::EnergyToTOTFEI4(G4double Energy, G4double threshold)
//{
//	double ToT=(1.0/(Energy+calib_C)+ calib_A)*calib_B;
//
//	return TMath::Floor(ToT);
//}





