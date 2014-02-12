/**
 *  Authors:
 *  Mathieu Benoit <benoit@lal.in2p3.fr>
 *  John Idarraga <idarraga@cern.ch>
 */

#include "AllPixFEI3StandardDigitizer.hh"
#include "AllPixTrackerHit.hh"
#include "G4DigiManager.hh"
#include "G4VDigitizerModule.hh"
#include "TMath.h"
#include "math.h"
#include "TF1.h"
#include "CLHEP/Random/RandGauss.h"
#include "CLHEP/Random/RandFlat.h"

//#include "SystemOfUnits.h"


#include "AllPixGeoDsc.hh"



using namespace TMath;

AllPixFEI3StandardDigitizer::AllPixFEI3StandardDigitizer(G4String modName, G4String hitsColName, G4String digitColName)
: AllPixDigitizerInterface (modName) {

	// Registration of digits collection name
	collectionName.push_back(digitColName);
	m_hitsColName.push_back(hitsColName);

	//// Unit for charge in FEIX average e/h pair creation energy in Silicon
	elec = 3.64*eV;


	// Example of detector description handle
	// provided by the interface
	AllPixGeoDsc * gD = GetDetectorGeoDscPtr();

	//////////////////////////
	// Bias and Temperature //
	//////////////////////////

	biasVoltage=100.0; //[V]
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

	//cout << TString::Format("Lx=%f Ly=%f npixX=%d npixY=%d",pitchX,pitchY,nPixX,nPixY) << endl;



	///////////////////////////////////////////////////
	// Silicon electron and hole transport constants //
	///////////////////////////////////////////////////

	// Default mobilities
	Default_Electron_Mobility=1562.0*cm2/s; // Electron mobility (cm2/Vs)
	Default_Hole_Mobility=450.0*cm2/s;// Hole mobility (cm2/Vs
	Default_Electron_D=40.43*cm2/s; // Electron mobility (cm2/s)
	Default_Hole_D=12*cm2/s;// Hole mobility (cm2/s



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

	G4cout << "!!!!!!!Radiation Damage Report !!!!!!!!" << endl
		 << TString::Format("depletionVoltage : %f depleted Depth : %f",depletionVoltage,depletedDepth/um) << endl
		 << TString::Format("Neff : %e Neff0 : %e ",Neff*cm3,Neff0*cm3);
	if(bulkType==false && Neff<0) cout << "Bulk inverted" << endl;






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
 	tup=	0.1*ns;
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

 	FEIX = 4;

 	switch (FEIX) {

 	case 3 :
 		MipTOT=60;
 		MipCharge=19200;
 		CounterDepth=255;
 		Lv1Unit = 25*ns;
 		chipNoise = 300;
 		m_digitIn.thl = 3500;
 		chargeSharingConstant = 0.0;
 		break;


 	case 4 :
 		MipTOT=10;
 		MipCharge=40000;
 		//MipCharge=20000;
 		CounterDepth=15;
 		Lv1Unit = 25*ns;
 		chipNoise = 125*elec;
 		m_digitIn.thl = 3200*elec;
 		chargeSharingConstant =0.0;
 		break;


 	case 666 :
 		MipTOT=5;
 		MipCharge=12000;
 		CounterDepth=32;
 		Lv1Unit = 25*ns;
 		chipNoise = 20*elec;
 		m_digitIn.thl = 800*elec;
 		chargeSharingConstant =0.01;
 		break;


 	default :
 		MipTOT=60;
 		MipCharge=22000;
 		CounterDepth=255;
 		Lv1Unit = 25*ns;
 		chipNoise = 300*elec;
 		m_digitIn.thl = 3500*elec;
 		chargeSharingConstant =0.01;

 	}


 	/////////////////////////////
 	// Sensor Selection //
 	/////////////////////////////
 	/*


 	 1=FEI3
 	 2,3= FEI4, SlimEdge or Conservative
 	 4=Omegapix2
 	 5= FEI3 Slim edge
 	 Default = FEI3

 	 */
 	Sensor = 1;

 	switch (Sensor) {

 	case 1 :
 	 	doSlimEdge=false;
 		GRShift =250*um;
 		break;
 	case 2 :
 	 	doSlimEdge=true;
 		GRShift =250*um;
 		break;
 	case 3 :
 	 	doSlimEdge=true;
 		GRShift =500*um;
 		break;
 	case 4 :
 	 	doSlimEdge=true;
 		GRShift =500*um;
 		break;
 	case 5 :
 	 	doSlimEdge=true;
 		GRShift =250*um;
 		break;

 	default :
 	 	doSlimEdge=false;
 		GRShift =0*um;


 		//fix me !!

 	}

		mobility = MobilityElectron(0,0,0);




}

AllPixFEI3StandardDigitizer::~AllPixFEI3StandardDigitizer(){


}

void AllPixFEI3StandardDigitizer::SetDetectorDigitInputs(G4double thl){

	// set digitization input values
	// thl
	m_digitIn.thl = thl; // <-- input !
}



G4double AllPixFEI3StandardDigitizer::MyErf(G4double x){
	TF1 * f1 = new TF1("expint","(2.0/TMath::Pi())*exp(-x**2)",-100.,100.);
	return f1->Integral(0,x);


}


vector<G4double>  AllPixFEI3StandardDigitizer::ComputeDriftTimeFullField(G4double x, G4double y, G4double z)
{
	G4double driftTime=0;
	G4double dt=dtIni;
	G4double xtemp=x;
	G4double ytemp=y;
	G4double ztemp=z;
	//cout << "z : " << z/um << endl;
	vector<G4double> step;
	int counter=0;
	while( ztemp>0 && ztemp<detectorThickness){

		driftTime+=dt;
		step=RKF5Integration(xtemp,ytemp,ztemp,dt);
		xtemp+=step[0];
		ytemp+=step[1];
		ztemp+=step[2];
		if(GetElectricFieldNorm(xtemp,ytemp,ztemp)==0)ztemp=detectorThickness;
		dt=SetDt(dt,step[3]);
		counter++;
		//cout << "!!!!!! Drifting youhou! !!!!! " << TString::Format("dt: %f z : %f drift : %f Ez : %f",dt/ns,ztemp/um,driftTime/ns,electricFieldZ*cm) << endl;
	}

	vector<G4double> output(4);
	output[0]=xtemp;
	output[1]=ytemp;
	output[2]=ztemp;
	output[3]=driftTime;

	return output;
}


void AllPixFEI3StandardDigitizer::ComputeElectricField(G4double /*x*/, G4double /*y*/, G4double z){

		Efield1D(z);
}


G4double AllPixFEI3StandardDigitizer::GetElectricFieldNorm(G4double /*x*/, G4double /*y*/, G4double /*z*/){

	return TMath::Sqrt(electricFieldX*electricFieldX +electricFieldY*electricFieldY +electricFieldZ*electricFieldZ );
}

G4double  AllPixFEI3StandardDigitizer::MobilityElectron(G4double /*x*/, G4double /*y*/, G4double /*z*/){

	//ComputeElectricField(x,y,z);
	G4double parElectricField = electricFieldZ;//GetElectricFieldNorm(x,y,z);
	G4double mobilite = Default_Electron_Mobility*
	TMath::Power((1.0/(1.+ TMath::Power((Default_Electron_Mobility*parElectricField)/
			Electron_Saturation_Velocity,Electron_Beta))),1.0/Electron_Beta);
	return mobilite;
	//cout << "!!!!!!!!!" << mobilite << " " << parElectricField*cm << endl;

}


G4double AllPixFEI3StandardDigitizer::MobilityHole(G4double x, G4double y, G4double z){

	ComputeElectricField(x,y,z);
	G4double parElectricField = GetElectricFieldNorm(x,y,z);
	G4double mobilite = Default_Hole_Mobility*
	TMath::Power((1.0/(1.+ TMath::Power((Default_Hole_Mobility*parElectricField)/
			Hole_Saturation_Velocity,Hole_Beta))),1.0/Hole_Beta);
	return mobilite;


}


G4double AllPixFEI3StandardDigitizer::ComputeDriftTimeUniformField(AllPixTrackerHit * /*a_hit*/)
{
	// drift in a uniform electric field from hit position to surface
	//G4cout << "Hit position with respect to pixel (um) "<< (a_hit->GetPosWithRespectToPixel().z() + detectorThickness/2.0)/um  <<endl;
	//return (a_hit->GetPosWithRespectToPixel().z() + detectorThickness/2.0 ) /(mobility*electricFieldY);
	double drift = ( CLHEP::RandGauss::shoot((depletedDepth/detectorThickness)*detectorThickness/2,(depletedDepth/detectorThickness)*detectorThickness/4)) /(mobility*electricFieldZ);
	//cout << "drift: " << mobility*s/cm2 << endl;
	return drift;
}

G4double AllPixFEI3StandardDigitizer::ComputeDiffusionRMS(G4double tDrift)
{
	return TMath::Sqrt(2.*Default_Electron_D*tDrift);

}

G4double AllPixFEI3StandardDigitizer::ApplyTrapping(G4double tDrift, G4double Energy)
{

	//cout << TString::Format("[Trapping!!!] Energy:%f after Trapping:%f",Energy, Energy*TMath::Exp(-tDrift/trappingTime)) << endl;
	return Energy*TMath::Exp(-tDrift/trappingTime);

}



G4double AllPixFEI3StandardDigitizer::IntegrateGaussian(G4double xhit,G4double yhit,G4double Sigma, G4double x1, G4double x2, G4double y1, G4double y2, G4double Energy )
{
	//G4cout << TString::Format("Integration borns x1=%f x2=%f y1=%f y2=%f",x1/um,x2/um,y1/um,y2/um) << endl;
	//G4cout << TString::Format("TMath::Erf test %f %f %f %f",TMath::Erf((x1 - xhit)/(Sqrt(2.)*Sigma)),TMath::Erf((x2 - xhit)/(TMath::Sqrt(2.)*Sigma)),TMath::Erf((y1 - yhit)/(TMath::Sqrt(2.)*Sigma)),TMath::Erf((y2 - yhit)/(TMath::Sqrt(2.)*Sigma)))<< endl;
	return (Energy*(-TMath::Erf((x1 - xhit)/(Sqrt(2.)*Sigma)) + TMath::Erf((x2 - xhit)/(TMath::Sqrt(2.)*Sigma)))
			*(-TMath::Erf((y1 - yhit)/(TMath::Sqrt(2.)*Sigma)) + TMath::Erf((y2 - yhit)/(TMath::Sqrt(2.0)*Sigma))))/4.;
}


void AllPixFEI3StandardDigitizer::Efield1D(G4double z){


	electricFieldZ=(detectorThickness-z)*biasVoltage/depletedDepth;
	if(z>depletedDepth){
		electricFieldZ=0;
		//cout << "bip" << endl;
	}

}




/**
 * Digitizer for FEI3Standard
 */
void AllPixFEI3StandardDigitizer::Digitize(){

	// Create the digits collection
	m_digitsCollection = new AllPixFEI3StandardDigitsCollection("AllPixFEI3StandardDigitizer", collectionName[0] );

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

	for(G4int itr  = 0 ; itr < nEntries ; itr++) {

		G4double eHitTotal = elec*CLHEP::RandGauss::shoot((*hitsCollection)[itr]->GetEdep()/elec,TMath::Sqrt((*hitsCollection)[itr]->GetEdep()/elec)*0.118);
		hitsETotal += eHitTotal;

		//under-depletion
		eHitTotal = eHitTotal*depletedDepth/detectorThickness;

		tempPixel.first  = (*hitsCollection)[itr]->GetPixelNbX();
		tempPixel.second = (*hitsCollection)[itr]->GetPixelNbY();

		//G4cout << "x : " << tempPixel.first << " ,  y : " << tempPixel.second << ", E = " << eHitTotal/keV << G4endl;


		G4double xpos = (*hitsCollection)[itr]->GetPosWithRespectToPixel().x();
		G4double ypos = (*hitsCollection)[itr]->GetPosWithRespectToPixel().y();
		// not used
		G4double zpos = (*hitsCollection)[itr]->GetPosWithRespectToPixel().z()+detectorThickness/2.;

		//cout << TString::Format("x,y,z : %5.5f %5.5f %5.5f",xpos/um,ypos/um,zpos/um)<<endl;


		for(G4int nQ  = 0 ; nQ < precision ; nQ++) {

		double eHit = double(eHitTotal)/precision;
		// Hit Info
		//Ugly Hack !!!
		//if(xpos>0)xpos=-pitchX/2+xpos;
		//else xpos=pitchX/2+xpos;

		//G4double zpos = CLHEP::RandGauss::shoot(zpos,10*um);



		if(doSlimEdge){
			if (isSlimEdge(tempPixel.first,tempPixel.second)) eHit = SlimEdgeEffect(tempPixel.first,xpos,eHit);
		}


		G4double driftTime;
		sigma = 0;
		if(doFullField==false){
			 driftTime = ComputeDriftTimeUniformField((*hitsCollection)[itr]);
			 sigma = ComputeDiffusionRMS(driftTime);
			 //cout << TString::Format("!!!!!!!!! vd/vdep : %f drift time : %f sigma : %f",depletedDepth/detectorThickness,driftTime,sigma) << endl;
		}
		else{

			// Until we get out position
			vector<G4double> data = ComputeDriftTimeFullField(xpos,ypos,zpos);
			driftTime = data[3];
			sigma = ComputeDiffusionRMS(driftTime);
			xpos=data[0];
			ypos=data[1];

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
		  for(int i=-1;i<=1;i++){
			for(int j=-1;j<=1;j++){
				extraPixel = tempPixel;
				extraPixel.first +=i;
				extraPixel.second+=j;
				if(extraPixel.first >= 0 && extraPixel.second>=0 && extraPixel.first < nPixX && extraPixel.second < nPixY)
					{
					//We compute contribution of the hit to each pixels

					double Etemp = IntegrateGaussian(xpos/nm,ypos/nm,sigma/nm,(-pitchX/2.0 + i*pitchX)/nm,(-pitchX/2.+(i+1)*pitchX)/nm,(-pitchY/2 + j*pitchY)/nm,(-pitchY/2 + (j+1)*pitchY)/nm, eHit );
					if(doTrapping==true) pixelsContent[extraPixel]+=ApplyTrapping(driftTime,Etemp);
					else pixelsContent[extraPixel] +=Etemp;

					//G4cout << TString::Format("[Digitizer] Pixel %i %i Energy=%f, Energy after Trapping=%f",extraPixel.first,extraPixel.second,pixelsContent[extraPixel]/elec,ApplyTrapping(driftTime,pixelsContent[extraPixel])/elec) << endl;
					//cout << TString::Format("Pixel %i %i, Energy collected = %f sigma=%f tdrift=%f",extraPixel.first,extraPixel.second,Etemp/keV,sigma/nm,driftTime) << endl;

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
		double threshold =CLHEP::RandGauss::shoot(m_digitIn.thl,chipNoise);

		//G4cout << "pixel : " << (*pCItr).first.first << " , " << (*pCItr).first.second
			//	<< ", E =  " << ((*pCItr).second)/keV << " keV | thl = "
			//	<< TString::Format("Threshold= %f keV", threshold/keV) << G4endl;
		//if((*pCItr).second > m_digitIn.thl)
		if((*pCItr).second > threshold)
		{
			// create one digit per pixel, I need to look at all the pixels first
			AllPixFEI3StandardDigit * digit = new AllPixFEI3StandardDigit;
			digit->SetPixelIDX((*pCItr).first.first);
			digit->SetPixelIDY((*pCItr).first.second);
			digit->SetPixelCounts(EnergyToTOT((*pCItr).second/elec,threshold));
			//digit->IncreasePixelCounts(); // Counting mode
			//cout << "TOT= "<< EnergyToTOT((*pCItr).second/elec,threshold/elec) << endl;
			//cout << "Energy= "<< (*pCItr).second/elec << endl;

			// MC only //
			// Replicating the same information in all pixels
			// FIXME !
			digit->SetPrimaryVertex(m_primaryVertex->GetPosition());
			digit->SetPixelEnergyDep((*pCItr).second);

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

}





vector<G4double>  AllPixFEI3StandardDigitizer::RKF5Integration(G4double x, G4double y, G4double z,G4double dt)
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




//______________________________________________________________________________
G4double  AllPixFEI3StandardDigitizer::SetDt(G4double dt, G4double ErreurMoy)
{

	double Dt=dt;
	if(isnan(ErreurMoy)){Dt=tup;}
	else if(ErreurMoy > Target){ Dt*=0.9;}
	else if(ErreurMoy < Target){ Dt*=1.1;};


	if(Dt<tlow) Dt=tlow;
	if(Dt>tup)  Dt=tup;

	return Dt;

}

G4int AllPixFEI3StandardDigitizer::EnergyToTOT(G4double Energy, G4double threshold)
{

	G4double dEdt = (MipCharge- threshold)/(MipTOT*Lv1Unit);
	//G4int TOT = TMath::FloorNint((Energy-threshold)/(dEdt*Lv1Unit*ns));
	G4int TOT = TMath::FloorNint((Energy-threshold)/(dEdt*Lv1Unit));

	if(Energy<threshold) TOT=0;
	if (TOT<0) TOT=0;
	if (TOT>CounterDepth) TOT=CounterDepth;

	return TOT;
}


//G4int AllPixFEI3StandardDigitizer::EnergyToTOTFEI4(G4double Energy, G4double threshold)
//{
//	double ToT=(1.0/(Energy+calib_C)+ calib_A)*calib_B;
//
//	return TMath::Floor(ToT);
//}

//______________________________________________________________________________
G4double AllPixFEI3StandardDigitizer::SlimEdgeEffect(G4int nX,G4double xpos,G4double eHit)
{

	G4double eEdge;
	G4double Etemp = 0.;

	double xpix= xpos + pitchX/2;
	//cout << "!!!!!!!!!SLIMEDGES!!!!!!!!!!! " << MipCharge << " " << xpix << endl;


	if(nX==0)
	{
		if((xpix < GRShift)){
			eEdge = eHit*(1-(GRShift-xpix)/(300*um));
			//cout << TString::Format("!!!!!!SLIM EDGES!!!!!!!! left nX:%d xpos:%f xpix:%f before:%f after:%f, R:%f",nX,xpos,xpix,eHit,Etemp,eHit*(1-(GRShift-xpix)/(GRShift))) << endl;

		}
		else{eEdge=eHit;};
		Etemp=eEdge;

		}

	else if(nX==(nPixX-1)){
		if(xpix > pitchX-GRShift){
			eEdge = eHit*(1-(xpix-(pitchX-GRShift))/(300*um));
			//cout << TString::Format("!!!!!!SLIM EDGES!!!!!!!! right NX:%d xpos:%f xpix:%f before:%f after:%f R:%f",nX,xpos,xpix,eHit,Etemp,eHit*(1-(xpix-(pitchX-GRShift))/(GRShift))) << endl;

		}
		else{eEdge=eHit;};
		Etemp=eEdge;


	}

	return Etemp;
}

G4bool AllPixFEI3StandardDigitizer::isSlimEdge(G4int nX, G4int /*nY*/)
{
	if(nX==0||nX==(nPixX-1)) return true;
	else return false;

}













