/**
 *  Author John Idarraga <idarraga@cern.ch>
 */

#include "AllPixMimosa26Digitizer.hh"
#include "AllPixTrackerHit.hh"
#include "G4DigiManager.hh"
#include "G4VDigitizerModule.hh"
#include "AllPixGeoDsc.hh"

#include "TRandom2.h"
#include "TMath.h"
#include "TF1.h"

#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGaussQ.h" // faster than RandGauss, less accurate

AllPixMimosa26Digitizer::AllPixMimosa26Digitizer(G4String modName, G4String hitsColName, G4String digitColName) 
: AllPixDigitizerInterface (modName) {

	// Registration of digits collection name
	collectionName.push_back(digitColName);
	m_hitsColName.push_back(hitsColName);

	// input parameters
	m_digitIn.thl = 0.2*keV;
	m_randomNoise = 0; // if 0, no random noise

	// Example of detector description handle
	// provided by the interface
	AllPixGeoDsc * gD = GetDetectorGeoDscPtr();
	////////////////////////////////
	// Geometry Related constants //
	////////////////////////////////
	pitchX = gD->GetPixelX();
	pitchY = gD->GetPixelY();

	hpitchX = gD->GetPixelX()/2.;
	hpitchY = gD->GetPixelY()/2.;

	thickness = gD->GetHalfSensorZ()*2;


	nPixX= gD->GetNPixelsX();
	nPixY= gD->GetNPixelsY();

	depletionDepth=5*um;
	crosstalk_prob=0.1;

	// This is defined in the stack and will free memory as we leave this scope
	G4double xs[5] = {0,hpitchX,0,hpitchX,pitchX/3.};
	G4double ys[5] = {0,0,hpitchY,hpitchY,pitchY/3.};
	// This remains in the heap
	xsig = new G4double[__SIG_SIZE/8];
	ysig = new G4double[__SIG_SIZE/8];
	// Copy the values
	memcpy(xsig, xs, __SIG_SIZE);
	memcpy(ysig, ys, __SIG_SIZE);

}

AllPixMimosa26Digitizer::~AllPixMimosa26Digitizer(){

	delete [] xsig;
	delete [] ysig;

}

int AllPixMimosa26Digitizer::indexofSmallestElement(double array[], int size)
{
    int index = 0;

    for(int i = 1; i < size; i++)
    {
        if(array[i] < array[index])
            index = i;
    }

    return index;
}



void AllPixMimosa26Digitizer::SetDetectorDigitInputs(G4double /*thl*/){

	// set digitization input values
	// thl

	//m_digitIn.thl = thl; // <-- input !
	//G4cout << "AllPixMimosa26Digitizer: threshhold = " << m_digitIn.thl/keV << " [keV]" << G4cout;

}

void AllPixMimosa26Digitizer::Digitize(){

	// create the digits collection
	// FIXME
	// should this be here ? inside the loop ?
	m_digitsCollection = new AllPixMimosa26DigitsCollection("AllPixMimosa26Digitizer", collectionName[0] );

	// get the digiManager
	G4DigiManager * digiMan = G4DigiManager::GetDMpointer();

	// BoxSD_0_HitsCollection
	G4int hcID = digiMan->GetHitsCollectionID(m_hitsColName[0]);

	AllPixTrackerHitsCollection * hitsCollection = 0;
	hitsCollection = (AllPixTrackerHitsCollection*)(digiMan->GetHitsCollection(hcID));

	// temporary data structure
	map<pair<G4int, G4int>, G4double > pixelsContent;
	pair<G4int, G4int> tempPixel;

	G4int nEntries = hitsCollection->entries();

	AllPixGeoDsc * gD = GetDetectorGeoDscPtr();


	for(G4int itr  = 0 ; itr < nEntries ; itr++) {




		// Position whitin pixel
		G4double xpos = (*hitsCollection)[itr]->GetPosWithRespectToPixel().x();
		G4double ypos = (*hitsCollection)[itr]->GetPosWithRespectToPixel().y();
		G4double zpos = (*hitsCollection)[itr]->GetPosWithRespectToPixel().z();

		//G4cout << "Hit position ib z : " << zpos << endl;

		G4double dist[5] = {0,0,0,0,0};
		for(int i =0;i<5;i++){
			dist[i] = TMath::Sqrt(pow(fabs(xpos)-xsig[i],2)+pow(fabs(ypos)-ysig[i],2));
		}

		int cluster_type = indexofSmallestElement(dist,5);

		if(zpos>(gD->GetSensorZ()/2.-depletionDepth)){
		switch (cluster_type) {

		case 0 :
			tempPixel.first  = (*hitsCollection)[itr]->GetPixelNbX();
			tempPixel.second = (*hitsCollection)[itr]->GetPixelNbY();
			pixelsContent[tempPixel] += (*hitsCollection)[itr]->GetEdep();
			break;


		case 1 :

			tempPixel.first  = (*hitsCollection)[itr]->GetPixelNbX();
			tempPixel.second = (*hitsCollection)[itr]->GetPixelNbY();
			pixelsContent[tempPixel] += (*hitsCollection)[itr]->GetEdep();

			if(xpos<0){
				tempPixel.first  = (*hitsCollection)[itr]->GetPixelNbX()-1;
				tempPixel.second = (*hitsCollection)[itr]->GetPixelNbY();
				if((tempPixel.first>=0 and tempPixel.second>=0) and (tempPixel.first<nPixX and tempPixel.second<nPixY)) pixelsContent[tempPixel] += (*hitsCollection)[itr]->GetEdep();

			}
			else {
				tempPixel.first  = (*hitsCollection)[itr]->GetPixelNbX()+1;
				tempPixel.second = (*hitsCollection)[itr]->GetPixelNbY();
				if((tempPixel.first>=0 and tempPixel.second>=0) and (tempPixel.first<nPixX and tempPixel.second<nPixY)) pixelsContent[tempPixel] += (*hitsCollection)[itr]->GetEdep();
			}

			break;

		case 2 :
			tempPixel.first  = (*hitsCollection)[itr]->GetPixelNbX();
			tempPixel.second = (*hitsCollection)[itr]->GetPixelNbY();
			pixelsContent[tempPixel] += (*hitsCollection)[itr]->GetEdep();

			if(ypos<0){
				tempPixel.first  = (*hitsCollection)[itr]->GetPixelNbX();
				tempPixel.second = (*hitsCollection)[itr]->GetPixelNbY()-1;
				if((tempPixel.first>=0 and tempPixel.second>=0) and (tempPixel.first<nPixX and tempPixel.second<nPixY)) pixelsContent[tempPixel] += (*hitsCollection)[itr]->GetEdep();

			}
			else {
				tempPixel.first  = (*hitsCollection)[itr]->GetPixelNbX();
				tempPixel.second = (*hitsCollection)[itr]->GetPixelNbY()+1;
				if((tempPixel.first>=0 and tempPixel.second>=0) and (tempPixel.first<nPixX and tempPixel.second<nPixY)) pixelsContent[tempPixel] += (*hitsCollection)[itr]->GetEdep();
			}
			break;

		case 3 :
			tempPixel.first  = (*hitsCollection)[itr]->GetPixelNbX();
			tempPixel.second = (*hitsCollection)[itr]->GetPixelNbY();
			pixelsContent[tempPixel] += (*hitsCollection)[itr]->GetEdep();

			if(xpos<0){

				if(ypos<0){

					tempPixel.first  = (*hitsCollection)[itr]->GetPixelNbX();
					tempPixel.second = (*hitsCollection)[itr]->GetPixelNbY()-1;
					if((tempPixel.first>=0 and tempPixel.second>=0) and (tempPixel.first<nPixX and tempPixel.second<nPixY)) pixelsContent[tempPixel] += (*hitsCollection)[itr]->GetEdep();

					tempPixel.first  = (*hitsCollection)[itr]->GetPixelNbX()-1;
					tempPixel.second = (*hitsCollection)[itr]->GetPixelNbY();
					if((tempPixel.first>=0 and tempPixel.second>=0) and (tempPixel.first<nPixX and tempPixel.second<nPixY)) pixelsContent[tempPixel] += (*hitsCollection)[itr]->GetEdep();

					tempPixel.first  = (*hitsCollection)[itr]->GetPixelNbX()-1;
					tempPixel.second = (*hitsCollection)[itr]->GetPixelNbY()-1;
					if((tempPixel.first>=0 and tempPixel.second>=0) and (tempPixel.first<nPixX and tempPixel.second<nPixY)) pixelsContent[tempPixel] += (*hitsCollection)[itr]->GetEdep();

				}
				else {
					tempPixel.first  = (*hitsCollection)[itr]->GetPixelNbX();
					tempPixel.second = (*hitsCollection)[itr]->GetPixelNbY()-1;
					if((tempPixel.first>=0 and tempPixel.second>=0) and (tempPixel.first<nPixX and tempPixel.second<nPixY)) pixelsContent[tempPixel] += (*hitsCollection)[itr]->GetEdep();

					tempPixel.first  = (*hitsCollection)[itr]->GetPixelNbX()-1;
					tempPixel.second = (*hitsCollection)[itr]->GetPixelNbY();
					if((tempPixel.first>=0 and tempPixel.second>=0) and (tempPixel.first<nPixX and tempPixel.second<nPixY)) pixelsContent[tempPixel] += (*hitsCollection)[itr]->GetEdep();

					tempPixel.first  = (*hitsCollection)[itr]->GetPixelNbX()-1;
					tempPixel.second = (*hitsCollection)[itr]->GetPixelNbY()+1;
					if((tempPixel.first>=0 and tempPixel.second>=0) and (tempPixel.first<nPixX and tempPixel.second<nPixY)) pixelsContent[tempPixel] += (*hitsCollection)[itr]->GetEdep();

				}


			}
			else{

				if(ypos<0){
					tempPixel.first  = (*hitsCollection)[itr]->GetPixelNbX();
					tempPixel.second = (*hitsCollection)[itr]->GetPixelNbY()-1;
					if((tempPixel.first>=0 and tempPixel.second>=0) and (tempPixel.first<nPixX and tempPixel.second<nPixY)) pixelsContent[tempPixel] += (*hitsCollection)[itr]->GetEdep();

					tempPixel.first  = (*hitsCollection)[itr]->GetPixelNbX()+1;
					tempPixel.second = (*hitsCollection)[itr]->GetPixelNbY();
					if((tempPixel.first>=0 and tempPixel.second>=0) and (tempPixel.first<nPixX and tempPixel.second<nPixY)) pixelsContent[tempPixel] += (*hitsCollection)[itr]->GetEdep();

					tempPixel.first  = (*hitsCollection)[itr]->GetPixelNbX()+1;
					tempPixel.second = (*hitsCollection)[itr]->GetPixelNbY()-1;
					if((tempPixel.first>=0 and tempPixel.second>=0) and (tempPixel.first<nPixX and tempPixel.second<nPixY)) pixelsContent[tempPixel] += (*hitsCollection)[itr]->GetEdep();

				}

				else {

					tempPixel.first  = (*hitsCollection)[itr]->GetPixelNbX();
					tempPixel.second = (*hitsCollection)[itr]->GetPixelNbY()+1;
					if((tempPixel.first>=0 and tempPixel.second>=0) and (tempPixel.first<nPixX and tempPixel.second<nPixY)) pixelsContent[tempPixel] += (*hitsCollection)[itr]->GetEdep();

					tempPixel.first  = (*hitsCollection)[itr]->GetPixelNbX()+1;
					tempPixel.second = (*hitsCollection)[itr]->GetPixelNbY();
					if((tempPixel.first>=0 and tempPixel.second>=0) and (tempPixel.first<nPixX and tempPixel.second<nPixY)) pixelsContent[tempPixel] += (*hitsCollection)[itr]->GetEdep();

					tempPixel.first  = (*hitsCollection)[itr]->GetPixelNbX()+1;
					tempPixel.second = (*hitsCollection)[itr]->GetPixelNbY()+1;
					if((tempPixel.first>=0 and tempPixel.second>=0) and (tempPixel.first<nPixX and tempPixel.second<nPixY)) pixelsContent[tempPixel] += (*hitsCollection)[itr]->GetEdep();

				}

			}

			break;

		case 4 :
			tempPixel.first  = (*hitsCollection)[itr]->GetPixelNbX();
			tempPixel.second = (*hitsCollection)[itr]->GetPixelNbY();
			pixelsContent[tempPixel] += (*hitsCollection)[itr]->GetEdep();

			if(xpos<0){

				if(ypos<0){

					tempPixel.first  = (*hitsCollection)[itr]->GetPixelNbX();
					tempPixel.second = (*hitsCollection)[itr]->GetPixelNbY()-1;
					if((tempPixel.first>=0 and tempPixel.second>=0) and (tempPixel.first<nPixX and tempPixel.second<nPixY)) pixelsContent[tempPixel] += (*hitsCollection)[itr]->GetEdep();

					tempPixel.first  = (*hitsCollection)[itr]->GetPixelNbX()-1;
					tempPixel.second = (*hitsCollection)[itr]->GetPixelNbY();
					if((tempPixel.first>=0 and tempPixel.second>=0) and (tempPixel.first<nPixX and tempPixel.second<nPixY)) pixelsContent[tempPixel] += (*hitsCollection)[itr]->GetEdep();

					if(CLHEP::RandFlat::shoot(1)<0.4){

						tempPixel.first  = (*hitsCollection)[itr]->GetPixelNbX()-1;
						tempPixel.second = (*hitsCollection)[itr]->GetPixelNbY()-1;
						if((tempPixel.first>=0 and tempPixel.second>=0) and (tempPixel.first<nPixX and tempPixel.second<nPixY)) pixelsContent[tempPixel] += (*hitsCollection)[itr]->GetEdep();
					}
				}

				else {
					tempPixel.first  = (*hitsCollection)[itr]->GetPixelNbX();
					tempPixel.second = (*hitsCollection)[itr]->GetPixelNbY()+1;
					if((tempPixel.first>=0 and tempPixel.second>=0) and (tempPixel.first<nPixX and tempPixel.second<nPixY)) pixelsContent[tempPixel] += (*hitsCollection)[itr]->GetEdep();

					tempPixel.first  = (*hitsCollection)[itr]->GetPixelNbX()-1;
					tempPixel.second = (*hitsCollection)[itr]->GetPixelNbY();
					if((tempPixel.first>=0 and tempPixel.second>=0) and (tempPixel.first<nPixX and tempPixel.second<nPixY)) pixelsContent[tempPixel] += (*hitsCollection)[itr]->GetEdep();

					if(CLHEP::RandFlat::shoot(1)<0.4){

						tempPixel.first  = (*hitsCollection)[itr]->GetPixelNbX()-1;
						tempPixel.second = (*hitsCollection)[itr]->GetPixelNbY()+1;
						if((tempPixel.first>=0 and tempPixel.second>=0) and (tempPixel.first<nPixX and tempPixel.second<nPixY)) pixelsContent[tempPixel] += (*hitsCollection)[itr]->GetEdep();
					}


				}


			}
			else{

				if(ypos<0){
					tempPixel.first  = (*hitsCollection)[itr]->GetPixelNbX();
					tempPixel.second = (*hitsCollection)[itr]->GetPixelNbY()-1;
					if((tempPixel.first>=0 and tempPixel.second>=0) and (tempPixel.first<nPixX and tempPixel.second<nPixY)) pixelsContent[tempPixel] += (*hitsCollection)[itr]->GetEdep();

					tempPixel.first  = (*hitsCollection)[itr]->GetPixelNbX()+1;
					tempPixel.second = (*hitsCollection)[itr]->GetPixelNbY();
					if((tempPixel.first>=0 and tempPixel.second>=0) and (tempPixel.first<nPixX and tempPixel.second<nPixY)) pixelsContent[tempPixel] += (*hitsCollection)[itr]->GetEdep();

					if(CLHEP::RandFlat::shoot(1)<0.4){

						tempPixel.first  = (*hitsCollection)[itr]->GetPixelNbX()+1;
						tempPixel.second = (*hitsCollection)[itr]->GetPixelNbY()-1;
						if((tempPixel.first>=0 and tempPixel.second>=0) and (tempPixel.first<nPixX and tempPixel.second<nPixY)) pixelsContent[tempPixel] += (*hitsCollection)[itr]->GetEdep();
					}
				}

				else {

					tempPixel.first  = (*hitsCollection)[itr]->GetPixelNbX();
					tempPixel.second = (*hitsCollection)[itr]->GetPixelNbY()+1;
					if((tempPixel.first>=0 and tempPixel.second>=0) and (tempPixel.first<nPixX and tempPixel.second<nPixY)) pixelsContent[tempPixel] += (*hitsCollection)[itr]->GetEdep();

					tempPixel.first  = (*hitsCollection)[itr]->GetPixelNbX()+1;
					tempPixel.second = (*hitsCollection)[itr]->GetPixelNbY();
					if((tempPixel.first>=0 and tempPixel.second>=0) and (tempPixel.first<nPixX and tempPixel.second<nPixY)) pixelsContent[tempPixel] += (*hitsCollection)[itr]->GetEdep();

					if(CLHEP::RandFlat::shoot(1)<0.4){

						tempPixel.first  = (*hitsCollection)[itr]->GetPixelNbX()+1;
						tempPixel.second = (*hitsCollection)[itr]->GetPixelNbY()+1;
						if((tempPixel.first>=0 and tempPixel.second>=0) and (tempPixel.first<nPixX and tempPixel.second<nPixY)) pixelsContent[tempPixel] += (*hitsCollection)[itr]->GetEdep();
					}
				}

			}
			break;
		}
		}

		}







/*		// track enters at Z positive and leave at Z negative, pixel surface a Z positive, so we take only hits within Z_surface and Z_surface - depletionDepth
		if(zpos>(thickness/2.-depletionDepth)){

		tempPixel.first  = (*hitsCollection)[itr]->GetPixelNbX();
		tempPixel.second = (*hitsCollection)[itr]->GetPixelNbY();
		pixelsContent[tempPixel] += (*hitsCollection)[itr]->GetEdep();

		//G4cout << "[Mimosa26] Hit Energy : " << (*hitsCollection)[itr]->GetEdep()/keV << " zpos : " << zpos << endl;

		//G4cout << "[Mimosa26] pitch X "  << hpitchX << " probX : " << fabs(xpos)/hpitchX << endl;

		G4double probX,probY,probXY;
		if(fabs(xpos)<hpitchX/2.){
			 probX = crosstalk_prob;
		}
		else{
			 probX = 1;
			//G4double probX = (fabs(xpos)/hpitchX);
		};

		if(fabs(ypos)<hpitchY/2.){
			 probY = crosstalk_prob;
		}
		else{

			 probY = 1;
			//G4double probY = (fabs(ypos)/hpitchY);
		};


		if(fabs(TMath::Sqrt(xpos*xpos+ypos*ypos))<TMath::Sqrt(hpitchX*hpitchX+hpitchY*hpitchY)/2.){
			 probY = crosstalk_prob;
		}
		else {
			 probXY = 1;
			//G4double probXY =(fabs(TMath::Sqrt(xpos*xpos+ypos*ypos))/TMath::Sqrt(hpitchX*hpitchX+hpitchY*hpitchY));
		};


		//G4cout << "[Mimosa26] "  << " xpos : " << xpos << " ypos :  " << ypos  <<  " probX : " << probX << " probY : " << probY << endl;



		if( CLHEP::RandFlat::shoot()<= probX){
			if(xpos>0){
			  tempPixel.first  += 1 ;
			  tempPixel.second = tempPixel.second;
			}
			else {
			  tempPixel.first  -= 1 ;
			  tempPixel.second = tempPixel.second;
			}
			//whitin matrix
			if((tempPixel.first>=0 and tempPixel.second>=0) and (tempPixel.first<nPixX and tempPixel.second<nPixY)) pixelsContent[tempPixel] += (*hitsCollection)[itr]->GetEdep();

		}



		if( CLHEP::RandFlat::shoot()<= probY){
			if(ypos>0){
			  tempPixel.first  = tempPixel.first ;
			  tempPixel.second +=1;
			}
			else {
			  tempPixel.first  = tempPixel.first ;
			  tempPixel.second -=1;
			}
			if((tempPixel.first>=0 and tempPixel.second>=0) and (tempPixel.first<nPixX and tempPixel.second<nPixY))	pixelsContent[tempPixel] += (*hitsCollection)[itr]->GetEdep();
		}



		if( CLHEP::RandFlat::shoot()<= probXY){
			if(xpos>0 and ypos>0){
			  tempPixel.first  +=1;
			  tempPixel.second +=1;
			}
			else if (xpos<0 and ypos<0) {
			  tempPixel.first  -=1;
			  tempPixel.second -=1;
			}

			else if (ypos<0 and xpos>0) {
				  tempPixel.first  +=1;
				  tempPixel.second -=1;
				}


		else {
			  tempPixel.first  -=1;
			  tempPixel.second +=1;
			}

			if((tempPixel.first>=0 and tempPixel.second>=0) and (tempPixel.first<nPixX and tempPixel.second<nPixY)) pixelsContent[tempPixel] += (*hitsCollection)[itr]->GetEdep();
			}

		}*/






	// now create digits, one per pixel // second entry in the map is the edep in the pixel
	map<pair<G4int, G4int>, G4double >::iterator pCItr = pixelsContent.begin();

	for( ; pCItr != pixelsContent.end() ; pCItr++)
	{

		//G4cout << "        " <<(*pCItr).first.first  << " --> " << (*pCItr).first.second << G4endl;

		if((*pCItr).second > m_digitIn.thl) // over threshold !
		{
			// create one digit per pixel, I need to look at all the pixels first
			AllPixMimosa26Digit * digit = new AllPixMimosa26Digit;
			digit->SetPixelIDX((*pCItr).first.first);
			digit->SetPixelIDY((*pCItr).first.second);
			//digit->SetPixelCounts(1);
			digit->IncreasePixelCounts(); // Couting mode

			// MC only //
			// Replicating the same information in all pixels
			// FIXME !
			digit->SetPrimaryVertex(m_primaryVertex->GetPosition());
			digit->SetPixelEnergyDep((*pCItr).second);

			m_digitsCollection->insert(digit);
		}

	}

	// Do random noise
	if ( m_randomNoise > 0 ) {

		int nNoisy = (int) CLHEP::RandGaussQ::shoot(double(m_randomNoise), double(m_randomNoise*0.1) );

		for (int ni = 0 ; ni < nNoisy ; ni++) {

			int xposn = CLHEP::RandFlat::shootInt(long(0),  long(nPixX) );
			int yposn = CLHEP::RandFlat::shootInt(long(0),  long(nPixY) );

			AllPixMimosa26Digit * digit = new AllPixMimosa26Digit;
			digit->SetPixelIDX(xposn);
			digit->SetPixelIDY(yposn);
			digit->IncreasePixelCounts(); // Couting mode

			digit->SetPrimaryVertex(m_primaryVertex->GetPosition());
			digit->SetPixelEnergyDep(0.); // FIXME !! edep in noise

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
