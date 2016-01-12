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
  map<pair<G4int, G4int>, MC_content> pixelsContent_MC; //contains information with only MC
  map<pair<G4int, G4int>, G4double > pixelsContent;
  pair<G4int, G4int> tempPixel;

  G4int nEntries = hitsCollection->entries();

  AllPixGeoDsc * gD = GetDetectorGeoDscPtr();


  //Nilou
  G4int noise_pos[7][2]=
    {
      {0, 0},
      {0, 0},
      {0, 0},
      {0, 0},
      {0, 0},
      {0, 0},
      {0, 0}
    };
  G4int pos1, pos2=0; //Nilou


  for(G4int itr  = 0 ; itr < nEntries ; itr++) {




    // Position whitin pixel
    G4double xpos = (*hitsCollection)[itr]->GetPosWithRespectToPixel().x();
    G4double ypos = (*hitsCollection)[itr]->GetPosWithRespectToPixel().y();
    G4double zpos = (*hitsCollection)[itr]->GetPosWithRespectToPixel().z();
    //G4cout << "nalipour: xpos=" << xpos << ", ypos=" << ypos << G4endl;

    //G4cout << "Hit position ib z : " << zpos << endl;

    G4double dist[5] = {0,0,0,0,0};
    for(int i =0;i<5;i++){
      dist[i] = TMath::Sqrt(pow(fabs(xpos)-xsig[i],2)+pow(fabs(ypos)-ysig[i],2));
    }

    int cluster_type = indexofSmallestElement(dist,5);

    tempPixel.first  = (*hitsCollection)[itr]->GetPixelNbX();
    tempPixel.second = (*hitsCollection)[itr]->GetPixelNbY();
    pixelsContent_MC[tempPixel].MC_energy += (*hitsCollection)[itr]->GetEdep(); //nalipour
    pixelsContent[tempPixel]=pixelsContent_MC[tempPixel].MC_energy; //nalipour

    if(pixelsContent_MC[tempPixel].posX_WithRespectToPixel == -11.0 || pixelsContent_MC[tempPixel].posY_WithRespectToPixel == -11.0)
      {
	pixelsContent_MC[tempPixel].posX_WithRespectToPixel=(*hitsCollection)[itr]->GetPosWithRespectToPixel().x();
	pixelsContent_MC[tempPixel].posY_WithRespectToPixel=(*hitsCollection)[itr]->GetPosWithRespectToPixel().y();
	pixelsContent_MC[tempPixel].posZ_WithRespectToPixel=(*hitsCollection)[itr]->GetPosWithRespectToPixel().z();
      }

    if(zpos>(gD->GetSensorZ()/2.-depletionDepth)){
      switch (cluster_type) {

      case 0 : //Center of pixel (cluster size=1)
	break;


      case 1 : //cluster size 2 (x=2, y=1)
	if(xpos<0)
    	  {
	    tempPixel.first  = (*hitsCollection)[itr]->GetPixelNbX()-1;
	    tempPixel.second = (*hitsCollection)[itr]->GetPixelNbY();
	    if((tempPixel.first>=0 and tempPixel.second>=0) and (tempPixel.first<nPixX and tempPixel.second<nPixY))
	      {
		pixelsContent[tempPixel] += (*hitsCollection)[itr]->GetEdep();
	      }
    	  }
	else
    	  {
	    tempPixel.first  = (*hitsCollection)[itr]->GetPixelNbX()+1;
	    tempPixel.second = (*hitsCollection)[itr]->GetPixelNbY();
	    if((tempPixel.first>=0 and tempPixel.second>=0) and (tempPixel.first<nPixX and tempPixel.second<nPixY))
	      {
		pixelsContent[tempPixel] += (*hitsCollection)[itr]->GetEdep();
	      }
    	  }
	break;

      case 2 : //cluster size 2 (x=1, y=2)

	if(ypos<0)
    	  {
	    tempPixel.first  = (*hitsCollection)[itr]->GetPixelNbX();
	    tempPixel.second = (*hitsCollection)[itr]->GetPixelNbY()-1;
	    if((tempPixel.first>=0 and tempPixel.second>=0) and (tempPixel.first<nPixX and tempPixel.second<nPixY))
	      {
		pixelsContent[tempPixel] += (*hitsCollection)[itr]->GetEdep();
	      }

    	  }
	else
    	  {
	    tempPixel.first  = (*hitsCollection)[itr]->GetPixelNbX();
	    tempPixel.second = (*hitsCollection)[itr]->GetPixelNbY()+1;
	    if((tempPixel.first>=0 and tempPixel.second>=0) and (tempPixel.first<nPixX and tempPixel.second<nPixY))
	      {
		pixelsContent[tempPixel] += (*hitsCollection)[itr]->GetEdep();
	      }
    	  }
	break;

      case 3 : //cluster size 4
	if(xpos<0)
    	  {
	    if(ypos<0)
	      {

		tempPixel.first  = (*hitsCollection)[itr]->GetPixelNbX();
		tempPixel.second = (*hitsCollection)[itr]->GetPixelNbY()-1;
		if((tempPixel.first>=0 and tempPixel.second>=0) and (tempPixel.first<nPixX and tempPixel.second<nPixY))
		  {
		    pixelsContent[tempPixel] += (*hitsCollection)[itr]->GetEdep();
		  }

		tempPixel.first  = (*hitsCollection)[itr]->GetPixelNbX()-1;
		tempPixel.second = (*hitsCollection)[itr]->GetPixelNbY();
		if((tempPixel.first>=0 and tempPixel.second>=0) and (tempPixel.first<nPixX and tempPixel.second<nPixY))
		  {
		    pixelsContent[tempPixel] += (*hitsCollection)[itr]->GetEdep();
		  }

		tempPixel.first  = (*hitsCollection)[itr]->GetPixelNbX()-1;
		tempPixel.second = (*hitsCollection)[itr]->GetPixelNbY()-1;
		if((tempPixel.first>=0 and tempPixel.second>=0) and (tempPixel.first<nPixX and tempPixel.second<nPixY))
		  {
		    pixelsContent[tempPixel] += (*hitsCollection)[itr]->GetEdep();
		  }

	      }
	    else
	      {
		tempPixel.first  = (*hitsCollection)[itr]->GetPixelNbX();
		tempPixel.second = (*hitsCollection)[itr]->GetPixelNbY()-1;
		if((tempPixel.first>=0 and tempPixel.second>=0) and (tempPixel.first<nPixX and tempPixel.second<nPixY))
		  {
		    pixelsContent[tempPixel] += (*hitsCollection)[itr]->GetEdep();
		  }

		tempPixel.first  = (*hitsCollection)[itr]->GetPixelNbX()-1;
		tempPixel.second = (*hitsCollection)[itr]->GetPixelNbY();
		if((tempPixel.first>=0 and tempPixel.second>=0) and (tempPixel.first<nPixX and tempPixel.second<nPixY))
		  {
		    pixelsContent[tempPixel] += (*hitsCollection)[itr]->GetEdep();
		  }

		tempPixel.first  = (*hitsCollection)[itr]->GetPixelNbX()-1;
		tempPixel.second = (*hitsCollection)[itr]->GetPixelNbY()+1;
		if((tempPixel.first>=0 and tempPixel.second>=0) and (tempPixel.first<nPixX and tempPixel.second<nPixY))
		  {
		    pixelsContent[tempPixel] += (*hitsCollection)[itr]->GetEdep();
		  }

	      }
    	  }
	else
    	  {
	    if(ypos<0)
	      {
		tempPixel.first  = (*hitsCollection)[itr]->GetPixelNbX();
		tempPixel.second = (*hitsCollection)[itr]->GetPixelNbY()-1;
		if((tempPixel.first>=0 and tempPixel.second>=0) and (tempPixel.first<nPixX and tempPixel.second<nPixY))
		  {
		    pixelsContent[tempPixel] += (*hitsCollection)[itr]->GetEdep();
		  }

		tempPixel.first  = (*hitsCollection)[itr]->GetPixelNbX()+1;
		tempPixel.second = (*hitsCollection)[itr]->GetPixelNbY();
		if((tempPixel.first>=0 and tempPixel.second>=0) and (tempPixel.first<nPixX and tempPixel.second<nPixY))
		  {
		    pixelsContent[tempPixel] += (*hitsCollection)[itr]->GetEdep();
		  }

		tempPixel.first  = (*hitsCollection)[itr]->GetPixelNbX()+1;
		tempPixel.second = (*hitsCollection)[itr]->GetPixelNbY()-1;
		if((tempPixel.first>=0 and tempPixel.second>=0) and (tempPixel.first<nPixX and tempPixel.second<nPixY))
		  {
		    pixelsContent[tempPixel] += (*hitsCollection)[itr]->GetEdep();
		  }

	      }
	    else
	      {
		tempPixel.first  = (*hitsCollection)[itr]->GetPixelNbX();
		tempPixel.second = (*hitsCollection)[itr]->GetPixelNbY()+1;
		if((tempPixel.first>=0 and tempPixel.second>=0) and (tempPixel.first<nPixX and tempPixel.second<nPixY))
		  {
		    pixelsContent[tempPixel] += (*hitsCollection)[itr]->GetEdep();
		  }

		tempPixel.first  = (*hitsCollection)[itr]->GetPixelNbX()+1;
		tempPixel.second = (*hitsCollection)[itr]->GetPixelNbY();
		if((tempPixel.first>=0 and tempPixel.second>=0) and (tempPixel.first<nPixX and tempPixel.second<nPixY))
		  {
		    pixelsContent[tempPixel] += (*hitsCollection)[itr]->GetEdep();
		  }

		tempPixel.first  = (*hitsCollection)[itr]->GetPixelNbX()+1;
		tempPixel.second = (*hitsCollection)[itr]->GetPixelNbY()+1;
		if((tempPixel.first>=0 and tempPixel.second>=0) and (tempPixel.first<nPixX and tempPixel.second<nPixY))
		  {
		    pixelsContent[tempPixel] += (*hitsCollection)[itr]->GetEdep();
		  }
	      }
    	  }
	// //nalipour start
	// if(CLHEP::RandFlat::shoot(1)<0.4)
	//   {
	//     G4int pos_x_add=(G4int) CLHEP::RandFlat::shoot(2);
	//     G4int pos_y_add=(G4int) CLHEP::RandFlat::shoot(2);
	//     G4cout <<"nalipour: case 3" << pos_x_add << G4endl;
	//     G4cout <<"nalipour: case 3" << pos_y_add << G4endl;
	    
	//     pos_x_add=TMath::Power(-1, pos_x_add);
	//     pos_y_add=TMath::Power(-1, pos_y_add);
	//     tempPixel.first  = (*hitsCollection)[itr]->GetPixelNbX()+pos_x_add;
	//     tempPixel.second = (*hitsCollection)[itr]->GetPixelNbY()+pos_y_add;
	//     if((tempPixel.first>=0 and tempPixel.second>=0) and (tempPixel.first<nPixX and tempPixel.second<nPixY))
	//       {
	// 	pixelsContent[tempPixel] += (*hitsCollection)[itr]->GetEdep();
	//       }
	//   }
	// //nalipour end
	break;

      case 4 : //cluster size 3

	if(xpos<0)
    	  {
	    if(ypos<0)
	      {
		tempPixel.first  = (*hitsCollection)[itr]->GetPixelNbX();
		tempPixel.second = (*hitsCollection)[itr]->GetPixelNbY()-1;
		if((tempPixel.first>=0 and tempPixel.second>=0) and (tempPixel.first<nPixX and tempPixel.second<nPixY))
		  {
		    pixelsContent[tempPixel] += (*hitsCollection)[itr]->GetEdep();
		  }

		tempPixel.first  = (*hitsCollection)[itr]->GetPixelNbX()-1;
		tempPixel.second = (*hitsCollection)[itr]->GetPixelNbY();
		if((tempPixel.first>=0 and tempPixel.second>=0) and (tempPixel.first<nPixX and tempPixel.second<nPixY))
		  {
		    pixelsContent[tempPixel] += (*hitsCollection)[itr]->GetEdep();
		  }

		if(CLHEP::RandFlat::shoot(1)<0.15)
		  {
		    tempPixel.first  = (*hitsCollection)[itr]->GetPixelNbX()-1;
		    tempPixel.second = (*hitsCollection)[itr]->GetPixelNbY()-1;
		    if((tempPixel.first>=0 and tempPixel.second>=0) and (tempPixel.first<nPixX and tempPixel.second<nPixY))
		      {
			pixelsContent[tempPixel] += (*hitsCollection)[itr]->GetEdep();
		      }
		  }
		//Nilou
		if(CLHEP::RandFlat::shoot(1)<0.15)
		  {
		    noise_pos[0][0]=-2;
		    noise_pos[0][1]=0;

		    noise_pos[1][0]=-1;
		    noise_pos[1][1]=1;

		    noise_pos[2][0]=0;
		    noise_pos[2][1]=1;

		    noise_pos[3][0]=1;
		    noise_pos[3][1]=0;

		    noise_pos[4][0]=1;
		    noise_pos[4][1]=-1;

		    noise_pos[5][0]=0;
		    noise_pos[5][1]=-2;
		
		    pos1 = (G4int) CLHEP::RandFlat::shoot(7);
		    pos2 = (G4int) CLHEP::RandFlat::shoot(7);
		    //G4cout << "Nilou: pos1 = " << pos1 << ", pos2=" << pos2 << G4endl;
		    if (pos1 != 6)
		      {
			tempPixel.first  = (*hitsCollection)[itr]->GetPixelNbX()+noise_pos[pos1][0];
			tempPixel.second = (*hitsCollection)[itr]->GetPixelNbY()+noise_pos[pos1][1];
			if((tempPixel.first>=0 and tempPixel.second>=0) and (tempPixel.first<nPixX and tempPixel.second<nPixY)) pixelsContent[tempPixel] += (*hitsCollection)[itr]->GetEdep();
		      }
		    else
		      {
			//G4cout << "Nilou: pos1 NO" << G4endl;
		      }
		    if (pos2 != 6)
		      {
			tempPixel.first  = (*hitsCollection)[itr]->GetPixelNbX()+noise_pos[pos2][0];
			tempPixel.second = (*hitsCollection)[itr]->GetPixelNbY()+noise_pos[pos2][1];
			if((tempPixel.first>=0 and tempPixel.second>=0) and (tempPixel.first<nPixX and tempPixel.second<nPixY)) pixelsContent[tempPixel] += (*hitsCollection)[itr]->GetEdep();
		      }
		    else
		      {
			//G4cout << "Nilou: pos2 NO" << G4endl;
		      }
		  }
		//End Nilou
	      }
	    else
	      {
		tempPixel.first  = (*hitsCollection)[itr]->GetPixelNbX();
		tempPixel.second = (*hitsCollection)[itr]->GetPixelNbY()+1;
		if((tempPixel.first>=0 and tempPixel.second>=0) and (tempPixel.first<nPixX and tempPixel.second<nPixY))
		  {
		    pixelsContent[tempPixel] += (*hitsCollection)[itr]->GetEdep();
		  }

		tempPixel.first  = (*hitsCollection)[itr]->GetPixelNbX()-1;
		tempPixel.second = (*hitsCollection)[itr]->GetPixelNbY();
		if((tempPixel.first>=0 and tempPixel.second>=0) and (tempPixel.first<nPixX and tempPixel.second<nPixY))
		  {
		    pixelsContent[tempPixel] += (*hitsCollection)[itr]->GetEdep();
		  }
		if(CLHEP::RandFlat::shoot(1)<0.15)
		  {
		    tempPixel.first  = (*hitsCollection)[itr]->GetPixelNbX()-1;
		    tempPixel.second = (*hitsCollection)[itr]->GetPixelNbY()+1;
		    if((tempPixel.first>=0 and tempPixel.second>=0) and (tempPixel.first<nPixX and tempPixel.second<nPixY)) pixelsContent[tempPixel] += (*hitsCollection)[itr]->GetEdep();
		  }

		//Nilou
		if(CLHEP::RandFlat::shoot(1)<0.15)
		  {
		    noise_pos[0][0]=-2;
		    noise_pos[0][1]=0;

		    noise_pos[1][0]=-1;
		    noise_pos[1][1]=-1;

		    noise_pos[2][0]=-1;
		    noise_pos[2][1]=0;

		    noise_pos[3][0]=1;
		    noise_pos[3][1]=0;

		    noise_pos[4][0]=1;
		    noise_pos[4][1]=1;

		    noise_pos[5][0]=0;
		    noise_pos[5][1]=2;
	      
		    pos1 = (G4int) CLHEP::RandFlat::shoot(7);
		    pos2 = (G4int) CLHEP::RandFlat::shoot(7);
		    //G4cout << "Nilou: pos1 = " << pos1 << ", pos2=" << pos2 << G4endl;
		    if (pos1 != 6)
		      {
			tempPixel.first  = (*hitsCollection)[itr]->GetPixelNbX()+noise_pos[pos1][0];
			tempPixel.second = (*hitsCollection)[itr]->GetPixelNbY()+noise_pos[pos1][1];
			if((tempPixel.first>=0 and tempPixel.second>=0) and (tempPixel.first<nPixX and tempPixel.second<nPixY)) pixelsContent[tempPixel] += (*hitsCollection)[itr]->GetEdep();
		      }
		    else
		      {
			//G4cout << "Nilou: pos1 NO" << G4endl;
		      }
		    if (pos2 != 6)
		      {
			tempPixel.first  = (*hitsCollection)[itr]->GetPixelNbX()+noise_pos[pos2][0];
			tempPixel.second = (*hitsCollection)[itr]->GetPixelNbY()+noise_pos[pos2][1];
			if((tempPixel.first>=0 and tempPixel.second>=0) and (tempPixel.first<nPixX and tempPixel.second<nPixY)) pixelsContent[tempPixel] += (*hitsCollection)[itr]->GetEdep();
		      }
		    else
		      {
			//G4cout << "Nilou: pos2 NO" << G4endl;
		      }
		  }
		//End Nilou

	      }
    	  }
	else
	  {
	    if(ypos<0)
	      {
		tempPixel.first  = (*hitsCollection)[itr]->GetPixelNbX();
		tempPixel.second = (*hitsCollection)[itr]->GetPixelNbY()-1;
		if((tempPixel.first>=0 and tempPixel.second>=0) and (tempPixel.first<nPixX and tempPixel.second<nPixY))
		  {
		    pixelsContent[tempPixel] += (*hitsCollection)[itr]->GetEdep();
		  }

		tempPixel.first  = (*hitsCollection)[itr]->GetPixelNbX()+1;
		tempPixel.second = (*hitsCollection)[itr]->GetPixelNbY();
		if((tempPixel.first>=0 and tempPixel.second>=0) and (tempPixel.first<nPixX and tempPixel.second<nPixY))
		  {
		    pixelsContent[tempPixel] += (*hitsCollection)[itr]->GetEdep();
		  }

		if(CLHEP::RandFlat::shoot(1)<0.15)
		  {
		    tempPixel.first  = (*hitsCollection)[itr]->GetPixelNbX()+1;
		    tempPixel.second = (*hitsCollection)[itr]->GetPixelNbY()-1;
		    if((tempPixel.first>=0 and tempPixel.second>=0) and (tempPixel.first<nPixX and tempPixel.second<nPixY)) pixelsContent[tempPixel] += (*hitsCollection)[itr]->GetEdep();
		  }

		//Nilou
		if(CLHEP::RandFlat::shoot(1)<0.15)
		  {
		    noise_pos[0][0]=2;
		    noise_pos[0][1]=0;

		    noise_pos[1][0]=1;
		    noise_pos[1][1]=1;

		    noise_pos[2][0]=0;
		    noise_pos[2][1]=1;

		    noise_pos[3][0]=-1;
		    noise_pos[3][1]=0;

		    noise_pos[4][0]=-1;
		    noise_pos[4][1]=-1;

		    noise_pos[5][0]=0;
		    noise_pos[5][1]=-2;
	      
		    pos1 = (G4int) CLHEP::RandFlat::shoot(7);
		    pos2 = (G4int) CLHEP::RandFlat::shoot(7);
		    //G4cout << "Nilou: pos1 = " << pos1 << ", pos2=" << pos2 << G4endl;
		    if (pos1 != 6)
		      {
			tempPixel.first  = (*hitsCollection)[itr]->GetPixelNbX()+noise_pos[pos1][0];
			tempPixel.second = (*hitsCollection)[itr]->GetPixelNbY()+noise_pos[pos1][1];
			if((tempPixel.first>=0 and tempPixel.second>=0) and (tempPixel.first<nPixX and tempPixel.second<nPixY)) pixelsContent[tempPixel] += (*hitsCollection)[itr]->GetEdep();
		      }
		    else
		      {
			//G4cout << "Nilou: pos1 NO" << G4endl;
		      }
		    if (pos2 != 6)
		      {
			tempPixel.first  = (*hitsCollection)[itr]->GetPixelNbX()+noise_pos[pos2][0];
			tempPixel.second = (*hitsCollection)[itr]->GetPixelNbY()+noise_pos[pos2][1];
			if((tempPixel.first>=0 and tempPixel.second>=0) and (tempPixel.first<nPixX and tempPixel.second<nPixY)) pixelsContent[tempPixel] += (*hitsCollection)[itr]->GetEdep();
		      }
		    else
		      {
			//G4cout << "Nilou: pos2 NO" << G4endl;
		      }
		  }
		//End Nilou



	      }

	    else
	      {
		tempPixel.first  = (*hitsCollection)[itr]->GetPixelNbX();
		tempPixel.second = (*hitsCollection)[itr]->GetPixelNbY()+1;
		if((tempPixel.first>=0 and tempPixel.second>=0) and (tempPixel.first<nPixX and tempPixel.second<nPixY))
		  {
		    pixelsContent[tempPixel] += (*hitsCollection)[itr]->GetEdep();
		  }

		tempPixel.first  = (*hitsCollection)[itr]->GetPixelNbX()+1;
		tempPixel.second = (*hitsCollection)[itr]->GetPixelNbY();
		if((tempPixel.first>=0 and tempPixel.second>=0) and (tempPixel.first<nPixX and tempPixel.second<nPixY))
		  {
		    pixelsContent[tempPixel] += (*hitsCollection)[itr]->GetEdep();
		  }

		if(CLHEP::RandFlat::shoot(1)<0.15)
		  {
		    tempPixel.first  = (*hitsCollection)[itr]->GetPixelNbX()+1;
		    tempPixel.second = (*hitsCollection)[itr]->GetPixelNbY()+1;
		    if((tempPixel.first>=0 and tempPixel.second>=0) and (tempPixel.first<nPixX and tempPixel.second<nPixY))
		      {
			pixelsContent[tempPixel] += (*hitsCollection)[itr]->GetEdep();
		      }
		  }



		//Nilou
		if(CLHEP::RandFlat::shoot(1)<0.15)
		  {
		    noise_pos[0][0]=2;
		    noise_pos[0][1]=0;

		    noise_pos[1][0]=-1;
		    noise_pos[1][1]=-1;

		    noise_pos[2][0]=0;
		    noise_pos[2][1]=-1;

		    noise_pos[3][0]=-1;
		    noise_pos[3][1]=0;

		    noise_pos[4][0]=-1;
		    noise_pos[4][1]=-1;

		    noise_pos[5][0]=0;
		    noise_pos[5][1]=2;
	      
		    pos1 = (G4int) CLHEP::RandFlat::shoot(7);
		    pos2 = (G4int) CLHEP::RandFlat::shoot(7);
		    //G4cout << "Nilou: pos1 = " << pos1 << ", pos2=" << pos2 << G4endl;
		    if (pos1 != 6)
		      {
			tempPixel.first  = (*hitsCollection)[itr]->GetPixelNbX()+noise_pos[pos1][0];
			tempPixel.second = (*hitsCollection)[itr]->GetPixelNbY()+noise_pos[pos1][1];
			if((tempPixel.first>=0 and tempPixel.second>=0) and (tempPixel.first<nPixX and tempPixel.second<nPixY)) pixelsContent[tempPixel] += (*hitsCollection)[itr]->GetEdep();
		      }
		    else
		      {
			//G4cout << "Nilou: pos1 NO" << G4endl;
		      }
		    if (pos2 != 6)
		      {
			tempPixel.first  = (*hitsCollection)[itr]->GetPixelNbX()+noise_pos[pos2][0];
			tempPixel.second = (*hitsCollection)[itr]->GetPixelNbY()+noise_pos[pos2][1];
			if((tempPixel.first>=0 and tempPixel.second>=0) and (tempPixel.first<nPixX and tempPixel.second<nPixY)) pixelsContent[tempPixel] += (*hitsCollection)[itr]->GetEdep();
		      }
		    else
		      {
			//G4cout << "Nilou: pos2 NO" << G4endl;
		      }
		  }
		//End Nilou
	      }

	  }


	break;
      }
    }

  }

  //------------------ RECORD DIGITS ------------------//
  // With charge sharing
  map<pair<G4int, G4int>, G4double >::iterator pCItr = pixelsContent.begin();
  for( ; pCItr != pixelsContent.end() ; pCItr++)
    {
      if((*pCItr).second > m_digitIn.thl) // over threshold !
	{
	  tempPixel.first=(*pCItr).first.first;
	  tempPixel.second=(*pCItr).first.second;

	  AllPixMimosa26Digit * digit = new AllPixMimosa26Digit;

	  digit->SetPixelIDX((*pCItr).first.first);
	  digit->SetPixelIDY((*pCItr).first.second);
	  digit->SetPixelEnergyDep((*pCItr).second/keV); //Energy with charge sharing
	  //digit->SetPixelCounts((*pCItr).second/eV); //TOT value
	  if (pixelsContent_MC.count(tempPixel))
	    {
	      digit->SetPixelEnergyMC(pixelsContent_MC[tempPixel].MC_energy/keV); //MC value
	      digit->Set_posX_WithRespectoToPixel(pixelsContent_MC[tempPixel].posX_WithRespectToPixel);
	      digit->Set_posY_WithRespectoToPixel(pixelsContent_MC[tempPixel].posY_WithRespectToPixel);
	      digit->Set_posZ_WithRespectoToPixel(pixelsContent_MC[tempPixel].posZ_WithRespectToPixel);
	      pixelsContent_MC.erase(tempPixel);
	    }
	  digit->IncreasePixelCounts(); // Counting mode

	  m_digitsCollection->insert(digit);
	}
    }
  //MC only
  map<pair<G4int, G4int>, MC_content >::iterator pCItr_MC = pixelsContent_MC.begin();
  for( ; pCItr_MC != pixelsContent_MC.end() ; pCItr_MC++)
    {
      tempPixel.first=(*pCItr_MC).first.first;
      tempPixel.second=(*pCItr_MC).first.second;

      AllPixMimosa26Digit * digit = new AllPixMimosa26Digit;

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


  /*

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
  digit->SetPixelEnergyDep((*pCItr).second/keV); //nalipour

  m_digitsCollection->insert(digit);
  }

  }
  */

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
      //digit->SetPixelEnergyDep(0.); // FIXME !! edep in noise //nalipour commented

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
