/*
 * AllpixPixelsParameterization.cpp
 *
 *  Created on: 6 April 2014
 *      Author: nalipour
 */

#include "AllPixPixelsParameterization.hh"

Allpix_PixelsParameterization::Allpix_PixelsParameterization(AllPixGeoDsc * geo) {
	this->fGeoPars = geo;
	hsensorX=fGeoPars->GetHalfSensorX();
	hsensorY=fGeoPars->GetHalfSensorY();
	hpixelX=fGeoPars->GetHalfPixelX();
	hpixelY=fGeoPars->GetHalfPixelY();

	npixelX = fGeoPars->GetNPixelsX();
	npixelY = fGeoPars->GetNPixelsY();
}



void Allpix_PixelsParameterization::ComputeTransformation(G4int copyId,
		G4VPhysicalVolume* Pixel) const {

  G4double XPos = posX(copyId);//+ fGeoPars->GetBumpOffsetX();//Nilou
  G4double YPos = posY(copyId);// + fGeoPars->GetBumpOffsetY();//Nilou
  G4double ZPos = 0;
  
  //G4cout << "[PArameterization] placing bump : " << copyId << endl;
  Pixel->SetTranslation(G4ThreeVector(XPos,YPos,ZPos));
  Pixel->SetRotation(0);
  
}

double Allpix_PixelsParameterization::posX(int id) const {
  
  G4int X =  id%npixelX;
  return X*hpixelX*2 + hpixelX - hsensorX;
  
  
}

double Allpix_PixelsParameterization::posY(int id) const {
  
  G4int Y = (id-(id%npixelX))/npixelY;
  return Y*hpixelY*2 + hpixelY - hsensorY;

  
}
