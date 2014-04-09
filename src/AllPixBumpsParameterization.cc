/*
 * AllpixBumpsParameterization.cpp
 *
 *  Created on: 31 janv. 2014
 *      Author: mbenoit
 */

#include "AllPixBumpsParameterization.hh"

Allpix_BumpsParameterization::Allpix_BumpsParameterization(AllPixGeoDsc * geo) {
	this->fGeoPars = geo;
	hsensorX=fGeoPars->GetHalfSensorX();
	hsensorY=fGeoPars->GetHalfSensorY();
	hpixelX=fGeoPars->GetHalfPixelX();
	hpixelY=fGeoPars->GetHalfPixelY();

	npixelX = fGeoPars->GetNPixelsX();
	npixelY = fGeoPars->GetNPixelsY();
}



void Allpix_BumpsParameterization::ComputeTransformation(G4int copyId,
		G4VPhysicalVolume* Bump) const {

	G4double XPos = posX(copyId) + fGeoPars->GetBumpOffsetX();
	G4double YPos = posY(copyId) + fGeoPars->GetBumpOffsetY();
	G4double ZPos = 0;

	//G4cout << "[PArameterization] placing bump : " << copyId << endl;
	Bump->SetTranslation(G4ThreeVector(XPos,YPos,ZPos));
	Bump->SetRotation(0);

}

double Allpix_BumpsParameterization::posX(int id) const {

	G4int X =  id%npixelX;
	return X*hpixelX*2 + hpixelX - hsensorX;


}

double Allpix_BumpsParameterization::posY(int id) const {

	G4int Y = (id-(id%npixelX))/npixelX;
	return Y*hpixelY*2 + hpixelY - hsensorY;


}
