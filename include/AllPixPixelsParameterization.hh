/*
 * AllpixPixelsParameterization.hh
 *
 *  Created on: 8 April 2014
 *      Author: nalipour
 */

#ifndef ALLPIXPIXELSPARAMETERIZATION_HH_
#define ALLPIXPIXELSPARAMETERIZATION_HH_

#include "G4VPVParameterisation.hh"
#include "G4ThreeVector.hh"
#include "G4VPhysicalVolume.hh"
#include "AllPixGeoDsc.hh"

class Allpix_PixelsParameterization : public G4VPVParameterisation {
public:
	Allpix_PixelsParameterization(AllPixGeoDsc * geo);
	void ComputeTransformation(G4int, G4VPhysicalVolume*) const;
	double posX(int id) const;
	double posY(int id) const;

 private:

	G4double hsensorX;
	G4double hsensorY;
	G4double hpixelX;
	G4double hpixelY;

	G4int npixelX;
	G4int npixelY;

	AllPixGeoDsc * fGeoPars;
};

#endif /* ALLPIXPIXELSPARAMETERIZATION_HH_ */
