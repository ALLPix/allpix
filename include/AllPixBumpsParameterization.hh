/*
 * AllpixBumpsParameterization.hh
 *
 *  Created on: 31 janv. 2014
 *      Author: mbenoit
 */

#ifndef ALLPIXBUMPSPARAMETERIZATION_HH_
#define ALLPIXBUMPSPARAMETERIZATION_HH_

#include "G4VPVParameterisation.hh"
#include "G4ThreeVector.hh"
#include "G4VPhysicalVolume.hh"
#include "AllPixGeoDsc.hh"

class Allpix_BumpsParameterization : public G4VPVParameterisation {
public:
	Allpix_BumpsParameterization(AllPixGeoDsc * geo);
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

#endif /* ALLPIXBUMPSPARAMETERIZATION_HH_ */
