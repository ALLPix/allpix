//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// John Idarraga <idarraga@cern.ch>
//
// -------------------------------------------------------------------

#include "MorourgoMagField.hh"
#include "globals.hh"
#include "TF1.h"
#include "TF3.h"
#include "TMath.h"

#define __morourgo_widh_at_half 1000.0 // 1000 mm
#define __sqrt2 1.414213562

TF1 * g_fieldX;
TF1 * g_fieldY;
TF1 * g_fieldZ;

double funcfieldx(double * , double * ){

	return 0;
}


double funcfieldy(double * , double * ){

	return 0;
}
double funcfieldz(double * xx, double * pars){

	double a = pars[0];
	double b = pars[1];
	double c = pars[2];
	double x = xx[0];

	double z = (x-b)/c;
	z = a * TMath::Exp( -0.5*( z*z ) );

	return z;
}


MorourgoMagField::MorourgoMagField(G4float max, G4float center)
{

	m_gradient = max;
	g_fieldZ = new TF1("morourgo", funcfieldz, -5000, 5000, 3);
	// a : max
	// b : center
	// c : width at half
	g_fieldZ->SetParameters(max, center, __morourgo_widh_at_half);

}

/////////////////////////////////////////////////////////////////////////  

MorourgoMagField::~MorourgoMagField()
{

}

/////////////////////////////////////////////////////////////////////////


void MorourgoMagField::GetFieldValue(const G4double yTrack[7],
		G4double B[3]     ) const
{

	/*
	// quadrupole approx
	B[0] = g_fieldZ->Eval(yTrack[2]) * yTrack[1];
	B[1] = g_fieldZ->Eval(yTrack[2]) * yTrack[0];
	B[2] = 0.;
	 */
	// Dipole, flat, no edges approximation
	//std::cout << g_fieldZ->Eval(yTrack[2]) << std::endl;
	B[0] = 0.;
	B[1] = -1. * g_fieldZ->Eval(yTrack[2]); // field goes down
	B[2] = 0.;
}

// -----------------------------------------------------------------
