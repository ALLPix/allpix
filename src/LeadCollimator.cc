/**
 *  John Idarraga <idarraga@cern.ch>, 2015
 *
 */

#include "AllPixDetectorConstruction.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Cons.hh"
#include "G4Trd.hh"
#include "G4Sphere.hh"
#include "G4Polyhedra.hh"
#include "G4SubtractionSolid.hh"
#include "G4NistManager.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

using namespace CLHEP;

G4LogicalVolume * PlaceAndHcc(G4Polyhedra * hcc, G4ThreeVector hccPos, G4LogicalVolume * mother, G4Material * Air);

G4LogicalVolume * AllPixDetectorConstruction::BuildCollimator(G4LogicalVolume * motherL, G4ThreeVector collPos) {

	// Material definition

	G4NistManager* nistManager = G4NistManager::Instance();
	//nistManager->ListMaterials("all");
	G4Material * Pb = nistManager->FindOrBuildMaterial("G4_Pb");

	// Print materials
	//G4cout << *(G4Material::GetMaterialTable()) << G4endl;


	double hc_rOuter_0 = 2*mm;
	double hc_rOuter_1 = hc_rOuter_0 * 0.5; // TODO ... calculate according to focus
	double hc_sep = 150*um;

	// Collimator Pb body
	G4Cons * collimatorBody = new G4Cons("CollimatorBody", 0, 3*1*cm, 0, 3*0.8*cm, 2*mm, pi, 2*pi);
	// Honeycomb core
	const G4double z[] = { 0, 4*cm };
	const G4double rInner[] = { 0, 0 };
	const G4double rOuter[] = { hc_rOuter_0, hc_rOuter_1 };
	G4Polyhedra * hcc = new G4Polyhedra(
			"oneTube",
			0.,		// initial phi starting angle
			2*pi,	// total phi angle
			6,		// number sides
			2,		// number of z planes
			z,
			rInner,
			rOuter
	);

	// Substract the honey comb structure
	// Loop in y
	G4VSolid * substractColl = 0x0;
	G4ThreeVector orgPos(0,0, -2*cm);
	G4ThreeVector runPos = orgPos;
	double inStep = 2*hc_rOuter_0 + hc_sep;
	double posx = 0., posy = 0.;

	for ( int i = 0 ; i <= 0 ; i++ ) { // loop in x

		for ( int j = -1 ; j <= -1 ; j++ ) { // loop in y

			posy = (double)j * inStep + ( (double)i*inStep/2. );
			posx = i * inStep;
			runPos.setY( posy );
			runPos.setX( posx );

			G4cout << "[COLL] Building collimator element " << i << "," << j << G4endl;

			// The first one
			if ( ! substractColl ) {
				substractColl = new G4SubtractionSolid("Coll", collimatorBody, hcc, 0x0, runPos );
			} else {
				substractColl = new G4SubtractionSolid("Coll", substractColl, hcc, 0x0, runPos );
			}
		}

	}


	G4LogicalVolume * CollimatorLogic = new G4LogicalVolume(
			substractColl,
			Pb,
			"CollimatorBodyLogic");

	new G4PVPlacement(
			0,                // no rotation
			collPos,
			CollimatorLogic, // its logical volume
			"CollimatorBody", // its name
			motherL,          // its mother  volume
			true,            // no boolean operations
			0,                // copy number
			true); 			  // checking overlaps



	vector<G4LogicalVolume *> hcc_logic_v;



	/*

	// Build first
	G4ThreeVector hccCenterPos(0, 0, -2*cm);
	G4ThreeVector hccPos(0, 0, -2*cm);
	hcc_logic_v.push_back( PlaceAndHcc(hcc, hccPos) );

	// Build sides

	hccPos.setY( 2*hc_rOuter_0 + hc_sep );
	hcc_logic_v.push_back( PlaceAndHcc(hcc, hccPos) );
	hccPos.setY( -1 * hccPos.y() );
	hcc_logic_v.push_back( PlaceAndHcc(hcc, hccPos) );

	// up
	hccPos = hccCenterPos;
	hccPos.setY( (2*hc_rOuter_0 + hc_sep) / 2. );
	hccPos.setX( 2*hc_rOuter_0 + hc_sep );
	hcc_logic_v.push_back( PlaceAndHcc(hcc, hccPos) );
	hccPos.setY( -1 * hccPos.y() );
	hcc_logic_v.push_back( PlaceAndHcc(hcc, hccPos) );

	// down
	hccPos = hccCenterPos;
	hccPos.setY( (2*hc_rOuter_0 + hc_sep) / 2. );
	hccPos.setX( -1 * (2*hc_rOuter_0 + hc_sep) );
	hcc_logic_v.push_back( PlaceAndHcc(hcc, hccPos) );
	hccPos.setY( -1 * hccPos.y() );
	hcc_logic_v.push_back( PlaceAndHcc(hcc, hccPos) );




	G4VisAttributes* hccVisAtt = new G4VisAttributes( G4Colour(0.1, 0.1, 0.9, 0.5) );
	hccVisAtt->SetForceSolid( true );
	//hccVisAtt->SetForceWireframe( true );
	hccVisAtt->SetForceAuxEdgeVisible( true );


	vector<G4LogicalVolume *>::iterator hccItr = hcc_logic_v.begin();
	vector<G4LogicalVolume *>::iterator hccItrE = hcc_logic_v.end();

	for ( ; hccItr != hccItrE ; hccItr++ ) {
		(*hccItr)->SetVisAttributes( hccVisAtt );
	}

	 */

	// Visualization attributes
	G4VisAttributes* CollVisAtt = new G4VisAttributes( G4Colour(1,1,1, 1) );
	CollVisAtt->SetForceSolid( true );
	//CollVisAtt->SetForceWireframe( true );
	CollVisAtt->SetForceAuxEdgeVisible( true );
	CollimatorLogic->SetVisAttributes( CollVisAtt );

	return CollimatorLogic;

}

G4LogicalVolume * PlaceAndHcc(G4Polyhedra * hcc, G4ThreeVector hccPos, G4LogicalVolume * mother, G4Material * Air) {

	G4LogicalVolume * hccLogic = new G4LogicalVolume(
			hcc,
			Air,
			"hccLogic");

	new G4PVPlacement(
			0,                // no rotation
			hccPos,
			hccLogic, 	      // its logical volume
			"hcc",            // its name
			mother, // its mother  volume
			false,            // no boolean operations
			0,                // copy number
			true); 			  // checking overlaps

	return hccLogic;
}
