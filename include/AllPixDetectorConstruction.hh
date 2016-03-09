/**
 * Allpix
 * Author: John Idarraga <idarraga@cern.ch> , 2010
 */

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
// $Id: AllPixDetectorConstruction.hh,v 1.5 2006/06/29 17:53:55 gunter Exp $
// GEANT4 tag $Name: geant4-09-02 $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef AllPixDetectorConstruction_h
#define AllPixDetectorConstruction_h 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4SDManager.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4PVDivision.hh"
#include "G4VSolid.hh"

#include "AllPixTrackerSD.hh"
#include "AllPixBumpsParameterization.hh"
#include "AllPixDetectorMessenger.hh"
#include "ReadGeoDescription.hh"

#include "G4ThreeVector.hh"

#include <vector>
#include <stdio.h>
#include <utility>
#include <string>
#include <map>

using namespace std;

class G4Box;
class G4LogicalVolume;
class G4VPhysicalVolume;
class AllPixGeoDsc;
class G4UniformMagField;
class G4QuadrupoleMagField;
class MorourgoMagField;

class AllPixDetectorConstruction : public G4VUserDetectorConstruction
{

public:

	AllPixDetectorConstruction();
	~AllPixDetectorConstruction();
	void DefineSensitiveDetector();
	void VolumesG4Properties();
	//void BuildMediPix(vector<G4ThreeVector>, vector<G4RotationMatrix *>);
	void BuildPixelDevices(map<int, AllPixGeoDsc *>);

	// check setup
	void CheckAllPixSetup();

	// test structures
	void SetBuildTestStructure(bool flg){m_buildTestStructureFlag = flg;};
	void SetBuildAppliances(bool flg){m_buildAppliancesFlag = flg;};
	G4LogicalVolume * BuildCollimator(G4LogicalVolume * motherL, G4ThreeVector collPos);

	void SetTestStructureType(int type){m_TestStructure_type = type;};
	void SetAppliancesType(int type){m_Appliances_type = type;};

	void SetTestStructurePosition(G4ThreeVector);
	void SetTestStructureRotation(G4ThreeVector);
	void SetTestStructureDetectorLink(G4int);
	void SetAppliancePosition(G4ThreeVector);
	void SetWrapperEnhancement(G4ThreeVector);
	void SetMaxStepLengthSensor(G4double);

	// mag field
	void SetPeakMagField(G4double fieldValue);

	// world volume from macro
	void SetWorldMaterial(G4String);

	// builders
	void BuildTestStructure(G4int);
	void BuildAppliances(G4int);

	void SetDetectorPosition(G4ThreeVector);
	void SetDetectorRotation(G4ThreeVector);
	void SetDetectorID(G4int);

	void SetLowTHL(G4double);
	void UpdateGeometry();

  // others
  void SetOutputFilePrefix(G4String);
  G4String GetOutputFilePrefix(){return m_outputFilePrefix;};

  G4VSolid * GetVSolidDetector (G4int detId) {
	  // Check first if the detId is good
	  // otherwise return null pointer
	  return m_Box_log[detId]->GetSolid();
  };

	// Specific EUTelescope
#ifdef _EUTELESCOPE
	void SetScintPos(G4ThreeVector);
#endif

public:

	G4VPhysicalVolume* Construct();

private:
	// flags
	bool m_clearanceToBuildGeometry;

	// pos rot detector
	vector<G4int>              m_detId;
	vector<G4int>::iterator    m_detIdItr;

	map<int, G4ThreeVector>      m_posVector; // position of medipix(es), key is detector Id
	map<int, G4RotationMatrix *> m_rotVector; // rotation
	G4int m_nPositions;
	G4int m_nRotations;
	G4int m_nIds;

	vector<G4double>           m_lowThlVector; // lowTHL
	// for user information.  Absolute position (center) of the Si wafers
	vector<G4ThreeVector>      m_absolutePosSiWafer;


	// pos,rot test structure
	map<int, G4ThreeVector> m_posVectorTestStructure; // position of test structure
	map<int, G4RotationMatrix *> m_rotVectorTestStructure; // rotation of test structure
	map<G4int, G4int> m_detectorLinkTestStructure;
	G4int m_nTestPositions;
	G4int m_nTestRotation;

	// pos appliance
	//vector<G4ThreeVector> m_posVectorAppliances;
	map<int, G4ThreeVector> m_posVectorAppliances;
	G4int m_nAppliancesPositions;
	//vector<G4ThreeVector> m_vectorWrapperEnhancement;
	map<int, G4ThreeVector> m_vectorWrapperEnhancement;
	G4int m_nWrapperEnhancement;

	// Geometry in G4
	AllPixGeoDsc * gD;   // geo bits for 1 detector
	ReadGeoDescription * m_geoDsc;

	// world
	G4LogicalVolume * expHall_log;
	G4VPhysicalVolume * expHall_phys;

	//Apliance support Timepix telescope
	
	G4LogicalVolume * supportLV;
	G4VPhysicalVolume * supportPV;	
			
	map<int, G4LogicalVolume *>    m_wrapper_log;
	map<int, G4VPhysicalVolume *>  m_wrapper_phys;

	map<int, G4LogicalVolume *>    m_PCB_log;
	map<int, G4VPhysicalVolume*>   m_PCB_phys;

	map<int, G4LogicalVolume *>    m_Box_log;
	map<int, G4VPhysicalVolume*>   m_Box_phys;

    map<int, G4LogicalVolume *>    m_Coverlayer_log;
    map<int, G4VPhysicalVolume*>   m_Coverlayer_phys;

	map<int, G4LogicalVolume *>    m_Chip_log;
	map<int, G4VPhysicalVolume*>   m_Chip_phys;

	map<int, G4LogicalVolume *>    m_Bumps_log;
	map<int, G4VPhysicalVolume*>   m_Bumps_phys;

	map<int, G4LogicalVolume *>    m_Bumps_Slice_log;
	map<int, G4LogicalVolume *>    m_Bumps_Cell_log;


	map<int, G4LogicalVolume *>    m_GuardRings_log;
	map<int, G4VPhysicalVolume*>   m_GuardRings_phys;

	map<int, G4LogicalVolume *>    m_Slice_log;
	map<int, G4LogicalVolume *>    m_Pixel_log;



	Allpix_BumpsParameterization * parameterization;

	/*
  //////////////////////////////////
  // wrapper
  G4LogicalVolume ** m_wrapper_log;
  G4VPhysicalVolume ** m_wrapper_phys;
  // PCB
  G4LogicalVolume ** m_PCB_log;
  G4VPhysicalVolume ** m_PCB_phys;
  // Box
  G4LogicalVolume ** m_Box_log;
  G4VPhysicalVolume ** m_Box_phys;
  G4LogicalVolume ** m_GuardRings_log;
  G4VPhysicalVolume ** m_GuardRings_phys;
  G4int m_detectorId;

  // slice and pixel
  G4LogicalVolume ** m_Slice_log;
  G4LogicalVolume ** m_Pixel_log;
  //////////////////////////////////
	 */

	// test Structure
	bool m_buildTestStructureFlag;
	G4LogicalVolume * m_TestStructure_log;
	G4VPhysicalVolume * m_TestStructure_phys;
	int m_TestStructure_type;
	
	G4LogicalVolume * m_TestStructure_log2;
	G4VPhysicalVolume * m_TestStructure_phys2;
	// appliances
	bool  m_buildAppliancesFlag;
	int m_Appliances_type;

	// messenger
	AllPixDetectorMessenger* m_detectorMessenger;

	// materials
	G4Material * m_Air;
	G4Material * m_Vacuum;
	G4Material * m_fillingWorldMaterial;
	bool m_userDefinedWorldMaterial;

	// mad field
	G4UniformMagField * m_magField;      // pointer to the magnetic field
	//MorourgoMagField * m_magField;

	// user limits
	G4UserLimits * m_ulim;
	G4double m_maxStepLengthSensor;

  // others
  G4String m_outputFilePrefix;

	// Eutelescope specifig
#ifdef _EUTELESCOPE
	vector<G4ThreeVector> m_scintPos;
	vector<G4ThreeVector> m_scintRot;
#endif

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#define _BUILD_MEDIPIX_MSG()						\
		G4cout << "  In order to place a detector you must include in your macro at least" << G4endl; \
		G4cout << "  the following lines: (setId must go first !)" << G4endl;			\
		G4cout << "    /allpix/det/setId 5" << G4endl;	\
		G4cout << "    /allpix/det/setPosition 0.0 0.0 0.0 mm" << G4endl;	\
		G4cout << "    /allpix/det/setRotation 0.0 0.0 0.0 deg" << G4endl;	\
		G4cout << "    /allpix/det/setLowTHL 13. keV" << G4endl;		\
		G4cout << "  see an example in \"allpix_vis.in\"" << G4endl;

#define _WRONG_CONFIG_ABORT_MSG()						\
		G4cout << "[ERROR] wrong configuration.  Aborting job." << G4endl;

#endif /*AllPixDetectorConstruction_h*/
