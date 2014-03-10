/**
 * Author: John Idarraga <idarraga@cern.ch> , 2010
 *
 */

#include "AllPixDetectorConstruction.hh"

#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4Element.hh"
#include "G4ElementTable.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Orb.hh"
#include "G4LogicalVolume.hh"
#include "G4RotationMatrix.hh"
#include "G4Transform3D.hh"
#include "G4PVPlacement.hh"
#include "G4OpBoundaryProcess.hh"
#include "G4VPhysicalVolume.hh"

#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4ThreeVector.hh"

#include "G4GDMLParser.hh"

#include "G4NistManager.hh"

#include "TString.h"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void AllPixDetectorConstruction::BuildTestStructure(int){

	switch (m_TestStructure_type) {
	case 0:
	{
		//Beampipe

		G4NistManager* nistman = G4NistManager::Instance();
		G4Material * Berylium = nistman->FindOrBuildMaterial("G4_Be");
		//G4cout << Berylium << G4endl;
		G4Tubs * beampipe = new G4Tubs("beampipe",
				//3.0*mm, //innerRadiusOfTheTube,
				27.4*mm, //innerRadiusOfTheTube,
				28.0*mm, //outerRadiusOfTheTube,
				130*mm, //hightOfTheTube,
				0.*deg, //startAngleOfTheTube,
				360.*deg); //spanningAngleOfTheTube)

		m_TestStructure_log = new G4LogicalVolume(beampipe,
				//Bone,
				Berylium,
				"Berylium",
				0,0,0);

		//Vacuum

		G4double pressure    = 1e-6*bar;
		G4double temperature = STP_Temperature;                      //from PhysicalConstants.h

		// PV=nRT -> P= (n/V) *RT -> n/V=P/RT
		G4double density     = ((pressure/(hep_pascal))/(286.9*temperature/kelvin))*(kg/m3);
		G4Material* beam = new G4Material("Beam ", density, 2,kStateGas,temperature,pressure);

		//G4cout << "[BeamPipe] Beampipe material with " << TString::Format("P=%3.5f bar rho= %3.3f mg/cm3",pressure/bar,density/(mg/cm3)) << endl;

		// air Material
		density = 1.290*mg/cm3;
		G4Material* Air = new G4Material("Air", density,2);

		//elements
		G4double a = 14.01*g/mole;
		G4Element* elN  = new G4Element("Nitrogen","N" , 7., a);
		a = 16.00*g/mole;
		G4Element* elO  = new G4Element("Oxygen"  ,"O" , 8., a);

		Air->AddElement(elN, 0.7);
		Air->AddElement(elO, 0.3);

		beam->AddElement(elN, 0.7);
		beam->AddElement(elO, 0.3);

		G4Tubs * innerbeampipe = new G4Tubs("innerbeampipe",
				//3.0*mm, //innerRadiusOfTheTube,
				0.0*mm, //innerRadiusOfTheTube,
				27.4*mm, //outerRadiusOfTheTube,
				150*mm, //heightOfTheTube,
				0.*deg, //startAngleOfTheTube,
				360.*deg); //spanningAngleOfTheTube)

		m_TestStructure_log2 = new G4LogicalVolume(innerbeampipe,
				//Bone,
				beam,
				"vacuum",
				0,0,0);

		//beampipe
		G4VisAttributes * visAtt_bp = new G4VisAttributes(G4Color(1,0,1,1));
		visAtt_bp->SetLineWidth(1);
		visAtt_bp->SetForceSolid(true);
		m_TestStructure_log->SetVisAttributes(visAtt_bp);
		G4RotationMatrix* matrix = new G4RotationMatrix();
		matrix->rotateX(0.*deg);

		m_TestStructure_phys = new G4PVPlacement(matrix,
				G4ThreeVector(0,0,0),
				m_TestStructure_log,
				"beampipe",
				expHall_log,
				false,
				0);
		//G4cout << "Beampipe !" << G4endl	;

		// Vacuum
		G4VisAttributes * visAtt_bp2 = new G4VisAttributes(G4Color::Black());
		visAtt_bp2->SetLineWidth(1);
		visAtt_bp2->SetForceSolid(true);
		m_TestStructure_log2->SetVisAttributes(G4VisAttributes::Invisible);
		G4RotationMatrix* matrix2 = new G4RotationMatrix();
		matrix2->rotateX(0.*deg);

		m_TestStructure_phys2 = new G4PVPlacement(matrix2,
				G4ThreeVector(0,0,0),
				m_TestStructure_log2,
				"vacuum",
				expHall_log,
				false,
				0);
		//G4cout << "Vacuum !" << G4endl;
		break;
	}
	case 1:
	{
#ifdef _EUTELESCOPE
		////////////////////////////////////////////////////////////////////////
		// Get GDML volume

		// Parsing
		G4GDMLParser parser;
		G4LogicalVolume * eud_log;

		//parser.Read("models/test_beam_telescope.gdml");
		parser.Read("share/GDML_EUDETAlHolder/EUDETAlHolder.gdml");

		//eud_log = parser.GetVolume("EUD0");
		eud_log = parser.GetVolume("Structure_105681480");

		G4VisAttributes * invisibleVisAtt = new G4VisAttributes(G4Color(1.0, 0.65, 0.0, 0.1));
		invisibleVisAtt->SetVisibility(false);

		int nDaugh = eud_log->GetNoDaughters();
		G4cout << "    Read volume has " << nDaugh << " daughters" << G4endl;

		G4VisAttributes * AlFoilVisAtt = new G4VisAttributes(G4Color(1, 1, 1, 1.0));
		AlFoilVisAtt->SetLineWidth(1);
		AlFoilVisAtt->SetVisibility(false);
		AlFoilVisAtt->SetForceSolid(true);
		eud_log->SetVisAttributes(AlFoilVisAtt);

		for(int dItr = 0 ; dItr < nDaugh ; dItr++){
			G4VPhysicalVolume * tempPhys = eud_log->GetDaughter(dItr);
			G4cout << tempPhys->GetName() << G4endl;
		}

		// I will extract the wrapper from this volume
		G4VSolid * eudAlPlaneSolid = eud_log->GetSolid();
		G4VSolid * wrapperSolid;
		G4SubtractionSolid * subtractionSolid;

		TString physName = "";

		map<int, G4ThreeVector>::iterator testStructItr = m_posVectorTestStructure.begin();

		for( ; testStructItr != m_posVectorTestStructure.end() ; testStructItr++){

			G4ThreeVector posRel = m_posVector[m_detectorLinkTestStructure[(*testStructItr).first]] - (*testStructItr).second;
			G4RotationMatrix * rotRel = m_rotVector[m_detectorLinkTestStructure[(*testStructItr).first]];
			wrapperSolid = m_wrapper_log[m_detectorLinkTestStructure[(*testStructItr).first]]->GetSolid(); // match to an specific detector
			subtractionSolid = new G4SubtractionSolid("EUDETAlPlane",
					eudAlPlaneSolid,
					wrapperSolid,
					rotRel,
					posRel); // match to an specific detector
			eud_log->SetSolid(subtractionSolid);

			physName = "test_phys_";
			physName += (*testStructItr).first;
			new G4PVPlacement(
					m_rotVectorTestStructure[(*testStructItr).first],
					(*testStructItr).second,
					eud_log,                // Logical volume
					physName.Data(),        // Name
					expHall_log,            // Mother volume logical
					false,                  // Unused boolean
					0,                      // copy number
					true);                  // overlap

			//eud_log->SetVisAttributes(AlFoilVisAtt);
		}
		
		
#endif
		break;
	}
	case 2:
	{
		//Timepix Telescope box

		//G4cout << "Structure type: " << m_TestStructure_type << G4endl;

		// materials
		G4NistManager* nistman = G4NistManager::Instance();

		G4Material * Alu = nistman->FindOrBuildMaterial("G4_Al");

		G4Box *box1=new G4Box("BoxTT1",75*mm,150*mm,216*mm);
		G4Box *box2=new G4Box("nBoxTT1",71*mm,146*mm,212*mm);
		G4Box *box3=new G4Box("window",30*mm,100*mm,215.960*mm);

		G4SubtractionSolid *boxTT1tmp = new G4SubtractionSolid("BoxTT1-nBoxTT1",box1,box2);
		G4SubtractionSolid *boxTT1 = new G4SubtractionSolid("BoxTT1tmp-window",boxTT1tmp,box3);

		m_TestStructure_log = new G4LogicalVolume(boxTT1,Alu,"box1_log",0,0,0);

		G4VisAttributes * visAtt_bp = new G4VisAttributes(G4Color(0.5, 0.5, 0.5,0.05));
		visAtt_bp->SetLineWidth(1);
		visAtt_bp->SetForceSolid(true);

		m_TestStructure_log->SetVisAttributes(visAtt_bp);

		m_TestStructure_phys = new G4PVPlacement(0,
				G4ThreeVector(0,5*cm,21.6*cm),
				m_TestStructure_log,
				"box1_phys",
				expHall_log,
				false,
				0);

		m_TestStructure_phys2 = new G4PVPlacement(0,
				G4ThreeVector(0,5*cm,82.7*cm),
				m_TestStructure_log,
				"box2_phys",
				expHall_log,
				false,
				1);

		///////////////////////////////////////
		// DUT Aluminum box

		// Get GDML volume
		G4GDMLParser parser;
		G4LogicalVolume * world_log;

		// Read structure
		parser.Read("models/clicpix_box.gdml");
		world_log = parser.GetVolume("Structure_144055344");

		G4VisAttributes * visAtt = new G4VisAttributes(G4Color(0.5, 0.5, 0.5,0.25));
		visAtt->SetLineWidth(1);
		visAtt->SetForceSolid(true);

		// Aluminum sheet for windows
		G4Box * aluSheet = new G4Box("aluSheet",199*mm/2,40*mm/2,20*um/2);
		G4LogicalVolume * aluSheet_log = new G4LogicalVolume(aluSheet,Alu,"aluSheet_log",0,0,0);
		aluSheet_log->SetVisAttributes(visAtt);

		// DUT Wrapper ThreeVector (detID==504)
		G4ThreeVector wrapper_pos = m_wrapper_phys[504]->GetObjectTranslation();
		//G4RotationMatrix * wrapper_rot = m_wrapper_phys[504]->GetObjectRotation();

		// Loop on daughter volumes i.e. various parts of the Alu box
		int nDaugh = world_log->GetNoDaughters();
		G4cout << "Read volume has " << nDaugh << " daughters" << G4endl;
		for(int dItr = 0 ; dItr < nDaugh ; dItr++)
		{
			G4VPhysicalVolume * temp_phys = world_log->GetDaughter(dItr);
			string name = temp_phys->GetName();
			G4ThreeVector pos = temp_phys->GetObjectTranslation();   // Volume position as in GDML
			G4RotationMatrix * rot = temp_phys->GetObjectRotation(); // Volume rotation as in GDML

			// Move the box down to align window with DUT based on macro value
			map<int, G4ThreeVector>::iterator testStructItr = m_posVectorTestStructure.begin();
			for( ; testStructItr != m_posVectorTestStructure.end() ; testStructItr++)
			{
				G4ThreeVector posRel = (*testStructItr).second;
				pos.setY(pos.getY()+posRel.getY());
			}

			pos.setZ(pos.getZ()+wrapper_pos.getZ()); // Translate w.r.t. wrapper Z position

			G4LogicalVolume * temp_log = temp_phys->GetLogicalVolume();
			new G4PVPlacement(
					rot,
					pos,                    // Position of physical volume
					temp_log,               // Logical volume
					name,                   // Name
					expHall_log,            // Mother volume logical
					false,                  // Unused boolean
					0,                      // copy number
					true);                  // overlap
			temp_log->SetVisAttributes(visAtt);

			// Alu sheets
			if (dItr == 5 || dItr == 11) // Side plates
			{
				pos.setY(pos.getY()-65);
				if (dItr == 5)  pos.setZ(pos.getZ()-1*mm); // move it inside
				if (dItr == 11) pos.setZ(pos.getZ()+1*mm); // move it inside

				new G4PVPlacement(
						rot,
						pos,
						aluSheet_log,
						"aluSheet_phys",
						expHall_log,
						false,
						1,
						true);
			}

		}
		
		break;
	}
	case 3:
	{
		//Am241 gamma Calibration source

		G4NistManager* nistman = G4NistManager::Instance();
		//G4Material * Alu = nistman->FindOrBuildMaterial("G4_PLEXIGLASS");
		G4Material * Pb = nistman->FindOrBuildMaterial("G4_Pb");

		// Define elements from NIST
		nistman->FindOrBuildElement("H");
		nistman->FindOrBuildElement("Be");
		G4Element* C  = nistman->FindOrBuildElement("C");
		nistman->FindOrBuildElement("N");
		nistman->FindOrBuildElement("O");
		nistman->FindOrBuildElement("Al");
		G4Element* Si = nistman->FindOrBuildElement("Si");
		nistman->FindOrBuildElement("Ti");
		G4Element* Cr = nistman->FindOrBuildElement("Cr");
		G4Element* Mn = nistman->FindOrBuildElement("Mn");
		G4Element* Fe = nistman->FindOrBuildElement("Fe");
		G4Element* Ni = nistman->FindOrBuildElement("Ni");
		nistman->FindOrBuildElement("W");
		nistman->FindOrBuildElement("Au");
		nistman->FindOrBuildElement("Pb");
		// Define pure NIST materials
		nistman->FindOrBuildMaterial("G4_Al");
		nistman->FindOrBuildMaterial("G4_Ti");
		nistman->FindOrBuildMaterial("G4_W");
		nistman->FindOrBuildMaterial("G4_Au");
		// Define other NIST materials
		nistman->FindOrBuildMaterial("G4_WATER");
		nistman->FindOrBuildMaterial("G4_KAPTON");
		//G4Material* Air = nistman->FindOrBuildMaterial("G4_AIR");
		// Define materials not in NIST
		G4double density;
		G4int ncomponents;
		G4double fractionmass;
		G4Material* StainlessSteel = new G4Material("StainlessSteel", density= 8.06*g/cm3, ncomponents=6);
		StainlessSteel->AddElement(C, fractionmass=0.001);
		StainlessSteel->AddElement(Si, fractionmass=0.007);
		StainlessSteel->AddElement(Cr, fractionmass=0.18);
		StainlessSteel->AddElement(Mn, fractionmass=0.01);
		StainlessSteel->AddElement(Fe, fractionmass=0.712);
		StainlessSteel->AddElement(Ni, fractionmass=0.09);

		G4Element * Cu = nistman->FindOrBuildElement("Cu");
		G4Element * Zn = nistman->FindOrBuildElement("Zn");
		G4Material * Brass = new G4Material("Brass", density=8.5*g/cm3, ncomponents=2);
		Brass->AddElement(Cu, fractionmass=0.7);
		Brass->AddElement(Zn, fractionmass=0.3);

		G4Tubs *HolderPrim = new G4Tubs("supportTop",
				0.0*mm, //innerRadiusOfTheTube,
				29.0*mm, //outerRadiusOfTheTube,
				27*mm/2, //heightOfTheTube,
				0.*deg, //startAngleOfTheTube,
				360.*deg); //spanningAngleOfTheTube)

		G4Tubs *hole = new G4Tubs("hole",
				0.0*mm, //innerRadiusOfTheTube,
				2.5*mm, //outerRadiusOfTheTube,
				10*mm/2, //heightOfTheTube,
				0.*deg, //startAngleOfTheTube,
				360.*deg); //spanningAngleOfTheTube)

		G4SubtractionSolid *holder = new G4SubtractionSolid("HolderPrim-hole",HolderPrim,hole,0,G4ThreeVector(0,-22.5*mm,-13*mm));

		G4Box *box1=new G4Box("BoxTT1",50*mm,100*mm,25*mm);

		m_TestStructure_log = new G4LogicalVolume(holder,
				Brass,
				"box1_log",
				0,0,0);

		m_TestStructure_log2 = new G4LogicalVolume(box1,
				Pb,
				"box2_log",
				0,0,0);

		G4VisAttributes * visAtt_bp = new G4VisAttributes(G4Color(0.5, 0.5, 0.5,0.25));
		visAtt_bp->SetLineWidth(1);
		visAtt_bp->SetForceSolid(true);
		m_TestStructure_log->SetVisAttributes(visAtt_bp);

		G4VisAttributes * visAtt_bp2 = new G4VisAttributes(G4Color(0.5, 0.1, 0.3,0.25));
		visAtt_bp2->SetLineWidth(1);
		visAtt_bp2->SetForceSolid(true);

		m_TestStructure_log2->SetVisAttributes(visAtt_bp2);

		m_TestStructure_phys = new G4PVPlacement(0,
				G4ThreeVector(0,22.5*mm,31*mm),
				m_TestStructure_log,
				"box1_phys",
				expHall_log,
				false,
				0);

		m_TestStructure_phys2 = new G4PVPlacement(0,
				G4ThreeVector(0,0*mm,-30.0*mm),
				//G4ThreeVector(0,-22.5*mm,14.0*mm),
				m_TestStructure_log2,
				"box2_phys",
				expHall_log,
				false,
				0);
		break;
	}
	case 4:
	{
		//Fe55 Calibration source

		G4NistManager* nistman = G4NistManager::Instance();
		G4Material * Plexi = nistman->FindOrBuildMaterial("G4_PLEXIGLASS");
		G4Material * Pb = nistman->FindOrBuildMaterial("G4_Pb");

		G4Tubs *HolderPrim = new G4Tubs("supportTop",
				0.0*mm, //innerRadiusOfTheTube,
				29.0*mm, //outerRadiusOfTheTube,
				27*mm/2, //heightOfTheTube,
				0.*deg, //startAngleOfTheTube,
				360.*deg); //spanningAngleOfTheTube)

		G4Tubs *hole = new G4Tubs("hole",
				0.0*mm, //innerRadiusOfTheTube,
				5.0*mm, //outerRadiusOfTheTube,
				10*mm/2, //heightOfTheTube,
				0.*deg, //startAngleOfTheTube,
				360.*deg); //spanningAngleOfTheTube)

		G4SubtractionSolid *holder = new G4SubtractionSolid("HolderPrim-hole",HolderPrim,hole,0,G4ThreeVector(0,-22.5*mm,-10*mm));

		G4Box *box1=new G4Box("BoxTT1",50*mm,100*mm,25*mm);

		m_TestStructure_log = new G4LogicalVolume(holder,
				Plexi,
				"box1_log",
				0,0,0);

		m_TestStructure_log2 = new G4LogicalVolume(box1,
				Pb,
				"box2_log",
				0,0,0);

		G4VisAttributes * visAtt_bp = new G4VisAttributes(G4Color(0.5, 0.5, 0.5,0.25));
		visAtt_bp->SetLineWidth(1);
		visAtt_bp->SetForceSolid(true);
		m_TestStructure_log->SetVisAttributes(visAtt_bp);

		G4VisAttributes * visAtt_bp2 = new G4VisAttributes(G4Color(0.5, 0.1, 0.3,0.25));
		visAtt_bp2->SetLineWidth(1);
		visAtt_bp2->SetForceSolid(true);

		m_TestStructure_log2->SetVisAttributes(visAtt_bp2);

		m_TestStructure_phys = new G4PVPlacement(0,
				G4ThreeVector(0,0*mm,31*mm),
				m_TestStructure_log,
				"box1_phys",
				expHall_log,
				false,
				0);

		m_TestStructure_phys2 = new G4PVPlacement(0,
				G4ThreeVector(0,22.5*mm,-30.0*mm),
				//G4ThreeVector(0,-22.5*mm,14.0*mm),
				m_TestStructure_log2,
				"box2_phys",
				expHall_log,
				false,
				0);
		break;
	}
	case 5:
	{
		//Cd109 Calibration source

		G4NistManager* nistman = G4NistManager::Instance();
		G4Material * Alu = nistman->FindOrBuildMaterial("G4_Al");
		G4Material * Pb = nistman->FindOrBuildMaterial("G4_Pb");

		G4Tubs *HolderPrim = new G4Tubs("supportTop",
				0.0*mm, //innerRadiusOfTheTube,
				29.0*mm, //outerRadiusOfTheTube,
				27*mm/2, //heightOfTheTube,
				0.*deg, //startAngleOfTheTube,
				360.*deg); //spanningAngleOfTheTube)

		G4Tubs *hole = new G4Tubs("hole",
				0.0*mm, //innerRadiusOfTheTube,
				2.5*mm, //outerRadiusOfTheTube,
				10*mm/2, //heightOfTheTube,
				0.*deg, //startAngleOfTheTube,
				360.*deg); //spanningAngleOfTheTube)

		G4SubtractionSolid *holder = new G4SubtractionSolid("HolderPrim-hole",HolderPrim,hole,0,G4ThreeVector(0,-22.5*mm,-13*mm));

		G4Box *box1=new G4Box("BoxTT1",50*mm,100*mm,25*mm);

		m_TestStructure_log = new G4LogicalVolume(holder,
				Alu,
				"box1_log",
				0,0,0);

		m_TestStructure_log2 = new G4LogicalVolume(box1,
				Pb,
				"box2_log",
				0,0,0);

		G4VisAttributes * visAtt_bp = new G4VisAttributes(G4Color(0.5, 0.5, 0.5,0.25));
		visAtt_bp->SetLineWidth(1);
		visAtt_bp->SetForceSolid(true);
		m_TestStructure_log->SetVisAttributes(visAtt_bp);

		G4VisAttributes * visAtt_bp2 = new G4VisAttributes(G4Color(0.5, 0.1, 0.3,0.25));
		visAtt_bp2->SetLineWidth(1);
		visAtt_bp2->SetForceSolid(true);

		m_TestStructure_log2->SetVisAttributes(visAtt_bp2);

		m_TestStructure_phys = new G4PVPlacement(0,
				G4ThreeVector(0,0*mm,31*mm),
				m_TestStructure_log,
				"box1_phys",
				expHall_log,
				false,
				0);

		m_TestStructure_phys2 = new G4PVPlacement(0,
				G4ThreeVector(0,0*mm,-30.0*mm),
				m_TestStructure_log2,
				"box2_phys",
				expHall_log,
				false,
				0);
		break;
	}
	case 6:
	{
		//Am241 alpha Calibration source

		G4NistManager* nistman = G4NistManager::Instance();
		G4Material * Pb = nistman->FindOrBuildMaterial("G4_Pb");

		// Define elements from NIST
		G4Element* C  = nistman->FindOrBuildElement("C");
		G4Element* Si = nistman->FindOrBuildElement("Si");
		G4Element* Cr = nistman->FindOrBuildElement("Cr");
		G4Element* Mn = nistman->FindOrBuildElement("Mn");
		G4Element* Fe = nistman->FindOrBuildElement("Fe");
		G4Element* Ni = nistman->FindOrBuildElement("Ni");
		// Define materials not in NIST
		G4double density;
		G4int ncomponents;
		G4double fractionmass;
		G4Material* StainlessSteel = new G4Material("StainlessSteel", density= 8.06*g/cm3, ncomponents=6);
		StainlessSteel->AddElement(C, fractionmass=0.001);
		StainlessSteel->AddElement(Si, fractionmass=0.007);
		StainlessSteel->AddElement(Cr, fractionmass=0.18);
		StainlessSteel->AddElement(Mn, fractionmass=0.01);
		StainlessSteel->AddElement(Fe, fractionmass=0.712);
		StainlessSteel->AddElement(Ni, fractionmass=0.09);

		G4Tubs *disk = new G4Tubs("disk",
				0.0*mm, //innerRadiusOfTheTube,
				11.0*mm, //outerRadiusOfTheTube,
				5*mm/2, //heightOfTheTube,
				0.*deg, //startAngleOfTheTube,
				360.*deg); //spanningAngleOfTheTube)

		G4Box *box1=new G4Box("BoxTT1",50*mm,100*mm,25*mm);

		m_TestStructure_log = new G4LogicalVolume(disk,
				StainlessSteel,
				"box1_log",
				0,0,0);

		m_TestStructure_log2 = new G4LogicalVolume(box1,
				Pb,
				"box2_log",
				0,0,0);

		G4VisAttributes * visAtt_bp = new G4VisAttributes(G4Color(0.5, 0.5, 0.5,0.25));
		visAtt_bp->SetLineWidth(1);
		visAtt_bp->SetForceSolid(true);
		m_TestStructure_log->SetVisAttributes(visAtt_bp);

		G4VisAttributes * visAtt_bp2 = new G4VisAttributes(G4Color(0.5, 0.1, 0.3,0.25));
		visAtt_bp2->SetLineWidth(1);
		visAtt_bp2->SetForceSolid(true);

		m_TestStructure_log2->SetVisAttributes(visAtt_bp2);

		m_TestStructure_phys = new G4PVPlacement(0,
				G4ThreeVector(0,0,24.5*mm),
				m_TestStructure_log,
				"box1_phys",
				expHall_log,
				false,
				0);

		m_TestStructure_phys2 = new G4PVPlacement(0,
				G4ThreeVector(0,0*mm,-30.0*mm),
				//G4ThreeVector(0,-22.5*mm,14.0*mm),
				m_TestStructure_log2,
				"box2_phys",
				expHall_log,
				false,
				0);
		break;
	}
	case 7:
	{
		// Test GDML

		// Get GDML volume
		G4GDMLParser parser;
		G4LogicalVolume * world_log;

		// Read structure
		parser.Read("models/clicpix_box.gdml");
		world_log = parser.GetVolume("Structure_11624736");

		G4VisAttributes * visAtt = new G4VisAttributes(G4Color(0.5, 0.5, 0.5,0.25));
		visAtt->SetLineWidth(1);
		visAtt->SetForceSolid(true);

		// Loop on daughter volumes i.e. various parts of the Alu box
		int nDaugh = world_log->GetNoDaughters();
		G4cout << "Read volume has " << nDaugh << " daughters" << G4endl;
		for(int dItr = 0 ; dItr < nDaugh ; dItr++){
			G4VPhysicalVolume * temp_phys = world_log->GetDaughter(dItr);
			string name = temp_phys->GetName();
			G4LogicalVolume * temp_log = temp_phys->GetLogicalVolume();
			new G4PVPlacement(0,
					  temp_phys->GetObjectTranslation(), // Position of physical volume
					  temp_log,               // Logical volume
					  name,                   // Name
					  expHall_log,            // Mother volume logical
					  false,                  // Unused boolean
					  0,                      // copy number
					  true);                  // overlap
			temp_log->SetVisAttributes(visAtt);
		}

		break;
	}

	case 8:
	{

		G4NistManager * nistman = G4NistManager::Instance();
		G4Material * mylar = nistman->FindOrBuildMaterial("G4_MYLAR");
		TString physName = "";
		TString fr_name = "" ;
		TString bk_name = "" ;
		TString fr_name_log = "" ;
		TString bk_name_log = "" ;



		G4VisAttributes * scintAtt = new G4VisAttributes(G4Color(1,1,1,0.5));
		scintAtt->SetLineWidth(1);
		scintAtt->SetForceSolid(true);



			map<int, G4ThreeVector>::iterator testStructItr = m_posVectorTestStructure.begin();

			for( ; testStructItr != m_posVectorTestStructure.end() ; testStructItr++){

				G4ThreeVector posRel = (*testStructItr).second;

				fr_name = "fr_box_";
				fr_name += (*testStructItr).first;

				G4Box* FrontWindow = new G4Box(fr_name.Data(),
						22.5*mm,
						32.5*mm,
						25*um);

				bk_name = "bk_box_";
				bk_name += (*testStructItr).first;

				G4Box* BackWindow = new G4Box(bk_name.Data(),
						22.5*mm,
						32.5*mm,
						25*um);



				// logical
				fr_name_log = "fr_log_";
				fr_name_log += (*testStructItr).first;

				G4LogicalVolume * window_log_fr = new G4LogicalVolume(
						FrontWindow,
						mylar,
						fr_name_log.Data());
				window_log_fr->SetVisAttributes(scintAtt);


				bk_name_log = "fr_log_";
				bk_name_log += (*testStructItr).first;

				G4LogicalVolume * window_log_bk = new G4LogicalVolume(
						BackWindow,
						mylar,
						bk_name_log.Data());
				window_log_bk->SetVisAttributes(scintAtt);


				//Phyisical

				physName = "test_fr_phys_";
				physName += (*testStructItr).first;

				new G4PVPlacement( 0,
						G4ThreeVector(posRel[0],posRel[1],posRel[2]-3.0125*mm),
						window_log_fr,
						physName.Data(),
						expHall_log,
						false,
						0,
						true);


				physName = "test_bk_phys_";
				physName += (*testStructItr).first;

				new G4PVPlacement( 0,
						G4ThreeVector(posRel[0],posRel[1],posRel[2]+12.0125*mm),
						window_log_bk,
						physName.Data(),
						expHall_log,
						false,
						0,
						true);
			}


		break;
	}




	default:
	{
		G4cout << "Unknown TestStructure Type" << G4endl;
		break;
	}
	}
	
#ifdef _EUTELESCOPE
	
		//////////////////////////////////////////////////////////
		// Scintillators for EUDET
		// Materials

		if( m_scintPos.empty() ) {
			G4cout << "[ERROR] no scintillators defined.  In the macro use the command"
					<< "       /allpix/eudet/scint1Pos 0.0  0.0  -24.0 mm"
					<< "       Can't recover ... giving up."
					<< G4endl;
			exit(1);
		}else{
			G4cout << "Building scintillators..." << G4endl;
		}

		// first plane at 0. mm
		//	G4double z1 = -24*mm;
		//	G4double z2 = -18*mm;
		// last plane at 490 mm
		//	G4double z3 = 523*mm;
		// G4double z4 = 529*mm;

		G4NistManager * nistman = G4NistManager::Instance();
		// Scintillator
		G4Material * scplastic = nistman->FindOrBuildMaterial("G4_POLYSTYRENE");
		// Scintillators
		G4Box* scintb = new G4Box("scintb",
				11.0*mm,
				5.4*mm,
				3*mm); // scintillators 6mm thick

		G4VisAttributes * scintAtt = new G4VisAttributes(G4Color(1,0,1,1));
		scintAtt->SetLineWidth(1);
		scintAtt->SetForceSolid(true);

		// Place scintillators
		vector<G4ThreeVector>::iterator scintItr = m_scintPos.begin();

		TString labelLog = "";
		TString labelPlacement = "";
		TString labelSD = "";
		Int_t cntr = 1;
		G4SDManager * SDman = G4SDManager::GetSDMpointer();

		for( ; scintItr != m_scintPos.end() ; scintItr++) {
			labelLog = "scint";
			labelLog += cntr;
			labelLog += "_log";

			labelPlacement = "Scint";
			labelPlacement += cntr;

			G4LogicalVolume * scint_log = new G4LogicalVolume(
					scintb,
					scplastic,
					labelLog.Data());
			scint_log->SetVisAttributes(scintAtt);
			G4RotationMatrix* matrix_s = new G4RotationMatrix();
			matrix_s->rotateX(0.*deg);

			new G4PVPlacement( matrix_s,
					(*scintItr),
					scint_log,
					labelPlacement.Data(),
					expHall_log,
					false,
					0,
					true);

			labelSD = "sdscint";
			labelSD += cntr;
			AllPixTrackerSD * scintTrack = new AllPixTrackerSD( labelSD.Data(), (*scintItr), 0);
			SDman->AddNewDetector( scintTrack );
			scint_log->SetSensitiveDetector( scintTrack );
			cntr++;
			}
#endif

}
