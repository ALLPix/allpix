/**
 * Author: John Idarraga <idarraga@cern.ch> , 2010
 *
 */
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
// $Id: AllPixDetectorConstruction.cc,v 1.15 2006/06/29 17:54:17 gunter Exp $
// GEANT4 tag $Name: geant4-08-01-patch-01 $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

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
#include "G4Sphere.hh"

#include "G4LogicalVolume.hh"
#include "G4RotationMatrix.hh"
#include "G4Transform3D.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4PVParameterised.hh"
#include "G4OpBoundaryProcess.hh"

#include "G4UserLimits.hh"

#include "G4SubtractionSolid.hh"
#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4Transform3D.hh"
#include "G4UniformMagField.hh"
#include "G4IntersectionSolid.hh"
#include "G4UnionSolid.hh"

#include "G4NistManager.hh"

#include "AllPixGeoDsc.hh"
#include "ReadGeoDescription.hh"

#include "G4DigiManager.hh"
#include "AllPixMimosa26Digitizer.hh"
#include "AllPixFEI3StandardDigitizer.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

AllPixDetectorConstruction::AllPixDetectorConstruction()
{
    //m_detConstruction = dc;
    //m_tgeoManager = gm;


    /*
    map<int, AllPixGeoDsc *>::iterator itr = m_geoMap.begin();
    for( ; itr != m_geoMap.end() ; itr++){
        cout << " --> " << (*itr).first << endl;
        (*itr).second->Dump();
    }
     */
    //gD = new AllPixGeoDsc();
    //gD->ReadGeoBits();

    //m_detectorId = 0;
    //m_detIdItr = m_detId.begin() - 1;

    m_detectorMessenger = new AllPixDetectorMessenger(this);
    m_nIds = 0;
    m_nPositions = 0;
    m_nRotations = 0;
    m_nTestPositions = 0;
    m_nTestRotation = 0;
    m_nAppliancesPositions = 0;
    m_nWrapperEnhancement = 0;

    m_buildAppliancesFlag = false;
    m_Appliances_type = 0;
    m_buildTestStructureFlag = false;
    m_TestStructure_type = 0;
    m_clearanceToBuildGeometry = false;

    m_Air = 0x0;
    m_Vacuum = 0x0;
    m_fillingWorldMaterial = 0x0;
    m_userDefinedWorldMaterial = false;
    gD = 0;
    m_maxStepLengthSensor = 10.0*um;
    m_ulim = 0x0;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

AllPixDetectorConstruction::~AllPixDetectorConstruction(){

    if(m_ulim) delete m_ulim;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void AllPixDetectorConstruction::VolumesG4Properties(){

}

void AllPixDetectorConstruction::DefineSensitiveDetector(){

}

G4VPhysicalVolume * AllPixDetectorConstruction::Construct()
{

    // Clean old geometry, if any
    G4GeometryManager::GetInstance()->OpenGeometry();
    G4PhysicalVolumeStore::GetInstance()->Clean();
    G4LogicalVolumeStore::GetInstance()->Clean();
    G4SolidStore::GetInstance()->Clean();

    // materials
    // Vacuum
    G4double z,a,density;
    m_Vacuum = new G4Material("Vacuum", z=1 , a=1.01*g/mole, density= 0.0001*g/cm3);

    // Air
    G4NistManager * nistman = G4NistManager::Instance();
    m_Air = nistman->FindOrBuildMaterial("G4_AIR");
    nistman->ListMaterials("all");

    // Air is the default.  Can be changed from the messenger using
    // /allpix/extras/setWorldMaterial
    // which calls AllPixDetectorConstruction::SetWorldMaterial(G4String mat)
    if(!m_userDefinedWorldMaterial) m_fillingWorldMaterial = m_Air;

    /////////////////////////////////////////
    // The experimental Hall.  World
    G4VisAttributes * invisibleVisAtt = new G4VisAttributes(G4Color(1.0, 0.65, 0.0, 0.1));
    //G4VisAttributes invisibleVisAtt(G4Color(1.0, 0.65, 0.0, 0.1));
    invisibleVisAtt->SetVisibility(false);
    invisibleVisAtt->SetForceSolid(false);

    G4Box* expHall_box = new G4Box("World",
                                   gD->GetHalfWorldDX(),
                                   gD->GetHalfWorldDY(),
                                   gD->GetHalfWorldDZ());

    expHall_log
            = new G4LogicalVolume(expHall_box,
                                  m_fillingWorldMaterial,
                                  "World",
                                  0,0,0);
    expHall_log->SetVisAttributes ( invisibleVisAtt );

    expHall_phys
            = new G4PVPlacement(0,
                                G4ThreeVector(0.,0.,0.),
                                expHall_log,
                                "World",
                                0x0,
                                false,
                                0);

    // There won't be buit-in detectors
    // user has to define position and rotation from macro

    // check if ready to build the rest of the geometry

    if ( m_clearanceToBuildGeometry ) {

        G4cout << "Clearance obtained to build the pixel devices." << G4endl;

        // Read the database
        m_geoDsc = new ReadGeoDescription("./models/pixeldetector.xml");

        // Don't keep in the database the detectors which are not used.
        // Searching the id's
        G4int nErased = m_geoDsc->UseTheseDetectorsOnly(m_detId);
        if(nErased) G4cout << "[INFO] " << nErased << " detectors have been erased from the data base"
                           <<   " because they have not been called from the macro !" << G4endl;

        // Get the detectors
        map<int, AllPixGeoDsc *> * geoMap = m_geoDsc->GetDetectorsMap();

        gD = (*(geoMap->begin())).second;

        CheckAllPixSetup();
        //BuildMediPix(m_posVector, m_rotVector);
        BuildPixelDevices(*geoMap);

        if(m_buildAppliancesFlag) // set in the DetectorMessenger
            BuildAppliances((*(geoMap->begin())).first);

        if(m_buildTestStructureFlag) // set in the DetectorMessenger
            BuildTestStructure(m_TestStructure_type);

    }

    // Report absolute center of Si wafers
    if(!m_absolutePosSiWafer.empty()) {
        G4cout << "Absolute position of the Si wafers. center(x,y,z) : " << G4endl;
        vector<G4ThreeVector>::iterator itr = m_absolutePosSiWafer.begin();
        G4int cntr = 0;
        for( ; itr != m_absolutePosSiWafer.end() ; itr++) {

            G4cout << "   device [" << cntr << "] : "
                   << (*itr).getX()/mm << " , "
                   << (*itr).getY()/mm << " , "
                   << (*itr).getZ()/mm << " [mm]"
                   << G4endl;
            cntr++;
        }
    }

    return expHall_phys;
}

void AllPixDetectorConstruction::CheckAllPixSetup(){

    G4cout << "*** Pixel device Geo configuration --> verifying setup " << G4endl;

    if(m_nPositions == 0){
        G4cout << "No detectors have been scheduled for construction." << G4endl;
        _BUILD_MEDIPIX_MSG();
        _WRONG_CONFIG_ABORT_MSG();
        exit(1);
    }

    if(m_nPositions != m_nRotations || m_nPositions != m_nIds){
        G4cout << "You must have at least a \"setRotation\", \"setPosition\" and \"setId\" statement" << G4endl;
        G4cout << "  for each detector." << G4endl;
        _BUILD_MEDIPIX_MSG();
        _WRONG_CONFIG_ABORT_MSG();
        exit(1);
    }

    G4cout << "*** device Geo configuration --> OK" << G4endl;

}

void AllPixDetectorConstruction::SetDetectorID(G4int id){

    m_detId.push_back(id);

    // get the iterator informed
    m_detIdItr = m_detId.end() - 1;

    m_nIds++;
}

void AllPixDetectorConstruction::SetDetectorPosition(G4ThreeVector pos){

    //m_posVector.push_back(pos);
    if(m_detId.empty()){
        _BUILD_MEDIPIX_MSG();
        exit(1);
    }

    m_posVector[*m_detIdItr] = pos; // last Id
    m_nPositions++;

}

void AllPixDetectorConstruction::SetDetectorRotation(G4ThreeVector rot){

    /*
    m_rotVector.push_back( new G4RotationMatrix() );
    m_rotVector[m_nRotations]->rotateX(rot.x());
    m_rotVector[m_nRotations]->rotateY(rot.y());
    m_rotVector[m_nRotations]->rotateZ(rot.z());
     */

    if(m_detId.empty()){
        _BUILD_MEDIPIX_MSG();
        exit(1);
    }

    m_rotVector[*m_detIdItr] = new G4RotationMatrix();
    m_rotVector[*m_detIdItr]->rotateX(rot.x());
    m_rotVector[*m_detIdItr]->rotateY(rot.y());
    m_rotVector[*m_detIdItr]->rotateZ(rot.z());

    //G4cout << m_rotVector[m_nRotations]->getPhi() << " "
    // << m_rotVector[m_nRotations]->getPsi() << " "
    // << m_rotVector[m_nRotations]->getTheta() << G4endl;

    m_nRotations++;
}

void AllPixDetectorConstruction::SetLowTHL(G4double lowTHL){
    m_lowThlVector.push_back(lowTHL);
}

/**
 * Postition of the test structure.
 * There could be many test structures,
 *  We use a vector also.
 */
void AllPixDetectorConstruction::SetTestStructurePosition(G4ThreeVector pos){
    m_posVectorTestStructure[m_nTestPositions] = pos;
    m_nTestPositions++;
}

void AllPixDetectorConstruction::SetTestStructureRotation(G4ThreeVector rot){
    m_rotVectorTestStructure[m_nTestRotation] = new G4RotationMatrix();
    m_rotVectorTestStructure[m_nTestRotation]->rotateX(rot.x());
    m_rotVectorTestStructure[m_nTestRotation]->rotateY(rot.y());
    m_rotVectorTestStructure[m_nTestRotation]->rotateZ(rot.z());
    m_nTestRotation++;
}

void AllPixDetectorConstruction::SetTestStructureDetectorLink(G4int id){
    m_detectorLinkTestStructure[m_nTestPositions - 1] = id;
}

/*  Apliances position.  Rotation is not needed since this appliances
 *   follow the medipix.  This positions will be used with respect to
 *   the wrapper.
 */
void AllPixDetectorConstruction::SetAppliancePosition(G4ThreeVector pos){
    m_posVectorAppliances[*m_detIdItr] = pos;
    m_nAppliancesPositions++;
}

void AllPixDetectorConstruction::SetWrapperEnhancement(G4ThreeVector vals){
    m_vectorWrapperEnhancement[*m_detIdItr] = vals;
    m_nWrapperEnhancement++;
}

/**
 *  The user does not pass the actual material name.  Just "Air" or "Vacuum".
 *  if none, Air is used as default
 */
void AllPixDetectorConstruction::SetWorldMaterial(G4String mat){

    if(mat == "Vacuum") {
        G4cout << "User action --> World volume material : Vacuum " << endl;
        m_fillingWorldMaterial = m_Vacuum;
        m_userDefinedWorldMaterial = true;
    }
    // else default Air

}

#ifdef _EUTELESCOPE
void AllPixDetectorConstruction::SetScintPos(G4ThreeVector pos){
    m_scintPos.push_back(pos);
}

#endif

//void AllPixDetectorConstruction::BuildMediPix(vector<G4ThreeVector> pos, vector<G4RotationMatrix *> rot) {
void AllPixDetectorConstruction::BuildPixelDevices(map<int, AllPixGeoDsc *> geoMap) {

    // materials
    G4NistManager * nistman = G4NistManager::Instance();

    //G4Element* elCarbon= man->FindOrBuildElement("C");

    //std::vector<G4String> mats = nistman->GetNistMaterialNames();
    //std::vector<G4String>::iterator itr = mats.begin();
    //for( ; itr != mats.end() ; itr++) G4cout << *itr << G4endl;

    // Si
    G4Material * Silicon = nistman->FindOrBuildMaterial("G4_Si");

    // this is supposed to be epoxy // FIXME
    G4Material * Epoxy = nistman->FindOrBuildMaterial("G4_PLEXIGLASS");

    // Elements
    G4Element* Sn = new G4Element("Tin", "Sn", 50., 118.710*g/mole);
    G4Element* Pb = new G4Element("Lead","Pb", 82., 207.2*g/mole);

    // Materials from Combination, SnPb eutectic
    G4Material* Solder = new G4Material("Solder", 8.4*g/cm3, 2);
    Solder->AddElement(Sn,63);
    Solder->AddElement(Pb,37);


    ///////////////////////////////////////////////////////////////////////
    // vis attributes
    G4VisAttributes * pixelVisAtt= new G4VisAttributes(G4Color::Blue());
    pixelVisAtt->SetLineWidth(1);
    pixelVisAtt->SetForceSolid(true);
    //pixelVisAtt->SetForceWireframe(true);

    G4VisAttributes * BoxVisAtt= new G4VisAttributes(G4Color(0,1,1,1));
    BoxVisAtt->SetLineWidth(2);
    BoxVisAtt->SetForceSolid(true);
    //BoxVisAtt->SetVisibility(false);

    G4VisAttributes * CoverlayerVisAtt= new G4VisAttributes(G4Color::White());
    CoverlayerVisAtt->SetLineWidth(2);
    CoverlayerVisAtt->SetForceSolid(true);

    G4VisAttributes * ChipVisAtt= new G4VisAttributes(G4Color::Gray());
    ChipVisAtt->SetLineWidth(2);
    ChipVisAtt->SetForceSolid(true);
    //BoxVisAtt->SetVisibility(false);

    G4VisAttributes * BumpBoxVisAtt = new G4VisAttributes(G4Color(0,1,0,1.0));
    BumpBoxVisAtt->SetLineWidth(1);
    BumpBoxVisAtt->SetForceSolid(false);
    BumpBoxVisAtt->SetVisibility(true);

    G4VisAttributes * BumpVisAtt = new G4VisAttributes(G4Color::Yellow());
    BumpVisAtt->SetLineWidth(2);
    BumpVisAtt->SetForceSolid(true);
    //BumpVisAtt->SetVisibility(true);
    //BumpVisAtt->SetForceAuxEdgeVisible(true);

    G4VisAttributes * pcbVisAtt = new G4VisAttributes(G4Color::Green());
    pcbVisAtt->SetLineWidth(1);
    pcbVisAtt->SetForceSolid(true);

    G4VisAttributes * guardRingsVisAtt = new G4VisAttributes(G4Color(0.5,0.5,0.5,1));
    guardRingsVisAtt->SetLineWidth(1);
    guardRingsVisAtt->SetForceSolid(true);

    G4VisAttributes * wrapperVisAtt = new G4VisAttributes(G4Color(1,0,0,0.9));
    wrapperVisAtt->SetLineWidth(1);
    wrapperVisAtt->SetForceSolid(false);
    wrapperVisAtt->SetVisibility(true);

    ///////////////////////////////////////////////////////////////////////


    //vector<G4ThreeVector>::iterator posItr = pos.begin();
    //map<int, G4ThreeVector>::iterator posItr = m_posVector.begin();

    G4int nOfDevices = (G4int) m_posVector.size();

    ///////////////////////////////////////////////////////////
    // The pixel detector !
    //
    pair<G4String, G4String> wrapperName = make_pair("wrapper", "");
    pair<G4String, G4String> PCBName = make_pair("PCB", "");
    pair<G4String, G4String> BoxName = make_pair("Box", "");
    pair<G4String, G4String> CoverlayerName = make_pair("Coverlayer", "");
    pair<G4String, G4String> SliceName = make_pair("Slice", "");
    pair<G4String, G4String> GuardRingsExtName = make_pair("GuardRingsExt", "");
    pair<G4String, G4String> GuardRingsName = make_pair("GuardRings", "");
    pair<G4String, G4String> PixelName = make_pair("Pixel", "");
    pair<G4String, G4String> ChipName = make_pair("Chip", "");
    pair<G4String, G4String> SDName = make_pair("BoxSD", "");
    pair<G4String, G4String> BumpName = make_pair("Bump", "");
    pair<G4String, G4String> BumpBoxName = make_pair("BumpBox", "");

    // log an phys
    G4cout << "Building " << nOfDevices << " device(s) ..." << G4endl;

    // SD manager
    G4SDManager * SDman = G4SDManager::GetSDMpointer();


    // User limits applied only to Si wafers.  Setting step.
    if ( !m_ulim ) m_ulim = new G4UserLimits(m_maxStepLengthSensor);

    //for ( ; posItr != pos.end() ; posItr++ )
    //for ( ; posItr != m_posVector.end() ; posItr++ )

    vector<int>::iterator detItr = m_detId.begin();

    for( ; detItr != m_detId.end() ; detItr++)
    {
        G4cout << "          start detector " << (*detItr) << " (" << geoMap[*detItr] << ")" << G4endl;

        /*------------------------------------ Solid definitions ----------------------------------------*/

        // replicated solids (same object/names everywhere)
        G4Box * Box_slice = new G4Box(SliceName.first,
                                      geoMap[*detItr]->GetHalfPixelX(),
                geoMap[*detItr]->GetHalfSensorY(),
                geoMap[*detItr]->GetHalfSensorZ());

        G4Box * Box_pixel = new G4Box(PixelName.first,
                                      geoMap[*detItr]->GetHalfPixelX(),
                geoMap[*detItr]->GetHalfPixelY(),
                geoMap[*detItr]->GetHalfPixelZ());

        //Bump itself
        G4UnionSolid * aBump = 0;
        G4Box * Bump_Box = 0;
        G4double bump_radius =0;
        G4double bump_dr =0;
        G4double bump_height=0;
        if ( geoMap[*detItr]->GetBumpHeight() != 0.0 and geoMap[*detItr]->GetHalfChipZ() != 0 ) {
            bump_radius =geoMap[*detItr]->GetBumpRadius();
            bump_height = geoMap[*detItr]->GetBumpHeight();
            bump_dr =  geoMap[*detItr]->GetBumpDr();
            G4Sphere * aBump_Sphere = new G4Sphere(BumpName.first+"sphere",0,bump_radius,0,360*deg,0,360*deg);
            G4Tubs * aBump_Tube= new G4Tubs(BumpName.first+"Tube", 0., bump_radius-bump_dr, bump_height/2., 0., 360 *deg);
            aBump = new G4UnionSolid(BumpName.first,aBump_Sphere,aBump_Tube);

            //bumps containing volume

            Bump_Box = new G4Box(BumpBoxName.first,
                                 geoMap[*detItr]->GetHalfSensorX(),
                    geoMap[*detItr]->GetHalfSensorY(),
                    bump_height/2.);
        }

        // Build names
        char temp[128];
        sprintf(temp, "%d", (*detItr));
        wrapperName.second = wrapperName.first + "_" + temp;
        PCBName.second = PCBName.first + "_" + temp;
        BoxName.second = BoxName.first + "_" + temp;
        CoverlayerName.second = CoverlayerName.first + "_" + temp;
        GuardRingsExtName.second = GuardRingsExtName.first  + "_" + temp;
        GuardRingsName.second = GuardRingsName.first  + "_" + temp;
        SliceName.second = BoxName.first + "_" + temp;
        PixelName.second = BoxName.first + "_" + temp;
        ChipName.second = ChipName.first + "_" + temp;
        SDName.second = SDName.first + "_" + temp; // BoxSD_XXX
        BumpName.second = BumpName.first + "_" + temp; // BoxSD_XXX
        BumpBoxName.second = BumpBoxName.first + "_" + temp; // BoxSD_XXX

        // Solids, I want different objects/names per medipix
        //  later on, they could be different
        G4Box * Box_box = new G4Box(BoxName.second,
                                    geoMap[*detItr]->GetHalfSensorX(),
                geoMap[*detItr]->GetHalfSensorY(),
                geoMap[*detItr]->GetHalfSensorZ());

        G4Box * Chip_box = 0;
        if ( geoMap[*detItr]->GetHalfChipZ() != 0 ) {
            // Chip box
            Chip_box = new G4Box(ChipName.second,
                                 geoMap[*detItr]->GetHalfChipX(),
                    geoMap[*detItr]->GetHalfChipY(),
                    geoMap[*detItr]->GetHalfChipZ());
        }

        // If the coverlayer is requested.  It is forced to fit the sensor in X,Y.
        // The user can pick the thickness.
        G4Box * Coverlayer_box = nullptr;
        if ( geoMap[*detItr]->IsCoverlayerON() ) {
            Coverlayer_box = new G4Box(
                        CoverlayerName.second,
                        geoMap[*detItr]->GetHalfSensorX(),
                    geoMap[*detItr]->GetHalfSensorY(),
                    geoMap[*detItr]->GetHalfCoverlayerZ()
                    );
        }

        // Guard rings will be GuardRingsExt - Box
        G4Box * Box_GuardRings_Ext = new G4Box(
                    GuardRingsExtName.second,
                    geoMap[*detItr]->GetHalfSensorX() + (geoMap[*detItr]->GetSensorExcessHRight() + geoMap[*detItr]->GetSensorExcessHLeft()),
                geoMap[*detItr]->GetHalfSensorY() + (geoMap[*detItr]->GetSensorExcessHTop() + geoMap[*detItr]->GetSensorExcessHBottom()),
                geoMap[*detItr]->GetHalfSensorZ()); // same depth as the sensor

        G4VSolid * Solid_GuardRings =  new G4SubtractionSolid(GuardRingsName.second,
                                                              Box_GuardRings_Ext,
                                                              Box_box);


        G4Box* PCB_box = 0;
        if ( geoMap[*detItr]->GetHalfPCBZ() != 0 ) {
            PCB_box = new G4Box(
                        PCBName.second,
                        geoMap[*detItr]->GetHalfPCBX(),
                    geoMap[*detItr]->GetHalfPCBY(),
                    geoMap[*detItr]->GetHalfPCBZ());
        }

        // The wrapper might be enhanced when the user set up
        //  Appliances to the detector (extra layers, etc).
        G4double wrapperHX = geoMap[*detItr]->GetHalfWrapperDX();
        G4double wrapperHY = geoMap[*detItr]->GetHalfWrapperDY();
        G4double wrapperHZ = geoMap[*detItr]->GetHalfWrapperDZ();

        // Apply the enhancement to the medipixes
        // We can have N medipixes and K enhancements, where K<=N.
        // For instance, for 2 medipixes.  We can have.
        // medipix 1 --> with enhancement
        // medipix 2 --> no enhancement
        if ( m_vectorWrapperEnhancement.find(*detItr) != m_vectorWrapperEnhancement.end() ) {
            wrapperHX += m_vectorWrapperEnhancement[*detItr].x()/2.; // half
            wrapperHY += m_vectorWrapperEnhancement[*detItr].y()/2.;
            wrapperHZ += m_vectorWrapperEnhancement[*detItr].z()/2.;
            // temp test
            //wrapperDZ += 10000*um;
        }

        G4cout << "Wrapper Dimensions [mm] : "
               << TString::Format("hX=%3.3f hY=%3.3f hZ=%3.3f",wrapperHX/mm,wrapperHY/mm,wrapperHZ/mm)
               << G4endl;

        G4Box* wrapper_box = new G4Box(wrapperName.second,
                                       wrapperHX,
                                       wrapperHY,
                                       wrapperHZ);

        //G4RotationMatrix yRot45deg;   // Rotates X and Z axes only
        //yRot45deg.rotateY(M_PI/4.*rad);
        //G4ThreeVector  translation(0, 0, 50*mm);
        //G4UnionSolid  wrapper_box_union("wrapper_box_union",
        //                                wrapper_box, wrapper_box, &yRot45deg, translation);

        /*------------------------------------ Logical and physical volumes definitions ----------------------------------------*/

        ///////////////////////////////////////////////////////////
        // wrapper
        m_wrapper_log[(*detItr)] = new G4LogicalVolume(wrapper_box,
                                                       m_fillingWorldMaterial,
                                                       wrapperName.second+"_log");
        m_wrapper_log[(*detItr)]->SetVisAttributes(wrapperVisAtt);


        G4double sensorOffsetX = geoMap[*detItr]->GetSensorXOffset();
        G4double sensorOffsetY = geoMap[*detItr]->GetSensorYOffset();

        G4ThreeVector posWrapper = m_posVector[(*detItr)];

        // Apply position Offset for the wrapper due to the enhancement
        if ( m_vectorWrapperEnhancement.find(*detItr) != m_vectorWrapperEnhancement.end() ) {
            posWrapper.setX(posWrapper.x() + m_vectorWrapperEnhancement[*detItr].x()/2.);
            posWrapper.setY(posWrapper.y() + m_vectorWrapperEnhancement[*detItr].y()/2.);
            posWrapper.setZ(posWrapper.z() + m_vectorWrapperEnhancement[*detItr].z()/2.);
        } else {
            posWrapper.setX(posWrapper.x() - sensorOffsetX );
            posWrapper.setY(posWrapper.y() - sensorOffsetY );
            posWrapper.setZ(posWrapper.z() );
        }


        //G4Transform3D transform( *m_rotVector[(*detItr)], posWrapper );

        // starting at user position --> vector pos
        m_wrapper_phys[(*detItr)] = new G4PVPlacement(
                    //transform,
                    m_rotVector[(*detItr)],
                    posWrapper,
                    m_wrapper_log[(*detItr)],
                wrapperName.second+"_phys",
                expHall_log,
                false,
                (*detItr), // copy number
                true);

        // Apply a translation to the wrapper first
        //G4Transform3D tA = G4TranslateX3D( -1*geoMap[*detItr]->GetSensorXOffset() );
        //G4Transform3D tB = G4TranslateY3D( -1*geoMap[*detItr]->GetSensorYOffset() );
        //G4Transform3D transform = tA * tB;
        //m_wrapper_phys[(*detItr)]->SetTranslation( transform.getTranslation() );

        ///////////////////////////////////////////////////////////
        // PCB
        // The PCB is placed respect to the wrapper.
        // Needs to be pushed -half Si wafer in z direction
        m_PCB_log[(*detItr)] = 0;
        if(geoMap[*detItr]->GetHalfPCBZ()!=0){

            m_PCB_log[(*detItr)] = new G4LogicalVolume(PCB_box,
                                                       Epoxy,
                                                       PCBName.second+"_log");
            m_PCB_log[(*detItr)]->SetVisAttributes(pcbVisAtt);

        }

        if ( geoMap[*detItr]->GetHalfChipZ()!=0 ) {

            ///////////////////////////////////////////////////////////
            // Chip
            // The Si wafer is placed respect to the wrapper.
            // Needs to be pushed -half Si wafer in z direction

            m_Chip_log[(*detItr)] = new G4LogicalVolume(Chip_box,
                                                        Silicon,
                                                        ChipName.second+"_log");
            m_Chip_log[(*detItr)]->SetVisAttributes(ChipVisAtt);

            if ( geoMap[*detItr]->GetBumpHeight()!=0.0 ) {
                ///////////////////////////////////////////////////////////
                // Bumps
                m_Bumps_log[(*detItr)] = new G4LogicalVolume(Bump_Box,
                                                             m_Air,
                                                             BumpBoxName.second+"_log");
                m_Bumps_log[(*detItr)]->SetVisAttributes(BumpBoxVisAtt);
            }


        }


        ///////////////////////////////////////////////////////////
        // Device
        // The Si wafer is placed respect to the wrapper.
        // Needs to be pushed -half Si wafer in z direction

        m_Box_log[(*detItr)] = new G4LogicalVolume(Box_box,
                                                   Silicon,
                                                   BoxName.second+"_log");

        m_Box_log[(*detItr)]->SetVisAttributes(BoxVisAtt);

        // positions
        G4ThreeVector posCoverlayer(0,0,0);
        G4ThreeVector posDevice(sensorOffsetX,sensorOffsetY,0);
        G4ThreeVector posBumps(0,0,0);
        G4ThreeVector posChip(0,0,0);
        G4ThreeVector posPCB(0,0,0);


        if ( geoMap[*detItr]->IsCoverlayerON() ) {

            posCoverlayer.setX( posDevice.x() );
            posCoverlayer.setY( posDevice.y() );
            posCoverlayer.setZ(
                        wrapperHZ
                        - geoMap[*detItr]->GetHalfCoverlayerZ()
                    );
        }

        // Apply position Offset for the detector due to the enhancement
        if ( m_vectorWrapperEnhancement.find(*detItr) != m_vectorWrapperEnhancement.end() ) {

            posDevice.setX(posDevice.x() - m_vectorWrapperEnhancement[*detItr].x()/2.);
            posDevice.setY(posDevice.y() - m_vectorWrapperEnhancement[*detItr].y()/2.);
            posDevice.setZ(
                        wrapperHZ
                        - 2.*geoMap[*detItr]->GetHalfCoverlayerZ()
                    - geoMap[*detItr]->GetHalfSensorZ()
                    - m_vectorWrapperEnhancement[*detItr].z()/2.
                    );
            //posDevice.z() - m_vectorWrapperEnhancement[*detItr].z()/2.);
        } else {
            posDevice.setX(posDevice.x() );
            posDevice.setY(posDevice.y() );
            posDevice.setZ(
                        wrapperHZ
                        - 2.*geoMap[*detItr]->GetHalfCoverlayerZ()
                    - geoMap[*detItr]->GetHalfSensorZ()
                    );
        }

        ///////////////////////////////////////////////////////////
        // Coverlayer if requested (typically made of Al, but user configurable)

        if ( geoMap[*detItr]->IsCoverlayerON() ) {

            // Find out about the material that the user requested
            G4Material * covermat = nistman->FindOrBuildMaterial( geoMap[*detItr]->GetCoverlayerMat() );

            // If this is an inexistent choice, them force Aluminum
            if ( covermat == nullptr ) {
                covermat = nistman->FindOrBuildMaterial( "G4_Al" );
            }

            m_Coverlayer_log[(*detItr)] = new G4LogicalVolume( Coverlayer_box,
                                                               covermat,
                                                               CoverlayerName.second+"_log");

            m_Coverlayer_log[(*detItr)]->SetVisAttributes(CoverlayerVisAtt);
        }

        // Calculation of position of the different physical volumes
        if ( geoMap[*detItr]->GetHalfChipZ() != 0 ) {

            posBumps.setX( posDevice.x() );
            posBumps.setY( posDevice.y() );
            posBumps.setZ(
                        wrapperHZ
                        - 2.*geoMap[*detItr]->GetHalfCoverlayerZ()
                    - 2.*geoMap[*detItr]->GetHalfSensorZ()
                    - (bump_height/2.)
                    );

            //posDevice.z() -
            //geoMap[*detItr]->GetHalfSensorZ() +
            //geoMap[*detItr]->GetChipZOffset() -
            //bump_height/2
            //);

            posChip.setX( posDevice.x() + geoMap[*detItr]->GetChipXOffset() );
            posChip.setY( posDevice.y() + geoMap[*detItr]->GetChipYOffset() );
            posChip.setZ(
                        wrapperHZ
                        - 2.*geoMap[*detItr]->GetHalfCoverlayerZ()
                    - 2.*geoMap[*detItr]->GetHalfSensorZ()
                    - bump_height
                    - geoMap[*detItr]->GetHalfChipZ()
                    );

            //posDevice.z() -
            //geoMap[*detItr]->GetHalfSensorZ() -
            //geoMap[*detItr]->GetHalfChipZ() +
            //geoMap[*detItr]->GetChipZOffset() -
            //bump_height
            //);

        } else {
            // Make sure no offset because of bumps for PCB is calculated if chip is not included
            bump_height = 0;
        }

        posPCB.setX( 0 ); //- 1.*geoMap[*detItr]->GetSensorXOffset() );
        posPCB.setY( 0 ); //- 1.*geoMap[*detItr]->GetSensorYOffset() );
        posPCB.setZ(
                    wrapperHZ
                    - 2.*geoMap[*detItr]->GetHalfCoverlayerZ()
                - 2.*geoMap[*detItr]->GetHalfSensorZ()
                - bump_height
                - 2.*geoMap[*detItr]->GetHalfChipZ()
                - geoMap[*detItr]->GetHalfPCBZ()
                );

        //posDevice.z() -
        //geoMap[*detItr]->GetHalfSensorZ() -
        //2*geoMap[*detItr]->GetHalfChipZ() -
        //geoMap[*detItr]->GetHalfPCBZ() -
        //bump_height
        //);

        G4cout << "- Coverlayer position  : " << posCoverlayer << G4endl;
        G4cout << "- Sensor position      : " << posDevice << G4endl;
        G4cout << "- Bumps position       : " << posBumps << G4endl;
        G4cout << "- Chip position        : " << posChip << G4endl;
        G4cout << "- PCB position         : " << posPCB << G4endl;


        /*------------------------------------ Physical placement  ----------------------------------------*/

        m_PCB_phys[(*detItr)] = 0;
        if(geoMap[*detItr]->GetHalfPCBZ()!=0){

            m_PCB_phys[(*detItr)] = new G4PVPlacement(
                        0,
                        posPCB,
                        m_PCB_log[(*detItr)],
                    PCBName.second+"_phys",
                    m_wrapper_log[(*detItr)], // mother log
                    false,
                    (*detItr),
                    true); // copy number
        }
        if(geoMap[*detItr]->GetHalfChipZ()!=0){

            m_Chip_phys[(*detItr)] = new G4PVPlacement(0,
                                                       posChip,
                                                       m_Chip_log[(*detItr)],
                    ChipName.second+"_phys",
                    m_wrapper_log[(*detItr)], // mother log
                    false,
                    (*detItr), // copy number
                    true); // check overlap
            if(geoMap[*detItr]->GetBumpHeight()!=0.0){
                m_Bumps_phys[(*detItr)] = new G4PVPlacement(0,
                                                            posBumps,
                                                            m_Bumps_log[(*detItr)],
                        BumpBoxName.second+"_phys",
                        m_wrapper_log[(*detItr)], // mother log
                        false,
                        (*detItr), // copy number
                        true); // check overlap
            };
        }

        m_Box_phys[(*detItr)] = new G4PVPlacement(
                    0,
                    posDevice,
                    m_Box_log[(*detItr)],
                BoxName.second+"_phys",
                m_wrapper_log[(*detItr)], // mother log
                false,
                (*detItr), // copy number
                true); // check overlap

        // coverlayer
        if ( geoMap[*detItr]->IsCoverlayerON() ) {

            m_Coverlayer_phys[(*detItr)] = new G4PVPlacement(
                        0,
                        posCoverlayer,
                        m_Coverlayer_log[(*detItr)],
                    CoverlayerName.second+"_phys",
                    m_wrapper_log[(*detItr)], // mother log
                    false,
                    (*detItr), // copy number
                    true); // check overlap

        }

        ///////////////////////////////////////////////////////////
        // bumps

        if ( geoMap[*detItr]->GetHalfChipZ() != 0 ) {

            m_Bumps_Cell_log[(*detItr)] = new G4LogicalVolume(aBump,Solder,BumpBoxName.second+"_log" );
            m_Bumps_Cell_log[(*detItr)]->SetVisAttributes(BumpVisAtt);

            //		m_Bumps_Slice_log[(*detItr)] = new G4LogicalVolume(Bump_Slice_Box,m_Air,BumpSliceName.second+"_log");
            //		m_Bumps_Slice_log[(*detItr)]->SetVisAttributes(BumpSliceVisAtt);


            parameterization = new Allpix_BumpsParameterization( geoMap[*detItr] );
            G4int NPixTot = geoMap[*detItr]->GetNPixelsX()*geoMap[*detItr]->GetNPixelsY();
            new G4PVParameterised(BumpName.second+"phys",
                                  m_Bumps_Cell_log[(*detItr)],        // logical volume
                    m_Bumps_log[(*detItr)],             // mother volume
                    kUndefined,                         // axis
                    NPixTot,                            // replicas
                    parameterization);                  // G4VPVParameterisation
        }


        ///////////////////////////////////////////////////////////
        // slices and pixels
        m_Slice_log[(*detItr)] = new G4LogicalVolume(Box_slice,
                                                     Silicon,
                                                     SliceName.second); // 0,0,0);
        //m_Slice_log[(*detItr)]->SetUserLimits(ulim);
        m_Pixel_log[(*detItr)] = new G4LogicalVolume(Box_pixel,
                                                     Silicon,
                                                     PixelName.second); // 0,0,0);

        if ( m_ulim ) m_Pixel_log[(*detItr)]->SetUserLimits(m_ulim);

        // divide in slices
        new G4PVDivision(
                    SliceName.second,
                    m_Slice_log[(*detItr)],
                m_Box_log[(*detItr)],
                kXAxis,
                //geoMap[*detItr]->GetPixelX(),
                geoMap[*detItr]->GetNPixelsX(),
                0); // offset

        new G4PVDivision(
                    PixelName.second,
                    m_Pixel_log[(*detItr)],
                m_Slice_log[(*detItr)],
                kYAxis,
                //geoMap[*detItr]->GetPixelY(),
                geoMap[*detItr]->GetNPixelsY(),
                0); // offset

        ///////////////////////////////////////////////////////////
        // Guard rings and excess area
        m_GuardRings_log[(*detItr)] = new G4LogicalVolume(Solid_GuardRings,
                                                          Silicon,
                                                          GuardRingsName.second+"_log");
        m_GuardRings_log[(*detItr)]->SetVisAttributes(guardRingsVisAtt);
        m_GuardRings_phys[(*detItr)] = new G4PVPlacement(0,
                                                         posDevice,
                                                         m_GuardRings_log[(*detItr)],
                GuardRingsName.second+"_phys",
                m_wrapper_log[(*detItr)], // mother log
                false,
                0,//(*detItr), // copy number
                true); // check overlap

        // Find out where the center of the Silicon has been placed
        // for user information

        m_absolutePosSiWafer.push_back(posWrapper + posDevice);

        //G4cout << posWrapper << " " << posDevice << endl;

        // SD --> pixels !
        // AllPixTrackerSD instance needs to know the absolute position of
        //  the device and rotation.  'm_absolutePosSiWafer' and 'rot' are
        //  the good figures.

        G4ThreeVector origin(0,0,0);

        AllPixTrackerSD * aTrackerSD = new AllPixTrackerSD( SDName.second,
                                                            (posWrapper + posDevice),
                                                            posDevice,
                                                            geoMap[*detItr],
                m_rotVector[(*detItr)] );

        SDman->AddNewDetector( aTrackerSD );
        m_Pixel_log[(*detItr)]->SetSensitiveDetector( aTrackerSD );

        // Store the hit Collection name in the geometry
        geoMap[*detItr]->SetHitsCollectionName( aTrackerSD->GetHitsCollectionName() );

        G4cout << "          detector " << (*detItr) << " ... done" << G4endl;

        //m_wrapper_phys[(*detItr)]->SetTranslation(posWrapper-posDevice);

        //m_wrapper_phys[(*detItr)]->SetRotation( m_rotVector[(*detItr)] );

    }


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4RunManager.hh"
#include "G4UserEventAction.hh"
#include "AllPixEventAction.hh"

void AllPixDetectorConstruction::UpdateGeometry()
{

    // Geometry building clearance.
    //  At this point the user had the chance to give the coordinates
    //  of the medipixes from the macro.
    m_clearanceToBuildGeometry = true;

    // build extra medipixes
    G4RunManager::GetRunManager()->DefineWorldVolume(Construct());

    // setup digitizers for new detectors
    AllPixEventAction * ea = (AllPixEventAction *) G4RunManager::GetRunManager()->GetUserEventAction();
    ea->SetupDigitizers();

    // setup all thl // FIXME
    /*
    G4int itr = 0;
    for( ; itr < (G4int)m_lowThlVector.size() ; itr++){
        ea->SetDetectorDigitInput(m_lowThlVector[itr], itr); // thl !
    }
     */
}

/**
 *  Someone get's informed throught the messenger about the prefix for the output file.
 */
void AllPixDetectorConstruction::SetOutputFilePrefix(G4String val){
    m_outputFilePrefix = val;
}
void AllPixDetectorConstruction::SetMaxStepLengthSensor(G4double val) {

    m_maxStepLengthSensor = val;

    if ( ! m_ulim ) {
        m_ulim = new G4UserLimits ( m_maxStepLengthSensor );
        G4cout << "Changing MaxStepLengthSensor to " << m_maxStepLengthSensor/um << "um" << G4endl;
    } else {
        G4cout << "[WARNING] Couldn't setup a new MaxStepLengthSensor. " << G4endl;
        G4cout << "          The command is probably called too soon in the macro." << G4endl;
        G4cout << "          Using default value of " << m_maxStepLengthSensor/um << "um" << G4endl;
    }

}

#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4QuadrupoleMagField.hh"
#include "G4UniformMagField.hh"
#include "MorourgoMagField.hh"
#include "G4PropagatorInField.hh"
void AllPixDetectorConstruction::SetPeakMagField(G4double fieldValue)
{
    //apply a global uniform magnetic field along Z axis
    G4FieldManager * fieldMgr
            = G4TransportationManager::GetTransportationManager()->GetFieldManager();
    G4TransportationManager* tmanager = G4TransportationManager::GetTransportationManager();
    tmanager->GetPropagatorInField()->SetLargestAcceptableStep(1*mm);
    if ( fieldValue != 0. )
    {

        // FIXME !!! --> feed this value from the macro
        //		m_magField = new MorourgoMagField(fieldValue, 252.5*mm);
        //		fieldMgr->SetDetectorField(m_magField);
        //		fieldMgr->CreateChordFinder(m_magField);
        m_magField = new G4UniformMagField ( G4ThreeVector(0.,0.,fieldValue) );
        fieldMgr->SetDetectorField(m_magField);
        fieldMgr->CreateChordFinder(m_magField);

        fieldMgr->SetMinimumEpsilonStep( 1e-7 );
        fieldMgr->SetMaximumEpsilonStep( 1e-6 );
        fieldMgr->SetDeltaOneStep( 0.05e-3 * mm );  // 0.5 micrometer


    } else {
        m_magField = 0x0;
        //fieldMgr->SetDetectorField(m_magField);
    }


    //if(m_magField) {

    // test at center of telescope (same point that center of magnet)
    //		G4double xyz[4] = { 1.*mm, 1.*mm, 0.*mm, 0.*mm };
    //		G4double fieldVal[3] = { 0, 0, 0 };
    //
    //		m_magField->GetFieldValue(xyz, fieldVal);
    //		std::cout << "y field (z=" << xyz[2]/mm << ") = " << fieldVal[1] << std::endl;
    //
    //		xyz[2] = 252.5*mm;
    //		m_magField->GetFieldValue(xyz, fieldVal);
    //		std::cout << "y field (z=" << xyz[2]/mm << ") = " << fieldVal[1] << std::endl;
    //
    //		xyz[2] = -3000*mm;
    //		m_magField->GetFieldValue(xyz, fieldVal);
    //		std::cout << "y field (z=" << xyz[2]/mm << ") = " << fieldVal[1] << std::endl;

    //}

}

