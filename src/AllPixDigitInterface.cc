/*
 * AllPixDigitInterface.cc
 *
 *  Created on: Jun 17, 2015
 *      Author: idarraga
 */


#include "AllPixDigitInterface.hh"

AllPixDigitInterface::AllPixDigitInterface() : G4VDigi() {

}

AllPixDigitInterface::~AllPixDigitInterface(){;}

void AllPixDigitInterface::SetPixelCountsMultiTHL(G4int, G4int){;}

void AllPixDigitInterface::SetPixelEnergyDepMultiTHL(G4double, G4double){;}

G4int AllPixDigitInterface::GetPixelCountsMultiTHL(G4int){return 0;}

void AllPixDigitInterface::IncreasePixelCountsMultiTHL(G4int){;}


void AllPixDigitInterface::Set_posX_WithRespectoToPixel(G4double /*pos*/){;}
void AllPixDigitInterface::Set_posY_WithRespectoToPixel(G4double /*pos*/){;}
void AllPixDigitInterface::Set_posZ_WithRespectoToPixel(G4double /*pos*/	){;}
G4double AllPixDigitInterface::Get_posX_WithRespectoToPixel(){return 0;}
G4double AllPixDigitInterface::Get_posY_WithRespectoToPixel(){return 0;}
G4double AllPixDigitInterface::Get_posZ_WithRespectoToPixel(){return 0;}
void  AllPixDigitInterface::SetPixelEnergyMC(G4double){;}
G4double  AllPixDigitInterface::GetPixelEnergyMC(){return 0;}
