/*
 * AllPixDigitInterface.cc
 *
 *  Created on: Jun 17, 2015
 *      Author: idarraga
 */


#include "AllPixDigitInterface.hh"

AllPixDigitInterface::AllPixDigitInterface() : G4VDigi() {

}

void AllPixDigitInterface::SetPixelCountsMultiTHL(G4int, G4int){;}

void AllPixDigitInterface::SetPixelEnergyDepMultiTHL(G4double, G4double){;}

G4int AllPixDigitInterface::GetPixelCountsMultiTHL(G4int){;}

void AllPixDigitInterface::IncreasePixelCountsMultiTHL(G4int){;}
