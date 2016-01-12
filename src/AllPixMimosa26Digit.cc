/**
 *  Author:  John Idarraga <idarraga@cern.ch>
 *  Digit for Mimosa26 detector
 */

#include "AllPixMimosa26Digit.hh"

G4Allocator<AllPixMimosa26Digit> AllPixMimosa26DigitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

AllPixMimosa26Digit::AllPixMimosa26Digit()
{
  pixelIDX = -1; 
  pixelIDY = -1; 
  pixelCounts = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

AllPixMimosa26Digit::~AllPixMimosa26Digit()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

AllPixMimosa26Digit::AllPixMimosa26Digit(const AllPixMimosa26Digit& right)
  : AllPixDigitInterface()
{
  pixelIDX = right.pixelIDX; 
  pixelIDY = right.pixelIDY; 
  pixelCounts = right.pixelCounts;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

const AllPixMimosa26Digit& AllPixMimosa26Digit::operator=(const AllPixMimosa26Digit& right)
{

  pixelIDX = right.pixelIDX; 
  pixelIDY = right.pixelIDY; 
  pixelCounts = right.pixelCounts;
  return *this;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

int AllPixMimosa26Digit::operator==(const AllPixMimosa26Digit& right) const
{ 
  return ((pixelIDX==right.pixelIDX)&&
	  (pixelIDY==right.pixelIDY)&&
	  (pixelCounts==right.pixelCounts));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void AllPixMimosa26Digit::Draw()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void AllPixMimosa26Digit::Print()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
