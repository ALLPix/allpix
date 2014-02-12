/**
 *  Author:  John Idarraga <idarraga@cern.ch>
 */

#include "AllPixFEI3StandardDigit.hh"

G4Allocator<AllPixFEI3StandardDigit> AllPixFEI3StandardDigitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

AllPixFEI3StandardDigit::AllPixFEI3StandardDigit()
{
  pixelIDX = -1; 
  pixelIDY = -1; 
  pixelCounts = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

AllPixFEI3StandardDigit::~AllPixFEI3StandardDigit()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

AllPixFEI3StandardDigit::AllPixFEI3StandardDigit(const AllPixFEI3StandardDigit& right)
  : AllPixDigitInterface()
{
  pixelIDX = right.pixelIDX; 
  pixelIDY = right.pixelIDY; 
  pixelCounts = right.pixelCounts;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

const AllPixFEI3StandardDigit& AllPixFEI3StandardDigit::operator=(const AllPixFEI3StandardDigit& right)
{

  pixelIDX = right.pixelIDX;
  pixelIDY = right.pixelIDY;
  pixelCounts = right.pixelCounts;
  return *this;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

int AllPixFEI3StandardDigit::operator==(const AllPixFEI3StandardDigit& right) const
{ 
  return ((pixelIDX==right.pixelIDX)&&
	  (pixelIDY==right.pixelIDY)&&
	  (pixelCounts==right.pixelCounts));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void AllPixFEI3StandardDigit::Draw()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void AllPixFEI3StandardDigit::Print()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
