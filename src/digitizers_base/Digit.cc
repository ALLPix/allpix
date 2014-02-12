/**
 *  Author:
 *    AAAuthor
 *
 *  allpix Authors:
 *   John Idarraga <idarraga@cern.ch>
 *   Mathieu Benoit <Mathieu.Benoit@cern.ch>
 */

#include "__Digit.hh"

G4Allocator<__Digit> __DigitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

__Digit::__Digit()
{
	pixelIDX = -1;
	pixelIDY = -1;
	pixelCounts = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

__Digit::~__Digit()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

__Digit::__Digit(const __Digit& right)
: AllPixDigitInterface()
{
	pixelIDX = right.pixelIDX;
	pixelIDY = right.pixelIDY;
	pixelCounts = right.pixelCounts;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

const __Digit& __Digit::operator=(const __Digit& right)
{

	pixelIDX = right.pixelIDX;
	pixelIDY = right.pixelIDY;
	pixelCounts = right.pixelCounts;
	return *this;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

int __Digit::operator==(const __Digit& right) const
		{
	return ((pixelIDX==right.pixelIDX)&&
			(pixelIDY==right.pixelIDY)&&
			(pixelCounts==right.pixelCounts));
		}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void __Digit::Draw()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void __Digit::Print()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
