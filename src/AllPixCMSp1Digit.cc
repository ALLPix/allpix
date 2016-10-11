/**
 *  Author:
 *    Paul Schuetze <paul.schuetze@desy.de>
 *
 *  allpix Authors:
 *   John Idarraga <idarraga@cern.ch>
 *   Mathieu Benoit <Mathieu.Benoit@cern.ch>
 */

#include "AllPixCMSp1Digit.hh"

G4Allocator<AllPixCMSp1Digit> AllPixCMSp1DigitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

AllPixCMSp1Digit::AllPixCMSp1Digit()
{
	pixelIDX = -1;
	pixelIDY = -1;
	pixelCounts = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

AllPixCMSp1Digit::~AllPixCMSp1Digit()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

AllPixCMSp1Digit::AllPixCMSp1Digit(const AllPixCMSp1Digit& right)
: AllPixDigitInterface()
{
	pixelIDX = right.pixelIDX;
	pixelIDY = right.pixelIDY;
	pixelCounts = right.pixelCounts;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

const AllPixCMSp1Digit& AllPixCMSp1Digit::operator=(const AllPixCMSp1Digit& right)
{

	pixelIDX = right.pixelIDX;
	pixelIDY = right.pixelIDY;
	pixelCounts = right.pixelCounts;
	return *this;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

int AllPixCMSp1Digit::operator==(const AllPixCMSp1Digit& right) const
		{
	return ((pixelIDX==right.pixelIDX)&&
			(pixelIDY==right.pixelIDY)&&
			(pixelCounts==right.pixelCounts));
		}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void AllPixCMSp1Digit::Draw()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void AllPixCMSp1Digit::Print()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
