/**
 *  Author:
 *    nalipour@cern.ch
 *
 *  allpix Authors:
 *   John Idarraga <idarraga@cern.ch>
 *   Mathieu Benoit <Mathieu.Benoit@cern.ch>
 */

#include "AllPixTMPXDigit.hh"

G4Allocator<AllPixTMPXDigit> AllPixTMPXDigitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

AllPixTMPXDigit::AllPixTMPXDigit()
{
	pixelIDX = -1;
	pixelIDY = -1;
	pixelCounts = -11;
	pixelEnergyDep = -11; // MC // Corrected MC charge (detector effects included, at Digi step)
	pixelEnergyMC = -11; //MC only //nalipour
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

AllPixTMPXDigit::~AllPixTMPXDigit()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

AllPixTMPXDigit::AllPixTMPXDigit(const AllPixTMPXDigit& right)
: AllPixDigitInterface()
{
	pixelIDX = right.pixelIDX;
	pixelIDY = right.pixelIDY;
	pixelCounts = right.pixelCounts;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

const AllPixTMPXDigit& AllPixTMPXDigit::operator=(const AllPixTMPXDigit& right)
{

	pixelIDX = right.pixelIDX;
	pixelIDY = right.pixelIDY;
	pixelCounts = right.pixelCounts;
	return *this;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

int AllPixTMPXDigit::operator==(const AllPixTMPXDigit& right) const
		{
	return ((pixelIDX==right.pixelIDX)&&
			(pixelIDY==right.pixelIDY)&&
			(pixelCounts==right.pixelCounts));
		}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void AllPixTMPXDigit::Draw()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void AllPixTMPXDigit::Print()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
