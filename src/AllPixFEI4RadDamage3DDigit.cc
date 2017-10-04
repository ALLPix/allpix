/**
 *  Author:
 *    Veronica Wallangen <veronica.wallangen@cern.ch>
 *
 *  allpix Authors:
 *   John Idarraga <idarraga@cern.ch>
 *   Mathieu Benoit <Mathieu.Benoit@cern.ch>
 */

#include "AllPixFEI4RadDamage3DDigit.hh"

G4Allocator<AllPixFEI4RadDamage3DDigit> AllPixFEI4RadDamage3DDigitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

AllPixFEI4RadDamage3DDigit::AllPixFEI4RadDamage3DDigit()
{
	pixelIDX = -1;
	pixelIDY = -1;
	pixelCounts = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

AllPixFEI4RadDamage3DDigit::~AllPixFEI4RadDamage3DDigit()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

AllPixFEI4RadDamage3DDigit::AllPixFEI4RadDamage3DDigit(const AllPixFEI4RadDamage3DDigit& right)
: AllPixDigitInterface()
{
	pixelIDX = right.pixelIDX;
	pixelIDY = right.pixelIDY;
	pixelCounts = right.pixelCounts;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

const AllPixFEI4RadDamage3DDigit& AllPixFEI4RadDamage3DDigit::operator=(const AllPixFEI4RadDamage3DDigit& right)
{

	pixelIDX = right.pixelIDX;
	pixelIDY = right.pixelIDY;
	pixelCounts = right.pixelCounts;
	return *this;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

int AllPixFEI4RadDamage3DDigit::operator==(const AllPixFEI4RadDamage3DDigit& right) const
		{
	return ((pixelIDX==right.pixelIDX)&&
			(pixelIDY==right.pixelIDY)&&
			(pixelCounts==right.pixelCounts));
		}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void AllPixFEI4RadDamage3DDigit::Draw()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void AllPixFEI4RadDamage3DDigit::Print()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
