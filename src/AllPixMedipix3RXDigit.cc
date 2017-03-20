/**
 *  Author:
 *    idarraga@cern.ch
 *
 *  allpix Authors:
 *   John Idarraga <idarraga@cern.ch>
 *   Mathieu Benoit <Mathieu.Benoit@cern.ch>
 */

#include "AllPixMedipix3RXDigit.hh"

G4Allocator<AllPixMedipix3RXDigit> AllPixMedipix3RXDigitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

AllPixMedipix3RXDigit::AllPixMedipix3RXDigit()
{
	pixelIDX = -1;
	pixelIDY = -1;
	nthresholds = 8;
	pixelCounts.clear();
	kinEParent = 0.;

	for(int i = 0 ; i < nthresholds ; i++) pixelCounts.push_back( 0 );


}

AllPixMedipix3RXDigit::AllPixMedipix3RXDigit(int nThresholds)
{
	pixelIDX = -1;
	pixelIDY = -1;
	nthresholds = nThresholds;
	pixelCounts.clear();

	for(int i = 0 ; i < nthresholds ; i++) pixelCounts.push_back( 0 );


}

void AllPixMedipix3RXDigit::IncreasePixelCountsMultiTHL(G4int thl){
	pixelCounts[thl]++;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

AllPixMedipix3RXDigit::~AllPixMedipix3RXDigit()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

AllPixMedipix3RXDigit::AllPixMedipix3RXDigit(const AllPixMedipix3RXDigit& right)
: AllPixDigitInterface()
{
	pixelIDX = right.pixelIDX;
	pixelIDY = right.pixelIDY;
	pixelCounts = right.pixelCounts;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

const AllPixMedipix3RXDigit& AllPixMedipix3RXDigit::operator=(const AllPixMedipix3RXDigit& right)
{

	pixelIDX = right.pixelIDX;
	pixelIDY = right.pixelIDY;
	pixelCounts = right.pixelCounts;
	return *this;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

int AllPixMedipix3RXDigit::operator==(const AllPixMedipix3RXDigit& right) const
		{
	return ((pixelIDX==right.pixelIDX)&&
			(pixelIDY==right.pixelIDY)&&
			(pixelCounts==right.pixelCounts));
		}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void AllPixMedipix3RXDigit::Draw()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void AllPixMedipix3RXDigit::Print()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
