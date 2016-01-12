/**
 *  Author:
 *    idarraga@cern.ch
 *
 *  allpix Authors:
 *   John Idarraga <idarraga@cern.ch>
 *   Mathieu Benoit <benoit@lal.in2p3.fr>
 */

#ifndef AllPixMedipix3RXDigit_h
#define AllPixMedipix3RXDigit_h 1

#include "AllPixDigitInterface.hh"
#include "G4TDigiCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"

#include <vector>
using namespace std;


/**
 *  Digit AllPixMedipix3RX class.
 */
class AllPixMedipix3RXDigit : public AllPixDigitInterface {

public:

	AllPixMedipix3RXDigit();
	AllPixMedipix3RXDigit(int nThresholds);
	~AllPixMedipix3RXDigit();

	AllPixMedipix3RXDigit(const AllPixMedipix3RXDigit&);
	const AllPixMedipix3RXDigit& operator=(const AllPixMedipix3RXDigit&);
	int operator==(const AllPixMedipix3RXDigit&) const;

	inline void* operator new(size_t);
	inline void  operator delete(void*);

	void Draw();
	void Print();

	void IncreasePixelCountsMultiTHL(G4int);

private:

	G4int pixelIDX;
	G4int pixelIDY;
	vector<G4int> pixelCounts;
	int nthresholds;
	G4double pixelEnergyDep; // MC // Corrected MC charge (detector effects included, at Digi step)
	G4ThreeVector primaryVertex;

public:

	inline void SetPixelIDX(G4int pidX)   {pixelIDX = pidX;};
	inline void SetPixelIDY(G4int pidY)   {pixelIDY = pidY;};
	inline void SetPixelCounts(G4int){};
	inline void SetPixelCounts(G4int thl, G4int pc) {
		if(thl >= 0 && thl <= nthresholds ) {
			pixelCounts[thl] = pc;
		}
	};
	inline void SetPixelEnergyDep(G4double ed)  {pixelEnergyDep = ed;}; // MC // Corrected MC charge (detector effects included, at Digi step)
	inline void SetPrimaryVertex(G4ThreeVector pv)  {primaryVertex = pv;}; // MC vertex //
	inline void IncreasePixelCounts(){};
	/*
	inline void IncreasePixelCounts(int thl)  {
		if( thl >= 0 && thl <= nthresholds ) {
			pixelCounts[thl]++;
		}
	};
*/
	inline G4int GetPixelIDX()   {return pixelIDX;};
	inline G4int GetPixelIDY()   {return pixelIDY;};
	inline G4int GetPixelCounts() {return 0;};
	/*
	inline G4int GetPixelCounts(G4int thl)  {
		if(thl >= 0 && thl <= nthresholds ) {
			return pixelCounts[thl];
		}
		return 0;
	};
	*/
	inline G4double GetPixelEnergyDep()  {return pixelEnergyDep;}; // MC //
	inline G4ThreeVector GetPrimaryVertex()  {return primaryVertex;}; // MC //

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

typedef G4TDigiCollection<AllPixMedipix3RXDigit> AllPixMedipix3RXDigitsCollection;

extern G4Allocator<AllPixMedipix3RXDigit> AllPixMedipix3RXDigitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void* AllPixMedipix3RXDigit::operator new(size_t)
{
	void* aDigi;
	aDigi = (void*) AllPixMedipix3RXDigitAllocator.MallocSingle();
	return aDigi;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void AllPixMedipix3RXDigit::operator delete(void* aDigi)
{
	AllPixMedipix3RXDigitAllocator.FreeSingle((AllPixMedipix3RXDigit*) aDigi);
}

#endif

