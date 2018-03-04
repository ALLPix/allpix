/**
 *  Author:  John Idarraga <idarraga@cern.ch>
 */

#ifndef AllPixDigitInterface_h
#define AllPixDigitInterface_h 1

#include "G4VDigi.hh"
#include "G4TDigiCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"

/**
 *  Interface Digit class.
 */

class AllPixDigitInterface : public G4VDigi {

public:

	AllPixDigitInterface();
	virtual ~AllPixDigitInterface();

	// This one can't be pure virtual
	AllPixDigitInterface(const AllPixDigitInterface &) : G4VDigi() { } ;

	// Don't need abstraction for the following
	//const AllPixDigitInterface & operator=(const AllPixDigitInterface &) { return *this; };
	//int operator==(const AllPixDigitInterface &) const;
	//virtual void* operator new(size_t);
	//virtual void  operator delete(void*);

	virtual void Draw() = 0;
	virtual void Print() = 0;

	virtual void SetPixelIDX(G4int){;};
	virtual void SetPixelIDY(G4int){;};
	virtual void SetPixelCounts(G4int){;};	//TOT (nalipour)
	virtual void SetPixelEnergyMC(G4double); //MC only
	virtual void SetPixelEnergyDep(G4double){;};     // MC + charge sharing (nalipour)
	virtual void SetPixelCountsMultiTHL(G4int, G4int);

	virtual void SetPixelEnergyDepMultiTHL(G4double, G4double);     // MC

	virtual G4int GetPixelCountsMultiTHL(G4int);

	virtual void IncreasePixelCounts() = 0;
	virtual void IncreasePixelCountsMultiTHL(G4int);

	virtual G4int GetPixelIDX(){return 0;} ;
	virtual G4int GetPixelIDY(){return 0;} ;
	virtual G4int GetPixelCounts(){return 0;} ; //TOT
	virtual G4double GetPixelEnergyDep(){return 0;} ;     // MC+charge sharing
	virtual G4double GetPixelEnergyMC() ;     //MC only (nalipour)
	virtual G4ThreeVector GetPrimaryVertex(){return G4ThreeVector() ;}; // MC //
  virtual G4double GetPOLangle() {return 0;};
  virtual G4double GetAZMangle() {return 0;};
  virtual G4double GetTruthEntryLocalX() {return 0;};
  virtual G4double GetTruthEntryLocalY() {return 0;};
  virtual G4double GetTruthExitLocalX() {return 0;};
  virtual G4double GetTruthExitLocalY() {return 0;};
  virtual G4double GetdeltaEfrac() {return 0;};
  virtual G4double Getpath_length_first_pixel() {return 0;};
  virtual G4double GetInitE() {return 0;};
  virtual G4int GetInitID() {return 0;};
	virtual void Set_posX_WithRespectoToPixel(G4double pos); //nalipour
	virtual void Set_posY_WithRespectoToPixel(G4double pos); //nalipour
	virtual void Set_posZ_WithRespectoToPixel(G4double pos); //nalipour
	virtual G4double Get_posX_WithRespectoToPixel(); //nalipour
	virtual G4double Get_posY_WithRespectoToPixel();//nalipour
	virtual G4double Get_posZ_WithRespectoToPixel();//nalipour

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

typedef G4TDigiCollection<AllPixDigitInterface> AllPixDigitsCollectionInterface;

extern G4Allocator<AllPixDigitInterface> AllPixDigitAllocatorInterface;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....



/*
inline void* AllPixMimosa26Digit::operator new(size_t)
{
	void* aDigi;
	aDigi = (void*) AllPixMimosa26DigitAllocator.MallocSingle();
	return aDigi;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void AllPixMimosa26Digit::operator delete(void* aDigi)
{
	AllPixMimosa26DigitAllocator.FreeSingle((AllPixMimosa26Digit*) aDigi);
}
*/
#endif

