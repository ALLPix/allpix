/**
 *  Author:
 *    Paul Schuetze <paul.schuetze@desy.de>
 *
 *  allpix Authors:
 *   John Idarraga <idarraga@cern.ch>
 *   Mathieu Benoit <benoit@lal.in2p3.fr>
 */

#ifndef AllPixCMSp1Digit_h
#define AllPixCMSp1Digit_h 1

#include "AllPixDigitInterface.hh"
#include "G4TDigiCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"

/**
 *  Digit AllPixCMSp1 class.
 */
class AllPixCMSp1Digit : public AllPixDigitInterface {

public:
  AllPixCMSp1Digit();
  ~AllPixCMSp1Digit();

  AllPixCMSp1Digit(const AllPixCMSp1Digit&);
  const AllPixCMSp1Digit& operator=(const AllPixCMSp1Digit&);
  int operator==(const AllPixCMSp1Digit&) const;
  
  inline void* operator new(size_t);
  inline void  operator delete(void*);
  
  void Draw();
  void Print();

private:
  
  G4int pixelIDX;    
  G4int pixelIDY;    
  G4int pixelCounts;
  G4double pixelEnergyDep; // MC // Corrected MC charge (detector effects included, at Digi step)
  G4ThreeVector primaryVertex;
  
public:
  
  inline void SetPixelIDX(G4int pidX)   {pixelIDX = pidX;};
  inline void SetPixelIDY(G4int pidY)   {pixelIDY = pidY;};
  inline void SetPixelCounts(G4int pc)  {pixelCounts = pc;};
  inline void SetPixelEnergyDep(G4double ed)  {pixelEnergyDep = ed;}; // MC // Corrected MC charge (detector effects included, at Digi step)
  inline void SetPrimaryVertex(G4ThreeVector pv)  {primaryVertex = pv;}; // MC vertex //
  inline void IncreasePixelCounts()     {pixelCounts++;};

  inline G4int GetPixelIDX()   {return pixelIDX;};
  inline G4int GetPixelIDY()   {return pixelIDY;};
  inline G4int GetPixelCounts()  {return pixelCounts;};
  inline G4double GetPixelEnergyDep()  {return pixelEnergyDep;}; // MC //
  inline G4ThreeVector GetPrimaryVertex()  {return primaryVertex;}; // MC //

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

typedef G4TDigiCollection<AllPixCMSp1Digit> AllPixCMSp1DigitsCollection;

extern G4Allocator<AllPixCMSp1Digit> AllPixCMSp1DigitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void* AllPixCMSp1Digit::operator new(size_t)
{
  void* aDigi;
  aDigi = (void*) AllPixCMSp1DigitAllocator.MallocSingle();
  return aDigi;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void AllPixCMSp1Digit::operator delete(void* aDigi)
{
  AllPixCMSp1DigitAllocator.FreeSingle((AllPixCMSp1Digit*) aDigi);
}

#endif

