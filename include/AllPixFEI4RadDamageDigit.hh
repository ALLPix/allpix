/**
 *  Author:
 *    bnachman@cern.ch
 *
 *  allpix Authors:
 *   John Idarraga <idarraga@cern.ch>
 *   Mathieu Benoit <benoit@lal.in2p3.fr>
 */

#ifndef AllPixFEI4RadDamageDigit_h
#define AllPixFEI4RadDamageDigit_h 1

#include "AllPixDigitInterface.hh"
#include "G4TDigiCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"

/**
 *  Digit AllPixFEI4RadDamage class.
 */
class AllPixFEI4RadDamageDigit : public AllPixDigitInterface {

public:
  AllPixFEI4RadDamageDigit();
  ~AllPixFEI4RadDamageDigit();

  AllPixFEI4RadDamageDigit(const AllPixFEI4RadDamageDigit&);
  const AllPixFEI4RadDamageDigit& operator=(const AllPixFEI4RadDamageDigit&);
  int operator==(const AllPixFEI4RadDamageDigit&) const;
  
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

typedef G4TDigiCollection<AllPixFEI4RadDamageDigit> AllPixFEI4RadDamageDigitsCollection;

extern G4Allocator<AllPixFEI4RadDamageDigit> AllPixFEI4RadDamageDigitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void* AllPixFEI4RadDamageDigit::operator new(size_t)
{
  void* aDigi;
  aDigi = (void*) AllPixFEI4RadDamageDigitAllocator.MallocSingle();
  return aDigi;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void AllPixFEI4RadDamageDigit::operator delete(void* aDigi)
{
  AllPixFEI4RadDamageDigitAllocator.FreeSingle((AllPixFEI4RadDamageDigit*) aDigi);
}

#endif

