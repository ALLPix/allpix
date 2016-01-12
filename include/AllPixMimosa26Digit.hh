/**
 *  Author:  John Idarraga <idarraga@cern.ch>
 */

#ifndef AllPixMimosa26Digit_h
#define AllPixMimosa26Digit_h 1

#include "AllPixDigitInterface.hh"
#include "G4TDigiCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"

/**
 *  Digit class for the Mimosa26 device.
 */

class AllPixMimosa26Digit : public AllPixDigitInterface {

public:
  AllPixMimosa26Digit();
  ~AllPixMimosa26Digit();

  AllPixMimosa26Digit(const AllPixMimosa26Digit&);
  const AllPixMimosa26Digit& operator=(const AllPixMimosa26Digit&);
  int operator==(const AllPixMimosa26Digit&) const;
  
  inline void* operator new(size_t);
  inline void  operator delete(void*);
  
  void Draw();
  void Print();

private:
  
  G4int pixelIDX;    
  G4int pixelIDY;    
  G4int pixelCounts;
  G4double pixelEnergyDep; // MC // Corrected MC charge (detector effects included, at Digi step)
  G4double pixelTruthEnergyDep; // MC // Truth energy (from Hits)
  G4ThreeVector primaryVertex;
  
public:
  
  inline void SetPixelIDX(G4int pidX)   {pixelIDX = pidX;};
  inline void SetPixelIDY(G4int pidY)   {pixelIDY = pidY;};
  inline void SetPixelCounts(G4int pc)  {pixelCounts = pc;};
  inline void SetPixelEnergyDep(G4double ed)  {pixelEnergyDep = ed;}; // MC // Corrected MC charge (detector effects included, at Digi step)
  inline void SetPixelTruthEnergyDep(G4double ed)  {pixelTruthEnergyDep = ed;}; // MC // Truth energy (from Hits)
  inline void SetPrimaryVertex(G4ThreeVector pv)  {primaryVertex = pv;}; // MC //

  inline void IncreasePixelCounts()     {pixelCounts++;};


  inline G4int GetPixelIDX()   {return pixelIDX;};
  inline G4int GetPixelIDY()   {return pixelIDY;};
  inline G4int GetPixelCounts()  {return pixelCounts;};
  inline G4double GetPixelEnergyDep()  {return pixelEnergyDep;}; // MC //
  inline G4ThreeVector GetPrimaryVertex()  {return primaryVertex;}; // MC //

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

typedef G4TDigiCollection<AllPixMimosa26Digit> AllPixMimosa26DigitsCollection;

extern G4Allocator<AllPixMimosa26Digit> AllPixMimosa26DigitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

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

#endif

