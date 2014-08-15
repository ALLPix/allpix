/**
 *  Author:  John Idarraga <idarraga@cern.ch>
 */

#ifndef AllPixFEI3StandardDigit_h
#define AllPixFEI3StandardDigit_h 1

#include "G4VDigi.hh"
#include "AllPixDigitInterface.hh"
#include "G4TDigiCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"

/**
 *  Digit class for the FEI3Standard device.
 */

class AllPixFEI3StandardDigit : public AllPixDigitInterface {

public:

  AllPixFEI3StandardDigit();
  virtual ~AllPixFEI3StandardDigit();

  AllPixFEI3StandardDigit(const AllPixFEI3StandardDigit &);
  const AllPixFEI3StandardDigit& operator=(const AllPixFEI3StandardDigit&);
  int operator==(const AllPixFEI3StandardDigit&) const;
  
  inline void* operator new(size_t);
  inline void  operator delete(void*);
  
  void Draw();
  void Print();

private:
  
  G4int pixelIDX;    
  G4int pixelIDY;    
  G4int pixelCounts;
  G4double pixelEnergyDep; // MC //
  G4ThreeVector primaryVertex;
  G4double pixelEnergyMC; //MC only
  
  G4double posX_WithRespectToPixel; //nalipour
  G4double posY_WithRespectToPixel; //nalipour
  G4double posZ_WithRespectToPixel; //nalipour

public:
  
  inline void SetPixelIDX(G4int pidX)   {pixelIDX = pidX;};
  inline void SetPixelIDY(G4int pidY)   {pixelIDY = pidY;};
  inline void SetPixelCounts(G4int pc)  {pixelCounts = pc;};
  inline void SetPixelEnergyDep(G4double ed)  {pixelEnergyDep = ed;}; // MC //
  inline void SetPrimaryVertex(G4ThreeVector pv)  {primaryVertex = pv;}; // MC //
  inline void SetPixelEnergyMC(G4double energyMC) {pixelEnergyMC = energyMC;};

  inline void IncreasePixelCounts()     {pixelCounts++;};


  inline G4int GetPixelIDX()   {return pixelIDX;};
  inline G4int GetPixelIDY()   {return pixelIDY;};
  inline G4int GetPixelCounts()  {return pixelCounts;};
  inline G4double GetPixelEnergyDep()  {return pixelEnergyDep;}; // MC //
  inline G4ThreeVector GetPrimaryVertex()  {return primaryVertex;}; // MC //
  inline G4double GetPixelEnergyMC() {return pixelEnergyMC;};



  inline void Set_posX_WithRespectoToPixel(G4double pos) {posX_WithRespectToPixel=pos;}; //nalipour
  inline void Set_posY_WithRespectoToPixel(G4double pos) {posY_WithRespectToPixel=pos;}; //nalipour
  inline void Set_posZ_WithRespectoToPixel(G4double pos) {posZ_WithRespectToPixel=pos;}; //nalipour
  inline G4double Get_posX_WithRespectoToPixel() {return posX_WithRespectToPixel;}; //nalipour
  inline G4double Get_posY_WithRespectoToPixel() {return posY_WithRespectToPixel;}; //nalipour
  inline G4double Get_posZ_WithRespectoToPixel() {return posZ_WithRespectToPixel;}; //nalipour

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

typedef G4TDigiCollection<AllPixFEI3StandardDigit> AllPixFEI3StandardDigitsCollection;

extern G4Allocator<AllPixFEI3StandardDigit> AllPixFEI3StandardDigitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void* AllPixFEI3StandardDigit::operator new(size_t)
{
  void* aDigi;
  aDigi = (void*) AllPixFEI3StandardDigitAllocator.MallocSingle();
  return aDigi;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void AllPixFEI3StandardDigit::operator delete(void* aDigi)
{
  AllPixFEI3StandardDigitAllocator.FreeSingle((AllPixFEI3StandardDigit*) aDigi);
}

#endif

