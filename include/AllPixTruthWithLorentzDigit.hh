/**
 *  Author:
 *    Ben Nachman <bnachman@cern.ch>
 *
 *  allpix Authors:
 *   John Idarraga <idarraga@cern.ch>
 *   Mathieu Benoit <benoit@lal.in2p3.fr>
 */

#ifndef AllPixTruthWithLorentzDigit_h
#define AllPixTruthWithLorentzDigit_h 1

#include "AllPixDigitInterface.hh"
#include "G4TDigiCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"

/**
 *  Digit AllPixTruthWithLorentz class.
 */
class AllPixTruthWithLorentzDigit : public AllPixDigitInterface {

public:
  AllPixTruthWithLorentzDigit();
  ~AllPixTruthWithLorentzDigit();

  AllPixTruthWithLorentzDigit(const AllPixTruthWithLorentzDigit&);
  const AllPixTruthWithLorentzDigit& operator=(const AllPixTruthWithLorentzDigit&);
  int operator==(const AllPixTruthWithLorentzDigit&) const;
  
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
  G4double POL;
  G4double AZM;
  G4double TruthEntryLocalX;
  G4double TruthEntryLocalY;
  G4double TruthExitLocalX;
  G4double TruthExitLocalY;
  G4double deltaEfrac;
  G4double path_length_first_pixel;
  G4double InitE;
  G4int InitID;

public:
  
  inline void SetPixelIDX(G4int pidX)   {pixelIDX = pidX;};
  inline void SetPixelIDY(G4int pidY)   {pixelIDY = pidY;};
  inline void SetPixelCounts(G4int pc)  {pixelCounts = pc;};
  inline void SetPOLangle(G4double pol) {POL = pol;};
  inline void SetAZMangle(G4double azm) {AZM = azm;};
  inline void SetPixelEnergyDep(G4double ed)  {pixelEnergyDep = ed;}; // MC // Corrected MC charge (detector effects included, at Digi step)
  inline void SetPrimaryVertex(G4ThreeVector pv)  {primaryVertex = pv;}; // MC vertex //
  inline void IncreasePixelCounts()     {pixelCounts++;};
  inline void SetTruthEntryLocalX(G4double pos) {TruthEntryLocalX = pos;};
  inline void SetTruthEntryLocalY(G4double pos) {TruthEntryLocalY = pos;};
  inline void SetTruthExitLocalX(G4double pos) {TruthExitLocalX = pos;};
  inline void SetTruthExitLocalY(G4double pos) {TruthExitLocalY = pos;};
  inline void SetdeltaEfrac(G4double ed) {deltaEfrac = ed;};
  inline void Setpath_length_first_pixel(G4double pos) {path_length_first_pixel = pos;};
  inline void SetInitE(G4double pos) {InitE = pos;};
  inline void SetInitID(G4int pos) {InitID = pos;};
  
  inline G4int GetPixelIDX()   {return pixelIDX;};
  inline G4int GetPixelIDY()   {return pixelIDY;};
  inline G4int GetPixelCounts()  {return pixelCounts;};
  inline G4double GetPixelEnergyDep()  {return pixelEnergyDep;}; // MC //
  inline G4ThreeVector GetPrimaryVertex()  {return primaryVertex;}; // MC //
  inline G4double GetPOLangle() {return POL;};
  inline G4double GetAZMangle() {return AZM;};
  inline G4double GetTruthEntryLocalX() {return TruthEntryLocalX;};
  inline G4double GetTruthEntryLocalY() {return TruthEntryLocalY;};
  inline G4double GetTruthExitLocalX() {return TruthExitLocalX;};
  inline G4double GetTruthExitLocalY() {return TruthExitLocalY;};
  inline G4double GetdeltaEfrac() {return deltaEfrac;};
  inline G4double Getpath_length_first_pixel() {return path_length_first_pixel;};
  inline G4double GetInitE() {return InitE;};
  inline G4int GetInitID() {return InitID;};

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

typedef G4TDigiCollection<AllPixTruthWithLorentzDigit> AllPixTruthWithLorentzDigitsCollection;

extern G4Allocator<AllPixTruthWithLorentzDigit> AllPixTruthWithLorentzDigitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void* AllPixTruthWithLorentzDigit::operator new(size_t)
{
  void* aDigi;
  aDigi = (void*) AllPixTruthWithLorentzDigitAllocator.MallocSingle();
  return aDigi;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void AllPixTruthWithLorentzDigit::operator delete(void* aDigi)
{
  AllPixTruthWithLorentzDigitAllocator.FreeSingle((AllPixTruthWithLorentzDigit*) aDigi);
}

#endif

