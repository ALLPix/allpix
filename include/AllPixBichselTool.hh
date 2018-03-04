#ifndef AllPixBichselTool_h
#define AllPixBichselTool_h 1
#include <vector>
#include "AllPixTrackerHit.hh"
//#include "AllPixDigitizerInterface.hh"

struct BichselData
{
  std::vector<double> Array_BetaGammaLog10;
  std::vector<std::vector<double> > Array_BetaGammaLog10_ColELog10;  // ColE = CollisionEnergy in eV
  std::vector<std::vector<double> > Array_BetaGammaLog10_IntXLog10;  // IntX = Integrated Xsection. The unit doesn't matter
  std::vector<double> Array_BetaGammaLog10_UpperBoundIntXLog10;      // upper bound of log10(IntX)
};


struct BichselHits
{
  G4ThreeVector pos;
  std::pair<G4int,G4int> PixelNb;
  G4double energy;
  G4int trackid;
};

class AllPixBichselTool
{

 public:
  AllPixBichselTool();
  virtual ~AllPixBichselTool();
  std::vector<std::pair<double,double> > BichselSim(double BetaGamma, int ParticleType, double TotalLength, double InciEnergy);
  double GetTotalBichselEnergy(){return TotalBichselEnergy;};
  int trfPDG(int pdgId);
  std::vector<std::pair<double,double> > ClusterHits(std::vector<std::pair<double,double> >& rawHitRecord, int n_pieces);
  void Initialize(G4int n_Col, G4double inpixelx, G4double inpixely, G4double inpixelz);
  std::vector<BichselHits> GetBichselhitsCollection(AllPixTrackerHitsCollection* hitsCollection,bool MatchNbofhits);


 private:
  std::pair<int,int> FastSearch(std::vector<double> vec, double item);
  G4double GetUpperBound(G4double BetaGammaLog10, BichselData& iData);
  double GetUpperBound(std::pair<int,int> indices_BetaGammaLog10, double BetaGammaLog10, BichselData& iData);
  double GetColE(std::pair<int,int> indices_BetaGammaLog10, double IntXLog10, BichselData& iData);
  double GetColE(double BetaGammaLog10, double IntXLog10, BichselData& iData) ;
  std::pair<int,int> GetBetaGammaIndices(G4double BetaGammaLog10, BichselData& iData);
  BichselData BichselDataTransform(G4int iParticleType);

  G4String tablepath;
  G4double pixelx,pixely,pixelz;
  G4double TotalBichselEnergy;//eV
  std::vector<BichselData> m_BichselData;   
  int                      m_nCols;            // number of collisions to simulate each time. This is mainly to save CPU time if necessary
 
};
#endif
