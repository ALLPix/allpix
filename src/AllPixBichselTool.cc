#include "AllPixBichselTool.hh"
#include <fstream>
#include <cmath>
#include <vector>
#include "TMath.h"
#include "TRandom.h"
#include <iostream>
#include "AllPixGeoDsc.hh"
#include "AllPixTrackerHit.hh"
/************************************************
The class AllpixBichselTool simulate the interactions between particles and silicon detector using Bichsel model.

Input: 
  1. n_Col: number of collisions that are grouped together as 1 collision. Large n_Col will save computation time
  2. AllpixTrackerHitsCollection from Geant4
  3. MatchNbofhits(bool): whether to keep the number of steps the same as Geant4
******************************************************/

AllPixBichselTool::AllPixBichselTool(){}

AllPixBichselTool::~AllPixBichselTool(){;}

void AllPixBichselTool::Initialize(G4int n_Col, G4double inpixelx, G4double inpixely, G4double inpixelz){
  //tablepath="/afs/cern.ch/user/f/fuwang/public/Allpix/data/Bichsel_%d%s.dat";
  tablepath="/afs/cern.ch/user/b/bnachman/work/public/RadDamage/newversion/allpix/share.r812211/Bichsel_%d%s.dat";
  //AllPixGeoDsc * gD = GetDetectorGeoDscPtr();
  pixelx=inpixelx;//gD->GetPixelX();//mm
  pixely=inpixely;//gD->GetPixelY();//mm
  pixelz=inpixelz;//gD->GetPixelZ();//mm
  m_nCols=n_Col;
  for(Int_t ParticleType=6;ParticleType<7;ParticleType++){
    m_BichselData.push_back(BichselDataTransform(ParticleType));
  }
}

std::vector<BichselHits> AllPixBichselTool::GetBichselhitsCollection(AllPixTrackerHitsCollection* hitsCollection,bool MatchNbofhits){
  /*
    Output: BichselHitsCollection ---- contains every hit simulated from Bichsel model 
    Input: 1. AllpixTrackerHitsCollection from Geant4
           2. MatchNbofhits: whether to keep the number of steps the same as Geant4
  */

  AllPixTrackerHit* hit0=(*hitsCollection)[0]; /* Initial hit information*/

  BichselHits Bichselhit;
  std::vector<BichselHits> BichselhitsCollection;
  G4ThreeVector endpos,dir;
  G4int nEntries = hitsCollection->entries();
  G4int nEntriestmp;//Number of primary hits, exclude secondaries

  for(G4int itr  = 0 ; itr < nEntries ; itr++) {
    if((*hitsCollection)[itr]->GetTrackID()<2) endpos=(*hitsCollection)[itr]->GetPos();//post position
    else{
      nEntriestmp--;
      //Secondary particle information: directly taking from Geant4
      Bichselhit.pos=(*hitsCollection)[itr]->GetPosWithRespectToPixel();
      Bichselhit.PixelNb.first=(*hitsCollection)[itr]->GetPixelNbX();
      Bichselhit.PixelNb.second=(*hitsCollection)[itr]->GetPixelNbY();
      Bichselhit.energy=(*hitsCollection)[itr]->GetEdep();
      Bichselhit.trackid=(*hitsCollection)[itr]->GetTrackID();
      BichselhitsCollection.push_back(Bichselhit);
    }
  }
  dir=endpos-hit0->GetPos();
  G4double TotalLength=dir.getR()*1e3;//um
  
  /********************* Simulation with Bichsel model *****************/

  G4double m=105.66;//MeV pion139.6,muon 105.66
  G4double InciEnergy=(*hitsCollection)[0]->GetKinEParent()/1e3+m;//MeV
  G4double BetaGamma=sqrt(pow(InciEnergy,2)-m*m)/m;//momentum/mass
  G4int ParticleType=trfPDG(hit0->GetTrackPdgId());

  std::vector<std::pair<G4double,G4double> > rawHitRecord;
  std::vector<std::pair<G4double,G4double> > trfHitRecord;
  rawHitRecord=BichselSim(BetaGamma,ParticleType,TotalLength,InciEnergy);//um,eV
 
  if(MatchNbofhits)//If true, keep the number of steps the same as Geant4
    trfHitRecord=ClusterHits(rawHitRecord,nEntriestmp);
  else trfHitRecord=rawHitRecord;
  /*********************************************************************/


  /********************* Simulation of hit position *****************/
  G4ThreeVector posStart_WithRespectToPixel =hit0->GetPosWithRespectToPixel();

  for(G4int hiti=0;hiti<trfHitRecord.size();hiti++){
    //z direction
    Bichselhit.pos.setZ(posStart_WithRespectToPixel.z()+trfHitRecord[hiti].first/TotalLength*dir.z());//mm
    //x direction
    G4double pos_respect_inipixelcenter=posStart_WithRespectToPixel.x()+trfHitRecord[hiti].first/TotalLength*dir.x();//With respect to the center of the initial pixel
    G4int deltaNb=floor((pos_respect_inipixelcenter-pixelx/2.)/pixelx)+1;
    Bichselhit.PixelNb.first=deltaNb+hit0->GetPixelNbX();
    Bichselhit.pos.setX(pos_respect_inipixelcenter-(G4double)deltaNb*pixelx);
    //y direction
    pos_respect_inipixelcenter=posStart_WithRespectToPixel.y()+trfHitRecord[hiti].first/TotalLength*dir.y();//With respect to the center of the initial pixel
    deltaNb=floor((pos_respect_inipixelcenter-pixely/2.)/pixely)+1;
    Bichselhit.PixelNb.second=deltaNb+hit0->GetPixelNbY();
    Bichselhit.pos.setY(pos_respect_inipixelcenter-(G4double)deltaNb*pixely);
    //energy
    Bichselhit.energy=trfHitRecord[hiti].second/1e6;//MeV
    //trackID
    Bichselhit.trackid=1;

    BichselhitsCollection.push_back(Bichselhit);
  }
  /*********************************************************************/
  return BichselhitsCollection;
}

std::vector<std::pair<G4double,G4double> > AllPixBichselTool::BichselSim(G4double BetaGamma, Int_t ParticleType, G4double TotalLength, G4double InciEnergy){ 
  TotalBichselEnergy=0;
  BichselData iData = m_BichselData[0];//[ParticleType-1]
  std::vector<std::pair<G4double,G4double> > rawHitRecord;  
  G4double accumLength = 0.;
  G4double BetaGammaLog10 = TMath::Log10(BetaGamma);
  std::pair<G4int,G4int> indices_BetaGammaLog10 = GetBetaGammaIndices(BetaGammaLog10, iData);

  // upper bound
  G4double IntXUpperBound = GetUpperBound(indices_BetaGammaLog10, BetaGammaLog10, iData);
  //check negative upper bound
  if(IntXUpperBound <= 0.){
    G4cout<<"Negative IntXUpperBound in BichselSimTool::BichselSim! (-1,-1) will be returned"<<G4endl;
    return rawHitRecord;
  }

  //mean free path
  G4double lambda = (1./IntXUpperBound) * 1.E4;   //um, unit of IntX is cm-1. It needs to be converted to micrometer-1
  // check nan lambda
  if(std::isnan(lambda)){
    G4cout<<"lambda is nan"<<G4endl;
    return rawHitRecord;
  } 


  //Set the maximum steps
  G4int LoopLimit=10000;
  if(fabs(1.0*TotalLength/lambda) > LoopLimit) LoopLimit=fabs(1.0*TotalLength/lambda); 

  // begin simulation
  G4int count = 0;
  gRandom->SetSeed();
  while(true){
    // infinite loop protection
    if(count >= (1.0*LoopLimit/m_nCols)){
      G4cout<<"Potential infinite loop in BichselSim. Exit Loop. A special flag will be returned (-1,-1). The total length is " << TotalLength << ". The lambda is " << lambda << "."<<G4endl;
      break;
    }

    // sample hit position -- exponential distribution
    G4double HitPosition = 0.;
    for(G4int iHit = 0; iHit < m_nCols; iHit++) HitPosition += gRandom->Exp(lambda); 
    // termination by hit position
    if(accumLength + HitPosition >= TotalLength) break;
   
    // sample single collision
    G4double TossEnergyLoss = -1.;
    while(TossEnergyLoss <= 0.){ // we have to do this because sometimes TossEnergyLoss will be negative due to too small TossIntX
      G4double TossIntX = gRandom->Uniform( 0., IntXUpperBound);
      TossEnergyLoss = GetColE(indices_BetaGammaLog10, TMath::Log10(TossIntX), iData);
    }
  
    bool fLastStep = true;
    // in case energy loss so far is larger than incident energy -- basically not possible ...
    if( ((TotalBichselEnergy + TossEnergyLoss)/1.E+6) > InciEnergy ){
      G4cout<<"Energy loss is larger than incident energy in AllpixBichselTool::BichselSim! This is usually delta-ray."<<G4endl;
      // then this is the last step
      TossEnergyLoss = InciEnergy*1.E+6 -TotalBichselEnergy;
      fLastStep = true;
    }

    // update
    accumLength += HitPosition;
    TotalBichselEnergy += TossEnergyLoss;
 
    // record this hit
    std::pair<G4double,G4double> oneHit;
    if(m_nCols == 1)  oneHit.first = accumLength; 
    else              oneHit.first = (accumLength - 1.0*HitPosition/2);     // as long as m_nCols is small enough (making sure lambda*m_nCols is withint resolution of a pixel), then taking middle point might still be reasonable
    oneHit.second = TossEnergyLoss;
    rawHitRecord.push_back(oneHit);

    count++;

    if(fLastStep)  break;
  }
   return rawHitRecord;
}

// assume vec is already sorted from small to large
std::pair<int,int> AllPixBichselTool::FastSearch(std::vector<double> vec, double item) {
  std::pair<int,int> output;

  int index_low = 0;
  int index_up = vec.size()-1;

  if((item < vec[index_low]) || (item > vec[index_up])){
    output.first = -1; output.second = -1;
    return output;
  }
  else if(item == vec[index_low]){
    output.first = index_low; output.second = index_low;
    return output;
  }
  else if(item == vec[index_up]){
    output.first = index_up; output.second = index_up;
    return output;
  }

  while( (index_up - index_low) != 1 ){
    int index_middle = int(1.0*(index_up + index_low)/2.);
    if(item < vec[index_middle])//
      index_up = index_middle;
    else if(item > vec[index_middle])
      index_low = index_middle;
    else{ // accurate hit. Though this is nearly impossible ...
      output.first = index_middle; output.second = index_middle;
      return output;
    }
  }

  output.first = index_low; output.second = index_up;
  return output;
}

double AllPixBichselTool::GetColE(std::pair<int,int> indices_BetaGammaLog10, double IntXLog10, BichselData& iData) {
  // std::pair<int,int> indices_BetaGammaLog10;
  // if(BetaGammaLog10 > iData.Array_BetaGammaLog10.back()){ // last one is used because when beta-gamma is very large, energy deposition behavior is very similar
  //   indices_BetaGammaLog10.first = iData.Array_BetaGammaLog10.size()-1;
  //   indices_BetaGammaLog10.second = iData.Array_BetaGammaLog10.size()-1;
  // }
  // else{
  //   indices_BetaGammaLog10 = FastSearch(iData.Array_BetaGammaLog10, BetaGammaLog10);
  // }

  if( (indices_BetaGammaLog10.first==-1) && (indices_BetaGammaLog10.second==-1) )
    return -1.;

  // BetaGammaLog10_2 then
  std::pair<int,int> indices_IntXLog10_x2 = FastSearch(iData.Array_BetaGammaLog10_IntXLog10[indices_BetaGammaLog10.second], IntXLog10);
  if (indices_IntXLog10_x2.first<0)  { return -1; }
  if (indices_IntXLog10_x2.second<0) { return -1; }
  double y21 = iData.Array_BetaGammaLog10_IntXLog10[indices_BetaGammaLog10.second][indices_IntXLog10_x2.first];
  double y22 = iData.Array_BetaGammaLog10_IntXLog10[indices_BetaGammaLog10.second][indices_IntXLog10_x2.second];
  double Est_x2 = ((y22 - IntXLog10)*iData.Array_BetaGammaLog10_ColELog10[indices_BetaGammaLog10.second][indices_IntXLog10_x2.first] + (IntXLog10 - y21)*iData.Array_BetaGammaLog10_ColELog10[indices_BetaGammaLog10.second][indices_IntXLog10_x2.second])/(y22-y21);

  // final estimation
  //double Est = ((BetaGammaLog10_2 - BetaGammaLog10)*Est_x1 + (BetaGammaLog10 - BetaGammaLog10_1)*Est_x2)/(BetaGammaLog10_2 - BetaGammaLog10_1);
  double Est = Est_x2;

  return TMath::Power(10., Est);
}



double AllPixBichselTool::GetColE(double BetaGammaLog10, double IntXLog10, BichselData& iData) {
  std::pair<int,int> indices_BetaGammaLog10 = GetBetaGammaIndices(BetaGammaLog10, iData);
  return GetColE(indices_BetaGammaLog10, IntXLog10, iData);
}

// IMPORTANT!! For this one, one should use interpolation, instead of fixed beta-gamma.
// Otherwise, dE/dx shape will get distorted again.
double AllPixBichselTool::GetUpperBound(std::pair<int,int> indices_BetaGammaLog10, double BetaGammaLog10, BichselData& iData) {
  // std::pair<int,int> indices_BetaGammaLog10;
  // if(BetaGammaLog10 > iData.Array_BetaGammaLog10.back()){
  //   indices_BetaGammaLog10.first = iData.Array_BetaGammaLog10.size()-1;
  //   indices_BetaGammaLog10.second = iData.Array_BetaGammaLog10.size()-1;
  // }
  // else{
  //   indices_BetaGammaLog10 = FastSearch(iData.Array_BetaGammaLog10, BetaGammaLog10);
  // }

  if (indices_BetaGammaLog10.first<0)  { return -1; }
  if (indices_BetaGammaLog10.second<0) { return -1; }
  double BetaGammaLog10_1 = iData.Array_BetaGammaLog10[indices_BetaGammaLog10.first];
  double BetaGammaLog10_2 = iData.Array_BetaGammaLog10[indices_BetaGammaLog10.second];

  // obtain estimation
  double Est_1 = iData.Array_BetaGammaLog10_UpperBoundIntXLog10[indices_BetaGammaLog10.first];
  double Est_2 = iData.Array_BetaGammaLog10_UpperBoundIntXLog10[indices_BetaGammaLog10.second];

  // final estimation
  double Est = ((BetaGammaLog10_2 - BetaGammaLog10)*Est_1 + (BetaGammaLog10 - BetaGammaLog10_1)*Est_2)/(BetaGammaLog10_2 - BetaGammaLog10_1);

  return TMath::Power(10., Est);
}


Double_t AllPixBichselTool::GetUpperBound(Double_t BetaGammaLog10, BichselData& iData) {
  std::pair<int,int> indices_BetaGammaLog10 = GetBetaGammaIndices(BetaGammaLog10, iData);
  return GetUpperBound(indices_BetaGammaLog10, BetaGammaLog10, iData);
}

std::pair<int,int> AllPixBichselTool::GetBetaGammaIndices(Double_t BetaGammaLog10, BichselData& iData) {
  std::pair<int,int> indices_BetaGammaLog10;
  if(BetaGammaLog10 > iData.Array_BetaGammaLog10.back()){ // last one is used because when beta-gamma is very large, energy deposition behavior is very similar
    indices_BetaGammaLog10.first = iData.Array_BetaGammaLog10.size()-1;
    indices_BetaGammaLog10.second = iData.Array_BetaGammaLog10.size()-1;
  }
  else{
    indices_BetaGammaLog10 = FastSearch(iData.Array_BetaGammaLog10, BetaGammaLog10);//找到iData.Array_BetaGammaLog10数组中等于BetaGammaLog10的那个的indices，iData.Array_BetaGammaLog10是一个每个元素都按从小到大排列的数组，
  }
  return indices_BetaGammaLog10;
}

//read from the Bichsel Model Table, store the data in the BichselData.
BichselData AllPixBichselTool::BichselDataTransform(G4int iParticleType){
  std::ifstream inputFile;
  TString inputFileName = TString::Format(tablepath, iParticleType, m_nCols == 1 ? "" : TString::Format("_%dsteps", m_nCols).Data());
  TString FullFileName=inputFileName;
  inputFile.open(FullFileName);
  BichselData iData;    
  G4double BetaGammaLog10 = 0; inputFile >> BetaGammaLog10;
  G4double ColELog10 = 0;      inputFile >> ColELog10;
  G4double IntXLog10 = 0;      inputFile >> IntXLog10;

  while(!inputFile.eof()){
    // check if this BetaGamma has already been stored
    if((iData.Array_BetaGammaLog10.size()==0)||(iData.Array_BetaGammaLog10.back()!=BetaGammaLog10)){ // a new BetaGamma
      if(iData.Array_BetaGammaLog10.size() != 0){//is not first betagamma
	iData.Array_BetaGammaLog10_UpperBoundIntXLog10.push_back(iData.Array_BetaGammaLog10_IntXLog10.back().back());
      }
      iData.Array_BetaGammaLog10.push_back(BetaGammaLog10);
      std::vector<G4double> new_ColELog10;  
      iData.Array_BetaGammaLog10_ColELog10.push_back(new_ColELog10);
      std::vector<G4double> new_IntXLog10;  
      iData.Array_BetaGammaLog10_IntXLog10.push_back(new_IntXLog10);
    }
    iData.Array_BetaGammaLog10_ColELog10.back().push_back(ColELog10);
    iData.Array_BetaGammaLog10_IntXLog10.back().push_back(IntXLog10);

    inputFile >> BetaGammaLog10;
    inputFile >> ColELog10;
    inputFile >> IntXLog10;
  }
  iData.Array_BetaGammaLog10_UpperBoundIntXLog10.push_back(iData.Array_BetaGammaLog10_IntXLog10.back().back());
  return iData;
}



G4int AllPixBichselTool::trfPDG(int pdgId) {
  if(std::fabs(pdgId) == 2212) return 1;   // proton
  if(std::fabs(pdgId) == 211)  return 2;   // pion
  // alpha is skipped -- 3
  if(std::fabs(pdgId) == 11)   return 4;   // electron
  if(std::fabs(pdgId) == 321)  return 5;   // kaon
  if(std::fabs(pdgId) == 13)   return 6;   // muon

  return -1;   // unsupported particle
}


//cluster thousand hits into ~20 groups
std::vector<std::pair<double,double> > AllPixBichselTool::ClusterHits(std::vector<std::pair<double,double> >& rawHitRecord, int n_pieces){
 
  std::vector<std::pair<double,double> > trfHitRecord;

  if((int)(rawHitRecord.size()) < n_pieces){ // each single collision is the most fundamental unit
    n_pieces = rawHitRecord.size();
  }

  int unitlength = int(1.0*rawHitRecord.size()/n_pieces);
  int index_start = 0;
  int index_end = unitlength-1;   // [index_start, index_end] are included
  while(true){
    // calculate weighted center of each slice
    double position = 0.;
    double energyloss = 0.;

    for(int index = index_start; index <= index_end; index++){
      position += (rawHitRecord[index].first * rawHitRecord[index].second);
      energyloss += rawHitRecord[index].second;
    }
    position = (energyloss == 0. ? 0. : position/energyloss);

    // store
    std::pair<double,double> oneHit;
    oneHit.first = position; oneHit.second = energyloss;
    trfHitRecord.push_back(oneHit);

    // procede to next slice
    index_start = index_end + 1;
    index_end = index_start + unitlength - 1;

    if(index_start > (int)(rawHitRecord.size()-1)){
      break;
    }

    if(index_end > (int)(rawHitRecord.size()-1)){
      index_end = rawHitRecord.size()-1;
    }
  }
  return trfHitRecord;
}
