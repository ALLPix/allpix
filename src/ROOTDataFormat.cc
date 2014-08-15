#include "ROOTDataFormat.hh"

ROOTDataFormat::ROOTDataFormat(Int_t detId)
{
  detectorID=detId;
}

ROOTDataFormat::~ROOTDataFormat()
{
}

/*
void ROOTDataFormat::set_posX_MC(vector<Int_t> vec)
{
  posX_MC=vec;
  nHits_MC=posX_MC.size();
}

void ROOTDataFormat::set_posY_MC(vector<Int_t> vec)
{
  posY_MC=vec;
}

void ROOTDataFormat::set_energy_MC(vector<Double_t> vec)
{
  energy_MC=vec;
}


vector<Int_t> ROOTDataFormat::get_posX_MC()
{
  return posX_MC;
}
vector<Int_t> ROOTDataFormat::get_posY_MC()
{
  return posY_MC;
}
vector<Double_t> ROOTDataFormat::get_energy_MC()
{
  return energy_MC;
}

void ROOTDataFormat::add_posX_MC(Int_t pos)
{
  posX_MC.push_back(pos);
}
void ROOTDataFormat::add_posY_MC(Int_t pos)
{
  posY_MC.push_back(pos);
}
void ROOTDataFormat::add_energy_MC(Double_t energy)
{
  energy_MC.push_back(energy);
}


void ROOTDataFormat::add_posX(Int_t pos)
{
	posX.push_back(pos);
}
void ROOTDataFormat::add_posY(Int_t pos)
{
	posY.push_back(pos);
}
void ROOTDataFormat::add_energy(Double_t e)
{
	energy.push_back(e);
}
void ROOTDataFormat::add_TOT(Int_t t)
{
	TOT.push_back(t);
}

vector<Int_t> ROOTDataFormat::get_posX()
{
  return posX;
}
vector<Int_t> ROOTDataFormat::get_posY()
{
  return posY;
}
vector<Double_t> ROOTDataFormat::get_energy()
{
  return energy;
}
vector<Int_t> ROOTDataFormat::get_TOT()
{
  return TOT;
}
*/
