/*
 * Generic output module for AllPix
 * Virtual class to inherit and overload functions from (see e.g. AllPixLCIOwriter)
 * 
 * Author: Paul Schuetze (paul.schuetze@desy.de)
 */

#pragma once

#ifndef AllPixWriter_h
#define AllPixWriter_h 1

#include "G4String.hh"

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>

using namespace std;


class AllPixWriter
{
public:

  virtual void Initialize(G4String folder, int runnr) = 0;

  virtual void WriteEvent(int runnr, int eventID, map<int,vector<vector<vector<int>>>> data) = 0;

  virtual void Close() = 0;

private:



};

#endif
