// @(#)codeTDv2Final/CTDR_Animation :$Name:  $:$Id: CTDR_Animation.h,v 1.7 2006/11/25 16:10:44 Benoit Exp $
// Author: Mathieu Benoit   24/11/06

//////////////////////////////////////////////////////////////////////////
// Source file : class CTDR_Animation					//
// This class is used to produce a root file containing the tracks of  	//
// the charge elements simulated. It creates a geometry, defined in the //
// the constructor, and everything needed to draw image and animation   //
// of the simulated event in the geometry . The root file can be read  	//
// using cint interpreter, and the geometry and the charge tracks can	//
// then be manipulated using functions provided with root's geometry  	//
// packages , see TGeoManager and TGeoTracks classes.	 		//
//////////////////////////////////////////////////////////////////////////

#ifndef ALLPIXDIGITANIMATION_H_DEF
#define ALLPIXDIGITANIMATION_H_DEF

#include <vector>
#include "TGeoManager.h"
#include "TVirtualGeoTrack.h"
#include "TFile.h"
#include "TGeoBBox.h"
#include "TString.h"
using namespace std;



class AllPixDigitAnimation : public TObject
{
      
      private :
             TGeoManager *Geo; // The geometry manager
             TFile *f; // The file where the objects are written
             int Ntracks; // The number of tracks
             TGeoMedium *medium; // the material
             TGeoVolume *top ; // the box
             TGeoBBox *PixBox; // The pixels
             TGeoVolume *PixVolume; //the pixels volume
             
             
             //double Lx,Ly,Lz;  // size of the box
	     	 int trackid;
	    	 int m_nx,m_ny;
	    	 double m_Lz;
	    	 double shiftx,shifty;
	   		 double m_pitchx,m_pitchy;
         	 Int_t MyPalette[100];
         	 
         	 double emax;
         	 
         	 double z_hit;
             
      public :
             AllPixDigitAnimation(int nx, int ny, double lz, double pitchx, double pitchy,int nHits, int eventid);
             void AddTrack(vector<double> x, vector<double> y, vector<double> z, vector<double> t,double sigma,double energy);
	     	 double modulo(double a , double b);
	    	 void FirePixel(int i , int j);	
		 void SubThresholdPixel(int i , int j);
	         void SetShift(double x, double y){
	     		shiftx=x;
	     		shifty=y;
	     };
	     
             ~AllPixDigitAnimation();
	     //ClassDef(AllPixDigitAnimation,1) // A geometry and event track visualization class
	     
      
};
      
#endif            
