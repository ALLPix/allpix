// @(#)codeTDv2Final/AllPixDigitAnimation :$Name:  $:$Id: AllPixDigitAnimation.h,v 1.7 init:2006/11/25 16:10:44 Benoit Exp $
// Author: Mathieu Benoit   rev. s26/04/2012

/*************************************************************************
 * . Class AllPixDigitAnimation					                 *
 *   Written by : Mathieu Benoit                                         *
 *   mathieu.benoit@cern.ch						 *
 *************************************************************************/

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <string>
#include "AllPixDigitAnimation.hh"
#include "G4SystemOfUnits.hh"
#include "TString.h"
#include "TGeoMatrix.h"
#include "TColor.h"
#include "TMath.h"
using namespace std;

//////////////////////////////////////////////////////////////////////////
// This class is used to produce a root file containing the tracks of  	//
// the charge elements simulated. It creates a geometry, defined in the //
// the constructor, and everything needed to draw image and animation   //
// of the simulated event in the geometry . The root file can be read  	//
// using cint interpreter, and the geometry and the charge tracks can	//
// then be manipulated using functions provided with root's geometry  	//
// packages , see TGeoManager and TGeoTracks classes.	 		//
//////////////////////////////////////////////////////////////////////////

ClassImp(AllPixDigitAnimation)

//______________________________________________________________________________
AllPixDigitAnimation::AllPixDigitAnimation(int nx, int ny, double Lz, double pitchx, double pitchy,int nHits, int eventid)
{
// This class builds a geometry and tracks object to put simulated data
// Lx,Ly,Lz are the dimensions of the box, inter, the interaction to be simulated, i and j ae ID for the
// root file name
             		

                        
            f= new TFile(TString::Format("AllPixDigitAnimation_%d.root",eventid),"RECREATE");
                        
            Geo= new TGeoManager(TString::Format("Tracks_%d",eventid),"Digitization graphical representation");
             
            emax = 2; 
            
            shiftx=0;
            shifty=0;
            
            z_hit=Lz/2+2.5;
            
   			Double_t r[]    = {0., 0.0, 1.0, 1.0, 1.0};
  			Double_t g[]    = {0., 0.0, 0.0, 1.0, 1.0};
   			Double_t b[]    = {0., 1.0, 0.0, 0.0, 1.0};
  			Double_t stop[] = {0., .25, .50, .75, 1.0};
   			Int_t FI = TColor::CreateGradientColorTable(5, stop, r, g, b, 100);
   			for (int i=0;i<100;i++) MyPalette[i] = FI+i;
              
            
              
                        
			this->nx=nx;
			this->ny=ny;
			this->Lz=Lz;
			this->pitchx=pitchx;
			this->pitchy=pitchy;
			
			TGeoMedium *medium=0; 
            top = Geo->MakeBox("TOP",medium,nx*pitchx/2,ny*pitchy/2,2*Lz/2);
            Geo->SetTopVolume(top);
            
            TGeoVolume *Sensor =  Geo->MakeBox("TOP",medium,nx*pitchx/2,ny*pitchy/2,Lz/2);
            Sensor->SetLineColor(kGray+1);
            top->AddNode(Sensor,1);
            
                        
			Double_t ori[3]={0,0,0};

                        trackid=0;
                        
            for(int i=0;i<nx;i++){
                    for(int j=0;j<ny;j++){
                            
                            
							ori[0]=-(nx-1)*pitchx/2 + (i)*pitchx; 
                            ori[1]=-(ny-1)*pitchy/2 + (j)*pitchy;
                            ori[2]=Lz/2 +0.5;
                            PixBox= new TGeoBBox(TString::Format("Pix_%i_%i",i,j),(pitchx-0.1*pitchx)/2,(pitchy-0.1*pitchy)/2,1,ori);
                            
                            PixVolume =new  TGeoVolume(TString::Format("Pix_%i_%i",i,j),PixBox);
                            PixVolume->SetLineColor(kYellow -3);
                            top->AddNode(PixVolume,1);
                            };};


             Geo->CloseGeometry();
             top->SetLineColor(kMagenta);
             top->SetVisibility(kTRUE);
    		 top->Draw("ogle");
    		 
             Geo->SetTopVisible();
             for(int i=0;i<nHits;i++){
                          Geo->AddTrack(i,100);
                          //Geo->GetCurrentTrack()->SetLineWidth(1);
                          //Geo->GetCurrentTrack()->SetMarkerSize(4);
                          };
		  
                                  
}

//______________________________________________________________________________
void AllPixDigitAnimation::AddTrack(vector<double> x, vector<double> y, vector<double> z, vector<double> t,double sigma, double energy)
{
// This method adds a step to the animation by reading the data from the simulated interaction

     if(x.size()>0){
     //cout  << "[animation] " << sigma << " " << t[t.size()-1] << endl;
     
     //cout << TString::Format("[Animation] Hit Position : x:%f y:%f z:%f",x[0],y[0],z[0]) << endl;

     
//     if(x[0]>nx*pitchx/(2)) {
//     		shiftx= (x[0] - modulo(x[0],nx*pitchx));
//		};
//	 
//     if(x[0]<-nx*pitchx/(2)){
//     		shiftx = -(x[0]-modulo(-x[0],nx*pitchx));
//		};
//     if(y[0]>ny*pitchy/(2)) {
//     		shifty= (y[0] - modulo(y[0],ny*pitchy));
//		};
//	 
//     if(y[0]<-ny*pitchy/(2)){
//     		shifty = -(y[0]-modulo(-y[0],ny*pitchy));
//		};		
     
			
    //cout << TString::Format("[Animation] shift x = %f shift y = %f",shiftx,shifty) << endl;		
		

     for(uint i=0;i<x.size();i++){
             
             Geo->GetTrack(trackid)->AddPoint(x[i]-shiftx,y[i]-shifty,-z[i],t[i]);
         
                  };
		
      
      double deltaX = x[x.size()-1]-x[0];
      double deltaY = y[y.size()-1]-y[0];
      double deltaZ = z[z.size()-1]-z[0];
      
      cout << "[Animation] Track angle = " << TMath::ATan(TMath::Sqrt(deltaX*deltaX+deltaY*deltaY)/deltaZ)*TMath::RadToDeg() << endl;
      
      
      
      
      TGeoVolume *tube ;
      if (sigma!=0) {
      	tube = Geo->MakeTube(TString::Format("hit%d",trackid),0,0,2*sigma,0.01);
      	}
      else {
      	tube = Geo->MakeTube(TString::Format("hit%d",trackid),0,0,0.2,0.01);
      }	
      
      int color = TMath::FloorNint(100*(energy/keV)/emax);
      if(color>99)color=99;
      
      //cout << "color = " << color << endl;
      //tube->SetFillColor(kBlue);
  	  tube->SetLineColor(8);
  
      
      TGeoHMatrix *matrix= new TGeoHMatrix();
      Double_t ori[3]={x[x.size()-1]-shiftx,y[x.size()-1]-shifty,z_hit};
      z_hit+=0.01;
      
      matrix->SetDx(ori[0]);
      matrix->SetDy(ori[1]);
      matrix->SetDz(ori[2]);
      top->AddNode(tube,1,matrix);
      }
      trackid++;
}

AllPixDigitAnimation::~AllPixDigitAnimation(){

    
	
	Geo->Write();
	f->Close();
	delete Geo;


}



double AllPixDigitAnimation::modulo(double a, double b){

	int result = static_cast<int>( a / b );
	return a - static_cast<double>( result ) * b;
}


void AllPixDigitAnimation::SubThresholdPixel(int i , int j){
	
	TString pixelName = TString::Format("Pix_%d_%d",i+nx/2-1,j+ny/2-1);
	
	cout << "[Fired]" << pixelName << endl;
	
	
	if (i+nx/2-1<nx && j+ny/2-1<ny &&i+nx/2-1 > 0 &&j+ny/2-1 > 0 ){

	
	Geo->GetVolume(pixelName)->SetLineColor(kRed);
    Geo->GetVolume(pixelName)->SetLineWidth(3);
	}
		
}


void AllPixDigitAnimation::FirePixel(int i , int j){
	
	TString pixelName = TString::Format("Pix_%d_%d",i+nx/2-1,j+ny/2-1);
	
	cout << "[Fired]" << pixelName << endl;
	
	
	if (i+nx/2-1<nx && j+ny/2-1<ny &&i+nx/2-1 > 0 &&j+ny/2-1 > 0 ){

	
	Geo->GetVolume(pixelName)->SetLineColor(kBlue);
    Geo->GetVolume(pixelName)->SetLineWidth(3);
	}
		
}
                          
