/**
 *  Author John Idarraga <idarraga@cern.ch>
 */

#ifndef allpix_dm_h
#define allpix_dm 1

#include <iostream>
#include <map>
#include <vector>
#include <stdlib.h>
#include "TROOT.h"
#include "TObject.h"
#include "TH2.h"
#include "TVector3.h"

#include "allpix_dm_consts.h"

/**
 * Contains only the pixels matrix
 *  and some members that may go to the
 *  ROOT file.
 */

class FrameContainer {

private:

	// Zero suppressed, X,counts or X,TOT
	std::map<int, int> m_frameXC;
	// Trigger
	std::map<int, int> m_lvl1;
	// FIXME -->  I am unable to deserialize here <int, double>
	// Truth energy (from Hits)
	std::map<int, int> m_frameXC_TruthE;
	// Corrected MC charge (detector effects included, at Digitization step)
	std::map<int, int> m_frameXC_E;


	//nalipour
	std::map<std::pair<Int_t, Int_t>, Double_t > RD53_E_MC; //True Monte Carlo energy
	std::map<std::pair<Int_t, Int_t>, Double_t > RD53_E;    //Energy taking into account the 

	// Counters
	// The number of pixels ON
	Int_t m_nEntriesPad;
	// The number of hits in the frame
	//  not necessarily the same number of entries
	//  e.x. if exposition time is large enough
	Int_t m_nHitsInPad;
	// The total TOT in the frame
	Int_t m_nChargeInPad;

	// If data is MC this flag is set to true
	Bool_t m_isMCData;

public:
	FrameContainer();
	virtual ~FrameContainer(){};
	void FillOneElement(Int_t, Int_t, Int_t, Int_t, Double_t, Double_t);
	void FillOneElement(Int_t, Int_t, Int_t, Int_t);
	void FillOneElementRD53_MC(Int_t xi, Int_t yi, Double_t truthE); //nalipour //MC
	void FillOneElementRD53(Int_t xi, Int_t yi, Double_t E); //nalipour //Considering charge sharing
	void SetLVL1(Int_t, Int_t, Int_t, Int_t);

	void ResetCountersPad();
	void CleanUpMatrix();
	void SetFrameAsMCData(){m_isMCData = true;};
	Int_t GetEntriesPad(){return m_nEntriesPad;};
	Int_t GetHitsInPad(){return m_nHitsInPad;};
	Int_t GetChargeInPad(){return m_nChargeInPad;};
	std::map<std::pair<Int_t, Int_t>, Double_t > return_RD53_E_MC(){return RD53_E_MC;}; //nalipour
	std::map<std::pair<Int_t, Int_t>, Double_t > return_RD53_E(){return RD53_E;}; //nalipour

	ClassDef(FrameContainer,3)
};


/** 
 * Contains the meta data and an object for the frame matrix.  It
 *  inherits the methods to fill up the frameMatrix from class
 *  FrameContainer
 */

class FrameStruct : public FrameContainer {

private:

	/* head info */
	Int_t    fFormat;
	Int_t    fWidth;
	Int_t    fHeight;

	/* all metadata */
	Int_t    fAcq_mode;
	Double_t fAcq_time;
	TString  fApplied_filters;
	Double_t fAuto_erase_interval;
	Int_t    fAutoerase_interval_counter;
	Bool_t   fBS_active;
	TString  fChipboardID;
	Double_t fCoinc_live_time;
	Byte_t   fCoincidence_delay;
	Byte_t   fCoincidence_mode;
	std::vector<Int_t> fCounters;
	std::vector<Int_t> fDACs;
	Double_t fHV;
	Int_t    fHw_timer;
	TString  fInterface;
	Double_t fMpx_clock;
	Int_t    fMpx_type;
	Int_t    fPolarity;
	Double_t fStart_time;
	TString  fStart_timeS;
	Byte_t   fTimepix_clock;
	Double_t fTrigger_time;

	/* MC only */
	std::vector<Double_t> m_primaryVertex_x;
	std::vector<Double_t> m_primaryVertex_y;
	std::vector<Double_t> m_primaryVertex_z;

	/* joint */
	Long_t   fFrameId;
	TString  fMPXDataSetNumber;

public:
	FrameStruct(TString);
	~FrameStruct(){};

	// filling Tree methods
	void FillMetaData(TString, Int_t);
	void IncreaseId(){ fFrameId++; }
	void SetId(Int_t id){ fFrameId=id; }
	void SetPrimaryVertex(Double_t x, Double_t y, Double_t z); // pushing back
	TString GetDataSet(){ return fMPXDataSetNumber; };
	Int_t GetFrameId(){return fFrameId;};
	void RewindMetaDataValues();
	void SetnX(int x){fWidth = x;};
	void SetnY(int y){fHeight = y;};

	ClassDef(FrameStruct,4)
};


/** 
 * Containst a FrameStruct and implements the methods needed to read
 *  the frame files, load up information to the frame.
 */

class FramesHandler {

private:

	Int_t m_nFrames;
	FrameStruct * m_aFrame;
	Int_t m_width; // replicate this information for speed purposes
	Int_t m_height; // replicate this information for speed purposes

	Bool_t getAFrameMatrix_flag;
	Bool_t getAFrameHist_flag;

	Int_t m_metaBit;
	Int_t m_metaCode;
	Bool_t m_ParseAndHold;

	Int_t nFrames256x256;
	Int_t nFramesXYC;

	// in case I am dealing with more than one detector
	Int_t m_detID;

public:

	FramesHandler(TString);
	~FramesHandler();
	Int_t IdentifyTypeOfInput(TString);

	/*****************************************************
	 *  Filling frame methods
	 */
	/* Reads from txt & dsc files.  Completely fills m_aFrame. */
	Bool_t readOneFrame(TString, TString);
	/* load a single frame pixel (X,Y,C), fills m_aFrame */
	Bool_t LoadFramePixel(Int_t, Int_t, Int_t, Double_t, Double_t);
	/* load a single frame pixel (X,Y,C), fills m_aFrame */
	Bool_t LoadFramePixel(Int_t, Int_t, Int_t);
	Bool_t LoadFramePixelRD53_MC(Int_t col, Int_t row, Double_t truthE);//nalipour
	Bool_t LoadFramePixelRD53(Int_t col, Int_t row, Double_t E);//nalipour
	/* load a lvl1 trigger */
	void SetLVL1(Int_t, Int_t, Int_t);
	/* load a single frame MetaData entry (METAString, Int_t metaCode)*/
	Bool_t LoadFrameMetaData(TString, Int_t);
	/* Load primary vertex info */
	Bool_t LoadPrimaryVertexInfo(Double_t vx, Double_t vy, Double_t vz){
		m_aFrame->SetPrimaryVertex(vx, vy, vz);
		return true;
	};
	void IncreaseCurrentFrameId();
	void SetCurrentFrameId(Int_t id){m_aFrame->SetId(id);};
	void SetAsMCData(){m_aFrame->SetFrameAsMCData();};
	void SetnX(int x){m_width = x; m_aFrame->SetnX(x);};
	void SetnY(int y){m_height = y; m_aFrame->SetnY(y);};

	/* *************************************************** */

	Int_t GetDetectorId(){return m_detID;};
	void  SetDetectorId(Int_t id){m_detID = id;};
	void RewindAll();
	//TH2I * getAFrameHist(TString, TString, TString);
	Int_t ** getAFrameMatrix(TString, TString);
	//TH2I * getHistFrame(Int_t, Int_t *);
	FrameStruct * getFrameStructObject(){return m_aFrame;};
	Bool_t frameSupervisor();
	void parseMetaLine(TString, Int_t);
	void push_back_nbytes(unsigned int *, char *, Int_t);

};

#endif
