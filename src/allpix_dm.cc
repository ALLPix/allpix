/**
 *  Author: John Idarraga <idarraga@cern.ch>
 *  Classes containing the frames
 */

#include <fstream>
#include <istream>
#include <iostream> //nalipour nilou for tests
#include <string>
#include "TH2.h"

#include "allpix_dm.h"

ClassImp(FrameContainer)
ClassImp(FrameStruct)

FrameContainer::FrameContainer(){

	m_nEntriesPad = 0;
	m_nHitsInPad = 0;
	m_nChargeInPad = 0;
	m_isMCData = false;

}

void FrameContainer::FillOneElement(Int_t xi, Int_t yi, Int_t width, Int_t counts) {

	// X,Y,C --> X,C : yi*width + xi
	Int_t X = yi*width + xi;

	m_frameXC[X] += counts;        // TOT or count(binary detector)

	// If the pixel didn't exist this is an extra entry
	if(m_frameXC.find(X) == m_frameXC.end()) m_nEntriesPad++;
	// But always an extra hit
	m_nHitsInPad++;
	// Increase the total counts
	m_nChargeInPad += counts;

}

void FrameContainer::FillOneElement(Int_t xi, Int_t yi, Int_t width, Int_t counts, Double_t truthE, Double_t E){

	// Fill pixel withouth MC info first
	FillOneElement(xi, yi, width, counts);

	// X,Y,C --> X,C : yi*width + xi
	Int_t X = yi*width + xi;

	m_frameXC_TruthE[X] += truthE; // Truth energy
	m_frameXC_E[X] += E;           // Energy with detector effects

	
//	//nalipour
//	std::cout << "Nilou: xi=" << xi << ", yi=" << yi << ", m_frameXC_TruthE[X] =" << m_frameXC_TruthE[X] << ", m_frameXC_E[X] =" << m_frameXC_E[X] << std::endl;
//	std::pair<int, int> tempPixel;
//	tempPixel.first=xi;
//	tempPixel.second=yi;
//
//	RD53_E_MC[tempPixel]+=truthE;
//	RD53_E[tempPixel]+=E;
}

//nalipour
void FrameContainer::FillOneElementRD53_MC(Int_t xi, Int_t yi, Double_t truthE)
{
	std::pair<int, int> tempPixel;
	tempPixel.first=xi;
	tempPixel.second=yi;

	RD53_E_MC[tempPixel]+=truthE;
}
//nalipour
void FrameContainer::FillOneElementRD53(Int_t xi, Int_t yi, Double_t E)
{std::pair<int, int> tempPixel;
	tempPixel.first=xi;
	tempPixel.second=yi;

	RD53_E[tempPixel]+=E;
}

void FrameContainer::SetLVL1(Int_t xi, Int_t yi, Int_t width, Int_t lvl1){

	// X,Y,C --> X,C : yi*width + xi
	Int_t X = yi*width + xi;

	m_lvl1[X] += lvl1;
}

void FrameContainer::ResetCountersPad(){

	m_nEntriesPad = 0;
	m_nHitsInPad = 0;
	m_nChargeInPad = 0;

}

FrameStruct::FrameStruct(TString dataset) 
: FrameContainer() {

	// clean matrix, every pixel to 0
	CleanUpMatrix();
	// reset counters m_nEntriesPad, m_nHitsInPad, m_nChargeInPad
	ResetCountersPad();
	// reset metadata
	RewindMetaDataValues();

	/* joint */
	fFrameId = -1;
	fMPXDataSetNumber = dataset;

}

void FrameStruct::RewindMetaDataValues(){

	/* head info */
	fFormat = 0;
	fWidth = 0;
	fHeight = 0;

	/* all metadata */
	fAcq_mode = 0;
	fAcq_time = 0.0;
	fApplied_filters = "";
	fAuto_erase_interval = 0.0;
	fAutoerase_interval_counter = 0;
	fBS_active = false;
	fChipboardID = "";
	fCoinc_live_time = 0.0;
	fCoincidence_delay = 0;
	fCoincidence_mode = 0;
	fCounters.clear();
	fDACs.clear();
	fHV = 0.0;
	fHw_timer = 0;
	fInterface = "";
	fMpx_clock = 0.0;
	fMpx_type = -1;
	fPolarity = -1;
	fStart_time = 0.0;
	fStart_timeS = "";
	fTimepix_clock = 0;
	fTrigger_time = 0.0;

	/* Joint data doesn't need to be rewinded
	 *  --> fFrameId
	 *  --> fMPXDataSetNumber
	 */
}

void FrameStruct::SetPrimaryVertex(Double_t vx, Double_t vy, Double_t vz){

	m_primaryVertex_x.push_back(vx);
	m_primaryVertex_y.push_back(vy);
	m_primaryVertex_z.push_back(vz);

}

void FrameStruct::FillMetaData(TString METAString, Int_t metaCode){

	TString tempEntry = "";

	switch (metaCode)
	{
	case ACQ_MODE_CODE:
		fAcq_mode = METAString.Atoi();
		break;
	case ACQ_TIME_CODE:
		fAcq_time = METAString.Atof();
		break;
	case CHIPBOARDID_CODE:
		fChipboardID = METAString;
		break;
	case DACS_CODE:
		for(Int_t i = 0 ; i < METAString.Length() ; i++)
		{
			if(METAString[i] != ' ')
			{
				tempEntry += METAString[i];
			}
			else if(METAString[i] == ' ' && tempEntry.Length() > 0)
			{
				//if(fDACsSize < MAX_DAQ_INT){
				//std::cout << tempEntry << std::endl;
				//fDACs[fDACsSize] = tempEntry.Atoi();
				//fDACsSize++;
				fDACs.push_back((Int_t)tempEntry.Atoi());
				//}
				tempEntry = "";
			}
		}
		tempEntry = "";
		break;
	case HW_TIMER_CODE:
		fHw_timer = METAString.Atoi();
		break;
	case INTERFACE_CODE:
		fInterface = METAString;
		break;
	case POLARITY_CODE:
		fPolarity = METAString.Atoi();
		break;
	case START_TIME_CODE:
		fStart_time = METAString.Atof();
		break;
	case START_TIME_S_CODE:
		fStart_timeS = METAString;
		break;
	case MPX_CLOCK_CODE:
		fMpx_clock = METAString.Atof();
		break;
	case HV_CODE:
		fHV = METAString.Atof();
		break;
	case TIMEPIX_CLOCK_CODE:
		fTimepix_clock = (Byte_t)METAString.Atoi();
		break;
	default:
		std::cout << "[WARNING] couldn't parse MetaData Info" << std::endl;
		break;
	}

	// FIXME ! .... have to decide how to give real dataSet numbers to data
	// PGP signature ?, md5 of the output file ? user decided ?
	// connect to db ?

	fMPXDataSetNumber = "1234ABCD";

}


FramesHandler::FramesHandler(TString dataset){

	//////////////////////////////////////////
	// Instanciate one frame
	// is going to be overwritten many times
	// but instanciated only once !
	m_aFrame = new FrameStruct(dataset);
	m_nFrames = 0;

	getAFrameMatrix_flag = false;
	getAFrameHist_flag = false;

	m_metaBit = -1;
	m_metaCode = -1;
	m_ParseAndHold = true;

	nFrames256x256 = 0;
	nFramesXYC = 0;

}

FramesHandler::~FramesHandler(){

	delete m_aFrame;
}

/*
TH2I * FramesHandler::getHistFrame(Int_t frameId, Int_t * frameMatrix){

	TString title = "frame ", hname = "frame_";
	title += frameId; hname += frameId;
	TH2I * histFrame = new TH2I(hname, title,
			MAX_FRAME_ROW, 0, MAX_FRAME_ROW-1,
			MAX_FRAME_COL, 0, MAX_FRAME_COL-1
	);

	Int_t iCntr = 0;
	Int_t jCntr = 0;
	for(Int_t i = 0 ; i < MAX_FRAME_ROW*MAX_FRAME_COL ; i++)
	{
		histFrame->Fill(jCntr, iCntr, frameMatrix[i]);
		iCntr++;
		if(iCntr >= MAX_FRAME_ROW)
		{
			iCntr = 0;
			jCntr++;
		}
		//std::cout << frameMatrix[i][j] << " ";
	}

	return histFrame;
}
 */
Int_t FramesHandler::IdentifyTypeOfInput(TString fullDSCFileName){

	fstream filestr;
	filestr.open(fullDSCFileName, fstream::in);

	Int_t lineCntr = 0;
	Char_t temp[2048];
	std::string tempS = "";
	std::string typeS = FRAME_TYPE_STRING;
	std::string binTypeS = BINARY_FRAME_TYPE_STRING;
	std::string fSave = "";
	size_t found = std::string::npos;

	while (filestr.good())// && lineCntr < 1)
	{

		filestr.getline(temp, 2048);
		tempS = temp;

		found = tempS.find(typeS);
		if (found != std::string::npos)
		{

			// found the format line, now pick up the format BITS
			for(int i = (int)typeS.length() ; i < (int)tempS.length() ; i++)
			{
				if(tempS[i] != ' ')
				{
					fSave.append(1,tempS[i]);
				}
				else
					break;
			}
		}

		/* check if binary */
		found = tempS.find(binTypeS);
		if (found != std::string::npos) {
			return FSAVE_BINARY;
		}

		lineCntr++;
	}

	filestr.close();

	return atoi(fSave.c_str());
}

void FrameContainer::CleanUpMatrix(){
	m_frameXC.clear();
	m_lvl1.clear();
	m_frameXC_TruthE.clear();
	m_frameXC_E.clear();
	
	//nalipour
	RD53_E_MC.clear();
	RD53_E.clear();
}

/* Rewind frame and metadata */
void FramesHandler::RewindAll(){

	// clean matrix, every pixel to 0
	m_aFrame->CleanUpMatrix();
	// reset counters m_nEntriesPad, m_nHitsInPad, m_nChargeInPad
	m_aFrame->ResetCountersPad();
	// reset metadata
	m_aFrame->RewindMetaDataValues();

}

Bool_t FramesHandler::readOneFrame(TString fullFileName, TString fullDSCFileName){

	RewindAll();
	m_aFrame->IncreaseId(); // the first time it'll come to 0 (initialized at -1)
	m_nFrames++;

	/* Identify the type of file, 256x256 or XYC */
	Int_t frameType = IdentifyTypeOfInput(fullDSCFileName);

	if(frameType == 0)
	{
		std::cout << "[ERROR] could not determine the frame type --> " << fullDSCFileName << std::endl;
		exit(1);
	}

	/* ******************* */
	/* open the frame file */
	TString typeS = "";
	fstream filestr;

	Int_t temp = 0;
	Long_t totalCounts = 0;
	Long_t pixelHits = 0;
	Int_t cntri = 0;
	Int_t cntrj = 0;
	Int_t wholePadCntr = 0;

	if(frameType == (FSAVE_ASCII | FSAVE_I16) ||
			frameType == (FSAVE_ASCII | FSAVE_U32) ||
			frameType == (FSAVE_ASCII | FSAVE_DOUBLE) )
	{
		filestr.open(fullFileName, fstream::in);

		typeS = TYPE_256x256_STRING;
		std::cout << "[INFO] opening: " << fullFileName << "  --|and|-->  " << fullDSCFileName
				<< " --> " << typeS << "(" << frameType << ")" << std::endl;

		while (filestr.good() && wholePadCntr < MAX_FRAME_OCC)
		{
			filestr >> temp;

			totalCounts += temp;
			pixelHits++;

			m_aFrame->FillOneElement(cntri, cntrj, m_width, temp);
			cntri++;

			if(cntri == MAX_FRAME_COL){
				cntrj++;
				cntri = 0;
			}
			wholePadCntr++;
		}
	}
	else if(frameType == (FSAVE_ASCII | FSAVE_I16 | FSAVE_SPARSEXY) ||
			frameType == (FSAVE_ASCII | FSAVE_U32 | FSAVE_SPARSEXY) ||
			frameType == (FSAVE_ASCII | FSAVE_DOUBLE | FSAVE_SPARSEXY) )
	{
		filestr.open(fullFileName, fstream::in);

		typeS = TYPE_XYC_STRING;
		std::cout << "[INFO] opening: " << fullFileName << "  --|and|-->  " << fullDSCFileName
				<< " --> " << typeS << "(" << frameType << ")" << std::endl;

		std::map<Int_t, Int_t> frameMap;
		for(Int_t itr = 0 ; itr < MAX_FRAME_OCC ; itr++)
		{
			frameMap[itr] = 0;
		}

		Int_t fillTime = 0;
		Int_t colr = 0, rowr = 0;
		while (filestr.good())
		{
			filestr >> temp;

			if(fillTime%3 == 0) /* got X */ /* There is a harmles extra
					     read at this point.  No risk.*/
			{
				colr = temp;
				//std::cout << "colr: " << colr << std::endl;
				fillTime++;
			}
			else if(fillTime%3 == 1) /* got Y */
			{
				rowr = temp;
				//std::cout << "rowr: " << rowr << std::endl;
				fillTime++;
			}
			else if(fillTime%3 == 2) /* got Counts */
			{
				//std::cout << "val: " << temp << std::endl;

				totalCounts += temp;
				pixelHits++;

				frameMap[rowr*256 + colr] = temp;
				fillTime = 0;
			}

			cntri++;

			if(cntri == MAX_FRAME_COL){
				cntrj++;
				cntri = 0;
			}

		}

		/* now read the dictionary and fill the matrix */
		std::map<Int_t, Int_t>::iterator it;
		cntri = 0;
		cntrj = 0;

		for ( it=frameMap.begin() ; it != frameMap.end(); it++ )
		{
			/*
	    if((*it).second > 0)
	    std::cout << (*it).first << " => " << (*it).second <<
	    std::endl;
			 */

			/* en el archivo : col, row */
			m_aFrame->FillOneElement(cntri, cntrj, m_width, (*it).second);
			cntri++;

			if(cntri == MAX_FRAME_COL){
				cntrj++;
				cntri = 0;
			}

			/*
	  if((*it).first == MAX_FRAME_OCC-1)
	    exit(1);
			 */
		}
	}
	else if(frameType == (FSAVE_ASCII | FSAVE_I16 | FSAVE_SPARSEX) ||
			frameType == (FSAVE_ASCII | FSAVE_U32 | FSAVE_SPARSEX) ||
			frameType == (FSAVE_ASCII | FSAVE_DOUBLE | FSAVE_SPARSEX) )
	{
		filestr.open(fullFileName, fstream::in);

		typeS = TYPE_XC_STRING;
		std::cout << "[INFO] opening: " << fullFileName << "  --|and|-->  " << fullDSCFileName
				<< " --> " << typeS << "(" << frameType << ")"<< std::endl;

		std::map<Int_t, Int_t> frameMap;
		for(Int_t itr = 0 ; itr < MAX_FRAME_OCC ; itr++)
		{
			frameMap[itr] = 0;
		}

		Int_t fillTime = 0;
		Int_t colr = 0, rowr = 0;
		while (filestr.good())
		{
			filestr >> temp;

			if(fillTime%2 == 0) /* got X */ /* There is a harmles extra
					     read at this point.  No risk.*/
			{
				colr = temp;
				//std::cout << "colr: " << colr << std::endl;
				fillTime++;
			}
			else if(fillTime%2 == 1) /* got Counts */
			{
				//std::cout << "val: " << temp << std::endl;

				totalCounts += temp;
				pixelHits++;

				frameMap[rowr*256 + colr] = temp;
				fillTime = 0;
			}

			cntri++;

			if(cntri == MAX_FRAME_COL){
				cntrj++;
				cntri = 0;
			}

		}

		/* now read the dictionary and fill the matrix */
		std::map<Int_t, Int_t>::iterator it;
		cntri = 0;
		cntrj = 0;

		for ( it=frameMap.begin() ; it != frameMap.end(); it++ )
		{
			/*
	    if((*it).second > 0)
	    std::cout << (*it).first << " => " << (*it).second <<
	    std::endl;
			 */

			/* en el archivo : col, row */
			m_aFrame->FillOneElement(cntri, cntrj, m_width, (*it).second);
			cntri++;

			if(cntri == MAX_FRAME_COL){
				cntrj++;
				cntri = 0;
			}

			/*
	  if((*it).first == MAX_FRAME_OCC-1)
	    exit(1);
			 */
		}
	}
	else if(frameType == FSAVE_BINARY) {

		typeS = TYPE_XYC_BIN_STRING;
		std::cout << "[INFO] opening: " << fullFileName << "  --|and|-->  " << fullDSCFileName
				<< " --> " << typeS << "(" << frameType << ")"<< std::endl;

		//filestr.open(fullFileName, fstream::binary);
		filestr.open(fullFileName, istream::in);

		unsigned int bytesRead = 0;
		char tempByte[4];
		unsigned int x = 0x0, y = 0x0, counts = 0x0;

		while (filestr.good())
		{

			// read X, 32 bits
			x = 0x0;
			filestr.get(tempByte[0]);
			filestr.get(tempByte[1]);
			filestr.get(tempByte[2]);
			filestr.get(tempByte[3]);
			push_back_nbytes(&x, tempByte, 4);

			// read Y, 32 bits
			y = 0x0;
			filestr.get(tempByte[0]);
			filestr.get(tempByte[1]);
			filestr.get(tempByte[2]);
			filestr.get(tempByte[3]);
			push_back_nbytes(&y, tempByte, 4);

			// Read Counts, 16 bits
			counts = 0x0;
			filestr.get(tempByte[0]);
			filestr.get(tempByte[1]);
			push_back_nbytes(&counts, tempByte, 2);

			//printf("%x, %x --> %x\n", x, y, counts);
			//printf("%d, %d --> %d\n", x, y, counts);

			bytesRead += 10;

			m_aFrame->FillOneElement(x, y, m_width, counts);

		}

	}
	else
	{
		/* Unknown case ... giving up here */
	}

	// FIXME ! // get this from the dsc file
	SetnX(256);
	SetnY(256);

	filestr.close();

	// Now fetching MetaData
	fstream META_filestr;
	Int_t nLineMetaFile = 0;
	META_filestr.exceptions ( ifstream::eofbit | ifstream::failbit | ifstream::badbit );

	// read metadata
	try {

		META_filestr.open(fullDSCFileName, fstream::in);

		Char_t * METAtemp;
		METAtemp = new char[META_DATA_LINE_SIZE];

		while (META_filestr.good()){
			META_filestr.getline(METAtemp, META_DATA_LINE_SIZE);
			if(m_ParseAndHold) parseMetaLine((TString)METAtemp, nLineMetaFile); // pass currentLine

			if(nLineMetaFile++ == m_metaBit){
				m_aFrame->FillMetaData((TString)METAtemp, m_metaCode);
				m_ParseAndHold = true;
			}
		}
	}
	catch (ifstream::failure &e){
		//std::cout << "[ERROR] Exception opening/reading file: " << fullDSCFileName << std::endl;
		//std::cout << "        probably file doesn't not exist.  giving up." << std::endl;
		//exit(1);
	}
	META_filestr.close();

	return true;
}

void FramesHandler::push_back_nbytes(unsigned int * val, char * bytes, Int_t nbytes) {

	// indexes go like this
	// 0 --> lower byte  0x......XX
	// 3 --> higher byte 0xXX......

	*val &= 0x00000000;
	unsigned int tempVal;

	for(int idx = nbytes-1 ; idx >= 0 ; idx--) {

		// Get the byte
		tempVal = 0x0;
		tempVal ^= bytes[idx];
		// Clean up to have info for only one byte
		tempVal &= 0x000000ff;
		// Switch it to the right place
		for(int sw = 0 ; sw < idx ; sw++){
			tempVal = tempVal << 8;
		}
		// XOR the value
		*val ^= tempVal;

	}

}

/* load a single frame pixel (X,Y,C) + truth info */
Bool_t FramesHandler::LoadFramePixel(Int_t col, Int_t row, Int_t counts, Double_t truthE, Double_t E){

	m_aFrame->FillOneElement(col, row, m_width, counts, truthE, E);

	return true;
}

//nalipour
Bool_t FramesHandler::LoadFramePixelRD53_MC(Int_t col, Int_t row, Double_t truthE){

	m_aFrame->FillOneElementRD53_MC(col, row, truthE);

	return true;
}
//nalipour
Bool_t FramesHandler::LoadFramePixelRD53(Int_t col, Int_t row, Double_t E){

	m_aFrame->FillOneElementRD53(col, row, E);

	return true;
}

/* load a single frame pixel (X,Y,C) */
Bool_t FramesHandler::LoadFramePixel(Int_t col, Int_t row, Int_t counts){

	m_aFrame->FillOneElement(col, row, m_width, counts);

	return true;
}


void FramesHandler::SetLVL1(Int_t col, Int_t row, Int_t lvl1){
	m_aFrame->SetLVL1(col, row, m_width, lvl1);
}

void FramesHandler::IncreaseCurrentFrameId(){
	m_aFrame->IncreaseId();
}

/* load a single frame MetaData entry (METAString, Int_t metaCode)*/
Bool_t FramesHandler::LoadFrameMetaData(TString metastring , Int_t metacode){

	m_aFrame->FillMetaData(metastring, metacode);

	return true;
}

void FramesHandler::parseMetaLine(TString aMetaLine, Int_t currentLine){

	// It is probably a good idea to change this to
	// "switch case" because it may be a bit faster
	if(aMetaLine.Contains(ACQ_MODE_STRING, TString::kExact)){
		//std::cout << "-> 1" << std::endl;
		m_metaBit = currentLine + ACTUAL_VALUE_OFFSET;
		m_metaCode = ACQ_MODE_CODE;
		m_ParseAndHold = false;
	}
	else if(aMetaLine.Contains(ACQ_TIME_STRING, TString::kExact)){
		//std::cout << "-> 2" << std::endl;
		m_metaBit = currentLine + ACTUAL_VALUE_OFFSET;
		m_metaCode = ACQ_TIME_CODE;
		m_ParseAndHold = false;
	}
	else if(aMetaLine.Contains(CHIPBOARDID_STRING, TString::kExact)){
		//std::cout << "-> 3" << std::endl;
		m_metaBit = currentLine + ACTUAL_VALUE_OFFSET;
		m_metaCode = CHIPBOARDID_CODE;
		m_ParseAndHold = false;
	}
	else if(aMetaLine.Contains(DACS_STRING, TString::kExact)){
		//std::cout << "-> 4" << std::endl;
		m_metaBit = currentLine + ACTUAL_VALUE_OFFSET;
		m_metaCode = DACS_CODE;
		m_ParseAndHold = false;
	}
	else if(aMetaLine.Contains(HW_TIMER_STRING, TString::kExact)){
		//std::cout << "-> 5" << std::endl;
		m_metaBit = currentLine + ACTUAL_VALUE_OFFSET;
		m_metaCode = HW_TIMER_CODE;
		m_ParseAndHold = false;
	}
	else if(aMetaLine.Contains(INTERFACE_STRING, TString::kExact)){
		//std::cout << "-> 6" << std::endl;
		m_metaBit = currentLine + ACTUAL_VALUE_OFFSET;
		m_metaCode = INTERFACE_CODE;
		m_ParseAndHold = false;
	}
	else if(aMetaLine.Contains(POLARITY_STRING, TString::kExact)){
		//std::cout << "-> 7" << std::endl;
		m_metaBit = currentLine + ACTUAL_VALUE_OFFSET;
		m_metaCode = POLARITY_CODE;
		m_ParseAndHold = false;
	}
	else if(aMetaLine.Contains(START_TIME_STRING, TString::kExact)){
		//std::cout << "-> 8" << std::endl;
		m_metaBit = currentLine + ACTUAL_VALUE_OFFSET;
		m_metaCode = START_TIME_CODE;
		m_ParseAndHold = false;
	}
	else if(aMetaLine.Contains(START_TIME_S_STRING, TString::kExact)){
		//std::cout << "-> 9" << std::endl;
		m_metaBit = currentLine + ACTUAL_VALUE_OFFSET;
		m_metaCode = START_TIME_S_CODE;
		m_ParseAndHold = false;
	}
	else if(aMetaLine.Contains(MPX_CLOCK_STRING, TString::kExact)){
		//std::cout << "-> 9" << std::endl;
		m_metaBit = currentLine + ACTUAL_VALUE_OFFSET;
		m_metaCode = MPX_CLOCK_CODE;
		m_ParseAndHold = false;
	}
	else if(aMetaLine.Contains(HV_STRING, TString::kExact)){
		m_metaBit = currentLine + ACTUAL_VALUE_OFFSET;
		m_metaCode = HV_CODE;
		m_ParseAndHold = false;
	}
	else if(aMetaLine.Contains(TIMEPIX_CLOCK_STRING, TString::kExact)){
		m_metaBit = currentLine + ACTUAL_VALUE_OFFSET;
		m_metaCode = TIMEPIX_CLOCK_CODE ;
		m_ParseAndHold = false;
	}
	else{
		//metaBit = -1;
	}

}
/*
TH2I * FramesHandler::getAFrameHist(TString fullFileName, TString fullDSCFileName, TString){

	getAFrameMatrix_flag = false; getAFrameHist_flag = true;
	Bool_t frameTrigger = readOneFrame(fullFileName, fullDSCFileName);

	// check the integrity of the frame
	if(!frameSupervisor())
	{
		std::cout << "[ERROR] frame " << fullFileName << "\n";
		std::cout << "        does not seem to be a well formatted frame\n";
		std::cout << "        giving up." << std::endl;
		exit(1);
	}

	TH2I * frameHist = 0;
	TH2I * h_emptyFrame = 0;

	Int_t * matrix_i;
	matrix_i = m_aFrame->GetPixelsMatrix();
	if (matrix_i != 0)
	{
		frameHist = getHistFrame(m_aFrame->GetFrameId(), matrix_i);
	}
	else
	{
		std::cout << "[ERROR] frame matrix not initialized ... givin up" << std::endl;
		exit(1);
	}


	if(frameTrigger)
		return frameHist;
	else
		return h_emptyFrame;
}
 */
Int_t ** FramesHandler::getAFrameMatrix(TString fullFileName, TString /*fullDSCFileName*/){

	getAFrameMatrix_flag = true; getAFrameHist_flag = false;
	//Bool_t frameTrigger = readOneFrame(fullFileName, fullDSCFileName);

	// check the integrity of the frame
	if(!frameSupervisor())
	{
		std::cout << "[ERROR] frame " << fullFileName << "\n";
		std::cout << "        does not seem to be a well formatted frame\n";
		std::cout << "        giving up." << std::endl;
		exit(1);
	}
	/*
    if (m_aFrame->GetPixelsMatrix() != 0)
    return m_aFrame->GetPixelsMatrix();
	 */

	// not used !
	return 0x0;
}


Bool_t FramesHandler::frameSupervisor(){

	/** check to perform when the matrix for one frame
	 *  has been completely loaded.
	 */

	if(m_aFrame->GetEntriesPad() != MAX_FRAME_OCC)
	{
		return false;
	}
	else
	{
		return true;
	}

}
