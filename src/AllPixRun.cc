/**
 *  Author John Idarraga <idarraga@cern.ch>
 */

#include "AllPixRun.hh"

#include "G4Event.hh"
#include "G4HCofThisEvent.hh"
#include "G4SDManager.hh"
#include "G4DigiManager.hh"
#include "G4Trajectory.hh"

// digits, frames
#include "AllPix_Frames_WriteToEntuple.h"
// hits
#include "AllPix_Hits_WriteToEntuple.h"
// dm
#include "allpix_dm.h"

// geometry
#include "ReadGeoDescription.hh"

//
#include "TString.h"
#include <map>
#include <time.h>


/**
 * This constructor is called once per run
 */
AllPixRun::AllPixRun(AllPixDetectorConstruction * det, TString ofp, TString dataset, TString tempDir, G4bool writeFlag, G4bool MCWriteFlag, G4bool eutelWriteFlag){

  m_writeTPixTelescopeFilesFlag = writeFlag;
  m_writeEUTelescopeFilesFlag = eutelWriteFlag;
  m_writeMCFilesFlag=MCWriteFlag; //nalipour: flag for MC

  m_detectorPtr = det;

  // Call for an instance to write.  Need to know how
  // many detectors do I have.  I have as many as
  // the number of digitizers.  At this point the digitizer
  // modules have already been created

  // Geo description
  extern ReadGeoDescription * g_GeoDsc; // already loaded ! :)
  map<int, AllPixGeoDsc *> * geoMap = g_GeoDsc->GetDetectorsMap();
  map<int, AllPixGeoDsc *>::iterator detItr;

  // digit collection info
  m_outputFilePrefix = ofp;
  m_datasetDigits = dataset;
  m_datasetDigits += "_det_";
  m_tempdir = tempDir;
  m_nOfDetectors = (G4int) geoMap->size();

  // create as many frame handlers as pixel detectors
  m_frames = new FramesHandler * [m_nOfDetectors];
  TString tempDataset = "";

  int cntr = 0;

  for( detItr = geoMap->begin() ; detItr != geoMap->end() ; detItr++){

    tempDataset = m_datasetDigits;
    tempDataset += (*detItr).first; // append detector id

    // FramesHandler constructor will communicate to FrameStruct the dataset.
    m_frames[cntr] = new FramesHandler(tempDataset);
    // parameters
    m_frames[cntr]->SetDetectorId((*detItr).first);
    m_frames[cntr]->SetnX( (*detItr).second->GetNPixelsX() );
    m_frames[cntr]->SetnY( (*detItr).second->GetNPixelsY() );

    // map detId to Index
    m_detIdToIndex[(*detItr).first] = cntr;

    //cout << tempDataset << endl;
    //cout << (*detItr).second->GetNPixelsX() << endl;
    //cout << (*detItr).second->GetNPixelsY() << endl;

    cntr++;
    MC_ROOT_data.push_back(new ROOTDataFormat((*detItr).first));  //nalipour: For each detector, initialize a ROOTDataFormat

  }

  /*

    G4DigiManager * fDM = G4DigiManager::GetDMpointer();
    G4int nDC = fDM->GetDCtable()->entries();
    m_nOfDetectors = nDC;
    m_datasetDigits = dataset;

    // append to dataset if there is more than 1 detector
    if(m_nOfDetectors > 1) m_datasetDigits += "_det_";
    m_tempdir = tempDir;

    // create as many frame handlers as pixel detectors
    m_frames = new FramesHandler * [m_nOfDetectors];

    TString tempDataset = m_datasetDigits;
    for(Int_t i = 0 ; i < m_nOfDetectors ; i++){

    if(m_nOfDetectors > 1) {
    tempDataset = m_datasetDigits;
    tempDataset += i; // append detector id
    }
    // FramesHandler constructor will communicate to FrameStruct
    // the dataset.
    m_frames[i] = new FramesHandler(tempDataset);
    // id control !
    m_frames[i]->SetDetectorId(i);
    // FIXME //
    m_frames[i]->SetnX(2048);
    m_frames[i]->SetnY(2048);
    }

  */

  //////////////////////////////////////////////////////
  // Simple hits
  // create as many SimpleHits as SD
  // SD manager
  G4SDManager * SDman = G4SDManager::GetSDMpointer();
  G4HCtable * SDTable = SDman->GetHCtable();
  G4int nSD = SDTable->entries();
  m_nOfSD = nSD;

  m_storableHits = new SimpleHits*[nSD];
  for (Int_t i = 0 ; i < nSD ; i++) {
    // FramesHandler constructor will communicate to FrameStruct
    // the dataset.
    m_storableHits[i] = new SimpleHits;
    m_storableHits[i]->Rewind();
  }

  // Append to name of digits file
  //  Will generate the name later on
  //m_datasetHits = dataset;

  // LCIO Bridge
  // coming from AllPixRunAction
  m_lciobridge_f = 0x0;
  m_lciobridge_dut_f = 0x0;
}

AllPixRun::~AllPixRun(){

  // erase frame handlers
  for (int i = 0 ; i < m_nOfDetectors ; i++) {
    delete m_frames[i]; // delete object using pointer
  }
  delete[] m_frames; // delete array

  // erase storable hits
  for (int i = 0 ; i < m_nOfSD ; i++) {
    delete m_storableHits[i];
  }
  delete[] m_storableHits;

}


void AllPixRun::FillROOTFiles(AllPixWriteROOTFile** rootFiles) //nalipour
{  
  for (uint itr=0; itr<MC_ROOT_data.size(); ++itr)
    {
	  (rootFiles[itr])->SetVectors(MC_ROOT_data[itr]);
	  (rootFiles[itr])->AllPixWriteROOTFillTree();
	  /*
      (rootFiles[itr])->posX_MC=MC_ROOT_data[itr]->get_posX_MC();
      (rootFiles[itr])->posY_MC=MC_ROOT_data[itr]->get_posY_MC();
      (rootFiles[itr])->energy_MC=MC_ROOT_data[itr]->get_energy_MC();

      (rootFiles[itr])->posX=MC_ROOT_data[itr]->get_posX();
      (rootFiles[itr])->posY=MC_ROOT_data[itr]->get_posY();
      (rootFiles[itr])->energy=MC_ROOT_data[itr]->get_energy();
      (rootFiles[itr])->TOT=MC_ROOT_data[itr]->get_TOT();
      (rootFiles[itr])->AllPixWriteROOTFillTree();
      */
    }
  MC_ROOT_data.clear();
}
/*
//nalipour
void AllPixRun::RecordHitsForROOTFiles(const G4Event* evt) // Fill MC_ROOT_data
{
  //MC hits
  G4HCofThisEvent* HCe = evt->GetHCofThisEvent();
  G4int nHC = HCe->GetNumberOfCollections();

  pair<G4int, G4int> tempPixel;
  for (G4int itrCol = 0 ; itrCol < nHC ; itrCol++)  //Iterate over the detectors
    {
      G4cout << "nalipour itrCol=" << itrCol << G4endl;
      m_hitsCollection = dynamic_cast<AllPixTrackerHitsCollection *> (HCe->GetHC(itrCol));

      map<pair<G4int, G4int>, G4double > pixelsContent;
      pair<G4int, G4int> tempPixel;
      G4int NbHits = m_hitsCollection->entries();
      for (G4int i = 0 ; i < NbHits ; i++) 
   	{
   	  tempPixel.first=(*m_hitsCollection)[i]->GetPixelNbX();
   	  tempPixel.second = (*m_hitsCollection)[i]->GetPixelNbY();
   	  pixelsContent[tempPixel]+=(*m_hitsCollection)[i]->GetEdep();
   	}
      
      map<pair<G4int, G4int>, G4double >::iterator pCItr = pixelsContent.begin();
      for( ; pCItr != pixelsContent.end() ; pCItr++)
      {
    	  int x=(*pCItr).first.first;
    	  int y=(*pCItr).first.second;
    	  double energy=(*pCItr).second;
    	  MC_ROOT_data[itrCol]->add_posX_MC(x);
    	  MC_ROOT_data[itrCol]->add_posY_MC(y);
    	  MC_ROOT_data[itrCol]->add_energy_MC(energy/keV);
      }
    }
}

//nalipour
void AllPixRun::RecordHitsForROOTFiles_withChargeSharing(const G4Event* evt) // Fill Data with charge sharing
{
  G4DCofThisEvent* DCe = evt->GetDCofThisEvent();
  if(!DCe)
    {
      G4cout << "No digits in this event !" << G4endl;
      return;
    }
  G4int nDC = DCe->GetNumberOfCollections();
  if(m_nOfDetectors != nDC)
    {
      G4cout << "[ERROR] the number of Digit Collections should not" << G4endl
	     << "        be found to be different than the number of" << G4endl
	     << "        detectors. nDetectors = " << m_nOfDetectors << G4endl
	     << "        nDC =  " << nDC <<  " ... Giving up." << G4endl;
      exit(1);
    }
  
  for (G4int i = 0 ; i < nDC ; i++) {
    
    // Get a hit in the Collection directly
    AllPixDigitsCollectionInterface * digitsCollection =
      static_cast<AllPixDigitsCollectionInterface *> (DCe->GetDC(i));
    G4int nDigits = digitsCollection->entries();    
    for (G4int itr  = 0 ; itr < nDigits ; itr++) 
      {
    	MC_ROOT_data[i]->add_posX((*digitsCollection)[itr]->GetPixelIDX());
    	MC_ROOT_data[i]->add_posY((*digitsCollection)[itr]->GetPixelIDY());
    	MC_ROOT_data[i]->add_energy((*digitsCollection)[itr]->GetPixelEnergyDep()/keV);
       	MC_ROOT_data[i]->add_TOT((*digitsCollection)[itr]->GetPixelCounts());
      }
  }
}*/

void AllPixRun::RecordDigits_all(const G4Event* evt) // Fill Data with charge sharing
{
  G4DCofThisEvent* DCe = evt->GetDCofThisEvent();
  if(!DCe)
    {
      G4cout << "No digits in this event !" << G4endl;
      return;
    }
  G4int nDC = DCe->GetNumberOfCollections();
  if(m_nOfDetectors != nDC)
    {
      G4cout << "[ERROR] the number of Digit Collections should not" << G4endl
	     << "        be found to be different than the number of" << G4endl
	     << "        detectors. nDetectors = " << m_nOfDetectors << G4endl
	     << "        nDC =  " << nDC <<  " ... Giving up." << G4endl;
      exit(1);
    }

  for (G4int i = 0 ; i < nDC ; i++) {

    // Get a hit in the Collection directly
    AllPixDigitsCollectionInterface * digitsCollection =
      static_cast<AllPixDigitsCollectionInterface *> (DCe->GetDC(i));
    G4int nDigits = digitsCollection->entries();
    for (G4int itr  = 0 ; itr < nDigits ; itr++)
      {
    	MC_ROOT_data[i]->add_posX((*digitsCollection)[itr]->GetPixelIDX());
    	MC_ROOT_data[i]->add_posY((*digitsCollection)[itr]->GetPixelIDY());
    	MC_ROOT_data[i]->add_energyTotal((*digitsCollection)[itr]->GetPixelEnergyDep());
       	MC_ROOT_data[i]->add_TOT((*digitsCollection)[itr]->GetPixelCounts());
       	MC_ROOT_data[i]->add_energyMC((*digitsCollection)[itr]->GetPixelEnergyMC());
       	MC_ROOT_data[i]->add_posX_WithRespectToPixel((*digitsCollection)[itr]->Get_posX_WithRespectoToPixel());
       	MC_ROOT_data[i]->add_posY_WithRespectToPixel((*digitsCollection)[itr]->Get_posY_WithRespectoToPixel());
       	MC_ROOT_data[i]->add_posZ_WithRespectToPixel((*digitsCollection)[itr]->Get_posZ_WithRespectoToPixel());
      }
  }
}

void AllPixRun::FillFramesNtuple(const G4Run* aRun){

  /*
  // geo description
  extern ReadGeoDescription * g_GeoDsc; // already loaded ! :)
  map<int, AllPixGeoDsc *> * geoMap = g_GeoDsc->GetDetectorsMap();
  map<int, AllPixGeoDsc *>::iterator detItr;
  int cntr;

  for( detItr = geoMap->begin() ; detItr != geoMap->end() ; detItr++){

  m_frames[cntr]->SetCurrentFrameId(aRun->GetRunID());
  m_frames[cntr]->SetAsMCData(); // <--- !! MC data !!
  WriteToNtuple::GetInstance(m_datasetDigits,
  m_tempdir,
  m_nOfDetectors,
  (*detItr).first)->fillVars(m_frames[cntr]);
  cntr++;
  }
  */

  for(G4int i = 0 ; i < m_nOfDetectors ; i++)
    {

      // fill one frame, set ID first
      m_frames[i]->SetCurrentFrameId(aRun->GetRunID());
      m_frames[i]->SetAsMCData(); // <--- !! MC data !!

      //cout << " AllPixRun::FillFramesNtuple " << m_frames[i]->GetDetectorId() << endl;

      WriteToNtuple::GetInstance(m_outputFilePrefix, m_datasetDigits,
				 m_tempdir,
				 m_nOfDetectors,
				 m_frames[i]->GetDetectorId())->fillVars(m_frames[i]);
    }

}

void AllPixRun::FillTelescopeFiles(const G4Run* aRun, G4String folderName, G4bool eventIDflag, G4bool sumTOTflag){

  /*
   * FILE FORMAT (header on a signle line)
   *
   * # Start time (string) : Sep 22 21:27:06.913 2012 # Start time : 1348342026.913 # Acq time : 0.000000
   * # ChipboardID : C10-W0108 # DACs : 5 100 255 127 127 0 405 7 130 128 80 62 128 128 # Mpx type : 3
   * # Timepix clock : 40 # Eventnr 305 # RelaxD 2 devs 4 TPX DAQ = 0x110402 = START_HW STOP_HW MPX_PAR_READ COMPRESS=1
   * <x> \t <y> \t <sumTOT> \t eventID
   *
   */

  runID = aRun->GetRunID();

  m_runTime = 1351163373; // reference time for Run 0: Thu Oct 25 10:09:33 2012 UTC

  time_t seconds;
  time(&seconds);

  struct tm * ptm;
  //ptm = gmtime(&m_runTime);
  ptm = gmtime(&seconds);
  ptm->tm_sec += runID; // increment by runID*seconds
  time_t newtime = mktime(ptm);
  //string tmp = ctime(&newtime);
  //string s_newtime = tmp.substr(0,tmp.size()-6); // remove newline character

  // format date to fit the filename requirements
  char filedate[80];
  strftime(filedate,sizeof(filedate),"%y%m%d-%H%M%S",gmtime(&newtime));

  char date[80];
  strftime(date,sizeof(date),"%b %d %H:%M:%S.000 %Y",gmtime(&newtime));

  // open Timepix Telescope output file for this run
  TString filename = TString::Format("%s/mpx-%s-0_%i.txt",folderName.data(),filedate,runID);
  
  G4cout << "Path to file: " << filename << G4endl;
  ofstream tpixtelescope_f(filename);

  // loop on data from this run, map filled by RecordTelescopeDigits
  map<int,vector<vector<vector<int> > > >::iterator itr;
  for (itr = m_data.begin(); itr != m_data.end(); itr++)
    {
      // for each telescope plane
      int detId = itr->first;

      // header line for each plane
      tpixtelescope_f << "# Start time (string) : " << date
		      << " # Start time : " << seconds+runID << ".000"
		      << " # Acq time : 0.000000 # ChipboardID : Chip_" << detId
		      << " # DACs : 5 100 255 127 127 0 405 7 130 128 80 62 128 128 # Mpx type : 3 # Timepix clock : 40 "
		      << " # Eventnr " << runID
		      << " # RelaxD 2 devs 4 TPX DAQ = 0x110402 = START_HW STOP_HW MPX_PAR_READ COMPRESS=1"
		      << endl;

      // maps to store multiple hits per pixel
      map<int,int> map_xyTOT;
      map<int,vector<int> > map_xyEvent;

      // loop on events
      vector<vector<vector<int> > > v_events = itr->second;
      for ( unsigned int i=0; i<v_events.size(); i++ )
	{
	  // loop on digits for this event
	  vector<vector<int> > v_digits = v_events[i];
	  for ( unsigned int j=0; j<v_digits.size(); j++ )
	    {
	      vector<int> v_data = v_digits[j];
	      int x    = v_data[0];
	      int y    = v_data[1];
	      int tot  = v_data[2];
	      int evID = v_data[3];

	      /*
	       * sumTOTflag set by /allpix/timepixtelescope/setSumTOT messenger
	       * in AllpixPrimaryGeneratorMesseneger
	       * if true, fill maps in order to sum up digits in same pixel
	       *
	       */
	      if ( sumTOTflag )
		{
		  map_xyTOT[x*1000+y] += tot;
		  map_xyEvent[x*1000+y].push_back(evID);
		}
	      else
		{
		  /*
		   * eventIDflag set by /allpix/timepixtelescope/setEventIDcolumn messenger
		   * if true, dump event number column
		   *
		   */
		  if ( eventIDflag )
		    {
		      tpixtelescope_f << x << "\t" << y << "\t" << tot << "\t" << evID << endl;
		    }
		  else
		    {
		      tpixtelescope_f << x << "\t" << y << "\t" << tot << endl;
		    }
		}
	    }
	}

      if ( sumTOTflag )
	{
	  // loop on maps in order to sum up digits in same pixels
	  map<int,int>::iterator map_itr;
	  for (map_itr = map_xyTOT.begin(); map_itr != map_xyTOT.end(); map_itr++)
	    {
	      int index = map_itr->first;
	      div_t division = div(index,1000);
	      int x = division.quot; // pixel x coordinate
	      int y = division.rem;  // pixel y coordinate
	      int tot = map_itr->second;

	      /*
	       * if more than on digit for a given pixel, evID is set to -1
	       * else print eventID
	       *
	       */
	      int evID = 0;
	      if (map_xyEvent[index].size() > 1)
		{
		  evID = -1;
		}
	      else
		{
		  evID = map_xyEvent[index][0];
		}

	      // dump info in text file
	      if(tot>=0) //nalipour not to write negative TOT 
		{
		  if ( eventIDflag )
		    {
		      tpixtelescope_f << x << "\t" << y << "\t" << tot << "\t" << evID << endl;
		    }
		  else
		    {
		      tpixtelescope_f << x << "\t" << y << "\t" << tot << endl;
		    }
		}
	    }
	}
    }

  // clear map and close file for next run
  // m_data.clear();
  tpixtelescope_f.close();
}



/**
 *  Called at the very end of the event
 *  I can get a handle on the HitsCollection here !
 *  and decide what I want to store.
 */
void AllPixRun::RecordEvent(const G4Event* evt) {

  RecordHits(evt);
  RecordDigits(evt);
  if (m_writeTPixTelescopeFilesFlag || m_writeEUTelescopeFilesFlag)
    RecordTelescopeDigits(evt);

  if(m_writeMCFilesFlag) //nalipour: Record MC hits
    {
	  RecordDigits_all(evt);
      //RecordHitsForROOTFiles(evt);
      //RecordHitsForROOTFiles_withChargeSharing(evt);
    }

}

void AllPixRun::RecordHits(const G4Event* evt) {

  // will need it later
  G4SDManager * SDman = G4SDManager::GetSDMpointer();

  // get the hit collection
  G4HCofThisEvent* HCe = evt->GetHCofThisEvent();

  /*
    if(!m_hitsCollection){
    G4cout << "No hits in this event !" << G4endl;
    return;
    }
  */

  G4int nHC = HCe->GetNumberOfCollections();

  for (G4int itrCol = 0 ; itrCol < nHC ; itrCol++) {
    // get a hit in the Collection directly
    m_hitsCollection = dynamic_cast<AllPixTrackerHitsCollection *> (HCe->GetHC(itrCol));

    G4int NbHits = m_hitsCollection->entries();
    //G4cout << "Recording hits : " << NbHits <<  G4endl;

    // rewind hits
    m_storableHits[itrCol]->Rewind();

    G4ThreeVector posTemp;
    TVector3 posTempStorable;

    for (G4int i = 0 ; i < NbHits ; i++) {

      //(*m_hitsCollection)[i]->Print();

      // Record total energy deposited in an event in all pixels.  SCALAR
      m_storableHits[itrCol]->edepTotal += ((*m_hitsCollection)[i]->GetEdep())/keV;

      // Interaction name.  VECTOR
      m_storableHits[itrCol]->interactions.push_back(
						     ((*m_hitsCollection)[i]->GetProcessName()).data()
						     );

      // Translate pos to storable format.  VECTOR
      posTemp = (*m_hitsCollection)[i]->GetPos();
      posTempStorable.SetXYZ( posTemp.x(), posTemp.y(), posTemp.z() );
      m_storableHits[itrCol]->pos.push_back(
					    posTempStorable
					    );

      //if(m_hitsCollection->GetName().contains("scint")) {
      //	G4cout << "--> " << m_hitsCollection->GetName() << " : posx = " << posTemp.x()/mm << " [mm]" << G4endl;
      //}

      // pdgId of parent particle
      m_storableHits[itrCol]->pdgId.push_back(
					      (*m_hitsCollection)[i]->GetTrackPdgId()
					      );

      // energy dep per hit
      m_storableHits[itrCol]->edep.push_back(
					     ((*m_hitsCollection)[i]->GetEdep())/keV
					     );

      // particle Id
      m_storableHits[itrCol]->trackId.push_back(
						(*m_hitsCollection)[i]->GetTrackID()
						);

      // parent Id
      m_storableHits[itrCol]->parentId.push_back(
						 (*m_hitsCollection)[i]->GetParentID()
						 );

      // Name of volume for this hit
      m_storableHits[itrCol]->trackVolumeName.push_back(
							(*m_hitsCollection)[i]->GetTrackVolumeName()
							);

      // Name of volume where particle was created
      m_storableHits[itrCol]->parentVolumeName.push_back(
							 (*m_hitsCollection)[i]->GetParentVolumeName()
							 );
    }

    // event id
    m_storableHits[itrCol]->event = evt->GetEventID();

    // run id
    m_storableHits[itrCol]->run = GetRunID();

    /** It is guaranteed by SD::ProcessHits that the hits contain
     *  energy deposit.  Thus we store every hit collection.
     */

    // fillVars rewinds values
    m_datasetHits = SDman->GetHCtable()->GetHCname(itrCol);
    Hits_WriteToNtuple::GetInstance(m_outputFilePrefix, m_datasetHits,
				    m_tempdir,
				    nHC, // here is the number of Hit Collections (SD), not detectors.
				    itrCol)->fillVars(m_storableHits[itrCol]);

  }
}

void AllPixRun::RecordDigits(const G4Event* evt){

  // Check if this event triggered for EUTELESCOPE only
#ifdef _EUTELESCOPE

  G4HCofThisEvent* HCe = evt->GetHCofThisEvent();

  G4int scintillatorsCntr = 0;
  G4int triggerCntr = 0;

  for(int i = 0 ; i < HCe->GetNumberOfCollections() ; i++){

    AllPixTrackerHitsCollection * hc = (AllPixTrackerHitsCollection *) HCe->GetHC(i);
    G4String hitName = hc->GetName();

    // look for scintillators
    if(hitName.contains("sdscint")) {

      // if containing hits !
      if((int)hc->GetSize() > 0) scintillatorsCntr++;

      AllPixTrackerHit * hit = 0x0;
      for( int j= 0 ; j < (int)hc->GetSize() ; j++ ) {

	hit = (AllPixTrackerHit *) hc->GetHit(j);
	// pdgId of primary particle.  Triggering on it.
	//G4int pdgPrim = evt->GetPrimaryVertex()->GetPrimary(0)->GetPDGcode();
	//if(hit->GetTrackPdgId() == pdgPrim) // pion for now
	triggerCntr++;

      }
    }
  }

  // trigger desicion
  if(scintillatorsCntr < __magic_trigger_cntr_ack){
    G4cout << "Trigger OFF --> " << scintillatorsCntr << " scintillators fired" << G4endl;
    return;
  } else {
    G4cout << "Trigger ON --> " << scintillatorsCntr << " scintillators fired" << G4endl;
  }
#endif

  // get digits in this event
  G4DCofThisEvent* DCe = evt->GetDCofThisEvent();
  if(!DCe){
    G4cout << "No digits in this event !" << G4endl;
    return;
  }

  G4int nDC = DCe->GetNumberOfCollections();

  // For now the number of DigitCollections should be
  // equal to the number of detectors
  if(m_nOfDetectors != nDC){
    G4cout << "[ERROR] the number of Digit Collections should not" << G4endl
	   << "        be found to be different than the number of" << G4endl
	   << "        detectors. nDetectors = " << m_nOfDetectors << G4endl
	   << "        nDC =  " << nDC <<  " ... Giving up." << G4endl;
    exit(1);
  }

  // lcio bridge
  // runId
  *m_lciobridge_f << "R " << GetRunID() << endl;
  *m_lciobridge_dut_f << "R " << GetRunID() << endl;

  /*
  // check event for information about the track
  //G4TrajectoryContainer * tC = evt->GetTrajectoryContainer();
  //G4Trajectory * tr;
  G4TrajectoryContainer * tC = evt->GetTrajectoryContainer();
  G4cout << "Trajectories : " << tC->GetVector()->size() << G4endl;
  TrajectoryVector * tV = tC->GetVector();
  TrajectoryVector::iterator tI = tV->begin();
  TrajectoryVector::iterator tIE = tV->end();

  G4Trajectory * tr;
  G4TrajectoryPoint* tP;
  for ( ; tI != tIE ; tI++ ) {
  tr = static_cast<G4Trajectory*>(*tI);
  //G4cout << "    " << tr->GetPDGEncoding() << G4endl;
  //G4cout << "       " << tr->GetPointEntries() << G4endl;
  tP = static_cast<G4TrajectoryPoint*>((*tI)->GetPoint(0));

  // Check which of the points is at the interior of the sensor
  // Use the G4VSolid::Inside(const G4ThreeVector& p)
  // Find first the Silicon G4VSolids (sensors).


  //G4VSolid * sol = m_detectorPtr->GetVSolidDetector(100);
  //G4cout << "       " << tP->GetPosition().z()/mm << G4endl;

  //if(sol->Inside(tP->GetPosition()) == kInside){
  //	G4cout << "       inside !!" << G4endl;
  //} else {
  //	G4cout << "       outside !!" << G4endl;
  //}

  }
  */

  // Ok, I have to match the right detector m_frames with the digit collection
  for (G4int i = 0 ; i < nDC ; i++) {

    // Get a hit in the Collection directly
    AllPixDigitsCollectionInterface * digitsCollection =
      static_cast<AllPixDigitsCollectionInterface *> (DCe->GetDC(i));
    G4int nDigits = digitsCollection->entries();

    // Pickup the right index
    TString theIndex_S = digitsCollection->GetName().c_str();
    int lastpos = theIndex_S.Index("_DigitCollection", 16, 0, TString::kExact);
    theIndex_S.Remove(lastpos, theIndex_S.Length()); // remove suffix
    theIndex_S.Remove(0,6); // remove BoxSD_
    int detId = atoi(theIndex_S.Data());

    // Clean up matrix
    //m_frames[m_detIdToIndex[detId]]->getFrameStructObject()->CleanUpMatrix();

    // lcio bridge
    // FIXME !
    if(detId >= 300) *m_lciobridge_f << detId << " ";
    if(detId < 300) *m_lciobridge_dut_f << detId << " ";

    for (G4int itr  = 0 ; itr < nDigits ; itr++) {

      // loading a frame
      m_frames[m_detIdToIndex[detId]]->LoadFramePixel(
						      (*digitsCollection)[itr]->GetPixelIDX(),
						      (*digitsCollection)[itr]->GetPixelIDY(),
						      (*digitsCollection)[itr]->GetPixelCounts(),
						      (*digitsCollection)[itr]->GetPixelEnergyDep()/keV,
						      0.
						      );

      m_frames[m_detIdToIndex[detId]]->LoadPrimaryVertexInfo(
							     (*digitsCollection)[itr]->GetPrimaryVertex().x()/mm,
							     (*digitsCollection)[itr]->GetPrimaryVertex().y()/mm,
							     (*digitsCollection)[itr]->GetPrimaryVertex().z()/mm
							     );

      // lcio bridge
      // FIXME !
      if (detId >= 300) {
	*m_lciobridge_f << (*digitsCollection)[itr]->GetPixelIDX() << " "
			<< (*digitsCollection)[itr]->GetPixelIDY() << " "
			<< (*digitsCollection)[itr]->GetPixelCounts() << " ";
      } else {
	*m_lciobridge_dut_f << (*digitsCollection)[itr]->GetPixelIDX() << " "
			    << (*digitsCollection)[itr]->GetPixelIDY() << " "
			    << (*digitsCollection)[itr]->GetPixelCounts() << " ";
      }
    }

    //if(detId >= 300 and nDigits == 0) *m_lciobridge_f << '*';

    // lcio bridge
    // FIXME !
    if(detId >= 300) *m_lciobridge_f << endl;
    if(detId < 300) *m_lciobridge_dut_f << endl;
  }
}

void AllPixRun::RecordTelescopeDigits(const G4Event* evt){

  // get digits in this event
  G4DCofThisEvent* DCe = evt->GetDCofThisEvent();
  if(!DCe){
    G4cout << "No digits in this event !" << G4endl;
    return;
  }

  G4int eventID = evt->GetEventID();
  //G4cout << "============================================================================> Event ID: " << eventID << G4endl;

  G4int nDC = DCe->GetNumberOfCollections();

  // For now the number of DigitCollections should be
  // equal to the number of detectors
  if(m_nOfDetectors != nDC){
    G4cout << "[ERROR] the number of Digit Collections should not" << G4endl
	   << "        be found to be different than the number of" << G4endl
	   << "        detectors. nDetectors = " << m_nOfDetectors << G4endl
	   << "        nDC =  " << nDC <<  " ... Giving up." << G4endl;
    exit(1);
  }

  // Ok, I have to match the right detector m_frames with the digit collection
  for (G4int i = 0 ; i < nDC ; i++) {

    // Get a hit in the Collection directly
    AllPixDigitsCollectionInterface * digitsCollection = static_cast<AllPixDigitsCollectionInterface *> (DCe->GetDC(i));
    G4int nDigits = digitsCollection->entries();

    // Pickup the right index
    TString theIndex_S = digitsCollection->GetName().c_str();
    int lastpos = theIndex_S.Index("_DigitCollection", 16, 0, TString::kExact);
    theIndex_S.Remove(lastpos, theIndex_S.Length()); // remove suffix
    theIndex_S.Remove(0,6); // remove BoxSD_
    int detId = atoi(theIndex_S.Data());

    vector<vector<int> > v_det;

    for (G4int itr  = 0 ; itr < nDigits ; itr++) {

      int x   = (*digitsCollection)[itr]->GetPixelIDX();
      int y   = (*digitsCollection)[itr]->GetPixelIDY();
      int tot = (*digitsCollection)[itr]->GetPixelCounts();
      //int tot = (*digitsCollection)[itr]->GetPixelEnergyDep();

      // store x,y,tot,eventID in vector
      vector<int> v_digit;
      v_digit.push_back(x);
      v_digit.push_back(y);
      v_digit.push_back(tot);
      v_digit.push_back(eventID);

      // append to vector for keeping
      v_det.push_back(v_digit);

    }

    // append to detector vector
    m_data[detId].push_back(v_det);

  }

}




















