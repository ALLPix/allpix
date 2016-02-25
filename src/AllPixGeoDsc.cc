/**
 * Allpix
 * Author: John Idarraga <idarraga@cern.ch> , 2010
 */

#include "AllPixGeoDsc.hh"



AllPixGeoDsc::AllPixGeoDsc() :
    m_coverlayer_hz(0),
    m_coverlayer_mat("G4Al"),
    m_coverlayer_ON(false)
{

	m_efieldfromfile = false;

}

AllPixGeoDsc::~AllPixGeoDsc() {

}

/* temporary ... this has to be changed by a db com or something like that */

void AllPixGeoDsc::Dump() {

    G4cout << "Dump geo description for object with Id : " << m_ID << G4endl;

    G4cout << "   Digitizer         : " << m_digitizer << G4endl;
    G4cout << "   npix_x            = " << m_npix_x << G4endl;
    G4cout << "   npix_y            = " << m_npix_y << G4endl;
    G4cout << "   npix_z            = " << m_npix_z << G4endl;
    G4cout << "   pixsize_x         = " << m_pixsize_x/mm << " [mm]" << G4endl;
    G4cout << "   pixsize_y         = " << m_pixsize_y/mm << G4endl;
    G4cout << "   pixsize_z         = " << m_pixsize_z/mm << G4endl;
    G4cout << "   sensor_hx         = " << m_sensor_hx/mm << ", posx "<< m_sensor_posx/mm << G4endl;
    G4cout << "   sensor_hy         = " << m_sensor_hy/mm << ", posy "<< m_sensor_posy/mm << G4endl;
    G4cout << "   sensor_hz         = " << m_sensor_hz/mm << ", posz "<< m_sensor_posz/mm << G4endl;
    G4cout << "   coverlayer_hz     = " << m_coverlayer_hz/mm << G4endl;
    G4cout << "   coverlayer_mat    = " << m_coverlayer_mat << G4endl;
    G4cout << "   chip_hx           = " << m_chip_hx/mm << ", posx "<< m_chip_posx/mm << G4endl;
    G4cout << "   chip_hy           = " << m_chip_hy/mm << ", posy "<< m_chip_posy/mm << G4endl;
    G4cout << "   chip_hz           = " << m_chip_hz/mm << ", posz "<< m_chip_posz/mm << G4endl;
    G4cout << "   pcb_hx            = " << m_pcb_hx/mm << G4endl;
    G4cout << "   pcb_hy            = " << m_pcb_hy/mm << G4endl;
    G4cout << "   pcb_hz            = " << m_pcb_hz/mm << G4endl;

}

void AllPixGeoDsc::SetEFieldMap(G4String valS){
	m_EFieldFile = valS;

	// Get EFieldFile from $allpix/valS
	struct stat buffer;
	if(!(stat (valS, &buffer))){
		
		m_efieldfromfile = true;

		ifstream efieldinput;
		efieldinput.open(valS);

		G4int nptsx, nptsy, nptsz;
		G4int currx, curry, currz;

		G4double ex, ey, ez;

		efieldinput >> nptsx >> nptsy >> nptsz;

		vector<G4ThreeVector> onedim;
		vector<vector<G4ThreeVector>> twodim;
		onedim.reserve(nptsx);
		twodim.reserve(nptsy);

		for (int k = 0; k < nptsz; k++) {
			for (int j = 0; j < nptsy; j++) {
				for (int i = 0; i < nptsx; i++) {
					efieldinput >> currx >> curry >> currz;
					efieldinput >> ex >> ey >> ez;
					onedim.push_back(G4ThreeVector(ex, ey, ez));
				}
				twodim.push_back(onedim);
				onedim.clear();
			}
			m_efieldmap.push_back(twodim);
			twodim.clear();
		}

		m_efieldmap_nx = nptsx;
		m_efieldmap_ny = nptsy;
		m_efieldmap_nz = nptsz;

		efieldinput.close();

	}else{
		if(!valS.isNull())
		{
			cout << "File for EField named, but not found: " << valS << endl;
			cout << "This cannot end well... Abort." << endl;
			m_efieldfromfile = false;
			exit(1);
		}
		cout << "Found no EFieldFile " << valS << endl;
		m_efieldfromfile = false;
	}

}

/*
void AllPixGeoDsc::operator=(AllPixGeoDsc & cp){

    cout << cp.GetSensorDigitizer() << " , " << m_digitizer << endl;
    m_digitizer = cp.GetSensorDigitizer();

}
*/

/*
void AllPixGeoDsc::ReadGeoBits(){

    ifstream ifs( "DetDscrDb/detgeom.txt" , ifstream::in );

#define NCHARS 128

    char temp[NCHARS];
    string tempS;
    double value = 0.;

    while (ifs.good()){

        ifs.getline(temp, NCHARS);
        tempS = temp;

        if(temp[0] == '#' || tempS.empty() ){ // comment or empty string
            continue;
        }

        for(int colItr = 0 ; colItr < 5 ; colItr++){
            value = GetColFromLine(tempS, colItr);
            switch (colItr) {
            case 0:
                m_pcb_hx= value*um;
                break;
            case 1:
                m_pcb_hy= value*um;
                break;
            case 2:
                m_pcb_hz = value*um;
                break;
            case 3:
                m_sensor_posx = value*um;
                break;
            case 4:
                m_sensor_posy = value*um;
                break;
            }
        }

    }

    ifs.close();

}

double AllPixGeoDsc::GetColFromLine(string tempS, int col){

    size_t pos;
    string cutS;

    int localcol = 0;
    double value = 0.;

    //pos = tempS.find(" ");

    while((pos=tempS.find(" ")) != string::npos)
    {
        cutS = tempS.substr(0, pos);
        value = atof(cutS.c_str());

        // in the desired col
        if(localcol == col){
            return value;
        }
        // shrink string
        tempS = tempS.substr(pos+1);
        localcol++;
    }

    // last entry has no extra space, checking last one
    value = atof(tempS.c_str());
    if(localcol == col)
        return value;

    // if reaching this point, the col couldn't be found
    G4cout << "[ERROR] ouch !, couldn't find that colum (AllPixGeoDsc::GetColFromLine)" << endl;
    G4cout << "        giving up, geometry depends on this." << endl;
    exit(1);

    return -1;
}
*/
