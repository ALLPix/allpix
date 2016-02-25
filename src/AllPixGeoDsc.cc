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

inline G4ThreeVector linear(const G4ThreeVector value0, const G4ThreeVector value1, const G4double p){

	G4ThreeVector result;

	for (size_t i = 0; i < 3; i++) {
		result[i] = (value0[i] + (value1[i] - value0[i])*(p));
	}

	return result;
}

inline G4ThreeVector bilinear(const G4ThreeVector* ecube, const G4double x, const G4double y, const G4int z){

	G4ThreeVector result;

	G4ThreeVector bil_y1 = linear(ecube[0+2*0+4*z], ecube[1+2*0+4*z], x);
  G4ThreeVector bil_y2 = linear(ecube[0+2*1+4*z], ecube[1+2*1+4*z], x);

  return linear(bil_y1, bil_y2, y);

}

inline G4ThreeVector trilinear(const G4ThreeVector* ecube, const G4ThreeVector pos){

	G4ThreeVector bil_z0, bil_z1;
	bil_z0 = bilinear(ecube, pos[0], pos[1], 0);
	bil_z1 = bilinear(ecube, pos[0], pos[1], 1);

	return linear(bil_z0, bil_z1, pos[2]);

}

inline int int_floor(double x)
{
  int i = (int)x;
  return i - ( i > x );
}

G4ThreeVector AllPixGeoDsc::GetEFieldFromMap(G4ThreeVector ppos){

	G4ThreeVector currentefield;
	
	G4double pixsize_x = GetPixelX();
	G4double pixsize_y = GetPixelY();
	G4double pixsize_z = GetPixelZ();
	
	// ppos is the position in mm inside one pixel cell
	
	ppos = G4ThreeVector(fmod(ppos[0],pixsize_x), fmod(ppos[1],pixsize_y), fmod(ppos[2],pixsize_z));
	
	// Assuming that point 1 and nx are basically at the same position. The "+1" takes account of the efieldmap starting at the iterator 1.
	// Get the position iside the grid

	G4ThreeVector pposgrid;
	pposgrid[0] = ppos[0]/pixsize_x*(G4double)(m_efieldmap_nx-1);
	pposgrid[1] = ppos[1]/pixsize_y*(G4double)(m_efieldmap_ny-1);
	pposgrid[2] = ppos[2]/pixsize_z*(G4double)(m_efieldmap_nz-1);
	
	// Got the position in units of the map coordinates. Do a 3D interpolation of the electric field.
	// Get the eight neighbors and do three linear interpolations.

	G4ThreeVector * ecube = new G4ThreeVector[8];

	for (size_t i = 0; i < 2; i++) {
		for (size_t j = 0; j < 2; j++) {
			for (size_t k = 0; k < 2; k++) {
				ecube[i+2*j+4*k] = m_efieldmap.at(int_floor(pposgrid[2]+k)).at(int_floor(pposgrid[1]+j)).at(int_floor(pposgrid[0]+i));
			}
		}
	}

	// Make pposgrid the position inside the cube

	for (size_t i = 0; i < 3; i++) {
		pposgrid[i] -= (double)(int_floor(pposgrid[i]));
	}
	
	currentefield = trilinear(ecube, pposgrid);
	// Try to just weight the eight edges with the 3D distance to the sampling point.

	delete[] ecube;

	return currentefield;

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
