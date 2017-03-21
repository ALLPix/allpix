
# AllPix     		    
## Generic simulation for pixel detectors	
                                                                                      
John Idarraga <idarraga@cern.ch>        
Mathieu Benoit <mbenoit@cern.ch>  
Samir Arfaoui  <sarfaoui@cern.ch>     

### Dependencies:

1) Geant4 v4.10 or newer must be installed with the following dependencies
satisfied

- GDML (xcerces)
- OpenGL (optional, only needed to use openGL vizualisation)
- Inventor  (optional, only needed to use the very convenient OpenInventor vizualisation)

2) ROOT 6.04 or newer .  A basic setup with xml parser is enough.

### Build:

Make sure you loaded your Geant4 and ROOT setups (check $G4LIB and $ROOTSYS vars for instance). For user installing allpix on lxplus with SLC6, a GEANT4 installation has been prepared for your convenience. The bash script setup_allpix_lxplus_geant4.10.3.1.sh may need however to be edited. In the file, change the folder assigned to G4WORKDIR to a folder where you have write access, objects and executable produced during allpix compilation will be put in there. Ex : 
```	  
    export G4WORKDIR=~/myG4WorkDirectory	
```	

The default is $HOME/Allpix/allpix-install, according to the following instructions.

Allpix is now compilable using Cmake. We suggest the following work folder structure. 	 
```	
|--- Allpix /  				# Mother folder containing source, build and install folder
	 |----- allpix 			# Source code folder to be checked out from Github 
	 |      |----- src 
	 |      |----- include
	 |      |----- models
	 |      |----- share
	 |      |----- macros 	 
	 |----- allpix-build 	# Build folder for cmake 
	 |----- allpix-install # Installation folder for cmake 
	 |		|----- bin		# allpix executable folder
```
First create the mother Allpix folder : 	
``` 
	mkdir AllpixMotherFolder 
	cd AllpixMotherFolder
``` 

Then checkout the github version of allpix : 
    git clone https://github.com/ALLPix/allpix.git

If running on lxplus, source the bash setup script (make sure you are using bash, if not, type ```exec bash```) to get access the the precompiled geant4, root and other dependencies : 

``` 
source allpix/setup_allpix_geant4.10.3.1.sh
```
Now create the build and install folder : 

	mkdir allpix-build allpix-install
	cd allpix-build 

Initialize cmake for compilation of allpix : 

	cmake ../allpix -DCMAKE_INSTALL_PREFIX=$G4WORKDIR	

if You wish to compile with the LCIO file writer option, use the following on LXPLUS :
```
	cmake ../allpix -DCMAKE_INSTALL_PREFIX=$G4WORKDIR -Dlcio=ON	
```

Please not that for a custom installation, the environnement variable LCIO should be set and poiting to your installation of LCIO. 
	
Once the Cmake environnement is set, to compile, simply execute : 

	make -jX install 

X is a number of processor 
	
If you follow the structure below, the allpix executable should be added to the PATH by the lxplus setup script. However, for custom installation, make sure allpix-install/bin is added to the PATH environnement variable using this command, you may want to add this command to your setup script : 

	export PATH=$PATH:PATHTOALLPIXINSTALLBINFOLDER
	
NOTE : the pixeldetector.xml file that is used for simulation is picked following a relative path to the allpix execution folder. So execute allpix executable from inside allpix folder, as before
NOTE : The makefile method is deprecated for this version, please use CMAKE instead 
	

### Folder Structure : 


the code on github follow this folder structure : 

```
	 |----- allpix 			
	 |      |----- src 
	 |      |----- include
	 |      |----- models
	 |      |----- share
	 |      |----- macros 	
```

- include and src : These contain the source code of the allpix framework, whenever you create a new digitizer, the .cc and  .h
files will be added to these directories 

- macros : This folder contains example of macro files to be used for different type of Simulation framework , see allpix Twiki
page for more details : https://twiki.cern.ch/twiki/bin/view/Main/AllPix

- models : This folder contains the pixel geometry database (pixeldetector.xml). If one imports geometry from a gdml file, it
should be put here. (ex : clicpix_box.gdml and clicpic_box_materials.xml)

- share :  This directory contains files and scripts indirectly related to allpix, for example , the python script to translate
allpix output to slcio file for processing in the EUTELESCOPE framework. See Twiki page for more details https://twiki.cern.ch/twiki/bin/view/Main/AllPix
	
	

### Preparing your Simulation : 

See : https://twiki.cern.ch/twiki/bin/view/Main/AllPix


	
### Running:

Setup your Geant4 and ROOT environments and run the following
way, from the source folder

Interactive run:
    allpix macros/telescope1_vis.in

Batch run : 
    allpix macros/telescope1.in batch
    
    
    
### Presentation and examples

Tutorial on simulation of silicon pixel detectors with Allpix
https://indico.desy.de/contributionDisplay.py?contribId=9&confId=16161

Measurement of the Lorentz angle in CMS pixel detector modules, by Paul Schuetze
https://indico.desy.de/contributionDisplay.py?contribId=11&confId=16161

First use of the Allpix framework and preliminary results, by Andreas Heggelund
https://indico.desy.de/getFile.py/access?contribId=56&sessionId=10&resId=0&materialId=slides&confId=16161

ALiBaVa Strip Sensor Analysis, by Thomas Eichhorn
https://indico.desy.de/getFile.py/access?contribId=21&sessionId=8&resId=0&materialId=slides&confId=16161

Full simulation of the LUCID experiment in the Low Earth Orbit radiation environment by T. Whyntiea,b and M.A. Harrisona
http://iopscience.iop.org/article/10.1088/1748-0221/10/03/C03043/meta

Calibration, simulation and test-beam characterisation of Timepix hybrid-pixel readout assemblies with ultra-thin sensors 
https://agenda.linearcollider.org/event/6000/contributions/27724/attachments/23014/35842/Samir_LCWS.pdf

Characterization of 3D Silicon Assemblies for ATLAS Pixel Upgrade, Thesis by Marcello Borri
https://www.escholar.manchester.ac.uk/api/datastream?publicationPid=uk-ac-man-scw:198996&datastreamId=FULL-TEXT.PDF

CLIC Telescope Optimization with ALLPIX Simulation, by Wu Qi
https://cds.cern.ch/record/2046145/files/CERN_WU_Qi.pdf

Recent results with HV-CMOS and planar sensors for the CLIC vertex detector, by Niloufar Alipour Tehrani
https://indico.cern.ch/event/391665/contributions/1827195/attachments/1230571/1803598/VCI2016_NiloufarAlipourTehrani.pdf

The FE-I4 telescope for particle tracking in testbeam experiments, by M. Benoit et al.
http://iopscience.iop.org/article/10.1088/1748-0221/11/07/P07003/meta








