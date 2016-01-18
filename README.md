
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

Make sure you loaded your Geant4 and ROOT setups (check $G4LIB and $ROOTSYS vars for instance). For user installing allpix on lxplus with SLC6, a GEANT4 installation has been prepared for you convenience. the bash script setup_allpix_lxplus_geant4.9.10.sh need to however edited. in the file, change the folder assigned to G4WORKDIR to a folder where you have write access, objects and executable produced during allpix compilation will be put in there. Ex : 
	
    export G4WORKDIR=~/myG4WorkDirectory	
	
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
First create the Allpix folder : 
	

	mkdir Allpix 

Then checkout the github version of allpix : 
    git clone https://github.com/ALLPix/allpix.git

By default, Allpix_v1.0 is cloned, to get the latest version, checkout the master branch : 

	cd allpix 
	git checkout master 
	cd ..

Now create the build and install folder : 

	mkdir allpix-build allpix-install
	cd allpix-build 

Initialize cmake for compilation of allpix : 

	cmake ../allpix -DCMAKE_CXX_COMPILER=/afs/cern.ch/sw/lcg/external/gcc/4.8/x86_64-slc6/bin/g++ -DCMAKE_C_COMPILER=/afs/cern.ch/sw/lcg/external/gcc/4.8/x86_64-slc6/bin/gcc -DCMAKE_INSTALL_PREFIX=../allpix-install

	
the flags -DCMAKE_CXX_COMPILER and -DCMAKE_C_COMPILER should point to a std-cxx-11 compatible compiler ie gcc4.8 
	
	
Once the Cmake environnement is set, to compile, simply execute : 

	make -jX install 

X is a number of processor 
	
make sure allpix-install/bin is added to the PATH environnement variable. 
	
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

