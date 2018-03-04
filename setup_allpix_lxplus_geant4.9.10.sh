

#gcc4.8
source /afs/cern.ch/sw/lcg/external/gcc/4.8/x86_64-slc6/setup.sh


#python2.7.4
export PATH="/afs/cern.ch/sw/lcg/external/Python/2.7.4/x86_64-slc6-gcc48-opt/bin:$PATH"
export LD_LIBRARY_PATH="/afs/cern.ch/sw/lcg/external/Python/2.7.4/x86_64-slc6-gcc48-opt/lib:$LD_LIBRARY_PATH" 


#numpy/scipy/sympy

export PYTHONPATH=$PYTHONPATH:/afs/cern.ch/work/m/mbenoit/public/AllPixSoftware/pytools/numpy/lib/python2.7/site-packages
export PYTHONPATH=$PYTHONPATH:/afs/cern.ch/work/m/mbenoit/public/AllPixSoftware/pytools/scipy/lib/python2.7/site-packages
export PYTHONPATH=$PYTHONPATH:/afs/cern.ch/work/m/mbenoit/public/AllPixSoftware/pytools/sympy/lib/python2.7/site-packages
export PYTHONPATH=$PYTHONPATH:/afs/cern.ch/work/m/mbenoit/public/AllPixSoftware/pytools/cython/lib/python2.7/site-packages
export PYTHONPATH=$PYTHONPATH:/afs/cern.ch/work/m/mbenoit/public/AllPixSoftware/pytools/fastcluster/lib/python2.7/site-packages


#root6
source /afs/cern.ch/work/m/mbenoit/public/AllPixSoftware/geant4/root-6.04.00/bin/thisroot.sh
export PATH=$ROOTSYS/bin:$PATH

#xerces 

export LD_LIBRARY_PATH=/afs/cern.ch/eng/clic/software/Pixel_TestBeam_Software/xerces/lib:$LD_LIBRARY_PATH


#open Inventor
export LD_RUN_PATH=/afs/cern.ch/eng/clic/software/Pixel_TestBeam_Software/OpenInventor/OIV_Install/lib
export PATH=$PATH:/afs/cern.ch/eng/clic/software/Pixel_TestBeam_Software/OpenInventor/OIV_Install/bin
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/afs/cern.ch/eng/clic/software/Pixel_TestBeam_Software/OpenInventor/OIV_Install/lib
export OIVHOME=/afs/cern.ch/eng/clic/software/Pixel_TestBeam_Software/OpenInventor/OIV_Install
export OIVLIBS=/afs/cern.ch/eng/clic/software/Pixel_TestBeam_Software/OpenInventor/OIV_Install/lib/libSoXt.so

export G4WORKDIR=~/myG4WorkDirectory 
#export G4WORKDIR=/home/mbenoit/Allpix/allpix-install/bin
export PATH=$PATH:$G4WORKDIR


source /afs/cern.ch/work/m/mbenoit/public/AllPixSoftware/geant4/4.9.10.1p02/share/Geant4-10.1.2/geant4make/geant4make.sh
export PATH=$PATH:/afs/cern.ch/user/b/bnachman/work/public/RadDamage/newversion/allpix-install/bin/