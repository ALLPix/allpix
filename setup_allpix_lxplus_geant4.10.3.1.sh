#cmake 3.3.2
export PATH=/afs/cern.ch/sw/lcg/contrib/CMake/3.3.2/Linux-x86_64/bin:${PATH}
#gcc4.8
#source /afs/cern.ch/sw/lcg/external/gcc/4.8/x86_64-slc6/setup.sh
# later for centos7
#source /cvmfs/sft.cern.ch/lcg/contrib/gcc/4.8.4/x86_64-centos7-gcc48-opt/setup.sh
source /cvmfs/sft.cern.ch/lcg/contrib/gcc/4.8.4/x86_64-slc6-gcc48-opt/setup.sh


#python2.7.4
export PATH="/afs/cern.ch/sw/lcg/external/Python/2.7.4/x86_64-slc6-gcc48-opt/bin:$PATH"
export LD_LIBRARY_PATH=/afs/cern.ch/sw/lcg/external/Python/2.7.4/x86_64-slc6-gcc48-opt/lib:$LD_LIBRARY_PATH

#LCIO
export LCIO="/afs/cern.ch/work/m/mbenoit/public/AllPixSoftware/geant4/LCIO/LCIO-install"
export ILCUTIL_DIR="/afs/cern.ch/work/m/mbenoit/public/AllPixSoftware/geant4/LCIO/iLCUtil-install"
export LD_LIBRARY_PATH=/afs/cern.ch/work/m/mbenoit/public/AllPixSoftware/geant4/LCIO/LCIO-install/lib:$LD_LIBRARY_PATH

#numpy/scipy/sympy
export PYTHONPATH=$PYTHONPATH:/afs/cern.ch/work/m/mbenoit/public/AllPixSoftware/pytools/numpy/lib/python2.7/site-packages
export PYTHONPATH=$PYTHONPATH:/afs/cern.ch/work/m/mbenoit/public/AllPixSoftware/pytools/scipy/lib/python2.7/site-packages
export PYTHONPATH=$PYTHONPATH:/afs/cern.ch/work/m/mbenoit/public/AllPixSoftware/pytools/sympy/lib/python2.7/site-packages
export PYTHONPATH=$PYTHONPATH:/afs/cern.ch/work/m/mbenoit/public/AllPixSoftware/pytools/cython/lib/python2.7/site-packages
export PYTHONPATH=$PYTHONPATH:/afs/cern.ch/work/m/mbenoit/public/AllPixSoftware/pytools/fastcluster/lib/python2.7/site-packages


#root6
#source /afs/cern.ch/work/m/mbenoit/public/AllPixSoftware/geant4/root-6.04.00/bin/thisroot.sh
source /cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.04.18/x86_64-slc6-gcc48-opt/root/bin/thisroot.sh
#export PATH=$ROOTSYS/bin:$PATH

#xerces 

export LD_LIBRARY_PATH=/afs/cern.ch/work/m/mbenoit/public/AllPixSoftware/xerces/lib:$LD_LIBRARY_PATH


#open Inventor
export LD_RUN_PATH=/afs/cern.ch/work/m/mbenoit/public/AllPixSoftware/OIV_Install/lib
export PATH=$PATH:/afs/cern.ch/work/m/mbenoit/public/AllPixSoftware/OIV_Install/bin
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/afs/cern.ch/work/m/mbenoit/public/AllPixSoftware/OIV_Install/lib
export OIVHOME=/afs/cern.ch/work/m/mbenoit/public/AllPixSoftware/OIV_Install
export OIVLIBS=/afs/cern.ch/work/m/mbenoit/public/AllPixSoftware/OIV_Install/lib/libSoXt.so


export G4WORKDIR=$HOME/Allpix/allpix-install
export PATH=$PATH:$G4WORKDIR/bin


source /afs/cern.ch/work/m/mbenoit/public/AllPixSoftware/geant4/4-10.3.1/share/Geant4-10.3.1/geant4make/geant4make.sh
export LD_LIBRARY_PATH=/afs/cern.ch/sw/lcg/external/gcc/4.8/x86_64-slc6/lib64:$LD_LIBRARY_PATH

