#!/bin/bash

source /afs/cern.ch/eng/clic/software/x86_64-slc6-gcc48/setup.sh

cd /afs/cern.ch/eng/clic/software/x86_64-slc6-gcc48/LCIO/v02-04 ; source setup.sh ; cd -

export PYTHONPATH="/afs/cern.ch/eng/clic/software/Pixel_TestBeam_Software/pytools/Cython-LCIO/lib/python2.7/site-packages:$PYTHONPATH"
