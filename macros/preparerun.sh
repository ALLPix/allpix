#!/bin/bash

if [ ${#@} -lt 4 ] ; then
    echo "use $0 nFrames(int) nPerFrame(int) fragment(filename) batchfile(newfilename)"
fi

c=0
max=$1
nPerFrame=$2
fragment=$3
newfile=$4
switch=0

cp $fragment $newfile

while [ $c -lt $max ]
  do
  echo "/run/beamOn $nPerFrame" >> $newfile
  let c=$c+1
done

echo "[DONE]"
