#!/bin/bash

# Author: John Idarraga <idarraga@cern.ch>
# LAL, Orsay

# Starting up a new digitizer for allpix

echo
echo "[-->?] Enter the name of your new Digitizer"
echo "       avoiding spaces ex: FEI7"

read digiName

if [ -e ./src/AllPix${digiName}Digitizer.cc ]; then
    echo "[WARNING] file AllPix${digiName}Digitizer.cc already exists !!!"
    echo "          hit any key to continue and OVERRIDE, or Ctrl+C to quit"
    read q1
fi

echo
echo "[-->?] Enter the author string"
echo "       ex(spaces are ok): someone <someone@cern.ch>"

read authorS

echo "[INFO] creating files for new digitizer --> \"${digiName}\""

nameFragment="AllPix${digiName}"
digiImp=${nameFragment}Digit.cc
digiHeader=${nameFragment}Digit.hh
digitizerImp=${nameFragment}Digitizer.cc
digitizerHeader=${nameFragment}Digitizer.hh

echo "[INFO] creating header and implementation files for Digits"
echo "       include/${digiHeader}   src/${digiImp}"
cp src/digitizers_base/Digit.hh include/${digiHeader}
cp src/digitizers_base/Digit.cc src/${digiImp}

echo "[INFO] creating header and implementation files for Digitizer"
echo "       include/${digitizerHeader}   src/${digitizerImp}"
cp src/digitizers_base/Digitizer.hh include/${digitizerHeader}
cp src/digitizers_base/Digitizer.cc src/${digitizerImp}

### 1
sed "s|__|${nameFragment}|g" < include/${digiHeader} > include/${digiHeader}_new
sed "s|AAAuthor|${authorS}|g" < include/${digiHeader}_new > include/${digiHeader}_new2
rm -f include/${digiHeader}_new
mv include/${digiHeader}_new2 include/${digiHeader}

### 2
sed "s|__|${nameFragment}|g" < src/${digiImp} > src/${digiImp}_new
sed "s|AAAuthor|${authorS}|g" < src/${digiImp}_new > src/${digiImp}_new2
rm -f src/${digiImp}_new
mv src/${digiImp}_new2 src/${digiImp}

### 3
sed "s|__|${nameFragment}|g" < include/${digitizerHeader} > include/${digitizerHeader}_new
sed "s|AAAuthor|${authorS}|g" < include/${digitizerHeader}_new > include/${digitizerHeader}_new2
rm -f include/${digitizerHeader}_new
mv include/${digitizerHeader}_new2 include/${digitizerHeader}

### 4
sed "s|__|${nameFragment}|g" < src/${digitizerImp} > src/${digitizerImp}_new
sed "s|AAAuthor|${authorS}|g" < src/${digitizerImp}_new > src/${digitizerImp}_new2
rm -f src/${digitizerImp}_new
mv src/${digitizerImp}_new2 src/${digitizerImp}

### Proceed with activation of new digitizer at AllPixEventAction::SetupDigitizers()
echo "[INFO] including new digitizer in the setup"
# save old src/AllPixSetupDigitizers.cc file
mv src/AllPixSetupDigitizers.cc src/AllPixSetupDigitizers.cc_save
# Create header crop from the existing AllPixSetupDigitizers.cc
sed -n '1,/__endofheader__/ p' < src/AllPixSetupDigitizers.cc_save > src/__header.cc_crop
# delete the key line
sed '/__endofheader__/d' < src/__header.cc_crop > src/digitizers_base/AllPixSetupDigitizers.cc_headers
rm -f src/__header.cc_crop
# copy headers and include the new one
cp src/digitizers_base/AllPixSetupDigitizers.cc_headers src/AllPixSetupDigitizers.cc
rm -f src/digitizers_base/AllPixSetupDigitizers.cc_headers
cat <<EOF>>src/AllPixSetupDigitizers.cc
// Included by newdigitizer.sh script --> ${digiName}
#include "${nameFragment}Digitizer.hh"
// __endofheader__
EOF

# copy the main code
cat src/digitizers_base/AllPixSetupDigitizers.cc_head >> src/AllPixSetupDigitizers.cc

# Keep the old digits !!!
# print section of file between two regular expressions (inclusive) 
sed -n '/__beginofdigitlist__/,/__endofdigitlist__/p'< src/AllPixSetupDigitizers.cc_save > src/__keepdigitlist.cc_crop
# delete line containing __beginofdigitlist__
sed '/__beginofdigitlist__/d' < src/__keepdigitlist.cc_crop > src/__keepdigitlist.cc_crop1
rm -f src/__keepdigitlist.cc_crop
# delete line containing __endofdigitlist__
sed '/__endofdigitlist__/d' < src/__keepdigitlist.cc_crop1 > src/__keepdigitlist.cc_crop2
rm -f src/__keepdigitlist.cc_crop1
# Append the good chunk into the new SetupDigitizers code
cat src/__keepdigitlist.cc_crop2 >> src/AllPixSetupDigitizers.cc
rm -f src/__keepdigitlist.cc_crop2

# include new digi code
cat <<EOF>>src/AllPixSetupDigitizers.cc
// Included by newdigitizer.sh script --> ${digiName}
else if (digitizerName == "${digiName}") {
			${nameFragment}Digitizer * dp = new ${nameFragment}Digitizer(digitizerModulesNames[itr] , hcName, digitColectionName);
			dmPtr = static_cast<AllPixDigitizerInterface *> (dp);
			cout << "    Setting up a " << digitizerName << " digitizer for det : " << detectorId << endl;
		}
EOF
# include the tail
cat src/digitizers_base/AllPixSetupDigitizers.cc_tail >> src/AllPixSetupDigitizers.cc

echo "*************************************************************************************"
echo "README !!!"
echo "1) please recompile allpix: gmake"
echo "2) If you want to use this digitizer you need to call it in the detector description."
echo "   The name of your new digitizer is: ${digiName}"
echo "   I can be used in the corresponding tag, i.e: <digitizer>${digiName}</digitizer>"
echo "    inside the file models/pixeldetector.xml which you should modify by hand."
echo "*************************************************************************************"

echo "[DONE]"
