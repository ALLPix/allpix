cd ..
rm -rf allpix-build
rm -rf allpix-install
mkdir allpix-build allpix-install
cd allpix-build 
cmake ../allpix -DCMAKE_INSTALL_PREFIX=$G4WORKDIR
make -j4 install