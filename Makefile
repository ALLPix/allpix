# --------------------------------------------------------------
# Makefile for examples module
# --------------------------------------------------------------

#export EUTELESCOPE=1

name = allpix
G4TARGET = $(name)
G4EXLIB = true


.PHONY: all
all: rootdict lib bin

rootdict:	./include/*.hh ./include/*.h
	@echo "Generating dictionary $@..."
	rootcint -v0 -f SelDict.cc -c -p -I./include \
		AllPix_Hits_WriteToEntuple.h \
		AllPix_Frames_WriteToEntuple.h allpix_dm.h \
		AllPixDigitAnimation.hh LinkDef.h
	@mv SelDict.cc ./src
	@mv SelDict.h ./include/
# the geant4 makefile will compile the dictionary

cleanroot:
	rm -f src/SelDict.* include/SelDict.* SelDict.* 

include $(G4INSTALL)/config/binmake.gmk
include ./flags.gmk

ifdef EUTELESCOPE
	INCFLAGS += -D_EUTELESCOPE
else
	INCFLAGS += -D_SINGLE
endif

