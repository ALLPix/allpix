'''
Writes the macro for the allpix simulation by giving the gear files corresponding to the placement of the telescope planes and the DUT and the alignment/prealignment files.
niloufar.alipour.tehrani@cern.ch
'''
#-------------------------------------------------#
import sys
from xml.dom import minidom    # To read XML file
import operator
#import numpy as np
#-------------------------------------------------#
#-------------------------------------------------#
def usage():
    print 'Usage: %s MacroName OutFolder GearFile Energy FrameNb'%(sys.argv[0])
    sys.exit( 2 )

#-------------------------------------------------#
def CreateFile(Energy, FrameNb, layerList, outputMacro, OutFolder) :
#     CreateFile(Energy, FrameNb, layerList, FinalAlignment, MacroName)
    file_general=open("GeneralMacro.in")		# open file GeneralMacro
    file_out=open(("%s")%outputMacro, "w")	# open file outputMacro
    
    lines = file_general.readlines()			# read line from GeneralMacro
    lines_control=[0]*2								# create strings with len=2 and fill with 0
    lines_DUT=[0]*2
    lines_telescope=[0]*2
    lines_buildDetectors=[0]*2


    for line in lines:															# read line in GeneralMacro 
        if line.find("#### Begin control ####") > -1:					# check for string segment
            lines_control[0]=lines.index(line)							# save index of string segment
        elif line.find("#### End control ####") > -1:
            lines_control[1]=lines.index(line)
        elif line.find("#### Begin DUT ####") > -1:
            lines_DUT[0]=lines.index(line)
        elif line.find("#### End DUT ####") > -1:
            lines_DUT[1]=lines.index(line)
        elif line.find("#### Begin telescope ####") > -1:
            lines_telescope[0]=lines.index(line)
        elif line.find("#### End telescope ####") > -1:
            lines_telescope[1]=lines.index(line)
        elif line.find("#### Begin build detectors ####") > -1:
            lines_buildDetectors[0]=lines.index(line)
        elif line.find("#### End build detectors ####") > -1:
            lines_buildDetectors[1]=lines.index(line)

    # --- Control ---
    file_out.write("".join(lines[lines_control[0]+1:lines_control[1]])) # write in outputfile, join lines with intermediate "" 
    																							# from start index to end index

  
    # --- Telescope ---
    file_out.write("#### Telescope planes #### \n")					# write "" in outputMacro
    for i in range(0, len(layerList)):									# loop over i in (0,length of layerlist)
        print "EUTel", i												
        ladder = layerList[i].getElementsByTagName("ladder")[0]
        sensitive = layerList[i].getElementsByTagName("sensitive")[0]
        planeNB=int(ladder.attributes['ID'].value)        
        print "Eutel plane:", planeNB
        file_out.write(("#### EUD%i #### \n")%planeNB)
        for line in lines[lines_telescope[0]+1:lines_telescope[1]]:
	         lineM=line
	         lineM=lineM.replace(("@ID@"), str(300+planeNB))				# replace "@attr@" with values from Alignment
	         lineM=lineM.replace(("@EUDpositionX@"), ladder.attributes['positionZ'].value)
	         lineM=lineM.replace(("@EUDpositionY@"), ladder.attributes['positionY'].value)
	         lineM=lineM.replace(("@EUDpositionZ@"), ladder.attributes['positionZ'].value)
	         lineM=lineM.replace(("@EUDalpha@"), ladder.attributes['rotationZY'].value)
	         lineM=lineM.replace(("@EUDbeta@"), ladder.attributes['rotationZX'].value)
	         lineM=lineM.replace(("@EUDgamma@"), ladder.attributes['rotationXY'].value)
	         lineM=lineM.replace(("@EUDalphaSensor@"), sensitive.attributes['rotation1'].value)
	         lineM=lineM.replace(("@EUDbetaSensor@"), sensitive.attributes['rotation1'].value) 
	         lineM=lineM.replace(("@EUDgammaSensor@"), sensitive.attributes['rotation1'].value)
	         file_out.write(lineM)

    # --- Build detectors ---
    for line in lines[lines_buildDetectors[0]+1:lines_buildDetectors[1]]:
        line=line.replace("@FolderPath@", OutFolder)
        line=line.replace("@BeamEnergy@", str(Energy))
        line=line.replace("@FrameNB@", str(FrameNb))
        file_out.write(line)

    file_out.close()


if __name__ == '__main__':					# script called as main function

    if len(sys.argv) < 5:					# print usage() if arguments incomplete
        usage()
    												
    MacroName=sys.argv[1]					# read in arguments
    OutFolder=sys.argv[2]
    GearFile=sys.argv[3]
    Energy=float(sys.argv[4])
    FrameNb=int(sys.argv[5])

   
#-------------------------------------------------------------------#
    #1. Create the macro
   
    #XML file 
    xmldoc = minidom.parse(GearFile)							# read layer from GearFile
    layerList = xmldoc.getElementsByTagName('layer') 
  
    CreateFile(Energy, FrameNb, layerList, MacroName, OutFolder)




