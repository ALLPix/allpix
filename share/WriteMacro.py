'''
Writes the macro for the allpix simulation by giving the gear files corresponding to the placement of the telescope planes and the DUT and the alignment/prealignment files.
niloufar.alipour.tehrani@cern.ch
'''
#-------------------------------------------------#
import sys
from xml.dom import minidom    # To read XML file
from pyLCIO import IOIMPL, EVENT
import operator
import numpy as np
#-------------------------------------------------#
#-------------------------------------------------#
def usage():
    print 'Usage: %s MacroName OutFolder preAlign Align DUTAlign GearFile Energy FrameNb'%(sys.argv[0])
    sys.exit( 2 )
#-------------------------------------------------#
def readAlignment(fileName): #From Samir #takes the .slcio files for the alignment and the prealignment
    #
    #planeData = defaultdict(list)
    alignmentVals=[]
    #
    reader = IOIMPL.LCFactory.getInstance().createLCReader()
    reader.open( fileName )
    nEvents = 0
    for event in reader:
        #
        # from pyLCIO import UTIL
        # UTIL.LCTOOLS.dumpEventDetailed( event )
        #
        for collectionName in event.getCollectionNames():
            collection = event.getCollection( collectionName )
            nElements = collection.getNumberOfElements()
            for obj in collection: 
                sensorID  = obj.getIntVal(0)
                xOff      = obj.getDoubleVal(0)
                yOff      = obj.getDoubleVal(1)
                zOff      = obj.getDoubleVal(2)
                alpha     = obj.getDoubleVal(3)
                beta      = obj.getDoubleVal(4)
                gamma     = obj.getDoubleVal(5)
                xOff_err  = obj.getDoubleVal(6)
                yOff_err  = obj.getDoubleVal(7)
                zOff_err  = obj.getDoubleVal(8)
                alpha_err = obj.getDoubleVal(9)
                beta_err  = obj.getDoubleVal(10)
                gamma_err = obj.getDoubleVal(11)
                alignmentVals.append([xOff, yOff, zOff, alpha, beta, gamma])
            nEvents += 1
    reader.close()

    return alignmentVals
#-------------------------------------------------#
def CreateFile(Energy, FrameNb, ladderList, FinalAlignment, outputMacro, OutFolder, DUTAlign) :
#     CreateFile(Energy, FrameNb, ladderList, FinalAlignment, MacroName)
    file_general=open("GeneralMacro.in")
    file_out=open(("%s")%outputMacro, "w")
    
    lines = file_general.readlines()
    lines_control=[0]*2
    lines_DUT=[0]*2
    lines_telescope=[0]*2
    lines_buildDetectors=[0]*2


    for line in lines:
        if line.find("#### Begin control ####") > -1:
            lines_control[0]=lines.index(line)
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
    file_out.write("".join(lines[lines_control[0]+1:lines_control[1]]))

    # --- DUT ---
    # Alignment for the DUT 
    f = open(DUTAlign)
    lineDUT = f.readlines()
    f.close()
    # if len(ladderList)>6:
    #     file_out.write("#### DUT #### \n")
    #     for line in lines[lines_DUT[0]+1:lines_DUT[1]]:
    #         lineM=line
    #         lineM=lineM.replace(("@DUTId@"), str(2042)) # SET the id of the DUT manually
    #         lineM=lineM.replace(("@DUTpositionX@"), lineDUT[0].split()[3])
    #         lineM=lineM.replace(("@DUTpositionY@"), lineDUT[0].split()[4])
    #         lineM=lineM.replace(("@DUTpositionZ@"), ladderList[6].attributes['positionZ'].value)
    #         lineM=lineM.replace(("@DUTalpha@"), lineDUT[0].split()[0])
    #         lineM=lineM.replace(("@DUTbeta@"), str(180+float(lineDUT[0].split()[1])))
    #         lineM=lineM.replace(("@DUTgamma@"), str(90+float(lineDUT[0].split()[2])))
            
    #         file_out.write(lineM)
    # --- Telescope ---
    file_out.write("#### Telescope planes #### \n")
    for i in range(0, len(ladderList)):
        planeNB=int(ladderList[i].attributes['ID'].value)
        if planeNB==6:
            print "DUT", i
            file_out.write("#### DUT #### \n")
            for line in lines[lines_DUT[0]+1:lines_DUT[1]]:
                lineM=line
                lineM=lineM.replace(("@DUTId@"), str(2042)) # SET the id of the DUT manually
                lineM=lineM.replace(("@DUTpositionX@"), lineDUT[0].split()[3])
                lineM=lineM.replace(("@DUTpositionY@"), lineDUT[0].split()[4])
                lineM=lineM.replace(("@DUTpositionZ@"), ladderList[i].attributes['positionZ'].value)
                lineM=lineM.replace(("@DUTalpha@"), lineDUT[0].split()[0])
                lineM=lineM.replace(("@DUTbeta@"), str(180+float(lineDUT[0].split()[1])))
                lineM=lineM.replace(("@DUTgamma@"), str(90+float(lineDUT[0].split()[2])))
                
                file_out.write(lineM)
        else:
            print "EUTel", i
            #planeNB=int(ladderList[i].attributes['ID'].value)
            print "Eutel plane:", planeNB
            file_out.write(("#### EUD%i #### \n")%planeNB)
            for line in lines[lines_telescope[0]+1:lines_telescope[1]]:
                lineM=line
                lineM=lineM.replace(("@ID@"), str(300+planeNB))
                lineM=lineM.replace(("@EUDpositionX@"), str(FinalAlignment[planeNB][0]))
                lineM=lineM.replace(("@EUDpositionY@"), str(FinalAlignment[planeNB][1]))
                lineM=lineM.replace(("@EUDpositionZ@"), ladderList[i].attributes['positionZ'].value)
                lineM=lineM.replace(("@EUDalpha@"), str(FinalAlignment[planeNB][3]))
                lineM=lineM.replace(("@EUDbeta@"), str(180+FinalAlignment[planeNB][4]))
                lineM=lineM.replace(("@EUDgamma@"), str(FinalAlignment[planeNB][5]))
                lineM=lineM.replace(("@EUDalphaSensor@"), str(FinalAlignment[planeNB][3]))
                lineM=lineM.replace(("@EUDbetaSensor@"), str(180+FinalAlignment[planeNB][4])) 
                lineM=lineM.replace(("@EUDgammaSensor@"), str(180+FinalAlignment[planeNB][5]))
            
                file_out.write(lineM)



    # --- Build detectors ---
    for line in lines[lines_buildDetectors[0]+1:lines_buildDetectors[1]]:
        line=line.replace("@FolderPath@", OutFolder)
        line=line.replace("@BeamEnergy@", str(Energy))
        line=line.replace("@FrameNB@", str(FrameNb))
        file_out.write(line)

    file_out.close()


if __name__ == '__main__':

    if len(sys.argv) < 8:
        usage()
    
    MacroName=sys.argv[1]
    OutFolder=sys.argv[2]
    preAlignFile=sys.argv[3]
    AlignFile=sys.argv[4]
    DUTAlign=sys.argv[5]
    GearFile=sys.argv[6]
    Energy=float(sys.argv[7])
    FrameNb=int(sys.argv[8])

    #The alignment file should contain: rotation1 rotation2 rotation3 translationX translationY
#-------------------------------------------------------------------#
    #1. Create the macro
    # READ ALIGNMENTS
    #pathAlign="/VertexScratch/workspace/nalipour/TestBeam/"+step+"/db"
    filesAlignment=[]
    filesAlignment.append(preAlignFile)
    filesAlignment.append(AlignFile)
 
    preAlign=readAlignment(filesAlignment[0])
    Align=readAlignment(filesAlignment[1])
    FinalAlignment=[[0 for i in range(len(preAlign))] for j in range(len(preAlign))]

    # Alignment for the telescope
    for i in range(len(preAlign)):
        for j in range(len(preAlign[i])):
            FinalAlignment[i][j]=(preAlign[i][j]+Align[i][j])
    
  

    #XML file 
    xmldoc = minidom.parse(GearFile)
    ladderList = xmldoc.getElementsByTagName('ladder') 


    CreateFile(Energy, FrameNb, ladderList, FinalAlignment, MacroName, OutFolder, DUTAlign)
