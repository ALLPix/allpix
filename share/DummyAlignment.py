'''
Created on July 19, 2013

Tool to create dummy alignment db file for EUTelescope reconstruction of allpix simulation.

@author: <a href="mailto:samir.arfaoui@cern.ch">Samir Arfaoui</a>
'''

from pyLCIO import EVENT, IMPL, IOIMPL, UTIL
from ROOT import TVector3, TLorentzVector, TRandom3, TMath, std
from time import time

import sys, math, glob
import collections
from collections import defaultdict


def ReadAlign(filename) :

    constants = {}
    aFile = open(filename)    
    lines = aFile.readlines()
    
    for line in lines : 
        words = line.split()
        constants[words[0]] = [words[1],words[2],words[3],words[4],words[5],words[6]]
    
    return collections.OrderedDict(sorted(constants.items()))
    
#########################################################################################
def doAlignmentDummy( outputFileName ):
    
    # create a writer and open the output file
    writer = IOIMPL.LCFactory.getInstance().createLCWriter()
    writer.open( outputFileName, EVENT.LCIO.WRITE_NEW )

    # create an event and set its parameters
    event = IMPL.LCEventImpl()
    event.setEventNumber( 0 )
    event.setRunNumber( 0 )
    event.setTimeStamp( int( time() * 1000000000. ) )
    
    alignmentColl = IMPL.LCCollectionVec( EVENT.LCIO.LCGENERICOBJECT )
            
    # collection parameters
    alignmentColl.parameters().setValue( 'DataDescription', 'sensorID, xOff, yOff, zOff, alpha, beta, gamma + 13 spare fields' )
    alignmentColl.parameters().setValue( 'TypeName', 'Alignment Constant' )
    
    # number of items
    nValues = 20
        
    for sensorID in xrange( 4, 10 ):

        alignObj = IMPL.LCGenericObjectImpl( nValues, 0, nValues )
        
        for i in xrange( nValues ):

            alignObj.setIntVal( i, 0 )
            alignObj.setDoubleVal( i, 0. )     
        
        alignmentColl.addElement( alignObj )

    event.addCollection( alignmentColl, 'alignment' )

    writer.writeEvent( event )
    
    writer.flush()
    writer.close()


#########################################################################################
def doAlignment( outputFileName,alignement_File):
    
    
    constants = ReadAlign(alignement_File)
    print constants
    
    # create a writer and open the output file
    writer = IOIMPL.LCFactory.getInstance().createLCWriter()
    writer.open( outputFileName, EVENT.LCIO.WRITE_NEW )

    # create an event and set its parameters
    event = IMPL.LCEventImpl()
    event.setEventNumber( 0 )
    event.setRunNumber( 0 )
    event.setTimeStamp( int( time() * 1000000000. ) )
    
    alignmentColl = IMPL.LCCollectionVec( EVENT.LCIO.LCGENERICOBJECT )
            
    # collection parameters
    alignmentColl.parameters().setValue( 'DataDescription', 'sensorID, xOff, yOff, zOff, alpha, beta, gamma + 13 spare fields' )
    alignmentColl.parameters().setValue( 'TypeName', 'Alignment Constant' )
    
    # number of items
    nValues = 20
        
    for sensor in constants:

        alignObj = IMPL.LCGenericObjectImpl( nValues, 0, nValues )
        
        for i in xrange( nValues ):

            if(i==0) : 
                alignObj.setIntVal( i, int(sensor) )
                alignObj.setDoubleVal( i, float(constants[sensor][i]) )             
            elif(i<6) : 
                alignObj.setIntVal( i, 0 )
                alignObj.setDoubleVal( i, float(constants[sensor][i]) )  
            else : 
                alignObj.setIntVal( i, 0 )
                alignObj.setDoubleVal( i, 0. )                 
        
        alignmentColl.addElement( alignObj )

    event.addCollection( alignmentColl, 'alignment' )

    writer.writeEvent( event )
    
    writer.flush()
    writer.close()

#########################################################################################
def usage():
    print 'Create dummy alignment DB LCIO file'
    print 'Usage:\n python %s <outputFileName> <alignement_File>' % ( sys.argv[0] )

#########################################################################################
if __name__ == '__main__':
    if len( sys.argv ) < 2 or len( sys.argv ) > 3:
        usage()
        sys.exit( 1 )
    
    if len( sys.argv ) == 3 :
        doAlignment( sys.argv[1],sys.argv[2])
    else : 
        doAlignmentDummy( sys.argv[1])    
    
