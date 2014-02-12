'''
Created on July 19, 2013

Tool to create dummy alignment db file for EUTelescope reconstruction of allpix simulation.

@author: <a href="mailto:samir.arfaoui@cern.ch">Samir Arfaoui</a>
'''

from pyLCIO import EVENT, IMPL, IOIMPL, UTIL
from ROOT import TVector3, TLorentzVector, TRandom3, TMath, std
from time import time

import sys, math, glob
from collections import defaultdict

#########################################################################################
def doAlignment( outputFileName ):
    
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
def usage():
    print 'Create dummy alignment DB LCIO file'
    print 'Usage:\n python %s <outputFileName>' % ( sys.argv[0] )

#########################################################################################
if __name__ == '__main__':
    if len( sys.argv ) < 2:
        usage()
        sys.exit( 1 )
    doAlignment( sys.argv[1] )
