import sys

from ROOT import *
import ROOT

def CompareResiduals( dataPath, simuPath, plane ):

    data = dataPath
    simu = simuPath
    
    #c1 = TCanvas()
    
    f_simu = TFile( simu )
    h_simu = f_simu.Get( 'Clustering/detector_%s/clusterSignal_d%s' % ( plane, plane ) )
    h_simu.SetDirectory( 0 )
    h_simu.SetLineColor( 2 )
    h_simu.SetTitle( 'Plane %i' % int( plane ) )
    h_simu.Scale( 1 / h_simu.Integral() )
    h_simu.Draw( '' )
    h_simu.GetXaxis().SetRangeUser( 0, 10 )
    gPad.Update()
    st_simu = h_simu.GetListOfFunctions().FindObject( 'stats' )
    st_simu.SetLineColor( 2 )
    st_simu.SetTextColor( 2 )
    st_simu.Draw()
    
    f_data = TFile( data )
    h_data = f_data.Get( 'Clustering/detector_%s/clusterSignal_d%s' % ( plane, plane ) )
    h_data.SetDirectory( 0 )
    h_data.SetLineColor( 4 )
    h_data.SetTitle( 'data' )
    h_data.Scale( 1 / h_data.Integral() )
    h_data.Draw( 'sames' )
    gPad.Update()
    st_data = h_data.GetListOfFunctions().FindObject( 'stats' )
    st_data.SetLineColor( 4 )
    st_data.SetTextColor( 4 )
    st_data.Draw()
    st_data.SetY1NDC( 0.75 )
    st_data.SetY2NDC( 0.59 )

    leg = TLegend( 0.55, 0.7, 0.75, 0.85 )
    leg.SetBorderSize( 1 )
    leg.AddEntry( h_simu, 'Simulation' )
    leg.AddEntry( h_data, 'Data' )
    leg.Draw()

    gPad.Update()

    return h_simu, h_data, leg
    #raw_input( '...' )

def usage():
    print 'Compare cluster size between data and allpix simulation'
    print 'Usage:\n python %s <data-clu.root> <simu-clu.root>' % sys.argv[0]


if __name__ == '__main__':
    if len( sys.argv ) < 3:
        usage()
        sys.exit( 1 )

    simu = []
    data = []
    legs = []

    c = TCanvas( 'c', 'c' )
    c.Divide( 3, 2 )
    for i in xrange( 1, 7 ):
        c.cd( i )
        h_simu, h_data, leg = CompareResiduals( sys.argv[1], sys.argv[2], i-1 )
        simu.append( h_simu )
        data.append( h_data )
        legs.append( leg )
        
    raw_input( '...' )
