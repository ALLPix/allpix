import sys

from ROOT import *
import ROOT

def CompareResiduals( dataPath, simuPath, plane, direction ):

    data = dataPath
    simu = simuPath
    
    #c1 = TCanvas()
    
    f_simu = TFile( simu )
    h_simu = f_simu.Get( 'Fitter/pl%s_residual%s' % ( plane, direction ) )

    h_simu.SetDirectory( 0 )
    h_simu.SetLineColor( 2 )
    h_simu.SetTitle( 'plane %s'%plane )
    h_simu.Scale( 1 / h_simu.Integral() )
    h_simu.Draw( '' )
    gPad.Update()
    st_simu = h_simu.GetListOfFunctions().FindObject( 'stats' )
    st_simu.SetLineColor( 2 )
    st_simu.SetTextColor( 2 )
    st_simu.Draw()
    
    f_data = TFile( data )
    h_data = f_data.Get( 'Fitter/pl%s_residual%s' % (plane, direction ) )
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
    st_data.SetY2NDC( 0.58 )
    st_data.SetY1NDC( 0.70 )

    leg = TLegend( 0.15, 0.7, 0.4, 0.85 )
    leg.SetBorderSize( 1 )
    leg.AddEntry( h_simu, 'Simulation' )
    leg.AddEntry( h_data, 'Data' )
    leg.Draw()

    

    gPad.SetLogy()
    h_simu.GetYaxis().SetRangeUser(0.0001,1.1*max([ h_simu.GetMaximum(),h_data.GetMaximum()]))
    gPad.Update()
    #raw_input( '...' )
    return h_simu,h_data,leg



def usage():
    print 'Compare residuals between data and allpix simulation'
    print 'Usage:\n python %s <data-track.root> <simu-track.root> <plane> <direction>' % sys.argv[0]


if __name__ == '__main__':
    if len( sys.argv ) < 5:
        usage()
        sys.exit( 1 )

    c1= TCanvas()

    if sys.argv[3]=='all' : 
	simu = []
	data = []
 	legends = []

	c1.Divide(3,2)
	for i in range(6) : 
		c1.cd(i+1)
		h_simu,h_data,leg = CompareResiduals( sys.argv[1], sys.argv[2], i , sys.argv[4] )
		simu.append(h_simu)
		data.append(h_data)
		legends.append(leg)

		#h_simu.SetTitle("plane %i"%i)		
		#h_simu.Draw()
		#h_data.Draw("same")
       		#leg.Draw()
	a=raw_input()

    else:
    	h_simu,h_data,leg = CompareResiduals( sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4] )

	h_simu.Draw()
	h_data.Draw("same")
        leg.Draw()
	a=raw_input()

    if len(sys.argv)>5 : 
	print "Saving plot as %s"%sys.argv[5]
	c1.SaveAs(sys.argv[5])







