import re
import os
from math import *




"""
Generate an allpix execution Script, using the oneCLICDetector.in_fragment template, using the provided parameters
"""

def BuildLadder(radius,angle,NChips,length,tilt):
    


	x=[]
	y=[]
	z=[]
	rotX=[]
	rotY=[]
	rotZ=[]


	for i in range(NChips) : 


	#center of the sensor position
		z.append((length/2.+(-length*NChips/2.)) + i*length)
		x.append(radius*cos(radians(angle)))
		y.append(radius*sin(radians(angle)))

	#rotations
		rotX.append(90.)
		rotY.append(angle+90+tilt)
		rotZ.append(0)


		#print "x=%f y=%f z=%f \n"%(x[i],y[i],z[i])

	return x,y,z,rotX,rotY,rotZ

	    


def BuildBarrel(radius,nLadder,nChipsPerLadder,ChipLength,tilt) :

	x=[]
	y=[]
	z=[]
	rotX=[]
	rotY=[]
	rotZ=[]
	
	for i in range(nLadder) :
		x.append([])
		y.append([])
		z.append([])
		rotX.append([])			
		rotY.append([])
		rotZ.append([])
		
	for i in range(nLadder) :
		x[i],y[i],z[i],rotX[i],rotY[i],rotZ[i]=BuildLadder(radius,i*360./nLadder,nChipsPerLadder,ChipLength,tilt)

	return x,y,z,rotX,rotY,rotZ
	    
def BuildDisk(innerRadius,zPosition,nModules,theta0,Length,tilt) :

	x=[]
	y=[]
	z=[]
	rotX=[]
	rotY=[]
	rotZ=[]
	
	xi=0
	yi=0
	zi=0
	thetai=theta0
	
	effRadius = innerRadius + 0.5*Length
	for i in range(nModules) :
		x.append(0)
		y.append(0)
		z.append(0)
		rotX.append(0)
		rotY.append(0)		
		rotZ.append(0)	
		
			
	for i in range(nModules) :
		
		thetai = (i)*360./nModules + theta0
		xi = effRadius*cos(radians(thetai))
		yi = effRadius*sin(radians(thetai))
		zi=  zPosition
		
		x[i],y[i],z[i],rotX[i],rotY[i],rotZ[i]= xi,yi,zi,0,tilt,90-thetai
		
		print "[Module] x=%3.3f y=%3.3f z=%3.3f Theta=%3.3f"%(xi,yi,zi,thetai)
		
	return x,y,z,rotX,rotY,rotZ


def BuildHeader(outfile) :
	outfile.write("#Script Generated using BuildTracker script \n")

def BuildFooter(outfile,HEPFile,i) :
	
	outfile.write("/allpix/extras/setTestStructurePosition 0. 0.  0. mm \n")
	outfile.write("/allpix/extras/setTestStructurePosition 0. 0.  0. mm \n")	
	outfile.write("/allpix/extras/setPeakField 4 T \n")
	
	
	
	outfile.write("/allpix/config/setOutputPrefixWithPath /afs/cern.ch/eng/clic/work/mbenoit/Tracker_allpix_results/Tracker%02d \n"%i)	
	
	outfile.write("/allpix/phys/Physics LIVERMORE_FTFP_BERT \n")
	outfile.write("/run/initialize \n")
	outfile.write("/allpix/config/HEPEvtFile %s \n"%HEPFile)
	
	
	outfile.write("/allpix/det/update\n")
	outfile.write("/vis/scene/create\n/vis/scene/add/axes 0. 0. 0. 10. cm\n\n#/vis/scene/add/volume World -1 2\n/vis/scene/add/volume World -1 2\n/vis/viewer/set/style s\n\n/run/verbose 0\n/control/verbose 0\n/control/saveHistory\n/tracking/verbose 0\n/allpix/phys/verbose 0\n\n#/vis/open OIX 1024x768-100+100\n#/vis/open OGLIXm 1024x768-100+100\n\n#/vis/open RayTracer\n#/vis/open OGLIQt\n#/vis/open OGLSQt\n#/vis/open OIX\n#/vis/open OGLIX 1024x768-100+100\n/vis/open OGLSXm\n#/vis/open DAWNFILE\n#/vis/open OGLSX\n#/vis/open OGL 600x600-0+0\n\n/vis/viewer/set/background 0.4 0.5 0.6\n/vis/viewer/set/viewpointThetaPhi 20 50\n#/vis/viewer/set/background 0 0 0\n/vis/viewer/zoom 2.0\n\n/vis/viewer/flush\n\n#\n# Draw trajectories at end of event, showing trajectory points as\n# markers of size 2 pixels\n/vis/scene/add/trajectories\n/vis/modeling/trajectories/create/drawByCharge\n/vis/modeling/trajectories/drawByCharge-0/default/setDrawStepPts true\n/vis/modeling/trajectories/drawByCharge-0/default/setStepPtsSize 2\n/vis/scene/endOfEventAction accumulate 1000 \n")
	#outfile.write("# GPS\n/gps/particle e-\n/gps/pos/type Plane\n#/gps/pos/rot1 0 0 1\n#/gps/pos/rot2 1 0 0\n/gps/pos/shape Rectangle\n/gps/pos/centre 0.0 0.0 0 mm\n/gps/pos/halfy <hY> um\n/gps/pos/halfx 1600. um\n/gps/direction  0 0 1\n/gps/ene/type User\n/gps/hist/type energy\n# spectra\n/gps/hist/point 100000 1\n/gps/source/list\n\n")

def BeamOn(outfile,nEvents) :
	for i in range(nEvents):
		outfile.write("/run/beamOn 1 \n")


def BuildGeometry(x,y,z,rotX,rotY,rotZ,ID0,outfile) :

	
	#outfile= open("tracker.in","w")
	template = "/allpix/det/setId        <ID>\n/allpix/det/setPosition <X> <Y> <Z> mm\n/allpix/det/setRotation  <rotX> <rotY> <rotZ> deg\n/allpix/det/setLowTHL 13. keV\n\n"


	ID=ID0
	ids=[]
		
	
	for i in range(len(x)):
		for j in range(len(x[i])) :
			det = template	
			
			det=det.replace("<X>","%f"%x[i][j])
			det=det.replace("<Y>","%f"%y[i][j])
			det=det.replace("<Z>","%f"%z[i][j])
			det=det.replace("<ID>","%d"%ID)
			ids.append(ID)
			ID+=1
			det=det.replace("<rotX>","%f"%rotX[i][j])
			det=det.replace("<rotY>","%f"%rotY[i][j])			
			det=det.replace("<rotZ>","%f"%rotZ[i][j])
			outfile.write(det)
			#print det
			
def BuildGeometryDisk(x,y,z,rotX,rotY,rotZ,ID0,outfile) :

	
	#outfile= open("tracker.in","w")
	template = "/allpix/det/setId        <ID>\n/allpix/det/setPosition <X> <Y> <Z> mm\n/allpix/det/setRotation  <rotX> <rotY> <rotZ> deg\n/allpix/det/setLowTHL 13. keV\n\n"


	ID=ID0
	ids=[]
	
	for i in range(len(x)):
		det = template	
		det=det.replace("<X>","%f"%x[i])
		det=det.replace("<Y>","%f"%y[i])
		det=det.replace("<Z>","%f"%z[i])
		det=det.replace("<ID>","%d"%ID)
		ids.append(ID)
		ID+=1
		det=det.replace("<rotX>","%f"%rotX[i])
		det=det.replace("<rotY>","%f"%rotY[i])			
		det=det.replace("<rotZ>","%f"%rotZ[i])
		outfile.write(det)
		#print det
		
def BuildGeometryDiskMirror(x,y,z,rotX,rotY,rotZ,ID0,outfile) :

	
	#outfile= open("tracker.in","w")
	template = "/allpix/det/setId        <ID>\n/allpix/det/setPosition <X> <Y> <Z> mm\n/allpix/det/setRotation  <rotX> <rotY> <rotZ> deg\n/allpix/det/setLowTHL 13. keV\n\n"


	ID=ID0
	ids=[]
	
	for i in range(len(x)):
		det = template	
		det=det.replace("<X>","%f"%x[i])
		det=det.replace("<Y>","%f"%y[i])
		det=det.replace("<Z>","%f"%-z[i])
		det=det.replace("<ID>","%d"%ID)
		ids.append(ID)
		ID+=1
		det=det.replace("<rotX>","%f"%rotX[i])
		det=det.replace("<rotY>","%f"%(180+rotY[i]))			
		det=det.replace("<rotZ>","%f"%(-rotZ[i]))
		outfile.write(det)
		#print det		

def GenerateScript(script,i):
        

        macro_path = "/afs/cern.ch/user/m/mbenoit/scratch0/Full_Tracker_Model/macros/"
	batchScriptName = "launch_run%d.sh"%(i)

	outfile= open(batchScriptName,"w")	
	outfile.write("source /afs/cern.ch/user/m/mbenoit/setup_lxplus.sh\n")
	#outfile.write("mkdir /tmp/mbenoit \n")	
	outfile.write("cd /afs/cern.ch/eng/clic/software/Pixel_TestBeam_Software/allpix  \n")                                                
	outfile.write("allpix "+macro_path+script+" batch \n")
	#outfile.write("cp Tracker* /afs/cern.ch/eng/clic/work/mbenoit/Tracker_allpix_results \n")  
	outfile.close()
	os.system("chmod a+rwx %s"%batchScriptName)
        
        return batchScriptName

"""
Launch the generated script to lxbatch
"""

def LaunchBatch(batchScripts,queue):
    
    for script in batchScripts:
        os.system("bsub -q %s %s \n"%(queue,script))             
                          








Scripts = []

filename = "disk_test.in"
print "Preaparing file %s"%filename
outfile= open(filename,"w")

BuildHeader(outfile)


x,y,z,rotX,rotY,rotZ=BuildDisk(33,160,15,0,61.44,2)
BuildGeometryDisk(x,y,z,rotX,rotY,rotZ,5000,outfile)
BuildGeometryDiskMirror(x,y,z,rotX,rotY,rotZ,5100,outfile)

x,y,z,rotX,rotY,rotZ=BuildDisk(33,161,15,12,61.44,2)
BuildGeometryDisk(x,y,z,rotX,rotY,rotZ,5200,outfile)
BuildGeometryDiskMirror(x,y,z,rotX,rotY,rotZ,5300,outfile)

BuildFooter(outfile,"/afs/cern.ch/user/m/mbenoit/scratch0/Full_Tracker_Model/HEPEVT_files/AllPairs3TeV8MeV6Deg.HEPEvt00",0)
outfile.close()


#for i in range(0,55):
#	filename = "tracker_layer0+1_occupancy_%d.in"%i
#	print "Preaparing file %s"%filename
#	outfile= open(filename,"w")
#	BuildHeader(outfile)
#	
#	x,y,z,rotX,rotY,rotZ=BuildBarrel(29,18,5,51.6,1.5)
#	BuildGeometry(x,y,z,rotX,rotY,rotZ,3000,outfile)
#	x,y,z,rotX,rotY,rotZ=BuildBarrel(30.87,18,5,51.6,1.5)
#	BuildGeometry(x,y,z,rotX,rotY,rotZ,3500,outfile)
#	
#	#ILD Layer 3+4 option 1
#	x,y,z,rotX,rotY,rotZ=BuildBarrel(41.65,13,5,51.6,1.5)
#	BuildGeometry(x,y,z,rotX,rotY,rotZ,4000,outfile)
#	x,y,z,rotX,rotY,rotZ=BuildBarrel(43.516,13,5,51.6,1.5)
#	BuildGeometry(x,y,z,rotX,rotY,rotZ,4200,outfile)
#
#	#ILD Layer 3+4 option 1
#	x,y,z,rotX,rotY,rotZ=BuildBarrel(54.91,17,5,51.6,1.5)
#	BuildGeometry(x,y,z,rotX,rotY,rotZ,4500,outfile)
#	x,y,z,rotX,rotY,rotZ=BuildBarrel(56.782,17,5,51.6,1.5)
#	BuildGeometry(x,y,z,rotX,rotY,rotZ,4800,outfile)
#
#	
#	BuildFooter(outfile,"/afs/cern.ch/user/m/mbenoit/scratch0/Full_Tracker_Model/HEPEVT_files/AllPairs3TeV8MeV6Deg.HEPEvt%02d"%i,i)
#	BeamOn(outfile,50000)
#	outfile.close()
#	Scripts.append(GenerateScript(filename,i))
#
#LaunchBatch(Scripts,"1nw")



#outfile= open("tracker_test2.in","w")
#
#BuildHeader(outfile)
#
##ILD Layer 1+2 option 1
#x,y,z,rotX,rotY,rotZ=BuildBarrel(29,18,5,51.6,1.5)
#BuildGeometry(x,y,z,rotX,rotY,rotZ,3000,outfile)
#x,y,z,rotX,rotY,rotZ=BuildBarrel(30.87,18,5,51.6,1.5)
#BuildGeometry(x,y,z,rotX,rotY,rotZ,3500,outfile)
#
###ILD Layer 3+4 option 1
##x,y,z,rotX,rotY,rotZ=BuildBarrel(41.65,13,5,51.6,1.5)
##BuildGeometry(x,y,z,rotX,rotY,rotZ,4000,outfile)
##x,y,z,rotX,rotY,rotZ=BuildBarrel(43.516,13,5,51.6,1.5)
##BuildGeometry(x,y,z,rotX,rotY,rotZ,4200,outfile)
##
###ILD Layer 3+4 option 1
##x,y,z,rotX,rotY,rotZ=BuildBarrel(54.91,17,5,51.6,1.5)
##BuildGeometry(x,y,z,rotX,rotY,rotZ,4500,outfile)
##x,y,z,rotX,rotY,rotZ=BuildBarrel(56.782,17,5,51.6,1.5)
##BuildGeometry(x,y,z,rotX,rotY,rotZ,4800,outfile)
#
#
#
#BuildFooter(outfile)
##BeamOn(outfile,7718124)
##BeamOn(outfile,2746845)


#outfile.close()


