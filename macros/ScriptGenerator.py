
import re
import os
from math import *

	
"""
Generate an allpix execution Script, using the oneCLICDetector.in_fragment template, using the provided parameters
"""

def ScriptGenerator(nTrigger):
    
    filename = "oneDetector_Am241gamma.in_fragment"
    f = open(filename, "r")
    text = f.read();
  
    for i in range(1,nTrigger) :
        text +="/run/beamOn 1000\n"  
#    print text
    
    outfileName = "Calibration_Am241_%d.in"%(nTrigger)
    outfile= open(outfileName,"w")
    outfile.write(text)
    outfile.close()
    return outfileName;



def GenerateScripts(liste):
        batchScripts = []
        macro_path = "/afs/cern.ch/user/m/mbenoit/scratch0/SmallPix_65nm_Frame/macros/"
	i=0
	
	for f in liste : 
	  batchScriptName = "launch_run%d.sh"%(i)
	  batchScripts.append(batchScriptName)
	  outfile= open(batchScriptName,"w")	
	  outfile.write("source /afs/cern.ch/user/m/mbenoit/setup_lxplus.sh\n")
	  outfile.write("mkdir /tmp/mbenoit \n")
	  outfile.write("rm /tmp/mbenoit/* \n")	  
	  outfile.write("cd /afs/cern.ch/eng/clic/software/Pixel_TestBeam_Software/allpix  \n")                                                
  	  ExecFileName = f
	  print f
	  outfile.write("allpix "+macro_path+ExecFileName+" batch \n")	
	  outfile.write("cp /tmp/mbenoit/* /afs/cern.ch/user/m/mbenoit/scratch0/SmallPix_65nm_Frame/allpix_output \n")
	  outfile.close()
	  os.system("chmod a+rwx %s"%batchScriptName)
          i+=1
	  
        return batchScripts

"""
Launch the generated script to lxbatch
"""

def LaunchBatch(batchScripts,queue):
    
    for script in batchScripts:
        os.system("bsub -o /dev/null -e /dev/null -q %s %s \n"%(queue,script))   




ScriptGenerator(25000)        
        
#batchs=GenerateScripts(liste)       
        
#LaunchBatch(batchs,"1nh")
