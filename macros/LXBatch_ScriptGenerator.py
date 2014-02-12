from ROOT import  *
import re
import os


""" 
Generate scripts to launch allpix on lxbatch for a list of X and Y rotation of the sensor. Execution script must be first generated using ScriptGenerator.py in macros folder
"""
def GenerateScripts(subruns,tr_per_frame,nTrigger):
        batchScripts = []
        macro_path = "/afs/cern.ch/user/m/mbenoit/scratch0/Timepix_Telescope_Simulation/macros/"

	i=0
	
	while i <subruns :

	  batchScriptName = "launch_run%d.sh"%(i)
	  batchScripts.append(batchScriptName)
	  outfile= open(batchScriptName,"w")
	  outfile.write("source /afs/cern.ch/user/m/mbenoit/setup_lxplus.sh\n")
	  outfile.write("mkdir /tmp/mbenoit \n")
	  outfile.write("cd /afs/cern.ch/eng/clic/software/Pixel_TestBeam_Software/allpix  \n")                                                
	  for j in range(1):
	  	ExecFileName = "TimepixTelescope_run%d_%dperFrame_%dFrames.in"%(i,tr_per_frame,nTrigger)
		if(i<subruns):
			outfile.write("allpix "+macro_path+ExecFileName+" batch \n")
		outfile.write("cp /tmp/mbenoit/Timepix_Telescope_run%d_* /afs/cern.ch/user/m/mbenoit/scratch0/Timepix_Telescope_Simulation/allpix_output \n"%i)
		i+=1
	  
	  
	  outfile.close()
	  os.system("chmod a+rwx %s"%batchScriptName)
        
        return batchScripts

"""
Launch the generated script to lxbatch
"""

def LaunchBatch(batchScripts,queue):
    
    for script in batchScripts:
        os.system("bsub -R \"rusage[mem=3500]\"-q %s %s \n"%(queue,script))             
                          

list=GenerateScripts(100,100,1000) 

LaunchBatch(list,"1nh")
