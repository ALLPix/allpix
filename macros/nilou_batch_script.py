from ROOT import  *
import re
import os

"""
nalipour: Generate scripts to run on the lxbatch
"""


def GenerateScript(run):
    allpixPath="/afs/cern.ch/work/n/nalipour/allpixCourse/svn_nalipour_allpix"
    macro="%s/macros/run%s.in"%(allpixPath, run)
    batchScriptName = "/afs/cern.ch/work/n/nalipour/testBeamAnalysis/Simulation/lxBatchScripts/launch_run%s.sh"%(run)
    

    outfile= open(batchScriptName,"w")
    outfile.write("echo Nilou \n")
    outfile.write("pwd \n")
    #outfile.write("ls run%s \n"%(run))
    #outfile.write("rm -rf run%s \n"%(run))
    outfile.write("source %s/setup_allpix.sh \n"%(allpixPath))
    outfile.write("cd %s \n"%(allpixPath))
    outfile.write("allpix "+macro+" batch \n")

    #outfile.write("pwd  \n")
    #outfile.write("cp -r /tmp/run%s /afs/cern.ch/work/n/nalipour/testBeamAnalysis/Simulation/ALLPix_results/. \n"%(run))
    outfile.write("cd /tmp \n")
    #outfile.write("pwd \n")
    outfile.write("find run%s -name \"*.txt\" -print > allfiles \n"%(run))
    outfile.write("tar czvf run%s.tar.gz --files-from allfiles \n"%(run))
    outfile.write("cp run%s.tar.gz /afs/cern.ch/work/n/nalipour/testBeamAnalysis/Simulation/ALLPix_results/. \n"%(run))



    outfile.close()
    os.system("chmod a+rwx %s"%batchScriptName)
    return batchScriptName



run="000007"
queue="1nd"
logname="/afs/cern.ch/work/n/nalipour/testBeamAnalysis/Simulation/lxBatchScripts/run%s.log"%(run)
batchScriptName=GenerateScript(run)
os.system("bsub -C0 -q %s -o %s -e %s -L /bin/bash %s \n"%(queue, logname, logname, batchScriptName))

