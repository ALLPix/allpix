from ROOT import  *
import re
import os

"""
nalipour: Generate scripts to run telescope converter on the batch
"""

def GenerateScript(run, dutID):
    converterPath="/afs/cern.ch/work/n/nalipour/allpixCourse/svn_nalipour_allpix/share"
    batchScriptName = "/afs/cern.ch/work/n/nalipour/testBeamAnalysis/Simulation/lxBatchScripts/converter_run%s.sh"%(run)

    outfile= open(batchScriptName,"w")
    outfile.write("source %s/setup_pyLCIO.sh \n"%(converterPath))
    outfile.write("cd %s \n"%(converterPath))
    outfile.write("python TelescopeConverter.py /afs/cern.ch/work/n/nalipour/testBeamAnalysis/Simulation/ALLPix_results/run%s.tar.gz /afs/cern.ch/work/n/nalipour/testBeamAnalysis/Simulation/ALLPix_results/run%s.slcio %s \n"%(run, run, dutID))

    outfile.close()
    os.system("chmod a+rwx %s"%batchScriptName)
    return batchScriptName



# run="000003"
# dutID="2044"
run="000007"
dutID="2048"
queue="1nd"
logname="/afs/cern.ch/work/n/nalipour/testBeamAnalysis/Simulation/lxBatchScripts/converter_run%s.log"%(run)
batchScriptName=GenerateScript(run, dutID)
os.system("bsub -C0 -q %s -o %s -e %s -L /bin/bash %s \n"%(queue, logname, logname, batchScriptName))
