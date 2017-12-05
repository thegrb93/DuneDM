import logging
import logging.config
import os
import re
import shutil
import subprocess
import sys
import time

def step_range(smin, smax, steps):
    if steps==1:
        return [smin]
    else:
        size = (smax-smin)/(steps-1)
        return [x * size + smin for x in range(steps)]

param_ranges = [
step_range(7.5,7.5,1), #vpmass
step_range(1,5,21), #chimass
step_range(0.01,0.01,1), #kappa
step_range(0.1,0.1,1) #alpha
]

jobcount = 700
nevents = 100000

iteration_state = [0, 0, 0, 0]
scriptdir = os.path.abspath(sys.argv[0]+"/..")+"/"
user = os.getenv("USER")
output_dir = "/pnfs/dune/scratch/users/"+user+"/darkmatter/"

if not os.path.exists(output_dir):
    os.makedirs(output_dir)
if not os.path.exists(scriptdir+"jobs"):
    os.makedirs(scriptdir+"jobs")
if os.path.exists(output_dir+"log.log"):
    os.remove(output_dir+"log.log")

failed = []

def submit_job(vpmass, chimass, kappa, alpha):
    output_filen = "DM_%.3f_%.3f_%.3f_%.3f" % (vpmass, chimass, kappa, alpha)
    script_filename = scriptdir+"jobs/"+output_filen+".sh"

    if os.path.exists(output_dir+output_filen):
        shutil.rmtree(output_dir+output_filen)
    os.makedirs(output_dir+output_filen)

    if not os.path.exists("jobs/"+output_filen+"_gridpack.tar.gz"):
        script = """generate_events %s
set ebeam1 0.93828
set ebeam2 80.93828
set gridpack 1
set vpmass %.3f
set chimass %.3f
set mguser 1 %.3f
set mguser 2 %.3f""" % (output_filen, vpmass, chimass, kappa, alpha)

        file_handle = open(scriptdir+"proc_card.dat", 'w')
        file_handle.write(script)
        file_handle.close()
    
        subprocess.call([sys.executable, '-O', "MG5_aMC_v2_6_0/darkphoton/bin/madevent", scriptdir+"proc_card.dat"])
        if os.path.exists("MG5_aMC_v2_6_0/darkphoton/"+output_filen+"_gridpack.tar.gz"):
            os.rename("MG5_aMC_v2_6_0/darkphoton/"+output_filen+"_gridpack.tar.gz", "jobs/"+output_filen+"_gridpack.tar.gz")
        else:
            failed.append(output_filen)

        output = open(script_filename, "w")
        output.write("""export jobname="""+output_filen+"""
source /grid/fermiapp/products/common/etc/setups.sh
source /grid/fermiapp/products/larsoft/setup
setup ifdhc
setup python v2_7_6
setup root v6_10_04d -q e14:nu:prof

cd $_CONDOR_SCRATCH_DIR
ifdh cp """+scriptdir+"""jobs/"$jobname"_gridpack.tar.gz ./"$jobname"_gridpack.tar.gz
ifdh cp """+scriptdir+"""ExRootAnalysis.tar.gz ./ExRootAnalysis.tar.gz
tar xzf "$jobname"_gridpack.tar.gz
tar xzf ExRootAnalysis.tar.gz
./run.sh """+str(nevents)+""" $PROCESS

ls .
gzip -d events.lhe.gz
./ExRootAnalysis/ExRootLHEFConverter events.lhe unweighted_events.root

if [ -f unweighted_events.root ]; then
    echo Copying to """+output_dir+output_filen+"""/"$PROCESS".root
    ifdh cp unweighted_events.root """+output_dir+output_filen+"""/"$PROCESS".root
fi

echo "ALL DONE!"
""")
        output.close()
        os.chmod(script_filename,0777)
    os.system("jobsub_submit.py -N "+str(jobcount)+" --resource-provides=usage_model=OPPORTUNISTIC --expected-lifetime 6h --OS=SL6 --group=dune -L "+output_dir+"log.log file:///"+script_filename)

for fail in failed:
    print("Job failed to generate: "+fail)

i = 0
index = 0
while i<4:
    vpmass = param_ranges[0][iteration_state[0]]
    chimass = param_ranges[1][iteration_state[1]]
    kappa = param_ranges[2][iteration_state[2]]
    alpha = param_ranges[3][iteration_state[3]]

    if vpmass>=chimass*2:
        submit_job(vpmass, chimass, kappa, alpha)
    
    index += 1
    i = 0
    while i<4:
        iteration_state[i] += 1
        if iteration_state[i]>=len(param_ranges[i]):
            iteration_state[i] = 0
            i += 1
        else:
            break

