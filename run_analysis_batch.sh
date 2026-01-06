#!/bin/bash
echo "Run script starting"
echo "Event range: $1 to $2"

arch=slc7_amd64_gcc10
rel=CMSSW_10_6_36

echo -e "------------------- START --------------------"
printf "Start time: "; TZ=CET /bin/date
printf "Job is running on node: "; /bin/hostname
printf "Job running as user: "; /usr/bin/id
printf "Job is running in directory: "; /bin/pwd -P

echo
echo -e "---------------- Environments ----------------"

echo -e "\n[0] source /cvmfs/cms.cern.ch/cmsset_default.sh"
source /cvmfs/cms.cern.ch/cmsset_default.sh

baseDir=`/bin/pwd -P`

echo -e "\n[1] export SCRAM_ARCH= $arch"
export SCRAM_ARCH=$arch

echo -e "\n[2] scramv1 project CMSSW $rel"
scramv1 project CMSSW $rel

# cd into CMSSW and do cmsenv
echo -e "\n[4] cd $rel/src/"
cd $rel/src/

echo -e "\n[5] cmsenv"
eval `scramv1 runtime -sh`

# go back to the base directory
cd ../../  

########Timber has been set up and is now running############

echo -e "\n------------------ Analyzer ------------------"

echo -e "current working directory is `pwd`"
workDir=`/bin/pwd -P`
printf "workDir ls -altr: " ; /bin/ls -altr

# copy the input file to be a local file with name input.root
rtfile=$3
xrdcp root://cms-xrd-global.cern.ch/$rtfile input.root

# run script with event range
echo "\n[1] root -l -q 'analyze_pmssm.C($1,$2)'"
root -b -l -q "analyze_pmssm.C($1,$2)"

echo "base directory has:"
/bin/ls -altr .

# clean up
rm input.root
rm triggerAndPresAndSelectionTight.root

echo -e "\n"
echo -e "-------------------- END ---------------------\n"
echo   "UnixTime-JobEnd: "$(date +%s)
