#!/bin/bash
echo "Run script starting"
# echo "All arguments passed:"
# echo $*
echo "Arg 1:"
echo $1
# echo "Arg 2:"
# echo $2
# echo "Arg 3:"
# echo $3

arch=slc7_amd64_gcc10
#rel=CMSSW_10_6_27
# rel=CMSSW_11_1_4
#rel=CMSSW_11_4_1
# rel=CMSSW_12_3_0
rel=CMSSW_10_6_36
#sandbox=$(ls tarball*.tgz)
# sandbox=$(ls TIMBER*.tgz)
#env=$(ls timber-env*.tgz)
#export X509_USER_PROXY=$1
#voms-proxy-info -all
#voms-proxy-info -all -file $1

echo -e "------------------- START --------------------"
printf "Start time: "; TZ=CET /bin/date
printf "Job is running on node: "; /bin/hostname
printf "Job running as user: "; /usr/bin/id
printf "Job is running in directory: "; /bin/pwd -P
# printf "Files in directory:"; /bin/ls -altr ./

echo
echo -e "---------------- Environments ----------------"

echo -e "\n[0] source /cvmfs/cms.cern.ch/cmsset_default.sh"
source /cvmfs/cms.cern.ch/cmsset_default.sh

baseDir=`/bin/pwd -P`

echo -e "\n[1] export SCRAM_ARCH= $arch"
# export SCRAM_ARCH=$arch
# export SCRAM_ARCH=slc7_amd64_gcc700
export SCRAM_ARCH=$arch

echo -e "\n[2] scramv1 project CMSSW $rel"
scramv1 project CMSSW $rel

# cd into CMSSW and do cmsenv
echo -e "\n[4] cd $rel/src/"
cd $rel/src/

echo -e "\n[5] cmsenv"
eval `scramv1 runtime -sh`

# echo -e "Root version first"
# root-config --version

# go back to the base directory
cd ../../  

# echo -e "Python3 version:"
# python3 --version

# echo -e "Voms info below"
# voms-proxy-info -all

########Timber has been set up and is now running############

########Coppied from pervious code##############

echo -e "\n------------------ Analyzer ------------------"

echo -e "current working directory is `pwd`"
workDir=`/bin/pwd -P`
printf "workDir ls -altr: " ; /bin/ls -altr

# copy the input file to be a local file with name input.root
rtfile=$1
xrdcp root://cms-xrd-global.cern.ch/$rtfile input.root
# ls -lthra
# example: rtfile=/store/mc/RunIISummer20UL18NanoAODv9/QCD_HT50to100_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/2820000/0029D547-3B36-7E4B-AE6D-DACA72AA0D9C.root

# copy nanoAOD file to workDir
# echo -e "copy nanoAOD file to $workDir"
# echo "\n[0] xrdcp root://cms-xrd-global.cern.ch/$rtfile $workDir"
# xrdcp "root://cms-xrd-global.cern.ch/$rtfile" $workDir

# file_basename=`basename $rtfile`
# echo -e "This is file_basename: $file_basename"

# simInfo=$(cut -d '/' -f5 <<< "$rtfile")  # QCD_HT50to100_TuneCP5_13TeV-madgraphMLM-pythia8
# id=$(cut -d '/' -f9 <<< "$rtfile")  # 0029D547-3B36-7E4B-AE6D-DACA72AA0D9C.root
# outfile=GenJet_$simInfo-$id

# echo "\n[1] ls -altr ./"
# /bin/ls -altr ./

# run script
echo "\n[1] root -l -q 'analyze_pmssm.C(30000,40000)'"
root -b -l -q 'analyze_pmssm.C(30000,40000)'

echo "base directory has:"
/bin/ls -altr .

# clean up
# rm $workDir/$file_basename
# rm -rf $baseDir/TIMBER
# rm *.html
# rm *.so
# rm -r build/
rm input.root
rm triggerAndPresAndSelectionTight.root
# rm *.c

# delete output file if too small
# This prevents late, parallel and failing condor jobs overwriting the otherwise good output
# if [ -f ${outfile} ]; then
#     if [ $(ls -l ${outfile} | awk '{ print $5 }' ) -lt 1000 ]; then
#         echo -e "\n[12] rm ${outfile}"
#         rm ${outfile}
#     fi
# fi

echo -e "\n"


#done
echo -e "-------------------- END ---------------------\n"
echo   "UnixTime-JobEnd: "$(date +%s)


# echo "Current working directory path: ---------------------------------------"
# pwd
# echo "Current working directory contents: -----------------------------------"
# # ls
# echo "-----------------------------------------------------------------------"
