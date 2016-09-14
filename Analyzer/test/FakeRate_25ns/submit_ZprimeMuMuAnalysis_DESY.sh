#!/bin/bash

echo "Processing on " `hostname` "at " `date` 

mkdir -p /nfs/dust2/cms/group/DAS2016/${USER}/jobdir
mkdir -p /nfs/dust2/cms/group/DAS2016/${USER}/histodir

cd /tmp
mkdir $$
workdir=${PWD}/$$

echo "Running ZprimeMuMu Analysis_FR with executables RunZprimeMuMuAnalysis_FR"
source /cvmfs/cms.cern.ch/cmsset_default.sh 
export SCRAM_ARCH=slc6_amd64_gcc530
exedir=`echo /afs/desy.de/user/s/school22/CMSSW_8_0_13/src/ZprimeDiLeptons/Analyzer/test/FakeRate_25ns`
cd ${exedir}
eval `scramv1 runtime -sh`

cd ${workdir};

savedir=`echo /nfs/dust2/cms/group/DAS2016/${USER}/histodir`

echo "Working dir is $workdir"
echo "Executable dir is $exedir"
echo "Saving dir is $savedir"

${exedir}/RunZprimeMuMuAnalysis_FR which ${exedir}/sig_input.txt 1 ${exedir}/bkg_input.txt 1 ${exedir}/data_input.txt 1 site year mc >& ${workdir}/RunZprimeMuMuAnalysis_FR.log 
cp -f ${workdir}/RunZprimeMuMuAnalysis_FR.log /nfs/dust2/cms/group/DAS2016/${USER}/jobdir/output.log
cp -f ${workdir}/CMSSW803-Analyse_ZprimeToMuMu_13TeV.root     ${savedir}/output.root
cp -f ${workdir}/CMSSW803-Analyse_ZprimeToMuMu_13TeV_cand.txt ${savedir}/output_cand.txt
# cleaning the worker node

#if [ -d "${workdir}" ]; then
#    rm -f -R *
#    rm -f *
#fi
