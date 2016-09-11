#!/bin/bash

echo "Processing on " `hostname` "at " `date` 

mkdir -p $$
mkdir -p /afs/desy.de/user/s/school22/CMSSW_8_0_13/src/ZprimeDiLeptons/Analyzer/test/macros_miniaod/jobdir
mkdir -p /afs/desy.de/user/s/school22/CMSSW_8_0_13/src/ZprimeDiLeptons/Analyzer/test/macros_miniaod/histodir

workdir=${PWD}
echo "Running ZprimeMuMu Analysis with executables RunZprimeMuMuAnalysis"
source /cvmfs/cms.cern.ch/cmsset_default.sh 
export SCRAM_ARCH=slc6_amd64_gcc530
exedir=`echo /afs/desy.de/user/s/school22/CMSSW_8_0_13/src/ZprimeDiLeptons/Analyzer/test/macros_miniaod`
cd ${exedir}
eval `scramv1 runtime -sh`

if [ -d "$$" ]; then
    workdir=`echo $$`;
    cd ${workdir};
fi

savedir=`echo /afs/desy.de/user/s/school22/CMSSW_8_0_13/src/ZprimeDiLeptons/Analyzer/test/macros_miniaod/histodir`

echo "Working dir is $workdir"
echo "Executable dir is $exedir"
echo "Saving dir is $savedir"

${exedir}/RunZprimeMuMuAnalysis which ${exedir}/sig_input.txt 1 ${exedir}/bkg_input.txt 1 ${exedir}/data_input.txt 1 site year mc >& ${workdir}/RunZprimeMuMuAnalysis.log 
cp -f ${workdir}/RunZprimeMuMuAnalysis.log /afs/desy.de/user/s/school22/CMSSW_8_0_13/src/ZprimeDiLeptons/Analyzer/test/macros_miniaod/jobdir/.
cp -f ${workdir}/CMSSW745-Analyse_ZprimeToMuMu_13TeV.root    ${savedir}/output.root
cp -f ${workdir}/CMSSW745-Analyse_ZprimeToMuMu_13TeV_cand.txt ${savedir}/output_cand.txt
# cleaning the worker node
if [ -d "$$" ]; then
    rm -f -R *
    rm -f *
fi
