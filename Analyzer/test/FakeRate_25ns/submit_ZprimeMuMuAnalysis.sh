#!/bin/bash

echo "Processing on " `hostname` "at " `date` 

mkdir -p /home/tmp/defilip/$$
mkdir -p /lustre/cms/store/user/defilip/ZprimeAnalysis76x/jobdir
mkdir -p /lustre/cms/store/user/defilip/ZprimeAnalysis76x/histodir

workdir=${PWD}
echo "Running ZprimeMuMu Fake Rate Analysis with executables RunZprimeMuMuAnalysis"
source /cmshome/nicola/slc6/logincms_cvmfs_slc6.sh
export SCRAM_ARCH=slc5_amd64_gcc481
exedir=`echo /cmshome/nicola/slc6/zprime/Analysis13TeV/CMSSW_7_6_3/src/ZprimeDiLeptons/Analyzer/test/FakeRate_25ns`
cd ${exedir}
eval `scramv1 runtime -sh`

if [ -d "/home/tmp/defilip/$$" ]; then
    workdir=`echo /home/tmp/defilip/$$`;
    cd ${workdir};
fi

savedir=`echo /lustre/cms/store/user/defilip/ZprimeAnalysis76x/histodir`

echo "Working dir is $workdir"
echo "Executable dir is $exedir"
echo "Saving dir is $savedir"

${exedir}/RunZprimeMuMuAnalysis_FR which ${exedir}/sig_input.txt 1 ${exedir}/bkg_input.txt 1 ${exedir}/data_input.txt 1 Bari year mc >& ${workdir}/RunZprimeMuMuAnalysis_FR.log 
cp -f ${workdir}/RunZprimeMuMuAnalysis_FR.log /lustre/cms/store/user/defilip/ZprimeAnalysis/jobdir/.
mv -f ${workdir}/CMSSW763-Analyse_ZprimeToMuMu_13TeV_FR.root     ${savedir}/output.root
mv -f ${workdir}/CMSSW763-Analyse_ZprimeToMuMu_13TeV_FR_cand.txt ${savedir}/output_cand.txt

# cleaning the worker node
if [ -d "/home/tmp/defilip/$$" ]; then
    rm -f -R *
    rm -f *
fi
