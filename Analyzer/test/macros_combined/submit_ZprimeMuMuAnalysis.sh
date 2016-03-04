#!/bin/bash

echo "Processing on " `hostname` "at " `date` 

mkdir -p /home/tmp/defilip/$$
mkdir -p /lustre/cms/store/user/defilip/ZprimeAnalysis76x/jobdir
mkdir -p /lustre/cms/store/user/defilip/ZprimeAnalysis76x/histodir

workdir=${PWD}
echo "Running ZprimeMuMu Analysis with executables RunZprimeMuMuAnalysis"
source /cmshome/nicola/slc6/logincms_cvmfs_slc6.sh
export SCRAM_ARCH=slc6_amd64_gcc491
exedir=`echo /cmshome/nicola/slc6/zprime/Analysis13TeV/CMSSW_7_6_3/src/ZprimeDiLeptons/Analyzer/test/macros_combined`
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

${exedir}/RunZprimeMuMuAnalysis which ${exedir}/sig_input.txt 1 ${exedir}/bkg_input.txt 1 ${exedir}/data_input.txt 1 Bari year mc >& ${workdir}/RunZprimeMuMuAnalysis.log 
cp -f ${workdir}/RunZprimeMuMuAnalysis.log /lustre/cms/store/user/defilip/ZprimeAnalysis76x/jobdir/.
cp -f ${workdir}/CMSSW763-Analyse_ZprimeToMuMu_13TeV.root    ${savedir}/output.root
cp -f ${workdir}/CMSSW763-Analyse_ZprimeToMuMu_13TeV_cand.txt ${savedir}/output_cand.txt
# cleaning the worker node
if [ -d "/home/tmp/defilip/$$" ]; then
    rm -f -R *
    rm -f *
fi
