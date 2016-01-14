#!/bin/bash

echo "Processing on " `hostname` "at " `date` 

mkdir -p /eos/uscms/store/user/ndefilip/ZprimeAnalysis/jobdir
mkdir -p /eos/uscms/store/user/ndefilip/ZprimeAnalysis/histodir

echo "Running ZprimeMuMu Analysis with executables RunZprimeMuMuAnalysis"
source /cvmfs/cms.cern.ch/cmsset_default.sh 
export LD_LIBRARY_PATH=lib:$LD_LIBRARY_PATH
export PATH=path:$PATH

if [ -d "$_CONDOR_SCRATCH_DIR/" ]; then
    workdir=`echo $_CONDOR_SCRATCH_DIR/`;
    cd ${workdir};
fi


savedir=`echo root://cmseos.fnal.gov///store/user/ndefilip/ZprimeAnalysis/histodir`

echo "Working dir is $workdir"
echo "Saving dir is $savedir"

echo "Compiling the macros"
bash compileZprimeMuMuAnalysis.sh

./RunZprimeMuMuAnalysis which sig_input.txt 1 bkg_input.txt 1 data_input.txt 1 site year mc >& ${workdir}/RunZprimeMuMuAnalysis.log 
xrdcp --force ${workdir}/RunZprimeMuMuAnalysis.log root://cmseos.fnal.gov///store/user/ndefilip/ZprimeAnalysis/jobdir/RunZprimeMuMuAnalysis.log
xrdcp --force ${workdir}/CMSSW745-Analyse_ZprimeToMuMu_13TeV.root     ${savedir}/output.root
xrdcp --force ${workdir}/CMSSW745-Analyse_ZprimeToMuMu_13TeV_cand.txt ${savedir}/output_cand.txt
