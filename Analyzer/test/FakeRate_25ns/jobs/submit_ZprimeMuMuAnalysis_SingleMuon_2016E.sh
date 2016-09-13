#!/bin/bash


mkdir -p /lustre/cms/store/user/defilip/ZprimeAnalysis/80X/jobs/jobsZprimeMuMu_FR
mkdir -p /lustre/cms/store/user/defilip/ZprimeAnalysis/80X/histos/histosZprimeMuMu_FR

echo "Running ZprimeMuMu Analysis_FR with executables RunZprimeMuMuAnalysis_FR"
source /cvmfs/cms.cern.ch/cmsset_default.sh

export LD_LIBRARY_PATH=/lustre/home/nicola/slc6/zprime/Analysis13TeV/CMSSW_8_0_13/biglib/slc6_amd64_gcc530:/lustre/home/nicola/slc6/zprime/Analysis13TeV/CMSSW_8_0_13/lib/slc6_amd64_gcc530:/lustre/home/nicola/slc6/zprime/Analysis13TeV/CMSSW_8_0_13/external/slc6_amd64_gcc530/lib:/cvmfs/cms.cern.ch/slc6_amd64_gcc530/cms/cmssw/CMSSW_8_0_13/biglib/slc6_amd64_gcc530:/cvmfs/cms.cern.ch/slc6_amd64_gcc530/cms/cmssw/CMSSW_8_0_13/lib/slc6_amd64_gcc530:/cvmfs/cms.cern.ch/slc6_amd64_gcc530/cms/cmssw/CMSSW_8_0_13/external/slc6_amd64_gcc530/lib:/cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/llvm/3.7.1-giojec/lib64:/cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/gcc/5.3.0/lib64:/cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/gcc/5.3.0/lib::/lustre/home/nicola/phantom_1_2_8/LHAPDF/lib:$LD_LIBRARY_PATH
export PATH=/cvmfs/cms.cern.ch/share/overrides/bin:/lustre/home/nicola/slc6/zprime/Analysis13TeV/CMSSW_8_0_13/bin/slc6_amd64_gcc530:/lustre/home/nicola/slc6/zprime/Analysis13TeV/CMSSW_8_0_13/external/slc6_amd64_gcc530/bin:/cvmfs/cms.cern.ch/slc6_amd64_gcc530/cms/cmssw/CMSSW_8_0_13/bin/slc6_amd64_gcc530:/cvmfs/cms.cern.ch/slc6_amd64_gcc530/cms/cmssw/CMSSW_8_0_13/external/slc6_amd64_gcc530/bin:/cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/llvm/3.7.1-giojec/bin:/cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/gcc/5.3.0/bin:/cvmfs/cms.cern.ch/common:/cvmfs/cms.cern.ch/bin:/usr//lustre/home/nicola/slc6/zprime/Analysis13TeV/CMSSW_8_0_13/biglib/slc6_amd64_gcc530:/lustre/home/nicola/slc6/zprime/Analysis13TeV/CMSSW_8_0_13/lib/slc6_amd64_gcc530:/lustre/home/nicola/slc6/zprime/Analysis13TeV/CMSSW_8_0_13/external/slc6_amd64_gcc530/lib:/cvmfs/cms.cern.ch/slc6_amd64_gcc530/cms/cmssw/CMSSW_8_0_13/biglib/slc6_amd64_gcc530:/cvmfs/cms.cern.ch/slc6_amd64_gcc530/cms/cmssw/CMSSW_8_0_13/lib/slc6_amd64_gcc530:/cvmfs/cms.cern.ch/slc6_amd64_gcc530/cms/cmssw/CMSSW_8_0_13/external/slc6_amd64_gcc530/lib:/cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/llvm/3.7.1-giojec/lib64:/cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/gcc/5.3.0/lib64:/cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/gcc/5.3.0/lib::/lustre/home/nicola/phantom_1_2_8/LHAPDF/lib64/qt-3.3/bin:/bin:/usr/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/sbin:/usr/lpp/mmfs/bin:/lustre/home/nicola/POWHEG-BOX_2015/FASTJET/bin/:/lustre/home/nicola/phantom_1_2_8/LHAPDF/bin:$PATH


if [ -d "$_CONDOR_SCRATCH_DIR" ]; then
    workdir=`echo $_CONDOR_SCRATCH_DIR`;
    cd ${workdir};
else 
    workdir=`echo $PWD`;
    cd ${workdir};
fi


savedir=`echo /lustre/cms/store/user/defilip/ZprimeAnalysis/80X/histos/histosZprimeMuMu_FR`

echo "Working dir is $workdir"
#echo "Executable dir is $exedir"
echo "Saving dir is $savedir"

echo "Compiling the macros"
bash compileZprimeMuMuAnalysis.sh

./RunZprimeMuMuAnalysis_FR data sig_input.txt 1 bkg_input.txt 1 data_input_4.txt 1 BARI 2016 NO >& ${workdir}/RunZprimeMuMuAnalysis_FR.log

mv -f ${workdir}/RunZprimeMuMuAnalysis_FR.log /lustre/cms/store/user/defilip/ZprimeAnalysis/80X/jobs/jobsZprimeMuMu_FR/output_SingleMuon_2016E.log
mv -f ${workdir}/CMSSW803-Analyse_ZprimeToMuMu_13TeV.root     ${savedir}/output_SingleMuon_2016E.root
mv -f ${workdir}/CMSSW803-Analyse_ZprimeToMuMu_13TeV_cand.txt ${savedir}/output_SingleMuon_2016E_cand.txt

