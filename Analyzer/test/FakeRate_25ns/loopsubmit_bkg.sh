#!/bin/bash

if [ "$1" == "" ]; then
  if [ "$2" == "" ]; then
      if [ "$3" == "" ]; then
          echo "Please provide arguments to the script: site configuration, data type and MC type"
          echo "Usage bash loopsubmit_bkg.sh <arg1> <arg2> <arg3>"
          exit
      fi
  fi
fi



echo "$1 configuration";
echo "$2 data"
echo "$3 simulation"

SCERN="CERN";
SFNAL="FNAL";
SDESY="DESY";

mkdir jobs;

###### Background
n=0;
m=0;

echo "Reading bkg_input_$3_AN.txt file"

cp -f bkg_input_$3_AN.txt bkg_input.txt
nlines=`wc -l bkg_input_$3_AN.txt | awk '{print $1}'`;

while [ $n -lt ${nlines} ]; do
  (( n = n + 1 ))
  (( m = ${nlines} - n ))
  echo $n $m
  mkdir -p BkgCards$3
  rm -f BkgCards$3/bkg_input_${n}.txt
  cat bkg_input.txt | head -1 > BkgCards$3/bkg_input_${n}.txt
  samplename=`cat BkgCards$3/bkg_input_${n}.txt | awk '{print $1}'`
  echo $samplename
  cat bkg_input.txt | tail -n $m >  bkg_input_tmp.txt
  mv  bkg_input_tmp.txt bkg_input.txt
  rm -f jobs/submit_ZprimeMuMuAnalysis_${samplename}.shx
  if [ $1 = ${SCERN} ]; then
      cat submit_ZprimeMuMuAnalysis_CERN.sh | sed "s?which?bkg?g" | sed "s?site?$1?g" | sed "s?mc?$3?g" |sed "s?year?$2?g" | sed "s?ZprimeMuMuAnalysis?RunZprimeMuMuAnalysis?g" | sed "s?jobdir?jobs/jobsZprimeMuMu_FR?g" | sed "s?histodir?histos/histosZprimeMuMu_FR?g" | sed "s?output?output_${samplename}?g" | sed "s?bkg_input.txt?BkgCards$3/bkg_input_${n}.txt?g" | sed "s?s_FR.log?s_FR_${samplename}.log?g" > jobs/submit_ZprimeMuMuAnalysis_${samplename}.sh
  elif  [ $1 = ${SFNAL} ]; then 
      cat submit_ZprimeMuMuAnalysis_FNAL.sh | sed "s?which?bkg?g" | sed "s?site?$1?g" | sed "s?mc?$3?g" |sed "s?year?$2?g" | sed "s?ZprimeMuMuAnalysis?RunZprimeMuMuAnalysis?g" | sed "s?jobdir?jobs/jobsZprimeMuMu?g" | sed "s?histodir?histos/histosZprimeMuMu?g" | sed "s?output?output_${samplename}?g" | sed "s?bkg_input.txt?BkgCards$3/bkg_input_${n}.txt?g" | sed "s?s.log?s_${samplename}.log?g" > jobs/submit_ZprimeMuMuAnalysis_${samplename}.sh
      cat condor_template.cfg  | sed "s?submit_ZprimeMuMuAnalysis_FNAL?submit_ZprimeMuMuAnalysis_${samplename}?g" | sed "s?sig_input_h150.txt?BkgCards$3/bkg_input_${n}.txt?g" | sed "s?mail?`whoami`?g" > jobs/condor_ZprimeMuMuAnalysis_${samplename}.cfg      
  elif  [ $1 = ${SDESY} ]; then
     cat submit_ZprimeMuMuAnalysis_DESY.sh | sed "s?which?bkg?g" | sed "s?site?$1?g" | sed "s?mc?$3?g" |sed "s?year?$2?g" | sed "s?ZprimeMuMuAnalysis?RunZprimeMuMuAnalysis?g" | sed "s?jobdir?jobs/jobsZprimeMuMu?g" | sed "s?histodir?histos/histosZprimeMuMu?g" | sed "s?output?output_${samplename}?g" | sed "s?bkg_input.txt?BkgCards$3/bkg_input_${n}.txt?g" | sed "s?s.log?s_${samplename}.log?g" > jobs/submit_ZprimeMuMuAnalysis_${samplename}.sh
  else
      cat submit_ZprimeMuMuAnalysis.sh | sed "s?which?bkg?g" | sed "s?mc?$3?g" |sed "s?year?$2?g" | sed "s?jobdir?jobs/jobsZprimeMuMu_FR_25ns?g" | sed "s?histodir?histos/histosZprimeMuMu_FR_25ns?g" | sed "s?output?output_${samplename}?g" | sed "s?bkg_input.txt?BkgCards$3/bkg_input_${n}.txt?g" | sed "s?s_FR.log?s_FR_25ns_${samplename}.log?g" > jobs/submit_ZprimeMuMuAnalysis_${samplename}.sh

  fi

  chmod u+xr jobs/submit_ZprimeMuMuAnalysis_${samplename}.sh

  cd jobs

  if [ $1 = ${SCERN} ]; then
      echo "Submitting jobs via LSF at CERN"
      bsub -q 8nh  submit_ZprimeMuMuAnalysis_${samplename}.sh
  elif  [ $1 = ${SFNAL} ]; then
      echo "Submitting jobs via CONDOR at FNAL"
      condor_submit  condor_ZprimeMuMuAnalysis_${samplename}.cfg
  elif  [ $1 = ${SDESY} ]; then
      echo "Submitting jobs via SGE"
      qsub submit_ZprimeMuMuAnalysis_${samplename}.sh   
  else
      echo "Submitting jobs via PBS"    
      qsub -q local submit_ZprimeMuMuAnalysis_${samplename}.sh
  fi
  cd ..
done 

