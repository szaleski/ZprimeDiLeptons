#!/bin/bash

if [ "$1" == "" ]; then
  if [ "$2" == "" ]; then
    if [ "$3" == "" ]; then
      if [ "$3" == "" ]; then
          echo "Please provide arguments to the script: site configuration, data type and MC type"
          echo "Usage bash loopsubmit_data.sh <arg1> <arg2> <arg3> <arg4>"
          exit
      fi
    fi
  fi
fi



echo "$1 configuration";
echo "$2 data"
echo "$3 simulation"
echo "$4 site"

SCERN="CERN";
SFNAL="FNAL";
SDESY="DESY";

mkdir -p jobs;

###### Background
n=0;
m=0;

echo "Reading data_input_$2_AN.txt file"

cp -f data_input_$2_AN.txt data_input.txt
nlines=`wc -l data_input_$2_AN.txt | awk '{print $1}'`;

while [ $n -lt ${nlines} ]; do
  (( n = n + 1 ))
  (( m = ${nlines} - n ))
  echo $n $m
  mkdir -p DataCards$2
  rm -f DataCards$2/data_input_${n}.txt
  cat data_input.txt | head -1 > DataCards$2/data_input_${n}.txt
  samplename=`cat DataCards$2/data_input_${n}.txt | awk '{print $1}'`
  echo $samplename
  cat data_input.txt | tail -n $m >  data_input_tmp.txt
  mv  data_input_tmp.txt data_input.txt
  rm -f jobs/submit_ZprimeMuMuAnalysis_${samplename}.shx
  if [ $4 = ${SCERN} ]; then
      cat submit_ZprimeMuMuAnalysis_CERN.sh | sed "s?which?data?g" | sed "s?site?$1?g" | sed "s?mc?$3?g" |sed "s?year?$2?g" | sed "s?ZprimeMuMuAnalysis?RunZprimeMuMuAnalysis?g" | sed "s?jobdir?jobs/jobsZprimeMuMu?g" | sed "s?histodir?histos/histosZprimeMuMu?g" | sed "s?output?output_${samplename}?g" | sed "s?data_input.txt?DataCards$2/data_input_${n}.txt?g" | sed "s?s.log?s_${samplename}.log?g" > jobs/submit_ZprimeMuMuAnalysis_${samplename}.sh
  elif  [ $4 = ${SFNAL} ]; then 
      cat submit_ZprimeMuMuAnalysis_FNAL.sh  | sed "s?path?$PATH?g"  | sed "s?lib?$LD_LIBRARY_PATH?g" | sed "s?which?data?g" | sed "s?site?$4?g" | sed "s?mc?$3?g" |sed "s?year?$2?g" | sed "s?jobdir?jobs/jobsZprimeMuMu?g" | sed "s?histodir?histos/histosZprimeMuMu?g" | sed "s?output?output_${samplename}?g" | sed "s?data_input.txt?data_input_${n}.txt?g" | sed "s?s.log?s_${samplename}.log?g" > jobs/submit_ZprimeMuMuAnalysis_${samplename}.sh
      cat condor_template.cfg  | sed "s?submit_ZprimeMuMuAnalysis_FNAL?submit_ZprimeMuMuAnalysis_${samplename}?g" | sed "s?sig_input_h150.txt?DataCards$2/data_input_${n}.txt?g" | sed "s?mail?`whoami`?g" > jobs/condor_ZprimeMuMuAnalysis_${samplename}.cfg      
  elif  [ $4 = ${SDESY} ]; then
     cat submit_ZprimeMuMuAnalysis_DESY.sh | sed "s?which?data?g" | sed "s?site?$1?g" | sed "s?mc?$3?g" |sed "s?year?$2?g" | sed "s?ZprimeMuMuAnalysis?RunZprimeMuMuAnalysis?g" | sed "s?jobdir?jobs/jobsZprimeMuMu?g" | sed "s?histodir?histos/histosZprimeMuMu?g" | sed "s?output?output_${samplename}?g" | sed "s?data_input.txt?DataCards$2/data_input_${n}.txt?g" | sed "s?s.log?s_${samplename}.log?g" > jobs/submit_ZprimeMuMuAnalysis_${samplename}.sh
  else
      cat submit_ZprimeMuMuAnalysis.sh | sed "s?which?data?g" | sed "s?mc?$3?g" |sed "s?year?$2?g" | sed "s?jobdir?jobs/jobsZprimeMuMu_combined?g" | sed "s?histodir?histos/histosZprimeMuMu_combined?g" | sed "s?output?output_${samplename}?g" | sed "s?data_input.txt?DataCards$2/data_input_${n}.txt?g" | sed "s?s.log?s_${samplename}.log?g" > jobs/submit_ZprimeMuMuAnalysis_${samplename}.sh

  fi

  chmod u+xr jobs/submit_ZprimeMuMuAnalysis_${samplename}.sh

  cd jobs

  if [ $4 = ${SCERN} ]; then
      echo "Submitting jobs via LSF at CERN"
      bsub -q 8nh  submit_ZprimeMuMuAnalysis_${samplename}.sh
  elif  [ $4 = ${SFNAL} ]; then
      echo "Submitting jobs via CONDOR at FNAL"
      # condor_submit  condor_ZprimeMuMuAnalysis_${samplename}.cfg
  elif  [ $4 = ${SDESY} ]; then
      echo "Submitting jobs via SGE"
      qsub submit_ZprimeMuMuAnalysis_${samplename}.sh   
  else
      echo "Submitting jobs via PBS"    
      qsub -q local submit_ZprimeMuMuAnalysis_${samplename}.sh
  fi
  cd ..
done 

