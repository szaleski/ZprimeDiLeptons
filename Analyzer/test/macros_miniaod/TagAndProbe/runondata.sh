#!/bin/bash


for file in `echo /lustre/cms/store/user/defilip/ZprimeAnalysis/Data2016_ZprimeMuMu_13TeV_merged_HLT/CMSSW_8_0_13_ZPRIMEMuMu_13TeV-DataG-V1-JSON-tree.root`; do
 echo $file
 basefile=`basename ${file}`;

 cat runondata.C | sed "s?filename?${file}?g" > tmp.C
 g++ -I $ROOTSYS/include tmp.C TagProbeMuon.C `root-config --glibs` `root-config --libs` `root-config --cflags` -o runondata
 ./runondata

done

