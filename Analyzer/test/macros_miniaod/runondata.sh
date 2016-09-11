#!/bin/bash


for file in `echo /nfs/dust2/cms/group/DAS2016/ZprimeDiLepton/Data2016_ZprimeMuMu_13TeV_merged_HLT/CMSSW_8_0_13_ZPRIMEMuMu_13TeV-DataD-V2-JSON-tree.root`; do
 echo $file
 basefile=`basename ${file}`;

 cat runondata.C | sed "s?filename?${file}?g" > tmp.C
 g++ -I $ROOTSYS/include tmp.C ZprimeMuMuPatMiniAodNewData.C `root-config --glibs` `root-config --libs` `root-config --cflags` -o runondata
 ./runondata

done

