#!/bin/bash


for file in `ls  /eos/uscms/store/user/cmsdas/2016/LONG_EXERCISES/ZprimeDiLeptons/Data2015_ZprimeMuMu_13TeV_merged/SingleMuon_Run2015B-16Oct2015-v1_Nov13_50ns.root`; do
 echo $file
 basefile=`basename ${file}`;

 cat runondata.C | sed "s?filename?${file}?g" > tmp.C
 g++ -I $ROOTSYS/include tmp.C ZprimeMuMuPat.C `root-config --glibs` `root-config --libs` `root-config --cflags` -o runondata
 ./runondata

done

