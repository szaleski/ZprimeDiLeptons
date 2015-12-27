#!/bin/bash

g++ -I $ROOTSYS/include checkentries.C `root-config --glibs` `root-config --libs` `root-config --cflags` -o checkentries

for file in `ls /lustre/cms/store/user/defilip/ZprimeAnalysis/Spring15_25ns_merged/*ZToMu*.root`; do
# echo $file
# cat checkentries.C | sed "s?filename?${file}?g" > tmp.C
# g++ -I $ROOTSYS/include tmp.C `root-config --glibs` `root-config --libs` `root-config --cflags` -o checkentries
./checkentries $file

done

