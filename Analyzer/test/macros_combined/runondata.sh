#!/bin/bash

# mkdir -p histo/ExpressPhysics_Run2015B/

# for file in `ls  Data2012/ExpressPhysics_Run2012C-PromptReco-v2/res/*.root`; do
for file in `cat a.txt`; do
 echo $file
 basefile=`basename ${file}`;

 cat runondata.C | sed "s?filename?${file}?g" > tmp.C
 g++ -I $ROOTSYS/include tmp.C ZprimeMuMuPat.C `root-config --glibs` `root-config --libs` `root-config --cflags` -o runondata
 ./runondata

# mv CMSSW745-Analyse_ZprimeToMuMu_13TeV.root histo/ExpressPhysics_Run2015B/histo_${basefile}

done

