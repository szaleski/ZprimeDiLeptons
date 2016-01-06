#!/bin/bash

g++ -I $ROOTSYS/include RunZprimeMuMuAnalysis_FR.C ZprimeMuMu_FR.C `root-config --glibs` `root-config --libs` `root-config --cflags` -o RunZprimeMuMuAnalysis_FR
