#!/bin/bash

g++ -I $ROOTSYS/include RunZprimeMuMuAnalysis.C ZprimeMuMuPat.C `root-config --glibs` `root-config --libs` `root-config --cflags` -o RunZprimeMuMuAnalysis
