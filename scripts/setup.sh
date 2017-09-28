#!/bin/sh

# usage() {
#     echo 'Usage : ./setup.sh -r<ROOT version> [-s<sub version for ROOT6>] -c<compiler>'
#     echo '         -r must be either 5 or 6'
#     echo '         -s must be either 6, 8, or 10 if -r is 6'
#     echo '         -c must be either gcc or clang'
#     kill -INT $$;
# }

# while getopts "r:c:s:h" o; do
#     case "${o}" in
#         r)
#             local rootver=${OPTARG}
#             if [ "${rootver}" != "5" ] && [ "${rootver}" != "6" ]
#             then
#                 echo "Invalid ROOT version specified\n"
#                 usage
#             fi
#             ;;
#         s)
#             local subver=${OPTARG}
#             ;;
#         c)
#             local compiler=${OPTARG}
#             if [ "${compiler}" != "gcc" ] && [ "${compiler}" != "clang" ]
#             then
#                 echo "Invalid compiler option\n"
#                 usage
#             fi
#             ;;
#         h)
#             echo "Printing help"
#             usage
#             ;;
#         *)
#             echo "Unknown option"
#             usage
#             ;;
#     esac
# done

# shift $((OPTIND-1))

echo cmsrel CMSSW_8_0_25
echo cd CMSSW_8_0_25/src
echo cmsenv
echo git cms-init
# brings in HEEP V70 into VID
echo git cms-merge-topic Sam-Harper:HEEPV70VID_8010_ReducedCheckout 
# for other E/gamma IDs in VID if you wish to have them
echo git cms-merge-topic ikrav:egm_id_80X_v3
# only necessary to run HEEP V70 on AOD (it will crash if this is not present looking for puppi candidates
echo git cms-merge-topic Sam-Harper:PackedCandNoPuppi
# we need this for the mva weights which runs in VID regardless if you need it or not
echo mkdir -p ../external/slc6_amd64_gcc530/data/RecoEgamma/ElectronIdentification/
echo git clone git@github.com:cms-data/RecoEgamma-ElectronIdentification ../external/slc6_amd64_gcc530/data/RecoEgamma/ElectronIdentification/data
# needed for HEEP modified electrons with value maps
echo git clone -b HEEPV70 git@github.com:Sam-Harper/HEEP.git
echo git cms-merge-topic Sam-Harper:NewEGModifiers_8010
# the zprime selection code
echo git clone https://github.com/cms-analysis/ZprimeDiLeptons
