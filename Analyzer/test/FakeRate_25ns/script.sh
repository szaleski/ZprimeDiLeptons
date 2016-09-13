

for file in `ls /lustre/cms/store/user/defilip/ZprimeAnalysis/Spring16_merged/QCD*`; do
 echo $file
 basefile=`basename $file`
 echo $basefile
 mv /lustre/cms/store/user/defilip/ZprimeAnalysis/Spring16_merged/$basefile /lustre/cms/store/user/defilip/ZprimeAnalysis/Spring16_merged/CMSSW803_MC_$basefile
done
