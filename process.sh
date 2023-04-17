#!/bin/bash
path=`pwd`

cd $path
Dauof92=("999999999" "888888888" "777777777")
pars=("qqh_bb" "qqh_cc" "qqh_gg" "ww_h0cuxx" "ww_h0uusd" "ww_h0ccds" "zz_h0dtdt" "zz_h0utut" "zz_h0uu_notd" "zz_h0cc_nots" "zzorww_h0cscs" "zzorww_h0udud")
#pars=("qqh_wwone" "qqh_wwtwo" "qqh_zzone" "qqh_zztwo" "qqh_ss" "qqh_bb" "qqh_cc" "qqh_gg" "ww_h0ccbs" "ww_h0cuxx" "ww_h0uusd" "ww_h0ccds" "zz_h0dtdt" "zz_h0utut" "zz_h0uu_notd" "zz_h0cc_nots")
#pars=("uusd" "ccds" "cuxx" "cc_nots" "dtdt" "utut" "uu_notd")
Rparas=("1" "2" "2.5")
Pparas=("-1" "0" "1" "2" "3")
#Gparas=("-1" "-0.5" "0" "0.5" "1" "2")
divPart=("2" "3" "4" "5" "6" "7" "8" "10" "12")

i92Dau=0
while [ "$i92Dau" -lt "1" ]
do
i92dau=${Dauof92[$i92Dau]}



ipar=3
while [ "$ipar" -lt "4" ]
do
par=${pars[$ipar]}

#jobPath=/afs/ihep.ac.cn/users/z/zhuyf/cefs/workspace/yudan/reconstruction/jetclustering/ww/${par}
#jobPath=/afs/ihep.ac.cn/users/z/zhuyf/cefs/workspace/yudan/reconstruction/inclusive/zz/${par}
#jobPath=/afs/ihep.ac.cn/users/z/zhuyf/cefs/workspace/yudan/simulation/softLink/${par}
jobPath=/afs/ihep.ac.cn/users/z/zhuyf/cefs/workspace/yudan/reconstruction/ZZ4Quarks/${par}
#jobPath=/afs/ihep.ac.cn/users/z/zhuyf/cefs/workspace/yudan/reconstruction/350GeV/${par}
#jobPath=/afs/ihep.ac.cn/users/z/zhuyf/cefs/workspace/yudan/reconstruction/normalization/${par}
#jobPath=/afs/ihep.ac.cn/users/z/zhuyf/cefs/workspace/yudan/reconstruction/ztoqq
#jobPath=/afs/ihep.ac.cn/users/z/zhuyf/cefs/workspace/yudan/reconstruction/E240.P${par}.e0.p0.whizard195
#jobPath=/cefs/data/FullSim/CEPC240/CEPC_v4/4fermions/E240.Pww_h.e0.p0.whizard195/
#numSlcio=`find ${jobPath} -name "*.slcio" |wc -l`
NumSlcio=`ls ${jobPath}/*.slcio |wc -l`
numSlcio=$((10#${NumSlcio}+1))

iDivPart=2
while [ "$iDivPart" -lt "3" ]
do
DivPart=${divPart[$iDivPart]}

iRpara=1
while [ "$iRpara" -lt "2" ]
do
Rpara=${Rparas[$iRpara]}

iPpara=2
while [ "$iPpara" -lt "3" ]
do
Ppara=${Pparas[${iPpara}]}

#iGpara=2
#while [ "$iGpara" -lt "3" ]
#do
#Gpara=${Gparas[${iGpara}]}

j=1
#while [ "$j" -lt 11 ]
while [ "$j" -lt "${numSlcio}" ]
do

if [ "$j" -lt 10 ]
then k="0000"$j
elif [ "$j" -ge 10 -a "$j" -lt 100 ]
then k="000"$j
elif [ "$j" -ge 100 -a "$j" -lt 1000 ]
then k="00"$j
fi

export RecoWorkDir=${jobPath}
#OUTPUTDATA=RecoJet/Out${Rpara}${Ppara}
OUTPUTDATA=test92Jet
#OUTPUTDATA=div${DivPart}/RecoJettwo
#echo $RecoWorkDir/$OUTPUTDATA/
#jobPathtwo=${jobPath}/VLC/GenJet/Out${Rpara}:${Ppara}:${Gpara}
#collectionUsed=MCPSIMUFSP
#collectionUsed=ArborPFOs
collectionCreated=GenJetsAll92Daus
#collectionCreated=MCPFastJets
mkdir -p $RecoWorkDir/$OUTPUTDATA/
#outlog=$RecoWorkDir/$OUTPUTDATA/recolog_${j}GeV.log
cp -fr /afs/ihep.ac.cn/users/z/zhuyf/cefs/workspace/testnew/ryuta/hig2zz/steer/knight.xml  $RecoWorkDir/$OUTPUTDATA/${k}.xml
#inputlcio=${jobPath}/div${DivPart}/RecoJet/Reco_${k}.slcio
#inputlcio=${jobPath}/GenJet/Out${Rpara}${Ppara}/Reco_${k}.slcio
#inputlcio=${jobPath}/Simu_ww_h0uusd_${k}.slcio
#inputlcio=${jobPath}/Reco_${par}_${k}.slcio
#inputlcio=${jobPath}/div5/GenJet/Reco_${k}.slcio
inputlcio=${jobPath}/${k}.slcio
echo $inputlcio
#inputlcio=${jobPath}/${par}.e0.p0.${k}_001000_rec.slcio
#inputlcio=${jobPath}/GenJet/Out${Rpara}${Ppara}/RECO_${i}.slcio

#M=$(${j}+50)
let "M=j+50"
#echo $M
if [ "$M" -lt 100 ]
then M="000"$M
elif [ "$M" -ge 100 -a "$M" -lt 1000 ]
then M="00"$M
fi

outputlcio=$RecoWorkDir/$OUTPUTDATA/${k}.slcio
#GEARFILE=/afs/ihep.ac.cn/users/z/zhuyf/cefs/workspace/yudan/simulation/nnhwwzz/ww/GearOutput.xml
GEARFILE=/afs/ihep.ac.cn/users/z/zhuyf/cefs/workspace/LCFIPLUS/Zpole_bb/AfterSimu/GearOutput.xml
sed -i "s#LCIOINPUT#${inputlcio}#g" $RecoWorkDir/$OUTPUTDATA/${k}.xml



#sed -i "s#OUTROOT#${outputroot}#g" $RecoWorkDir/$OUTPUTDATA/${i}.xml
sed -i "s#LCIOOUTPUT#${outputlcio}#g" $RecoWorkDir/$OUTPUTDATA/${k}.xml
sed -i "s#GEARFILE#$GEARFILE#g" $RecoWorkDir/$OUTPUTDATA/${k}.xml
#sed -i "s/ColectionUsed/${collectionUsed}/g" $RecoWorkDir/$OUTPUTDATA/${k}.xml
sed -i "s/ColectionCreated/${collectionCreated}/g" $RecoWorkDir/$OUTPUTDATA/${k}.xml
sed -i "s/NUMJets/${DivPart}/g" $RecoWorkDir/$OUTPUTDATA/${k}.xml
sed -i "s/Rparameter/${Rpara}/g" $RecoWorkDir/$OUTPUTDATA/${k}.xml
sed -i "s/Pparameter/${Ppara}/g" $RecoWorkDir/$OUTPUTDATA/${k}.xml
#sed -i "s/Gparameter/${Gpara}/g" $RecoWorkDir/$OUTPUTDATA/${k}.xml
sed -i "s/REPLACE/${i92dau}/g" $RecoWorkDir/$OUTPUTDATA/${k}.xml

echo "$RecoWorkDir/$OUTPUTDATA/${k}.xml"

# the following code used to write a sh file to be hep_sub
echo \
"#! /bin/bash
#source /afs/ihep.ac.cn/users/z/zhuyf/cefs/workspace/testnew/ryuta/hig2zz/setup.sh
unset MARLIN_DLL
source /afs/ihep.ac.cn/users/z/zhuyf/cefs/workspace/testnew/ryuta/hig2zz/env.sh
#echo $LD_LIBRARY_PATH
#ldd /afs/ihep.ac.cn/users/z/zhuyf/cefs/workspace/testnew/ryuta/hig2zz/lib/libHiggs2zz.so
Marlin $RecoWorkDir/$OUTPUTDATA/${k}.xml
" > $RecoWorkDir/$OUTPUTDATA/submit_${k}.sh
echo "$RecoWorkDir/$OUTPUTDATA/submit_${k}.sh"

#export PATH=/afs/ihep.ac.cn/soft/common/sysgroup/hep_job/bin:$PATH
chmod +x $RecoWorkDir/$OUTPUTDATA/submit_${k}.sh
#hep_sub $RecoWorkDir/$OUTPUTDATA/submit_${k}.sh
sh $RecoWorkDir/$OUTPUTDATA/submit_${k}.sh
#sh $SimuWorkDir/$OUTPUTDATA/job_E${nlay}L_E${ecalsize}mm.sh
let "j+=1"
done

#let "iGpara+=1"
#done

let "iPpara+=1"
done

let "iRpara+=1"
done

let "iDivPart+=1"
done

let "ipar+=1"
done

let "i92Dau+=1"
done
