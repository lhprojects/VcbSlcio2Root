#!/bin/bash
path=`pwd`

cd $path
#pars=("sw_l0mu" "sw_l0tau")
#pars=("sw_sl0qq")
#pars=("sze_sl0dd" "sze_sl0uu")
#pars=("sznu_sl0nu_down" "sznu_sl0nu_up")
#pars=("ww_h0ccbs" "ww_h0ccds" "ww_h0cuxx" "ww_h0uubd" "ww_h0uusd")
#pars=("ww_sl0muq" "ww_sl0tauq")
#pars=("zz_h0cc_nots" "zz_h0dtdt" "zz_h0utut" "zz_h0uu_notd")
#pars=("zz_sl0mu_down" "zz_sl0mu_up" "zz_sl0tau_down" "zz_sl0tau_up" "zz_sl0nu_down" "zz_sl0nu_up")
#pars=("zzorww_h0cscs" "zzorww_h0udud")
#pars=("qqh_X" "e1e1h_X" "e2e2h_X" "e3e3h_X" "nnh_X")
#pars=("bhabha" "e2e2" "e3e3" "qq")
pars=("bb" "cc" "uu" "dd" "ss")

ipar=0
while [ "$ipar" -lt "5" ]
do
par=${pars[$ipar]}

jobPath=/cefs/higgs/liangh/liangh

#inputDir=/cefs/data/DstData/CEPC240/CEPC_v4/4fermions/E240.Pww_sl.e0.p0.whizard195
#inputDir=/cefs/data/DstData/CEPC240/CEPC_v4/4fermions/E240.Psw_l.e0.p0.whizard195
#inputDir=/cefs/data/DstData/CEPC240/CEPC_v4/4fermions/E240.Psw_sl.e0.p0.whizard195
#inputDir=/cefs/data/DstData/CEPC240/CEPC_v4/4fermions/E240.Psze_sl.e0.p0.whizard195
#inputDir=/cefs/data/DstData/CEPC240/CEPC_v4/4fermions/E240.Psznu_sl.e0.p0.whizard195
#inputDir=/cefs/data/DstData/CEPC240/CEPC_v4/4fermions/E240.Pww_h.e0.p0.whizard195
#inputDir=/cefs/data/DstData/CEPC240/CEPC_v4/4fermions/E240.Pzz_h.e0.p0.whizard195
#inputDir=/cefs/data/DstData/CEPC240/CEPC_v4/4fermions/E240.Pzz_sl.e0.p0.whizard195
#inputDir=/cefs/data/DstData/CEPC240/CEPC_v4/4fermions/E240.Pzzorww_h.e0.p0.whizard195
#inputDir=/cefs/data/DstData/CEPC240/CEPC_v4/2fermions/E240.P${par}.e0.p0.whizard195
#inputDir=/cefs/data/DstData/CEPC240/CEPC_v4/higgs/E240.P${par}.e0.p0.whizard195
#inputDir=/cefs/data/DstData/CEPC240/CEPC_v4/higgs/E240.Pe1e1h_X.e0.p0.whizard195
#inputDir=/cefs/data/DstData/CEPC240/CEPC_v4/higgs/E240.Pe2e2h_X.e0.p0.whizard195
#inputDir=/cefs/data/DstData/CEPC240/CEPC_v4/higgs/E240.Pe3e3h_X.e0.p0.whizard195
#inputDir=/cefs/data/DstData/CEPC240/CEPC_v4/higgs/E240.Pnnh_X.e0.p0.whizard195
#inputDir=/cefs/data/DstData/CEPC240/CEPC_v4/2fermions/E240.Pbhabha.e0.p0.whizard195
inputDir=/cefs/data/DstData/Zpole/CEPC_v4/E91.2.P${par}.e0.p0.whizard195


RecoWorkDir=${jobPath}
OUTPUTDATA=${par}
echo $RecoWorkDir/$OUTPUTDATA/

mkdir -p $RecoWorkDir/$OUTPUTDATA/



Numjob=`ls ${inputDir}/${par}.e0.p0.*_000000_dst.slcio |wc -l`
numJob=$((10#${Numjob}))
echo ${numJob}

if [ "${numJob}" -lt "100" ]; then
want="${numJob}"
elif [ "${numJob}" -gt "100" ]; then
want=100
fi


j=6
while [ "$j" -lt "${want}" ]
#while [ "$j" -lt "${numJob}" ]
do



if [ "$j" -lt 10 ]
then k="0000"$j
elif [ "$j" -ge 10 -a "$j" -lt 100 ]
then k="000"$j
elif [ "$j" -ge 100 -a "$j" -lt 1000 ]
then k="00"$j
fi




have=`ls $RecoWorkDir/$OUTPUTDATA/submit_${par}_*.sh |wc -l`
Have=$((10#${have}))


outjob=${Have}

ojob=`printf "%05d" $outjob`

while [ -f "$RecoWorkDir/$OUTPUTDATA/submit_${par}_${ojob}.sh" ]
do
let "outjob+=1"
ojob=`printf "%05d" $outjob`
done


if [ ! -f "$RecoWorkDir/$OUTPUTDATA/submit_${par}_${ojob}.sh" ]; then
outputroot=$RecoWorkDir/$OUTPUTDATA/${ojob}.root
fi



cp -fr /cefs/higgs/liangh/liangh/hig2zz/steer/knight.xml  $RecoWorkDir/$OUTPUTDATA/${par}_${ojob}.xml

#inputlcio=${inputDir}/${par}.e0.p0.${k}_${class}_dst.slcio
#echo $inputlcio


outputlcio=$RecoWorkDir/$OUTPUTDATA/${ojob}.slcio
GEARFILE=/afs/ihep.ac.cn/users/z/zhuyf/cefs/workspace/LCFIPLUS/GearOutput.xml

Numclass=`ls ${inputDir}/${par}.e0.p0.${k}_*_dst.slcio |wc -l`
numClass=$((10#${Numclass}))
echo ${numClass}


times=`expr ${numClass} \* 1000`

iclass=0000
#while [ "$iclass" -lt "${times}" ]
while [ "$iclass" -lt "2000" ]
do

class=`printf "%06d" $iclass`


if [ -f "${inputDir}/${par}.e0.p0.${k}_${class}_dst.slcio" ]; then
line=${inputDir}/${par}.e0.p0.${k}_${class}_dst.slcio
echo $line

sed -i "\#LCIOInputFiles#a${line}" $RecoWorkDir/$OUTPUTDATA/${par}_${ojob}.xml
fi

let "iclass+=1000"
done


#sed -i "s#LCIOINPUT#${inputlcio}#g" $RecoWorkDir/$OUTPUTDATA/${par}_${ojob}.xml










sed -i "s#OUTROOT#${outputroot}#g" $RecoWorkDir/$OUTPUTDATA/${par}_${ojob}.xml
#sed -i "s#LCIOOUTPUT#${outputlcio}#g" $RecoWorkDir/$OUTPUTDATA/${k}.xml
sed -i "s#GEARFILE#$GEARFILE#g" $RecoWorkDir/$OUTPUTDATA/${par}_${ojob}.xml
#sed -i "s#LCIOOUTPUT#${outputlcio}#g" $RecoWorkDir/$OUTPUTDATA/${k}.xml

# the following code used to write a sh file to be hep_sub
echo \
"#! /bin/bash
unset MARLIN_DLL
source /cvmfs/cepc.ihep.ac.cn/software/cepcenv/setup.sh
cepcenv -r /cvmfs/cepc.ihep.ac.cn/software/cepcsoft use 0.1.0-rc9
export MARLIN_DLL=/cefs/higgs/liangh/liangh/hig2zz/lib/libHiggs2zz.so.0.1.0:$MARLIN_DLL
Marlin $RecoWorkDir/$OUTPUTDATA/${par}_${ojob}.xml
" > $RecoWorkDir/$OUTPUTDATA/submit_${par}_${ojob}.sh
echo "$RecoWorkDir/$OUTPUTDATA/submit_${par}_${ojob}.sh"

chmod +x $RecoWorkDir/$OUTPUTDATA/submit_${par}_${ojob}.sh
hep_sub $RecoWorkDir/$OUTPUTDATA/submit_${par}_${ojob}.sh
#sh $RecoWorkDir/$OUTPUTDATA/submit_${par}_${ojob}.sh




let "j+=1"
done



let "ipar+=1"
done

