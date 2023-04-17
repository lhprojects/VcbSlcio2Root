#!/bin/bash

#source /besfs/groups/higgs/Software/v01-17-05/root/5.34.07/bin/thisroot.sh

#export INSDIR="/afs/ihep.ac.cn/soft/common/gcc/whizard227/install"
#export GCCPATH=${INSDIR}
#export PATH_TO_GFORTRAN=$GCCPATH/bin
#alias gcc=$GCCPATH/bin/gcc
#alias g++=$GCCPATH/bin/g++
#export LD_GFORTRAN=${GCCDIR}/lib

#export LCIO=/afs/ihep.ac.cn/soft/common/gcc/v01-17-05/lcio/v02-04-03/
#export DRUIDDIR=/afs/ihep.ac.cn/users/m/manqi/Software/Druid_SL6
#export PATH=$DRUIDDIR/bin:/afs/ihep.ac.cn/users/m/manqi/ArborV3Deve/SL6/Marlin/v01-05/bin:/afs/ihep.ac.cn/users/m/manqi/Software/ilcsoft/v01-16/CMake/2.8.5/bin:/usr/local/bin:/usr/bin:/bin:/usr/local/sbin:/usr/sbin:/sbin:$ROOTSYS/bin:/usr/X11R6/bin:$LCIO/bin:/opt/glite/bin/:/opt/lcg/bin/:/afs/ihep.ac.cn/soft/common/gcc/v01-17-05/CMake/2.8.5/bin
#export LD_LIBRARY_PATH=$ROOTSYS/lib:/home/manqi/Softwares/ilcsoft/v01-08-01/CLHEP/2.0.4.2/lib:$LCIO/lib:/usr/lib:/lib:${GCCPATH}/lib64:${GCCPATH}/lib:


source /cvmfs/cepc.ihep.ac.cn/software/cepcenv/setup.sh
cepcenv -r /cvmfs/cepc.ihep.ac.cn/software/cepcsoft use 0.1.0-rc9

unset MARLIN_DLL

export MARLIN_DLL=/afs/ihep.ac.cn/users/z/zhuyf/cefs/workspace/testnew/ryuta/hig2zz/lib/libHiggs2zz.so.0.1.0


