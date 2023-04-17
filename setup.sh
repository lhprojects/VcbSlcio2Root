#!/usr/bin/env bash

# Unset MARLIN_DLL
unset MARLIN_DLL

# Set Variables for MARLIN Execution
shopt -s expand_aliases
source /besfs/groups/higgs/Software/v01-17-05_slc6/init_ilcsoft.sh

# Add MARLIN Library Path 
#export LD_LIBRARY_PATH=$PWD/lib:$LD_LIBRARY_PATH
#export MARLIN_DLL=$PWD/lib/libHiggs2zz.so:$MARLIN_DLL

#export PATH=/afs/ihep.ac.cn/users/z/zhuyf/cefs/workspace/FastJet/fastjet-install/bin:$PATH
#export LD_LIBRARY_PATH=/afs/ihep.ac.cn/users/z/zhuyf/cefs/workspace/FastJet/fastjet-install/lib:$LD_LIBRARY_PATH

# For Condor Job Submit
export PATH=/afs/ihep.ac.cn/soft/common/sysgroup/hep_job/bin:$PATH


# PyROOT 
export PYTHONPATH=$ROOTSYS/lib:$PYTHONPATH
