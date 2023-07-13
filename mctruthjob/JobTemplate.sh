#! /bin/bash
source /cvmfs/cepc.ihep.ac.cn/software/cepcenv/setup.sh
cepcenv -r /cvmfs/cepc.ihep.ac.cn/software/cepcsoft use 0.1.0-rc9

export MARLIN_DLL=/cefs/higgs/liangh/VTX/VcbSlcio2Root/lib/libVcbCreateRoot.so:$MARLIN_DLL

Marlin PLACEHOLDER_XML

