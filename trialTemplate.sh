#!/bin/bash

TRIAL=$1

export SCRAM_ARCH=slc6_amd64_gcc530
export VO_CMS_SW_DIR=/cms/cvmfs/cms.cern.ch
source ${VO_CMS_SW_DIR}/cmsset_default.sh
export CMS_PATH=${VO_CMS_SW_DIR}
cd /nfs/fanae/user/vischia/workarea/cmssw/combine/CMSSW_7_4_7/src/
    #source /afs/cern.ch/cms/LCG/LCG-2/UI/cms_ui_env.csh 
eval `scramv1 r -sh`
cmsenv
    #source /afs/cern.ch/cms/ccs/wm/scripts/Crab/crab.csh 
source /cms/cvmfs/cms.cern.ch/crab3/crab.sh
cmsenv

cd /nfs/fanae/user/vischia/workarea/cmssw/tthMultilepton/Clusterization/
root -l -b runTrials.C\(${TRIAL}\)
