### Instructions

## Updating a git repo to a new branch 

Commit all outstanding change to a new branch.
```
git fetch
git checkout NEW_TAG_NAME
```

## Production Code

The following installation instructions assume the user wants to process Run2016 data or Spring16 MC.

```
cmsrel CMSSW_9_4_10
cd CMSSW_9_4_10/src
cmsenv
git cms-init
git cms-merge-topic -u pastika:AddAxis1_946p1 #needed for AXIS1
#https://twiki.cern.ch/twiki/bin/view/CMS/MissingETUncertaintyPrescription#Instructions_for_9_4_X_X_9_for_2
git cms-merge-topic cms-met:METFixEE2017_949_v2 #NEED for METFix in 2017

git clone git@github.com:cms-jet/JetToolbox.git JMEAnalysis/JetToolbox -b jetToolbox_94X_v1 #needed for JEtToolBox
#DeepAK8
git clone ssh://git@gitlab.cern.ch:7999/TreeMaker/NNKit.git -b cmssw-improvements-4
scram setup /cvmfs/cms.cern.ch/slc6_amd64_gcc700/cms/cmssw/CMSSW_10_3_0_pre4/config/toolbox/slc6_amd64_gcc700/tools/selected/mxnet-predict.xml
#94x V2 Egamma ID: https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2#Recipe_for_regular_users_formats
git cms-merge-topic guitargeek:EgammaID_9_4_X
#LPC ntuple code
git clone -b master git@github.com:susy2015/StopTupleMaker.git

scram b -j9


you also need Fall17_17Nov2017_V8_MC.db You can get this from the Jet db link below
Cert_314472-317080_13TeV_PromptReco_Collisions18_JSON.txt 
Cert_314472-318876_13TeV_PromptReco_Collisions18_JSON.txt  
Fall17_17Nov2017BCDEF_V6_DATA.db                          
Fall17_17Nov2017_V8_MC.db

Jet db: https://twiki.cern.ch/twiki/bin/view/CMS/JECDataMC
Global Tag: https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideFrontierConditions#Global_Tags_for_2018_data_taking Need for running on crab.
```

To produce ntuples with crab submission (google doc for production signup is https://docs.google.com/spreadsheets/d/17Hpp5S-UkiKvuugKxqbW0-3aLhiJrJP8MEEyHce_Lzw/edit#gid=0):
```
cd SusyAnaTools/SkimsAUX/workdir/prod/80X_crab_example
# Edit the line in MultiCrab3.py
# selSubmitKey = 'TEST ALL'
# to the sample keys you'd like to submit. Note we use space to seperate the sample keys.
# selSubmitKey = 'TTJets HTMHT'
# Will let you submit all the samples matching "*TTJets*" and "*HTMHT*"
# ----
# After submission, one can change the selSubmitKey to the following to monitor and re-submit failed jobs for you automatically
# selSubmitKey = 'TEST STATUS TTJets HTMHT'
```

To test the ntuple production interactively:
```
cd SusyAnaTools/SkimsAUX/workdir/prod/80X_crab_example
# do "cmsRun treeMaker_stopRA2.py" but with commandline options extracted from the MultiCrab3.py.
# For instance, for 2016 MC, do the following from
# https://github.com/susy2015/SusyAnaTools/blob/master/SkimsAUX/workdir/prod/80X_crab_example/MultiCrab3.py#L271:
# cmsRun treeMaker_stopRA2.py mcInfo=1 GlobalTag=80X_mcRun2_asymptotic_2016_miniAODv2 specialFix=JEC jecDBname=Spring16_25nsV1_MC maxEvents=1000
```

