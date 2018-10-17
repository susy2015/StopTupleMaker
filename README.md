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
cmsrel CMSSW_9_4_6_patch1
cd CMSSW_9_4_6_patch1/src
cmsenv
git cms-init
git clone -b jetToolbox_94X git@github.com:cms-jet/JetToolbox.git
or
git clone git@github.com:cms-jet/JetToolbox.git JMEAnalysis/JetToolbox -b jetToolbox_94X_v1
git cms-merge-topic -u pastika:AddAxis1_946p1

git clone -b master git@github.com:susy2015/StopTupleMaker.git

scram b -j9


you also need Fall17_17Nov2017_V8_MC.db
Cert_314472-317080_13TeV_PromptReco_Collisions18_JSON.txt 
Cert_314472-318876_13TeV_PromptReco_Collisions18_JSON.txt  
Fall17_17Nov2017BCDEF_V6_DATA.db                          
Fall17_17Nov2017_V8_MC.db
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

