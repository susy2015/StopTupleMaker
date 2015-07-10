CommonTools

  This repository holds tools which are common to all sub-analysis groups.  


## Code Recipe for Fat Top Tagger

```
cd CMSSW_7_3_0/src/
cmsenv
git cms-merge-topic gkasieczka:htt-v2 # V2 HEP
git clone git@github.com:gkasieczka/ShowerDeconstruction.git ## Private code
git clone -b TestMiniAOD git@github.com:susy2015/recipeAUX.git
git clone -b align_leptDef_better_search_bins git@github.com:susy2015/SusyAnaTools.git
```

### PUPPI
```
git cms-addpkg CommonTools/PileupAlgos
git remote add nhan-remote https://github.com/nhanvtran/cmssw.git
git fetch nhan-remote puppi-bugfix-for-miniaod
git cherry-pick 0585bf21ae098f14f144b9a3d361178e6cc830e6
scram b
```
https://twiki.cern.ch/twiki/bin/viewauth/CMS/PUPPI


### CMS Top Tagger
https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD2014#Jets

### HEP V2 Top Tagger
https://twiki.cern.ch/twiki/bin/viewauth/CMS/JME14004Organization#HEPTopTagge

### Soft Drop
git clone https://github.com/cms-jet/JetToolbox JMEAnalysis/JetToolbox -b jetToolbox_73X

### Shower Decontruction
Code remain private
