import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras

import os
import sys
import re
import tarfile


## --------------------------
## -- Command line options --
## --------------------------

import FWCore.ParameterSet.VarParsing as VarParsing
options = VarParsing.VarParsing ('standard')
#92X_dataRun2_Prompt_v7
options.register('era', "Run2_25ns", VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.string, "Run2_25ns or Run2_50ns")
options.register('ntpVersion', "Ntp_BLAH", VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.string, "ntpVersion: to be same as the tag of the release. But can be used to produce 72X ntuple as well!")
options.register('GlobalTag', "94X_mc2017_realistic_v12", VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.string, "74X PromptReco: 74X_dataRun2_Prompt_v0")
options.register('cmsswVersion', '101X', VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.string, "'36X' for example. Used for specific MC fix")
options.register('specialFix', 'JEC', VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.string, "special fixes ==>   JEC : use external JEC; IVF : fix IVF; BADMUON : bad muon filters")
options.register('jecDBname', "Fall17_17Nov2017_V8_MC", VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.string, "Summer15_25nsV6_DATA for data")
options.register('hltName', 'HLT', VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.string, "HLT menu to use for trigger matching")

options.register('mcInfo', True, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.bool, "process MonteCarlo data, default is data")

options.register('doPDFs', False, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.bool, "switch to enable the production of PDF weights for NNPDF3.0, CT10, MMHT2014, n.b. you need to do `scram setup lhapdf` before running (default=False)")

options.register('addJetsForZinv', True, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.bool, "switch to add top projected jets for Zinv")

options.register('externalFilterList', '', VarParsing.VarParsing.multiplicity.list, VarParsing.VarParsing.varType.string, "event list for filters")

options.register('jetCorrections', 'L2Relative', VarParsing.VarParsing.multiplicity.list, VarParsing.VarParsing.varType.string, "Level of jet corrections to use: Note the factors are read from DB via GlobalTag")
options.jetCorrections.append('L3Absolute')

options.register('mcVersion', '', VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.string, "'36X' for example. Used for specific MC fix")
options.register('dataVersion', '', VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.string, "'36X' for example. Used for specific DATA fix")

options.register('jetTypes', 'AK4PF', VarParsing.VarParsing.multiplicity.list, VarParsing.VarParsing.varType.string, "Additional jet types that will be produced (AK4Calo and AK4PF, cross cleaned in PF2PAT, are included anyway)")
options.register('hltSelection', '*', VarParsing.VarParsing.multiplicity.list, VarParsing.VarParsing.varType.string, "hlTriggers (OR) used to filter events. for data: ''HLT_Mu9', 'HLT_IsoMu9', 'HLT_IsoMu13_v*''; for MC, HLT_Mu9")
options.register('addKeep', '', VarParsing.VarParsing.multiplicity.list, VarParsing.VarParsing.varType.string, "Additional keep and drop statements to trim the event content")

options.register('debug', False, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.bool, "switch on/off debug mode")
options.register('verbose', 0, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.int, "verbose of debug")

options.register('test', False, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.bool, "switch on/off debug mode")

options.register('doPtHatWeighting', False, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.bool, "PtHat weighting for QCD flat samples, default is False")

options.register('fileslist', '', VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.string, "name of a file with input source list")

options.register('fastsim', False, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.bool, "fastsim sample or not, default is False")

#options.register('doTopTagger', False, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.bool, "do top tagger or not, default is True")

options.register('usePhiCorrMET', False, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.bool, "use phi corrected MET or not, default is False")

options.register('reducedfilterTags', True, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.bool, "use phi corrected MET or not, default is True")

options.register('smsModel', 'T1tttt', VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.string, "SMS model name")
options.register('smsMotherMass',  -1, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.int, "SMS mother mass")
options.register('smsDaughterMass',  -1, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.int, "SMS daughter mass")
options.register('selSMSpts', False, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.bool, "select model pobools")

options.register('pythia8', True, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.bool, "pythi8 or not, default True")

options.parseArguments()
options._tagOrder =[]

print options

## -------------------------
## -- Check CMSSW version --
## -------------------------

procCMSSWver = os.environ['CMSSW_RELEASE_BASE'].split("/")[-1]
print "procCMSSWver : ", procCMSSWver, "\n"

#if not "CMSSW_8_0" in procCMSSWver:
#   print "You should be using CMSSW 80X!! Please change your release area"
#   sys.exit("ERROR: Not using 80X release")

#if not options.cmsswVersion == "80X":
#   print "You should be using CMSSW 80X as option!! Please change to be consistent with the release area"
#   sys.exit("ERROR: Not using 80X option")


## ------------------------
## -- Define the process --
## ------------------------

process = cms.Process("SUSY")

if options.era == "Run2_25ns":
   process = cms.Process("SUSY", eras.Run2_25ns)
elif options.era == "Run2_50ns":
   process = cms.Process("SUSY", eras.Run2_50ns)

process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load("Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')

## -- MessageLogger --
process.MessageLogger.cerr.FwkReport.reportEvery = 100
if options.debug:
   process.MessageLogger.cerr.FwkReport.reportEvery = 1
process.MessageLogger.suppressWarning = cms.untracked.vstring('ecalLaserCorrFilter','manystripclus53X','toomanystripclus53X')

## -- Options and Output Report --
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.options.allowUnscheduled = cms.untracked.bool(True)

## -- Maximal Number of Events --
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )#cms.untracked.PSet(input = cms.untracked.int32(-1))#options.maxEvents) )

if options.debug and options.verbose ==1:
   process.SimpleMemoryCheck = cms.Service('SimpleMemoryCheck',
                                        ignoreTotal=cms.untracked.int32(0),
                                        oncePerEventMode = cms.untracked.bool(False)
                                        )
   process.Timing = cms.Service("Timing")

## -- Input Source --
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
      #'/store/relval/CMSSW_3_8_0_pre8/RelValTTbar/GEN-SIM-RECO/START38_V6-v1/0004/847D00B0-608E-DF11-A37D-003048678FA0.root'
      #'/store/data/Run2017C/SingleMuon/MINIAOD/17Nov2017-v1/40000/0015635A-0BD9-E711-A76C-02163E0133BB.root'
    #'/store/relval/CMSSW_9_4_5_cand1/RelValTTbar_13/MINIAODSIM/94X_mc2017_realistic_v14_RelVal_rmaod-v1/10000/A8356B71-6E2E-E811-8A63-0CC47A7C3424.root'      
   '/store/mc/RunIIFall17MiniAOD/TTJets_SingleLeptFromT_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/94X_mc2017_realistic_v10-v1/00000/00257B91-1808-E811-BD39-0242AC130002.root'
   #'/store/mc/RunIIFall17MiniAOD/WWTo2L2Nu_NNPDF31_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/94X_mc2017_realistic_v10-v1/00000/2C881587-4CF6-E711-B958-0CC47AA989C0.root'
      #'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAOD/WWTo2L2Nu_NNPDF31_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/94X_mc2017_realistic_v10-v1/40000/04291F1E-A501-E811-8EC9-6CC2173D4980.root'
    )
)

if options.files:
   process.source.fileNames = options.files
elif options.fileslist:
   inputfiles = cms.untracked.vstring()
   if os.path.exists(options.fileslist) == False or os.path.isfile(options.fileslist) == False:
      sys.exit(5)
   else:
      ifile = open(options.fileslist, 'r')
      for line in ifile.readlines():
         inputfiles.append(line)
   print "inputfiles : \n", inputfiles, "\n"
   process.source.fileNames = inputfiles
else:
   process.source.fileNames = [
      #'file:/uscms_data/d3/mkilpatr/CMSSW_9_4_2/src/AnalysisBase/Analyzer/test/9EE984CF-39E7-E711-A918-001E67DFF67C.root'
     #'/store/mc/RunIIFall17MiniAOD/QCD_Pt_1400to1800_TuneCP5_13TeV_pythia8/MINIAODSIM/94X_mc2017_realistic_v10-v1/60000/06634AB3-14E6-E711-8DA6-0CC47A4D7632.root'
       #'/store/mc/RunIIFall17MiniAOD/TTJets_SingleLeptFromT_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/94X_mc2017_realistic_v10-v1/00000/00257B91-1808-E811-BD39-0242AC130002.root'
      #'/store/data/Run2018A/MET/MINIAOD/PromptReco-v3/000/316/569/00000/08838869-3266-E811-B1E2-02163E019F28.root '
      #'/store/relval/CMSSW_10_1_7/RelValNuGun/MINIAODSIM/PU25ns_101X_upgrade2018_realistic_HEmiss_v1-v1/10000/40CC6E9C-0A80-E811-B01B-0CC47A78A2EC.root '
      #'/store/relval/CMSSW_10_1_7/RelValTTbar_13/MINIAODSIM/PU25ns_101X_upgrade2018_realistic_HEmiss_v1-v1/10000/6002343D-1780-E811-906C-0CC47A7C3422.root'
      '/store/user/benwu/Stop18/NtupleSyncMiniAOD/00257B91-1808-E811-BD39-0242AC130002.root'
      #'/store/mc/RunIIFall17MiniAODv2/SMS-T2tt_mStop-850_mLSP-100_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/9C042E6F-797C-E811-A553-FA163E7ABBFC.root'
      #'/store/mc/RunIIFall17MiniAOD/TTJets_SingleLeptFromT_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/94X_mc2017_realistic_v10-v1/00000/00257B91-1808-E811-BD39-0242AC130002.root'
      #'/store/relval/CMSSW_10_1_7/EGamma/MINIAOD/101X_dataRun2_Prompt_HEmiss_v1_RelVal_EGamma2018B-v1/10000/00BCDF45-3680-E811-9772-0CC47A7C3638.root'
      #'/store/relval/CMSSW_10_1_7/EGamma/MINIAOD/101X_dataRun2_Prompt_v11_RelVal_EGamma2018B-v1/10000/169FB560-2E80-E811-BB39-0CC47A7C360E.root' 
      ]

## ---------------------
## -- Calibration tag --
## ---------------------

from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag

if options.GlobalTag:
   process.GlobalTag = GlobalTag(process.GlobalTag, options.GlobalTag, '')

if options.mcInfo == False: 
   options.jetCorrections.append('L2L3Residual')
options.jetCorrections.insert(0, 'L1FastJet')

########################################################################################################################################################

## get leptons
process.load("StopTupleMaker.SkimsAUX.prodMuons_cfi")
process.load("StopTupleMaker.SkimsAUX.prodElectrons_cfi")

########################################################################################################################################################

## ---------------------------------------
## -- Create all needed jet collections --
## ---------------------------------------

## define the JECs JET Corrections
jetCorrLevelLists = ['L1FastJet', 'L2Relative', 'L3Absolute']
jetCorrectionLevels = ('AK4PFchs', cms.vstring(jetCorrLevelLists), 'None')
if options.mcInfo == False:
      jetCorrLevelLists = ['L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual']
      jetCorrectionLevels = ('AK4PFchs', cms.vstring(jetCorrLevelLists), 'None')

AK4_btagDiscriminators = [ 'pfCombinedInclusiveSecondaryVertexV2BJetTags',
                           'softPFMuonBJetTags',
                           'softPFElectronBJetTags',
                           'pfJetBProbabilityBJetTags',
                           'pfJetProbabilityBJetTags',
                           'pfCombinedCvsLJetTags',
                           'pfCombinedCvsBJetTags',
                           'pfCombinedSecondaryVertexV2BJetTags',
                           'pfDeepCSVJetTags:probudsg',
                           'pfDeepCSVJetTags:probb',
                           'pfDeepCSVJetTags:probc',
                           'pfDeepCSVJetTags:probbb',
                           'pfDeepCSVDiscriminatorsJetTags:BvsAll',
                           ]


from RecoJets.JetProducers.ak4PFJets_cfi import ak4PFJets

#Clean well ID and isolated e/mu out of PF candidates
process.pfPVCleanedCHS = cms.EDFilter("CandPtrSelector",
                                      src = cms.InputTag("packedPFCandidates"),
                                      cut = cms.string("fromPV") )
process.pfNoMuonCHSNoMu =  cms.EDProducer("CandPtrProjector", 
                                          src = cms.InputTag("pfPVCleanedCHS"), 
                                          veto = cms.InputTag("prodMuons", "mu2Clean"))
process.pfNoElectronCHSNoEle = cms.EDProducer("CandPtrProjector", 
                                              src = cms.InputTag("pfNoMuonCHSNoMu"), 
                                              veto = cms.InputTag("prodElectrons", "ele2Clean"))

#recluster AK4 jets from leptonclened PF candidates 
process.ak4PFJetsCHSNoLep = ak4PFJets.clone(src = 'pfNoElectronCHSNoEle', doAreaFastjet = True) # no idea while doArea is false by default, but it's True in RECO so we have to set it

##Override default JEC source 
process.jec = cms.ESSource("PoolDBESSource",
    DBParameters = cms.PSet(
      messageLevel = cms.untracked.int32(0)
    ),
    timetype = cms.string('runnumber'),
    toGet = cms.VPSet(
      cms.PSet(
        record = cms.string('JetCorrectionsRecord'),
        tag    = cms.string('JetCorrectorParametersCollection_'+options.jecDBname+"_AK4PFchs"),
        label  = cms.untracked.string('AK4PFchs')
      ),
      ## here you add as many jet types as you need
      ## note that the tag name is specific for the particular sqlite file 
    ),
    # from page 19 on slides https://indico.cern.ch/event/405326/contribution/2/attachments/811719/1112498/Pythia8.pdf
    # connect = cms.string('sqlite:PY8_RunIISpring15DR74_bx25_MC.db')
    connect = cms.string('sqlite:'+options.jecDBname+'.db')
  )
## add an es_prefer statement to resolve a possible conflict from simultaneous connection to a global tag
process.es_prefer_jec = cms.ESPrefer('PoolDBESSource','jec')

#Calculate QG variables for AK4 jet collections 
process.load('RecoJets.JetProducers.QGTagger_cfi')
process.QGTagger.srcJets   = cms.InputTag('slimmedJets') 

process.QGTaggerNoLep = process.QGTagger.clone()
process.QGTaggerNoLep.srcJets = cms.InputTag("ak4PFJetsCHSNoLep")

from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection
updateJetCollection(
      process,
      postfix = "",
      labelName = 'AK4PFCHS',
      jetSource = cms.InputTag('slimmedJets'),
      pvSource = cms.InputTag('offlineSlimmedPrimaryVertices'),
      pfCandidates = cms.InputTag('packedPFCandidates'),
      svSource = cms.InputTag('slimmedSecondaryVertices'),
      elSource = cms.InputTag('slimmedElectrons'),
      muSource = cms.InputTag('slimmedMuons'),
      jetCorrections = jetCorrectionLevels,
      btagDiscriminators = AK4_btagDiscriminators,
)
process.updatedPatJetsAK4PFCHS.userData.userFloats.src += ['QGTagger:qgLikelihood','QGTagger:ptD', 'QGTagger:axis2', 'QGTagger:axis1']
process.updatedPatJetsAK4PFCHS.userData.userInts.src += ['QGTagger:mult']

from PhysicsTools.PatAlgos.tools.jetTools import addJetCollection
addJetCollection(
      process,
      postfix = "",
      labelName = 'AK4PFCHSNoLep',
      jetSource = cms.InputTag('ak4PFJetsCHSNoLep'),
      pvSource = cms.InputTag('offlineSlimmedPrimaryVertices'),
      pfCandidates = cms.InputTag('pfNoElectronCHSNoEle'),
      svSource = cms.InputTag('slimmedSecondaryVertices'),
      elSource = cms.InputTag('slimmedElectrons'),
      muSource = cms.InputTag('slimmedMuons'),
      jetCorrections = jetCorrectionLevels,
      btagDiscriminators = AK4_btagDiscriminators,
      genParticles = cms.InputTag('prunedGenParticles'),
      genJetCollection = cms.InputTag('slimmedGenJets'),
)
process.patJetsAK4PFCHSNoLep.userData.userFloats.src += ['QGTaggerNoLep:qgLikelihood','QGTaggerNoLep:ptD', 'QGTaggerNoLep:axis2', 'QGTaggerNoLep:axis1']
process.patJetsAK4PFCHSNoLep.userData.userInts.src += ['QGTaggerNoLep:mult']

#if "JEC" in options.specialFix:
#  print ("\nApplying fix to JEC issues in %s ...\n" %(options.cmsswVersion))
#  # JEC can be downloaded from https://twiki.cern.ch/twiki/bin/viewauth/CMS/JECDataMC
#  #inputDB = "sqlite_file:" + os.environ['CMSSW_BASE'] + "/src/StopTupleMaker/SkimsAUX/data/PY8_RunIISpring15DR74_bx25_MC.db"
#  #print inputDB
#   
#
#  #if options.cmsswVersion == "80X":
#  updateJetCollection(
#      process,
#      jetSource = cms.InputTag('slimmedJets'),
#      postfix = 'UpdatedJEC',
#      jetCorrections = ('AK4PFchs', jetCorrLevelLists, 'None')
#  )
#  process.updatedPatJetsUpdatedJEC.userData.userFloats.src += ['QGTagger:qgLikelihood','QGTagger:ptD', 'QGTagger:axis2', 'QGTagger:axis1']
#  process.updatedPatJetsUpdatedJEC.userData.userInts.src += ['QGTagger:mult']
#  # update the MET to account for the new JECs
#  from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD
#  runMetCorAndUncFromMiniAOD(
#      process,
#      isData=not options.mcInfo, # controls gen met
#      #jetCollUnskimmed="updatedPatJetsUpdatedJEC",
#      #jetColl="updatedPatJetsUpdatedJEC",
#      #postfix="Update"
#  )
#  #process.fix80XJEC = cms.Sequence( process.patJetCorrFactorsUpdatedJEC + process.updatedPatJetsUpdatedJEC ) ## NS: What is this patJetCorrFactorsUpdatedJEC ??
    

from JMEAnalysis.JetToolbox.jetToolbox_cff import jetToolbox

outputModuleName = "out"
jetToolbox( process, 'ak8', 'ak8JetSubsNoLep', outputModuleName, 
            runOnMC = options.mcInfo, 
            PUMethod='Puppi', 
            newPFCollection=True,
            nameNewPFCollection='pfNoElectronCHSNoEle',
            addSoftDropSubjets = True, 
            addSoftDrop = True, 
            addNsub = True, 
            bTagDiscriminators = ['pfCombinedInclusiveSecondaryVertexV2BJetTags',
                                  'pfCombinedSecondaryVertexV2BJetTags',
                                  'pfDeepCSVJetTags:probudsg',
                                  'pfDeepCSVJetTags:probb',
                                  'pfDeepCSVJetTags:probc',
                                  'pfDeepCSVJetTags:probbb',
                                  ], 
            addCMSTopTagger = False,
            postFix="NoLep")

# Keep this behind the cleaned version for now, otherwise everything will be lepton cleaned
jetToolbox( process, 'ak8', 'ak8JetSubs', outputModuleName, 
            runOnMC = options.mcInfo, 
            PUMethod='Puppi', 
            addSoftDropSubjets = True, 
            addSoftDrop = True, 
            addNsub = True, 
            bTagDiscriminators = ['pfCombinedInclusiveSecondaryVertexV2BJetTags',
                                  'pfCombinedSecondaryVertexV2BJetTags',
                                  'pfDeepCSVJetTags:probudsg',
                                  'pfDeepCSVJetTags:probb',
                                  'pfDeepCSVJetTags:probc',
                                  'pfDeepCSVJetTags:probbb',
                                  ], 
            addCMSTopTagger = False)

# Hack to remove the stupid output module created by jetToolbox
# which makes this infuriating "jettoolbox.root" file
if hasattr(process, outputModuleName):
   delattr(process, outputModuleName)

process.load("StopTupleMaker.SkimsAUX.prodJets_cfi")

process.prodJets.jetSrc = cms.InputTag('updatedPatJetsAK4PFCHS')
process.prodJets.puppiJetsSrc = cms.InputTag('packedPatJetsAK8PFPuppiSoftDrop')
process.prodJets.debug = cms.bool(options.debug)

process.prodJetsNoLep = process.prodJets.clone()
process.prodJetsNoLep.jetSrc = cms.InputTag('patJetsAK4PFCHSNoLep')
process.prodJetsNoLep.ak8PFJetsPuppi_label = cms.string('ak8PFJetsPuppiNoLep')
process.prodJetsNoLep.qgTaggerKey = cms.string('QGTaggerNoLep')
process.prodJetsNoLep.puppiJetsSrc = cms.InputTag('packedPatJetsAK8PFPuppiNoLepSoftDrop')
process.prodJetsNoLep.NjettinessAK8Puppi_label = cms.string('NjettinessAK8PuppiNoLep')

########################################################################################################################################################

process.load('StopTupleMaker.SkimsAUX.prodJetIDEventFilter_cfi')
process.prodJetIDEventFilter.JetSource = cms.InputTag("updatedPatJetsAK4PFCHS")
process.prodJetIDEventFilter.MinJetPt  = cms.double(10.0)
process.prodJetIDEventFilter.MaxJetEta = cms.double(999.0)

process.prodJetIDEventFilterNoLep = process.prodJetIDEventFilter.clone()
process.prodJetIDEventFilterNoLep.JetSource = cms.InputTag("packedPatJetsAK8PFPuppiNoLepSoftDrop")

########################################################################################################################################################

process.load("StopTupleMaker.SkimsAUX.ISRJetProducer_cfi")

########################################################################################################################################################

#recalculate MET because of updated JEC
#https://twiki.cern.ch/twiki/bin/view/CMS/MissingETUncertaintyPrescription#PF_MET
from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD
runMetCorAndUncFromMiniAOD(
   process,
   isData=not options.mcInfo, # controls gen met
#   fixEE2017 = True,
#   fixEE2017Params = {'userawPt': True, 'PtThreshold':50.0, 'MinEtaThreshold':2.65, 'MaxEtaThreshold': 3.139} ,
#   postfix = "ModifiedMET"
)

process.load("StopTupleMaker.SkimsAUX.prodMET_cfi")
#make sure we pick up the MET produced above
process.prodMET.metSrc = cms.InputTag("slimmedMETs", "", process.name_())

########################################################################################################################################################

## -- Other Analysis related configuration --
import HLTrigger.HLTfilters.hltHighLevel_cfi as hlt
if options.hltSelection:
   process.hltFilter = hlt.hltHighLevel.clone(
      TriggerResultsTag = cms.InputTag("TriggerResults","",options.hltName),
      HLTPaths = cms.vstring(options.hltSelection),
      throw = True, # Don't throw?!
      andOr = True
   )

# HT 
process.load("StopTupleMaker.Skims.htProducer_cfi")
process.load("StopTupleMaker.Skims.htFilter_cfi")
process.htPFchs = process.ht.clone()
process.htPFchs.JetCollection = cms.InputTag("ak4patJetsPFchsPt50Eta25")
process.htPFchsFilter = process.htFilter.clone()
process.htPFchsFilter.HTSource = cms.InputTag("htPFchs")

# MHT
process.load("StopTupleMaker.Skims.mhtProducer_cfi")
process.load("StopTupleMaker.Skims.mhtFilter_cfi")
process.mhtPFchs = process.mht.clone()
process.mhtPFchs.JetCollection = cms.InputTag("ak4patJetsPFchsPt30")
process.mhtPFchsFilter = process.mhtFilter.clone()
process.mhtPFchsFilter.MHTSource = cms.InputTag("mhtPFchs")

#MT2
#process.load("StopTupleMaker.SkimsAUX.mt2Producer_cfi")
##process.load("StopTupleMaker.SkimsAUX.mt2Filter_cfi")
#process.mt2PFchs = process.mt2.clone()
#process.mt2PFchs.JetTag = cms.InputTag("ak4patJetsPFchsPt30")
#process.mt2PFchs.METTag = cms.InputTag("slimmedMETs")


# Delta Phi
process.load("StopTupleMaker.Skims.jetMHTDPhiFilter_cfi")
process.ak4jetMHTPFchsDPhiFilter = process.jetMHTDPhiFilter.clone()
process.ak4jetMHTPFchsDPhiFilter.JetSource = cms.InputTag("ak4patJetsPFchsPt30")
process.ak4jetMHTPFchsDPhiFilter.MHTSource = cms.InputTag("mhtPFchs")

#process.ra2Objects_task = cms.Task(ra2PFchsJets_task, process.htPFchs, process.mhtPFchs)

#process.ra2Objects = cms.Sequence( 
#                                 process.ra2Objects_task
#                                 )

process.load("PhysicsTools.PatAlgos.selectionLayer1.muonCountFilter_cfi")
process.load("PhysicsTools.PatAlgos.selectionLayer1.electronCountFilter_cfi")

############################# START SUSYPAT specifics ####################################
process.prefilterCounter        = cms.EDProducer("EventCountProducer")
process.postStdCleaningCounter  = cms.EDProducer("EventCountProducer")

process.cleanpatseq_task  = cms.Task(process.postStdCleaningCounter)

# Standard Event cleaning 
process.load("StopTupleMaker.SkimsAUX.prodFilterOutScraping_cfi")
process.load("StopTupleMaker.SkimsAUX.prodGoodVertices_cfi")
process.load("StopTupleMaker.SkimsAUX.prodSecondaryVertex_cfi")

# an example sequence to create skimmed susypat-tuples
process.cleanpatseq = cms.Sequence(
#                      process.ra2StdCleaning          *
                      #process.postStdCleaningCounter  #*
                      process.cleanpatseq_task
                      )
############################# EDN SUSYPAT specifics ####################################

process.dummyCounter = cms.EDProducer("EventCountProducer")

process.load('StopTupleMaker.SkimsAUX.weightProducer_cfi')
process.weightProducer.inputPUfileMC   = cms.untracked.string("")
process.weightProducer.inputPUfileData = cms.untracked.string("")
if options.doPtHatWeighting:
  process.weightProducer.Method     = cms.string("PtHat")
  process.weightProducer.Exponent   = cms.double(-4.5)
  process.weightProducer.XS         = cms.double(1.0)
  process.weightProducer.NumberEvts = cms.double(1.0)
  process.weightProducer.Lumi       = cms.double(1.0)
  process.weightProducer.weightWARNingUpThreshold  = cms.double(2.0)

from JetMETCorrections.Configuration.DefaultJEC_cff import *
process.ak4PFJetschsL1FastL2L3Residual = ak4PFJetsL1FastL2L3Residual.clone( algorithm = cms.string('AK4PFchs'), src = 'slimmedJets' )
process.ak4PFJetschsL1FastL2L3 = ak4PFJetsL1FastL2L3.clone( algorithm = cms.string('AK4PFchs'), src = 'slimmedJets' )

# Default is dR = 0.3, dz < 0.05, pt > 10, reliso < 0.1
process.load("StopTupleMaker.Skims.trackIsolationMaker_cfi")
process.trackIsolation = process.trackIsolationFilter.clone()
#process.trackIsolation.pfCandidatesTag = cms.InputTag("packedPFCandidates")
process.trackIsolation.doTrkIsoVeto = cms.bool(False)

process.loosetrackIsolation = process.trackIsolation.clone()
#process.loosetrackIsolation.minPt_PFCandidate = cms.double(5.0)
process.loosetrackIsolation.isoCut            = cms.double(0.5)

#process.refalltrackIsolation = process.trackIsolation.clone()
#process.refalltrackIsolation.mintPt_PFCandidate = cms.double (-1.0)
#process.refalltrackIsolation.isoCut           = cms.double(9999.0)

process.load('StopTupleMaker.Skims.StopJets_drt_from_AOD_cff')
process.load("StopTupleMaker.SkimsAUX.nJetsForSkimsRA2_cfi")
process.load("StopTupleMaker.SkimsAUX.jetMHTDPhiForSkimsRA2_cfi")

# ak4 jets
process.ak4stopJetsPFchsPt30 = process.stopJetsPFchsPt30.clone(jetSrc = "slimmedJets")
process.ak4stopJetsPFchsPt50Eta24 = process.stopJetsPFchsPt50Eta24.clone(jetSrc = "slimmedJets")

process.ak4nJetsForSkimsStop = process.nJetsForSkimsRA2.clone()
process.ak4nJetsForSkimsStop.JetSource = cms.InputTag("ak4stopJetsPFchsPt30")
process.ak4jetMHTDPhiForSkimsStop = process.jetMHTDPhiForSkimsRA2.clone()
process.ak4jetMHTDPhiForSkimsStop.MHTSource = cms.InputTag("slimmedMETs")
process.ak4jetMHTDPhiForSkimsStop.JetSource = cms.InputTag("ak4stopJetsPFchsPt30")

process.ak4stophtPFchs = process.htPFchs.clone()
process.ak4stophtPFchs.JetCollection = cms.InputTag("ak4stopJetsPFchsPt50Eta24")

process.ak4stopmhtPFchs = process.mhtPFchs.clone()
process.ak4stopmhtPFchs.JetCollection = cms.InputTag("ak4stopJetsPFchsPt30")
#

process.prepareCutvars_task = cms.Task(process.ak4stopJetsPFchsPt30, process.ak4stopJetsPFchsPt50Eta24, process.ak4nJetsForSkimsStop, process.ak4jetMHTDPhiForSkimsStop, process.ak4stophtPFchs, process.ak4stopmhtPFchs )

process.prepareCutvars_seq = cms.Sequence( process.prepareCutvars_task )

#############################
# Joe just put these here, maybe the go better slightlt higher up
############################

process.load("StopTupleMaker.SkimsAUX.prodGenJets_cfi")
process.load("StopTupleMaker.SkimsAUX.prodGenInfo_cfi")
process.load("StopTupleMaker.SkimsAUX.prodIsoTrks_cfi")
process.load("StopTupleMaker.SkimsAUX.prodEventInfo_cfi")
process.load("StopTupleMaker.SkimsAUX.PhotonIDisoProducer_cfi")

process.load("StopTupleMaker.Skims.StopBTagJets_cff")
process.stopBJets.JetSrc = cms.InputTag("stopJetsPFchsPt30")

#if hasattr(process, 'goodPhotons'):
#  print blahAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA

process.load("StopTupleMaker.Skims.StopDPhiSelection_cff")
process.jetsMETDPhiFilter.jetSrc = cms.InputTag("stopJetsPFchsPt30")
if options.usePhiCorrMET == True:
   process.jetsMETDPhiFilter.metSrc = cms.InputTag("slimmedMETs")
else:
   process.jetsMETDPhiFilter.metSrc = cms.InputTag("slimmedMETs")
process.jetsMETDPhiFilter.dPhiCuts = cms.untracked.vdouble(0.5, 0.5, 0.3)

process.stopCount1BJets = process.stopCountBJets.clone()
process.stopCount1BJets.minNumber = cms.uint32(1)

process.stopCount2BJets = process.stopCountBJets.clone()
process.stopCount2BJets.minNumber = cms.uint32(2)

#process.load("StopTupleMaker.Skims.StopType3TopTagger_cff")
#if options.usePhiCorrMET == True:
#   process.type3topTagger.metSrc = cms.InputTag("slimmedMETs")
#else:
#   process.type3topTagger.metSrc = cms.InputTag("slimmedMETs")
#process.type3topTagger.taggingMode = cms.untracked.bool(True)
#process.type3topTagger.jetSrc = cms.InputTag("stopJetsPFchsPt30")

process.metPFchsFilter = process.mhtPFchsFilter.clone()
if options.usePhiCorrMET == True:
   process.metPFchsFilter.MHTSource = cms.InputTag("slimmedMETs")
else:
   process.metPFchsFilter.MHTSource = cms.InputTag("slimmedMETs")

process.met175PFchsFilter = process.metPFchsFilter.clone()
process.met175PFchsFilter.MinMHT = cms.double(175)

process.met200PFchsFilter = process.metPFchsFilter.clone()
process.met200PFchsFilter.MinMHT = cms.double(200)

process.met350PFchsFilter = process.metPFchsFilter.clone()
process.met350PFchsFilter.MinMHT = cms.double(350)

process.TFileService = cms.Service("TFileService",
   fileName = cms.string('stopFlatNtuples.root')
)

process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

process.load("StopTupleMaker.SkimsAUX.genDecayStringMakerPythia8_cfi")
process.printDecayPythia8.src = cms.InputTag("prunedGenParticles")
process.printDecayPythia8.keyDecayStrs = cms.vstring("t", "tbar", "~chi_1+", "~chi_1-")
process.printDecayPythia8.printDecay = cms.untracked.bool(options.debug)


#process.load("StopTupleMaker.TopTagger.groomProd_cfi")
#process.groomProdak4 = process.groomProd.clone()
#process.groomProdak4.jetSrc = cms.InputTag("ak4patJetsPFchsPt10")
#process.groomProdak4.groomingOpt = cms.untracked.int32(1)
#process.groomProdak4.debug = cms.untracked.bool(options.debug)


# See https://twiki.cern.ch/twiki/bin/view/CMSPublic/ReMiniAOD03Feb2017Notes#MET_Recipes
# Also https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2
# This is special treatment for reMINIAOD DATA...
if "BADMUON" in options.specialFix and options.mcInfo == False:
#   # Note that calo MET is stored only in the slimmedMETs collection, therefore addcalomet is True only in the slimmedMETsDefault (with slimmedMETs as source)
   process.prodMETslimmedMETsDefault     = process.prodMET.clone(metSrc = cms.InputTag("slimmedMETs"),            addcalomet = cms.bool(True) ) # This collection is the muon cleaned MET only.
   process.prodMETslimmedMETsUncorrected = process.prodMET.clone(metSrc = cms.InputTag("slimmedMETsUncorrected"), addcalomet = cms.bool(False)) # This is the uncleaned MET collection.
   process.prodMETslimmedMETsEGClean     = process.prodMET.clone(metSrc = cms.InputTag("slimmedMETsEGClean"),     addcalomet = cms.bool(False)) # This collection is the e/gamma cleaned MET only.
   # The most correct MET collection is slimmedMETsMuEGClean, this is the collection corrected by both e/gamma and muon effects
   process.prodMET.metSrc = cms.InputTag("slimmedMETsMuEGClean") # This collection is the muon and e/gamma cleaned MET.
   process.prodMET.addcalomet = cms.bool(False)

#Addition of Filter Decision Bits and Trigger Results
process.load("StopTupleMaker.SkimsAUX.prodTriggerResults_cfi")
process.load("StopTupleMaker.SkimsAUX.prodFilterFlags_cfi")

#rerun HBHE noise filter manually
process.load('CommonTools.RecoAlgos.HBHENoiseFilterResultProducer_cfi')
process.HBHENoiseFilterResultProducer.minZeros = cms.int32(99999)
process.HBHENoiseFilterResultProducer.IgnoreTS4TS5ifJetInLowBVRegion=cms.bool(False)
process.HBHENoiseFilterResultProducer.defaultDecision = cms.string("HBHENoiseFilterResultRun2Loose")

#New met filters for 2016
process.load('RecoMET.METFilters.BadChargedCandidateFilter_cfi')
process.load('RecoMET.METFilters.BadPFMuonFilter_cfi')
process.BadChargedCandidateFilter.taggingMode = cms.bool(True)
process.BadChargedCandidateFilter.muons = cms.InputTag("slimmedMuons")
process.BadChargedCandidateFilter.PFCandidates = cms.InputTag("packedPFCandidates")
process.BadPFMuonFilter.taggingMode = cms.bool(True)
process.BadPFMuonFilter.muons = cms.InputTag("slimmedMuons")
process.BadPFMuonFilter.PFCandidates = cms.InputTag("packedPFCandidates")

#Bad muon filger for Moriond 2017 https://hypernews.cern.ch/HyperNews/CMS/get/physics-validation/2786.html
process.load('RecoMET.METFilters.badGlobalMuonTaggersMiniAOD_cff')
process.badGlobalMuonTaggerMAOD.taggingMode = cms.bool(True)
process.cloneGlobalMuonTaggerMAOD.taggingMode = cms.bool(True)

process.triggerProducer.trigTagSrc = cms.InputTag("TriggerResults","",options.hltName)
process.METFilters = process.filterDecisionProducer.clone( filterName  =   cms.string("Flag_METFilters") , trigTagSrc = cms.InputTag("TriggerResults",processName=cms.InputTag.skipCurrentProcess()) )
process.CSCTightHaloFilter = process.filterDecisionProducer.clone( filterName  =   cms.string("Flag_CSCTightHaloFilter") , trigTagSrc = cms.InputTag("TriggerResults",processName=cms.InputTag.skipCurrentProcess()))
#process.globalTightHalo2016Filter = process.filterDecisionProducer.clone( filterName  =   cms.string("Flag_globalTightHalo2016Filter") , trigTagSrc = cms.InputTag("TriggerResults",processName=cms.InputTag.skipCurrentProcess()))
process.globalSuperTightHalo2016Filter = process.filterDecisionProducer.clone( filterName  =   cms.string("Flag_globalSuperTightHalo2016Filter") ,      trigTagSrc = cms.InputTag("TriggerResults",processName=cms.InputTag.skipCurrentProcess()))
process.goodVerticesFilter = process.filterDecisionProducer.clone( filterName  =   cms.string("Flag_goodVertices") , trigTagSrc = cms.InputTag("TriggerResults",processName=cms.InputTag.skipCurrentProcess()))
process.ecalBadCalibFilter = process.filterDecisionProducer.clone( filterName  =   cms.string("Flag_ecalBadCalibFilter"), trigTagSrc = cms.InputTag("TriggerResults",processName=cms.InputTag.skipCurrentProcess()) )
process.HBHENoiseIsoFilter = process.filterDecisionProducer.clone( filterName  =   cms.string("Flag_HBHENoiseIsoFilter"), trigTagSrc = cms.InputTag("TriggerResults",processName=cms.InputTag.skipCurrentProcess()) )
process.HBHENoiseFilter = process.filterDecisionProducer.clone( filterName  =   cms.string("Flag_HBHENoiseFilter") )
process.EcalDeadCellTriggerPrimitiveFilter = process.filterDecisionProducer.clone( filterName  =   cms.string("Flag_EcalDeadCellTriggerPrimitiveFilter"),trigTagSrc = cms.InputTag("TriggerResults",processName=cms.InputTag.skipCurrentProcess()) )


process.filterDecisionProducerPAT = process.filterDecisionProducer.clone()
process.filterDecisionProducerPAT.trigTagSrc = cms.InputTag("TriggerResults","","PAT")
process.noBadMuonsFilter = process.filterDecisionProducerPAT.clone( filterName  =   cms.string("Flag_noBadMuons"),trigTagSrc = cms.InputTag("TriggerResults",processName=cms.InputTag.skipCurrentProcess()) )
process.badMuonsFilter = process.filterDecisionProducerPAT.clone( filterName = cms.string("Flag_badMuons"),trigTagSrc = cms.InputTag("TriggerResults",processName=cms.InputTag.skipCurrentProcess()) )
process.duplicateMuonsFilter = process.filterDecisionProducerPAT.clone( filterName = cms.string("Flag_duplicateMuons"),trigTagSrc = cms.InputTag("TriggerResults",processName=cms.InputTag.skipCurrentProcess()) )

process.prodMuonsNoIso = process.prodMuons.clone()
process.prodMuonsNoIso.DoMuonIsolation = cms.int32(0)

process.prodElectronsNoIso = process.prodElectrons.clone()
process.prodElectronsNoIso.DoElectronIsolation = cms.int32(0)

process.load("StopTupleMaker.StopTreeMaker.stopTreeMaker_cfi")
process.stopTreeMaker.debug = cms.bool(options.debug)
process.stopTreeMaker.TreeName = cms.string("AUX")

process.ntpVersion = cms.EDFilter(
   "prodNtupleVersionString", 
   inputStr = cms.vstring(options.ntpVersion, options.GlobalTag, options.cmsswVersion, options.specialFix, options.hltName, options.era, options.jecDBname, procCMSSWver)
)

#process.stopTreeMaker.vectorString.append(cms.InputTag("ntpVersion"))

process.stopTreeMaker.vectorInt.append(cms.InputTag("triggerProducer", "PassTrigger"))
process.stopTreeMaker.vectorInt.append(cms.InputTag("triggerProducer", "TriggerPrescales"))
process.stopTreeMaker.vectorString.append(cms.InputTag("triggerProducer", "TriggerNames"))

# prodGoodVertices has the same as vtxSize in prodEventInfo...
#process.stopTreeMaker.varsInt.append(cms.InputTag("prodGoodVertices"))
#process.stopTreeMaker.varsInt.append(cms.InputTag("prodFilterOutScraping"))
process.stopTreeMaker.varsBool.append(cms.InputTag("prodJetIDEventFilter", "AK4looseJetID"))
process.stopTreeMaker.varsBool.append(cms.InputTag("prodJetIDEventFilter", "AK4tightJetID"))
process.stopTreeMaker.varsBool.append(cms.InputTag("prodJetIDEventFilter", "AK4tightlepvetoJetID"))
process.stopTreeMaker.varsBool.append(cms.InputTag("prodJetIDEventFilterNoLep", "AK4looseJetID"))
process.stopTreeMaker.varsBoolNamesInTree.append("prodJetIDEventFilterNoLep:AK4looseJetID|prodJetsNoLep_AK4looseJetID")
process.stopTreeMaker.varsBool.append(cms.InputTag("prodJetIDEventFilterNoLep", "AK4tightJetID"))
process.stopTreeMaker.varsBoolNamesInTree.append("prodJetIDEventFilterNoLep:AK4tightJetID|prodJetsNoLep_AK4tightJetID_NoLep")
process.stopTreeMaker.varsBool.append(cms.InputTag("prodJetIDEventFilterNoLep", "AK4tightlepvetoJetID"))
process.stopTreeMaker.varsBoolNamesInTree.append("prodJetIDEventFilterNoLep:AK4tightlepvetoJetID|prodJetsNoLep_AK4tightlepvetoJetID_NoLep")
process.stopTreeMaker.varsInt.append(cms.InputTag("METFilters"))
process.stopTreeMaker.varsInt.append(cms.InputTag("CSCTightHaloFilter")) # 74X txt files are ready for the 2015 working point, use this and not the flag in miniAOD 
process.stopTreeMaker.varsInt.append(cms.InputTag("globalSuperTightHalo2016Filter"))
process.stopTreeMaker.varsInt.append(cms.InputTag("goodVerticesFilter"))
process.stopTreeMaker.varsInt.append(cms.InputTag("ecalBadCalibFilter"))
process.stopTreeMaker.varsInt.append(cms.InputTag("HBHENoiseIsoFilter"))
process.stopTreeMaker.varsInt.append(cms.InputTag("EcalDeadCellTriggerPrimitiveFilter"))
process.stopTreeMaker.varsBool.append(cms.InputTag("BadChargedCandidateFilter"))
process.stopTreeMaker.varsBool.append(cms.InputTag("BadPFMuonFilter"))

process.stopTreeMaker.varsInt.append(cms.InputTag("noBadMuonsFilter"))
process.stopTreeMaker.varsInt.append(cms.InputTag("badMuonsFilter"))
process.stopTreeMaker.varsInt.append(cms.InputTag("duplicateMuonsFilter"))

if options.fastsim == False:
   process.stopTreeMaker.varsBool.append(cms.InputTag("HBHENoiseFilterResultProducer", "HBHENoiseFilterResult"))
   process.stopTreeMaker.varsBoolNamesInTree.append("HBHENoiseFilterResultProducer:HBHENoiseFilterResult|HBHENoiseFilter")

   process.stopTreeMaker.varsBool.append(cms.InputTag("HBHENoiseFilterResultProducer", "HBHEIsoNoiseFilterResult"))
   process.stopTreeMaker.varsBoolNamesInTree.append("HBHENoiseFilterResultProducer:HBHEIsoNoiseFilterResult|HBHEIsoNoiseFilter")

process.stopTreeMaker.varsInt.append(cms.InputTag("prodMuons", "nMuons"))
process.stopTreeMaker.varsIntNamesInTree.append("prodMuons:nMuons|nMuons_CUT")
process.stopTreeMaker.varsInt.append(cms.InputTag("prodMuonsNoIso", "nMuons"))
process.stopTreeMaker.vectorTLorentzVector.append(cms.InputTag("prodMuonsNoIso", "muonsLVec"))

process.stopTreeMaker.vectorFloat.extend([cms.InputTag("prodMuonsNoIso", "muonsCharge"), cms.InputTag("prodMuonsNoIso", "muonsMtw"), cms.InputTag("prodMuonsNoIso", "muonsRelIso"), cms.InputTag("prodMuonsNoIso", "muonsMiniIso"), cms.InputTag("prodMuonsNoIso", "muonspfActivity")])
process.stopTreeMaker.vectorInt.extend([cms.InputTag("prodMuonsNoIso","muonsFlagLoose"),cms.InputTag("prodMuonsNoIso", "muonsFlagMedium"), cms.InputTag("prodMuonsNoIso", "muonsFlagTight")])

#ANDRES Gamma Var  
from PhysicsTools.SelectorUtils.tools.vid_id_tools import *

switchOnVIDPhotonIdProducer(process, DataFormat.MiniAOD)

# Define which IDs we want to produce
my_photon_id_modules = ['RecoEgamma.PhotonIdentification.Identification.cutBasedPhotonID_Fall17_94X_V1_cff']

# Add them to the VID producer
for idmod in my_photon_id_modules:
   setupAllVIDIdsInModule(process, idmod, setupVIDPhotonSelection)

from RecoEgamma.PhotonIdentification.egmPhotonIDs_cff import * 
loadEgmIdSequence(process, DataFormat.MiniAOD)

# Set ID tags
process.goodPhotons.loosePhotonID = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-Fall17-94X-V1-loose")
process.goodPhotons.mediumPhotonID = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-Fall17-94X-V1-medium")
process.goodPhotons.tightPhotonID = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-Fall17-94X-V1-tight")
            
process.stopTreeMaker.vectorBool.append(cms.InputTag("goodPhotons", "loosePhotonID"))
process.stopTreeMaker.vectorBool.append(cms.InputTag("goodPhotons", "mediumPhotonID"))
process.stopTreeMaker.vectorBool.append(cms.InputTag("goodPhotons", "tightPhotonID"))
                                                  
process.stopTreeMaker.vectorFloat.append(cms.InputTag("goodPhotons", "photonPFGammaIso"))
process.stopTreeMaker.vectorFloat.append(cms.InputTag("goodPhotons", "photonISEB"))
process.stopTreeMaker.vectorFloat.append(cms.InputTag("goodPhotons", "photonGenMatched"))
process.stopTreeMaker.vectorFloat.append(cms.InputTag("goodPhotons", "photonHadTowOverEM"))
process.stopTreeMaker.vectorFloat.append(cms.InputTag("goodPhotons", "photonSigmaIetaIeta"))
process.stopTreeMaker.vectorFloat.append(cms.InputTag("goodPhotons", "photonPFChargedIso"))
process.stopTreeMaker.vectorFloat.append(cms.InputTag("goodPhotons", "photonPFNeutralIso"))
process.stopTreeMaker.vectorFloat.append(cms.InputTag("goodPhotons", "photonPFChargedIsoRhoCorr"))
process.stopTreeMaker.vectorFloat.append(cms.InputTag("goodPhotons", "photonPFNeutralIsoRhoCorr"))
process.stopTreeMaker.vectorFloat.append(cms.InputTag("goodPhotons", "photonPFGammaIsoRhoCorr"))
process.stopTreeMaker.vectorFloat.append(cms.InputTag("goodPhotons", "photonHasPixelSeed"))
process.stopTreeMaker.vectorBool.append(cms.InputTag("goodPhotons", "photonNonPrompt"))
process.stopTreeMaker.vectorTLorentzVector.append(cms.InputTag("goodPhotons", "photonLVec"))
process.stopTreeMaker.vectorTLorentzVector.append(cms.InputTag("goodPhotons", "photonLVecGen"))
process.stopTreeMaker.vectorTLorentzVector.append(cms.InputTag("goodPhotons", "genPartonLVec"))

process.stopTreeMaker.vectorFloat.append(cms.InputTag("prodSecondaryVertex", "svPT"))
process.stopTreeMaker.vectorFloat.append(cms.InputTag("prodSecondaryVertex", "svETA"))
process.stopTreeMaker.vectorFloat.append(cms.InputTag("prodSecondaryVertex", "svPhi"))
process.stopTreeMaker.vectorFloat.append(cms.InputTag("prodSecondaryVertex", "svMass"))
process.stopTreeMaker.vectorFloat.append(cms.InputTag("prodSecondaryVertex", "svNTracks"))
process.stopTreeMaker.vectorFloat.append(cms.InputTag("prodSecondaryVertex", "svChi2"))
process.stopTreeMaker.vectorFloat.append(cms.InputTag("prodSecondaryVertex", "svNDF"))
process.stopTreeMaker.vectorFloat.append(cms.InputTag("prodSecondaryVertex", "svDXY"))
process.stopTreeMaker.vectorFloat.append(cms.InputTag("prodSecondaryVertex", "svDXYerr"))
process.stopTreeMaker.vectorFloat.append(cms.InputTag("prodSecondaryVertex", "svD3D"))
process.stopTreeMaker.vectorFloat.append(cms.InputTag("prodSecondaryVertex", "svD3Derr"))
process.stopTreeMaker.vectorFloat.append(cms.InputTag("prodSecondaryVertex", "svCosThetaSVPS"))
process.stopTreeMaker.vectorTLorentzVector.append(cms.InputTag("prodSecondaryVertex", "svSoftLVec"))
process.stopTreeMaker.vectorTLorentzVector.append(cms.InputTag("prodSecondaryVertex", "svLVec"))

"""
if "BADMUON" in options.specialFix:
   print ("\nAdding bad muon special filter information in prodMuon & prodMuonNoIso ...\n")
   process.prodMuons.specialFix      = cms.bool(True)
   process.prodMuons.badGlobalMuonTaggerSrc = cms.InputTag("badGlobalMuonTaggerMAOD", "bad")
   process.prodMuons.cloneGlobalMuonTaggerSrc = cms.InputTag("cloneGlobalMuonTaggerMAOD", "bad")
   process.prodMuonsNoIso.specialFix = cms.bool(True)
   process.prodMuonsNoIso.badGlobalMuonTaggerSrc = cms.InputTag("badGlobalMuonTaggerMAOD", "bad")
   process.prodMuonsNoIso.cloneGlobalMuonTaggerSrc = cms.InputTag("cloneGlobalMuonTaggerMAOD", "bad")
   process.stopTreeMaker.vectorInt.append(cms.InputTag("prodMuonsNoIso", "specialFixtype"))
   process.stopTreeMaker.vectorTLorentzVector.append(cms.InputTag("prodMuonsNoIso", "specialFixMuonsLVec"))
   process.stopTreeMaker.vectorDouble.append(cms.InputTag("prodMuonsNoIso", "specialFixMuonsCharge"))
"""

process.stopTreeMaker.varsInt.append(cms.InputTag("prodElectrons", "nElectrons"))
process.stopTreeMaker.varsIntNamesInTree.append("prodElectrons:nElectrons|nElectrons_CUT")
process.stopTreeMaker.varsInt.append(cms.InputTag("prodElectronsNoIso", "nElectrons"))
process.stopTreeMaker.vectorTLorentzVector.append(cms.InputTag("prodElectronsNoIso", "elesLVec"))
process.stopTreeMaker.vectorFloat.extend([cms.InputTag("prodElectronsNoIso", "elesCharge"), cms.InputTag("prodElectronsNoIso", "elesMtw"), cms.InputTag("prodElectronsNoIso", "elesRelIso"), cms.InputTag("prodElectronsNoIso", "elesMiniIso"), cms.InputTag("prodElectronsNoIso", "elespfActivity")])
process.stopTreeMaker.vectorBool.extend([cms.InputTag("prodElectronsNoIso", "elesisEB")])
process.stopTreeMaker.vectorInt.extend([cms.InputTag("prodElectronsNoIso", "elesFlagTight"), cms.InputTag("prodElectronsNoIso", "elesFlagMedium"), cms.InputTag("prodElectronsNoIso", "elesFlagLoose"), cms.InputTag("prodElectronsNoIso", "elesFlagVeto")])

my_electron_id_modules = ['RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Fall17_94X_V1_cff']

switchOnVIDElectronIdProducer(process, DataFormat.MiniAOD)
for idmod in my_electron_id_modules:
   setupAllVIDIdsInModule(process, idmod, setupVIDElectronSelection)

from RecoEgamma.ElectronIdentification.egmGsfElectronIDs_cff import *
process.load('RecoEgamma.ElectronIdentification.egmGsfElectronIDs_cff')

process.prodElectronsNoIso.vetoElectronID = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Fall17-94X-V1-veto")
process.prodElectronsNoIso.looseElectronID = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Fall17-94X-V1-loose")
process.prodElectronsNoIso.mediumElectronID = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Fall17-94X-V1-medium")
process.prodElectronsNoIso.tightElectronID = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Fall17-94X-V1-tight")

process.stopTreeMaker.vectorBool.append(cms.InputTag("prodElectronsNoIso","vetoElectronID"))
process.stopTreeMaker.vectorBool.append(cms.InputTag("prodElectronsNoIso","looseElectronID"))
process.stopTreeMaker.vectorBool.append(cms.InputTag("prodElectronsNoIso","mediumElectronID"))
process.stopTreeMaker.vectorBool.append(cms.InputTag("prodElectronsNoIso","tightElectronID"))

#process.stopTreeMaker.varsInt.append(cms.InputTag("prodJets", "nJets"))
process.stopTreeMaker.vectorTLorentzVector.append(cms.InputTag("prodJets", "jetsLVec"))
#process.stopTreeMaker.vectorTLorentzVector.append(cms.InputTag("prodGenJets", "genjetsLVec"))
process.stopTreeMaker.vectorInt.append(cms.InputTag("prodJets", "recoJetsFlavor"))

process.stopTreeMaker.vectorFloat.append(cms.InputTag("prodJets", "recoJetsJecUnc"))
process.stopTreeMaker.vectorFloat.append(cms.InputTag("prodJets", "recoJetsJecScaleRawToFull"))

process.stopTreeMaker.vectorFloat.append(cms.InputTag("prodJets", "qgLikelihood"))
process.stopTreeMaker.vectorFloat.append(cms.InputTag("prodJets", "qgPtD"))
process.stopTreeMaker.vectorFloat.append(cms.InputTag("prodJets", "qgAxis2"))
process.stopTreeMaker.vectorFloat.append(cms.InputTag("prodJets", "qgMult"))
#process.stopTreeMaker.vectorDouble.append(cms.InputTag("prodJets", "qgPtDrLog"))
process.stopTreeMaker.vectorFloat.append(cms.InputTag("prodJets", "qgAxis1"))
#process.stopTreeMaker.vectorInt.append(cms.InputTag("prodJets", "qgnMult"))
#process.stopTreeMaker.vectorInt.append(cms.InputTag("prodJets", "qgcMult"))

process.stopTreeMaker.vectorFloat.append(cms.InputTag("prodJets", "recoJetschargedHadronEnergyFraction"))
process.stopTreeMaker.vectorFloat.append(cms.InputTag("prodJets", "recoJetschargedEmEnergyFraction"))
process.stopTreeMaker.vectorFloat.append(cms.InputTag("prodJets", "recoJetsneutralEmEnergyFraction"))
process.stopTreeMaker.vectorFloat.append(cms.InputTag("prodJets", "recoJetsmuonEnergyFraction"))
process.stopTreeMaker.vectorFloat.append(cms.InputTag("prodJets", "recoJetsneutralEnergyFraction"))
process.stopTreeMaker.vectorFloat.append(cms.InputTag("prodJets", "recoJetsHFEMEnergyFraction"))
process.stopTreeMaker.vectorFloat.append(cms.InputTag("prodJets", "recoJetsHFHadronEnergyFraction"))

process.stopTreeMaker.vectorFloat.append(cms.InputTag("prodJets", "PhotonEnergyFraction"))
process.stopTreeMaker.vectorFloat.append(cms.InputTag("prodJets", "ElectronEnergyFraction"))

process.stopTreeMaker.vectorFloat.append(cms.InputTag("prodJets", "ChargedHadronMultiplicity"))
process.stopTreeMaker.vectorFloat.append(cms.InputTag("prodJets", "NeutralHadronMultiplicity"))
process.stopTreeMaker.vectorFloat.append(cms.InputTag("prodJets", "PhotonMultiplicity"))
process.stopTreeMaker.vectorFloat.append(cms.InputTag("prodJets", "ElectronMultiplicity"))
process.stopTreeMaker.vectorFloat.append(cms.InputTag("prodJets", "MuonMultiplicity"))

process.stopTreeMaker.vectorFloat.append(cms.InputTag("prodJets", "recoJetsCSVv2"))
process.stopTreeMaker.vectorFloatNamesInTree.append("prodJets:recoJetsCSVv2|recoJetsCSVv2")
process.stopTreeMaker.vectorFloat.append(cms.InputTag("prodJets", "recoJetsCharge"))
process.stopTreeMaker.vectorFloatNamesInTree.append("prodJets:recoJetsCharge|recoJetsCharge")

process.stopTreeMaker.vectorFloat.append(cms.InputTag("prodJets", "DeepCSVb"))
process.stopTreeMaker.vectorFloat.append(cms.InputTag("prodJets", "DeepCSVc"))
process.stopTreeMaker.vectorFloat.append(cms.InputTag("prodJets", "DeepCSVl"))
process.stopTreeMaker.vectorFloat.append(cms.InputTag("prodJets", "DeepCSVbb"))
process.stopTreeMaker.vectorFloat.append(cms.InputTag("prodJets", "DeepCSVcc"))

process.stopTreeMaker.vectorFloat.append(cms.InputTag("prodJets", "CversusL"));
process.stopTreeMaker.vectorFloat.append(cms.InputTag("prodJets", "CversusB"));

process.stopTreeMaker.vectorInt.append(cms.InputTag("prodJets", "muMatchedJetIdx"))
process.stopTreeMaker.vectorInt.append(cms.InputTag("prodJets", "eleMatchedJetIdx"))
process.stopTreeMaker.vectorInt.append(cms.InputTag("prodJets", "looseisoTrksMatchedJetIdx"))
process.stopTreeMaker.vectorInt.append(cms.InputTag("prodJets", "trksForIsoVetoMatchedJetIdx"))

process.stopTreeMaker.vectorTLorentzVector.append(cms.InputTag("prodJets", "puppiJetsLVec"))
process.stopTreeMaker.vectorFloat.append(cms.InputTag("prodJets", "puppiTau1"))
process.stopTreeMaker.vectorFloat.append(cms.InputTag("prodJets", "puppiTau2"))
process.stopTreeMaker.vectorFloat.append(cms.InputTag("prodJets", "puppiTau3"))

process.stopTreeMaker.vectorVectorTLorentzVector.append(cms.InputTag("prodJets", "puppiAK8SubjetLVec"))
process.stopTreeMaker.vectorVectorFloat.append(cms.InputTag("prodJets", "puppiAK8SubjetMult"))
process.stopTreeMaker.vectorVectorFloat.append(cms.InputTag("prodJets", "puppiAK8SubjetPtD"))
process.stopTreeMaker.vectorVectorFloat.append(cms.InputTag("prodJets", "puppiAK8SubjetAxis1"))
process.stopTreeMaker.vectorVectorFloat.append(cms.InputTag("prodJets", "puppiAK8SubjetAxis2"))
process.stopTreeMaker.vectorVectorFloat.append(cms.InputTag("prodJets", "puppiAK8SubjetBDisc"))

if options.addJetsForZinv == True: 
   process.stopTreeMaker.vectorTLorentzVector.append(cms.InputTag("prodJetsNoLep", "jetsLVec"))
   process.stopTreeMaker.vectorTLorentzVectorNamesInTree.append("prodJetsNoLep:jetsLVec|prodJetsNoLep_jetsLVec")

   process.stopTreeMaker.vectorFloat.append(cms.InputTag("prodJetsNoLep", "recoJetsJecUnc"))
   process.stopTreeMaker.vectorFloatNamesInTree.append("prodJetsNoLep:recoJetsJecUnc|prodJetsNoLep_recoJetsJecUnc")

   process.stopTreeMaker.vectorFloat.append(cms.InputTag("prodJetsNoLep", "qgLikelihood"))
   process.stopTreeMaker.vectorFloat.append(cms.InputTag("prodJetsNoLep", "qgPtD"))
   process.stopTreeMaker.vectorFloat.append(cms.InputTag("prodJetsNoLep", "qgAxis2"))
   process.stopTreeMaker.vectorFloat.append(cms.InputTag("prodJetsNoLep", "qgAxis1"))
   process.stopTreeMaker.vectorFloat.append(cms.InputTag("prodJetsNoLep", "qgMult"))

   process.stopTreeMaker.vectorFloat.append(cms.InputTag("prodJetsNoLep", "recoJetschargedHadronEnergyFraction"))
   process.stopTreeMaker.vectorFloatNamesInTree.append("prodJetsNoLep:recoJetschargedHadronEnergyFraction|prodJetsNoLep_recoJetschargedHadronEnergyFraction")
   process.stopTreeMaker.vectorFloat.append(cms.InputTag("prodJetsNoLep", "recoJetsneutralEmEnergyFraction"))
   process.stopTreeMaker.vectorFloatNamesInTree.append("prodJetsNoLep:recoJetsneutralEmEnergyFraction|prodJetsNoLep_recoJetsneutralEmEnergyFraction")
   process.stopTreeMaker.vectorFloat.append(cms.InputTag("prodJetsNoLep", "recoJetschargedEmEnergyFraction"))
   process.stopTreeMaker.vectorFloatNamesInTree.append("prodJetsNoLep:recoJetschargedEmEnergyFraction|prodJetsNoLep_recoJetschargedEmEnergyFraction")
   process.stopTreeMaker.vectorFloat.append(cms.InputTag("prodJetsNoLep", "recoJetsmuonEnergyFraction"))
   process.stopTreeMaker.vectorFloatNamesInTree.append("prodJetsNoLep:recoJetsmuonEnergyFraction|prodJetsNoLep_recoJetsmuonEnergyFraction")
    
   process.stopTreeMaker.vectorFloat.append(cms.InputTag("prodJetsNoLep", "ChargedHadronMultiplicity"))
   process.stopTreeMaker.vectorFloatNamesInTree.append("prodJetsNoLep:ChargedHadronMultiplicity|prodJetsNoLep_ChargedHadronMultiplicity")
  
   process.stopTreeMaker.vectorFloat.append(cms.InputTag("prodJetsNoLep", "recoJetsneutralEnergyFraction"))
   process.stopTreeMaker.vectorFloat.append(cms.InputTag("prodJetsNoLep", "recoJetsHFEMEnergyFraction"));
   process.stopTreeMaker.vectorFloat.append(cms.InputTag("prodJetsNoLep", "recoJetsHFHadronEnergyFraction"))

   process.stopTreeMaker.vectorFloat.append(cms.InputTag("prodJetsNoLep", "PhotonEnergyFraction"))
   process.stopTreeMaker.vectorFloat.append(cms.InputTag("prodJetsNoLep", "ElectronEnergyFraction"))

   process.stopTreeMaker.vectorFloat.append(cms.InputTag("prodJetsNoLep", "NeutralHadronMultiplicity"))
   process.stopTreeMaker.vectorFloat.append(cms.InputTag("prodJetsNoLep", "PhotonMultiplicity"))
   process.stopTreeMaker.vectorFloat.append(cms.InputTag("prodJetsNoLep", "ElectronMultiplicity"))
   process.stopTreeMaker.vectorFloat.append(cms.InputTag("prodJetsNoLep", "MuonMultiplicity"))
   
   process.prodJetsNoLep.jetPBJetTags = cms.string('pfJetBProbabilityBJetTags')
   process.prodJetsNoLep.jetPNegBJetTags= cms.string('pfNegativeOnlyJetBProbabilityBJetTags')
   process.prodJetsNoLep.jetPPosBJetTags= cms.string('pfPositiveOnlyJetBProbabilityBJetTags')
   process.prodJetsNoLep.jetBPBJetTags= cms.string('jetBPBJetTags')
   process.prodJetsNoLep.jetBPNegBJetTags= cms.string('jetBPNegBJetTags')
   process.prodJetsNoLep.jetBPPosBJetTags= cms.string('jetBPPosBJetTags')
    
   process.prodJetsNoLep.jetPBJetTags = cms.string('jetPBJetTags')
   process.prodJetsNoLep.jetPNegBJetTags= cms.string('jetPNegBJetTags')
   process.prodJetsNoLep.jetPPosBJetTags= cms.string('jetPPosBJetTags')
   
   process.prodJetsNoLep.jetBPBJetTags= cms.string('jetBPBJetTags')
   process.prodJetsNoLep.jetBPNegBJetTags= cms.string('jetBPNegBJetTags')
   process.prodJetsNoLep.jetBPPosBJetTags= cms.string('jetBPPosBJetTags')
           
   #process.prodJetsNoLep.deepCSVBJetTags= cms.string('pfDeepCSVJetTags')
   #process.prodJetsNoLep.deepCSVNegBJetTags= cms.string('deepCSVNegBJetTags')
   #process.prodJetsNoLep.deepCSVPosBJetTags= cms.string('deepCSVPosBJetTags')

   process.stopTreeMaker.vectorInt.append(cms.InputTag("prodJetsNoLep", "recoJetsFlavor"))

   process.stopTreeMaker.vectorFloat.append(cms.InputTag("prodJetsNoLep", "DeepCSVb"))
   process.stopTreeMaker.vectorFloat.append(cms.InputTag("prodJetsNoLep", "DeepCSVc"))
   process.stopTreeMaker.vectorFloat.append(cms.InputTag("prodJetsNoLep", "DeepCSVl"))
   process.stopTreeMaker.vectorFloat.append(cms.InputTag("prodJetsNoLep", "DeepCSVbb"))
   process.stopTreeMaker.vectorFloat.append(cms.InputTag("prodJetsNoLep", "DeepCSVcc"))
      
       
   process.stopTreeMaker.vectorFloat.append(cms.InputTag("prodJets", "CversusL"));
   process.stopTreeMaker.vectorFloat.append(cms.InputTag("prodJets", "CversusB"));

   process.stopTreeMaker.vectorFloat.append(cms.InputTag("prodJetsNoLep", "recoJetsCSVv2"))
   process.stopTreeMaker.vectorFloatNamesInTree.append("prodJetsNoLep:recoJetsCSVv2|prodJetsNoLep_recoJetsCSVv2")
   process.stopTreeMaker.vectorFloat.append(cms.InputTag("prodJetsNoLep", "recoJetsCharge"))
   process.stopTreeMaker.vectorFloatNamesInTree.append("prodJetsNoLep:recoJetsCharge|prodJetsNoLep_recoJetsCharge")

   process.stopTreeMaker.vectorFloat.append(cms.InputTag("prodJetsNoLep", "recoJetsJecScaleRawToFull"))
   process.stopTreeMaker.vectorFloatNamesInTree.append("prodJetsNoLep:recoJetsJecScaleRawToFull|prodJetsNoLep_recoJetsJecScaleRawToFull")

   process.stopTreeMaker.vectorTLorentzVector.append(cms.InputTag("prodJetsNoLep", "puppiJetsLVec"))
   process.stopTreeMaker.vectorFloat.append(cms.InputTag("prodJetsNoLep", "puppiTau1"))
   process.stopTreeMaker.vectorFloat.append(cms.InputTag("prodJetsNoLep", "puppiTau2"))
   process.stopTreeMaker.vectorFloat.append(cms.InputTag("prodJetsNoLep", "puppiTau3"))
    
   process.stopTreeMaker.vectorVectorTLorentzVector.append(cms.InputTag("prodJetsNoLep", "puppiAK8SubjetLVec"))
   process.stopTreeMaker.vectorVectorFloat.append(cms.InputTag("prodJetsNoLep", "puppiAK8SubjetMult"))
   process.stopTreeMaker.vectorVectorFloat.append(cms.InputTag("prodJetsNoLep", "puppiAK8SubjetPtD"))
   process.stopTreeMaker.vectorVectorFloat.append(cms.InputTag("prodJetsNoLep", "puppiAK8SubjetAxis1"))
   process.stopTreeMaker.vectorVectorFloat.append(cms.InputTag("prodJetsNoLep", "puppiAK8SubjetAxis2"))
   process.stopTreeMaker.vectorVectorFloat.append(cms.InputTag("prodJetsNoLep", "puppiAK8SubjetBDisc"))

if options.mcInfo == True:
   process.prodGenInfo.debug = cms.bool(options.debug)
   process.stopTreeMaker.vectorString.append(cms.InputTag("prodGenInfo", "genDecayStrVec"))
   process.stopTreeMaker.vectorInt.extend([cms.InputTag("prodGenInfo", "genDecayIdxVec"), cms.InputTag("prodGenInfo", "genDecayPdgIdVec"), cms.InputTag("prodGenInfo", "genDecayMomIdxVec"), cms.InputTag("prodGenInfo", "genDecayMomRefVec"), cms.InputTag("prodGenInfo", "WemuVec"), cms.InputTag("prodGenInfo", "WtauVec"), cms.InputTag("prodGenInfo", "WtauemuVec"), cms.InputTag("prodGenInfo", "WtauprongsVec"), cms.InputTag("prodGenInfo", "WtaunuVec"),  cms.InputTag("prodGenInfo","selPDGid")])
   process.stopTreeMaker.vectorIntNamesInTree.extend(["prodGenInfo:WemuVec|W_emuVec", "prodGenInfo:WtauVec|W_tauVec", "prodGenInfo:WtauemuVec|W_tau_emuVec", "prodGenInfo:WtauprongsVec|W_tau_prongsVec", "prodGenInfo:WtaunuVec|W_tau_nuVec"])
   process.stopTreeMaker.vectorFloat.extend([cms.InputTag("prodGenInfo", "WemupfActivityVec"), cms.InputTag("prodGenInfo", "WtauemupfActivityVec"), cms.InputTag("prodGenInfo", "WtauprongspfActivityVec")])
   process.stopTreeMaker.vectorFloatNamesInTree.extend(["prodGenInfo:WemupfActivityVec|W_emu_pfActivityVec", "prodGenInfo:WtauemupfActivityVec|W_tau_emu_pfActivityVec", "prodGenInfo:WtauprongspfActivityVec|W_tau_prongs_pfActivityVec"])
   process.stopTreeMaker.vectorTLorentzVector.append(cms.InputTag("prodGenInfo", "genDecayLVec"))
   process.stopTreeMaker.vectorTLorentzVector.append(cms.InputTag("prodGenInfo", "selGenParticle"))

   #isrJets
   process.ISRJetProducer.debug = cms.bool(options.debug)
   process.stopTreeMaker.varsInt.append(cms.InputTag("ISRJetProducer", "NJetsISR"))

   process.genHT = cms.EDProducer('GenHTProducer')
   #process.stopTreeMaker.varsFloat.append(cms.InputTag("genHT", "genHT"))
   process.stopTreeMaker.vectorTLorentzVector.append(cms.InputTag("prodGenJets", "genjetsLVec"))
   process.stopTreeMaker.varsFloat.extend([cms.InputTag("prodMET:genmet"), cms.InputTag("prodMET:genmetphi")])

   #process.PDFWeights = cms.EDProducer('PDFWeightProducer')
   #process.stopTreeMaker.varsFloat.append(cms.InputTag("PDFWeights", "x1"))
   #process.stopTreeMaker.varsFloat.append(cms.InputTag("PDFWeights", "x2"))
   #process.stopTreeMaker.varsFloat.append(cms.InputTag("PDFWeights", "q"))
   #process.stopTreeMaker.varsInt.append(cms.InputTag("PDFWeights", "id1"))
   #process.stopTreeMaker.varsInt.append(cms.InputTag("PDFWeights", "id2"))
   #process.stopTreeMaker.vectorFloat.append(cms.InputTag("PDFWeights", "ScaleWeightsMiniAOD"))
   '''
   if options.doPDFs == True:
      process.stopTreeMaker.vectorFloat.append(cms.InputTag("PDFWeights", "PDFweights"))
      process.stopTreeMaker.vectorInt.append(cms.InputTag("PDFWeights", "PDFids"))
      process.stopTreeMaker.vectorFloat.append(cms.InputTag("PDFWeights", "PDFweightsMiniAOD"))
      process.stopTreeMaker.vectorInt.append(cms.InputTag("PDFWeights", "PDFidsMiniAOD"))
   '''
   if options.fastsim == True:
      process.load("StopTupleMaker.SkimsAUX.susyscan_cfi")
      process.SusyScanProducer.shouldScan = cms.bool(options.fastsim)
      process.stopTreeMaker.varsFloat.extend([cms.InputTag("SusyScanProducer", "SusyMotherMass"), cms.InputTag("SusyScanProducer", "SusyLSPMass")])

      process.prodJets.jetOtherSrc = cms.InputTag('slimmedJets')

#process.stopTreeMaker.vectorTLorentzVector.append(cms.InputTag("prodIsoTrks:trksForIsoVetoLVec"))

process.stopTreeMaker.vectorFloat.extend([cms.InputTag("prodIsoTrks:trksForIsoVetocharge"), cms.InputTag("prodIsoTrks:trksForIsoVetodz"), cms.InputTag("prodIsoTrks:trksForIsoVetoiso"), cms.InputTag("prodIsoTrks:trksForIsoVetopfActivity"), cms.InputTag("prodIsoTrks:looseisoTrkscharge"), cms.InputTag("prodIsoTrks:looseisoTrksdz"), cms.InputTag("prodIsoTrks:looseisoTrksiso"), cms.InputTag("prodIsoTrks:looseisoTrksmtw"), cms.InputTag("prodIsoTrks:looseisoTrkspfActivity")])
process.stopTreeMaker.vectorFloatNamesInTree.extend([  "prodIsoTrks:looseisoTrkscharge|loose_isoTrks_charge", "prodIsoTrks:looseisoTrksdz|loose_isoTrks_dz", "prodIsoTrks:looseisoTrksiso|loose_isoTrks_iso", "prodIsoTrks:looseisoTrksmtw|loose_isoTrks_mtw", "prodIsoTrks:looseisoTrkspfActivity|loose_isoTrks_pfActivity"])

process.stopTreeMaker.vectorInt.extend([cms.InputTag("prodIsoTrks:trksForIsoVetopdgId"), cms.InputTag("prodIsoTrks:trksForIsoVetoidx"), cms.InputTag("prodIsoTrks:looseisoTrkspdgId"), cms.InputTag("prodIsoTrks:looseisoTrksidx"), cms.InputTag("prodIsoTrks:forVetoIsoTrksidx")])
process.stopTreeMaker.vectorIntNamesInTree.extend([ "prodIsoTrks:looseisoTrkspdgId|loose_isoTrks_pdgId", "prodIsoTrks:looseisoTrksidx|loose_isoTrks_idx"])

process.stopTreeMaker.vectorTLorentzVector.append(cms.InputTag("prodIsoTrks:looseisoTrksLVec"))
process.stopTreeMaker.vectorTLorentzVectorNamesInTree.append("prodIsoTrks:looseisoTrksLVec|loose_isoTrksLVec")

process.stopTreeMaker.varsInt.extend([cms.InputTag("prodIsoTrks:loosenIsoTrks"), cms.InputTag("prodIsoTrks:nIsoTrksForVeto")])
process.stopTreeMaker.varsIntNamesInTree.extend(["prodIsoTrks:loosenIsoTrks|loose_nIsoTrks", "prodIsoTrks:nIsoTrksForVeto|nIsoTrks_CUT"])

process.stopTreeMaker.varsFloat.extend([cms.InputTag("ak4stopmhtPFchs:mht"), cms.InputTag("ak4stopmhtPFchs:mhtphi")])

#process.stopTreeMaker.varsDouble.append(cms.InputTag("mt2PFchs:mt2"))

process.stopTreeMaker.varsFloat.append(cms.InputTag("ak4stophtPFchs"))
process.stopTreeMaker.varsFloatNamesInTree.append("ak4stophtPFchs|ht")

process.stopTreeMaker.varsInt.append(cms.InputTag("ak4nJetsForSkimsStop:nJets"))
process.stopTreeMaker.varsIntNamesInTree.append("ak4nJetsForSkimsStop:nJets|nJets_CUT")

#process.stopTreeMaker.varsDouble.extend([cms.InputTag("prodMET:met"), cms.InputTag("prodMET:metphi")])

if "BADMUON" in options.specialFix and options.mcInfo == False:
 #if options.cmsswVersion == "80X":
#   # Note that this default met from prodMET is both e/gamma and muon corrected which is the recommended one
#   process.stopTreeMaker.varsDouble.extend([cms.InputTag("prodMET:met"), cms.InputTag("prodMET:metphi")])
   # Note that calo MET is stored only in the slimmedMETs collection, therefore addcalomet is True only in the slimmedMETsDefault (with slimmedMETs as source)
   process.stopTreeMaker.varsFloat.extend([cms.InputTag("prodMETslimmedMETsDefault:calomet"), cms.InputTag("prodMETslimmedMETsDefault:calometphi")])
#   process.stopTreeMaker.vectorDouble.extend([cms.InputTag("prodMET:metMagUp"), cms.InputTag("prodMET:metMagDown"), cms.InputTag("prodMET:metPhiUp"), cms.InputTag("prodMET:metPhiDown")])

   # Store other different met
   process.stopTreeMaker.varsFloat.extend([cms.InputTag("prodMETslimmedMETsDefault:met"), cms.InputTag("prodMETslimmedMETsDefault:metphi")])
   process.stopTreeMaker.varsFloatNamesInTree.extend(["prodMETslimmedMETsDefault:met|metMuCleanOnly", "prodMETslimmedMETsDefault:metphi|metphiMuCleanOnly"])

   process.stopTreeMaker.varsFloat.extend([cms.InputTag("prodMETslimmedMETsUncorrected:met"), cms.InputTag("prodMETslimmedMETsUncorrected:metphi")])
   process.stopTreeMaker.varsFloatNamesInTree.extend(["prodMETslimmedMETsUncorrected:met|metNoClean", "prodMETslimmedMETsUncorrected:metphi|metphiNoClean"])

   process.stopTreeMaker.varsFloat.extend([cms.InputTag("prodMETslimmedMETsEGClean:met"), cms.InputTag("prodMETslimmedMETsEGClean:metphi")])
   process.stopTreeMaker.varsFloatNamesInTree.extend(["prodMETslimmedMETsEGClean:met|metEGCleanOnly", "prodMETslimmedMETsEGClean:metphi|metphiEGCleanOnly"])
else:
   process.stopTreeMaker.varsFloat.extend([cms.InputTag("prodMET:met"), cms.InputTag("prodMET:metphi")])
   process.stopTreeMaker.varsFloat.extend([cms.InputTag("prodMET:calomet"), cms.InputTag("prodMET:calometphi")])
   process.stopTreeMaker.vectorFloat.extend([cms.InputTag("prodMET:metMagUp"), cms.InputTag("prodMET:metMagDown"), cms.InputTag("prodMET:metPhiUp"), cms.InputTag("prodMET:metPhiDown")])

process.stopTreeMaker.varsFloat.extend([cms.InputTag("ak4jetMHTDPhiForSkimsStop:dPhi0"), cms.InputTag("ak4jetMHTDPhiForSkimsStop:dPhi1"), cms.InputTag("ak4jetMHTDPhiForSkimsStop:dPhi2")])
process.stopTreeMaker.varsFloatNamesInTree.extend(["ak4jetMHTDPhiForSkimsStop:dPhi0|dPhi0_CUT", "ak4jetMHTDPhiForSkimsStop:dPhi1|dPhi1_CUT", "ak4jetMHTDPhiForSkimsStop:dPhi2|dPhi2_CUT"])

process.stopTreeMaker.varsInt.extend([cms.InputTag("prodEventInfo:vtxSize"), cms.InputTag("prodEventInfo:npv"), cms.InputTag("prodEventInfo:nm1"), cms.InputTag("prodEventInfo:n0"), cms.InputTag("prodEventInfo:np1")])
process.stopTreeMaker.varsFloat.extend([cms.InputTag("prodEventInfo:trunpv"), cms.InputTag("prodEventInfo:avgnpv"), cms.InputTag("prodEventInfo:storedWeight")])
process.stopTreeMaker.varsFloatNamesInTree.extend(["prodEventInfo:trunpv|tru_npv", "prodEventInfo:avgnpv|avg_npv", "prodEventInfo:storedWeight|stored_weight"])


process.stopTreeMaker.varsFloat.append(cms.InputTag("weightProducer:weight"))
process.stopTreeMaker.varsFloatNamesInTree.append("weightProducer:weight|evtWeight")

if options.fastsim == False:
   process.trig_filter_task = cms.Task( process.HBHENoiseFilterResultProducer, process.triggerProducer, process.CSCTightHaloFilter, process.globalSuperTightHalo2016Filter, process.goodVerticesFilter, process.ecalBadCalibFilter, process.EcalDeadCellTriggerPrimitiveFilter, process.BadChargedCandidateFilter, process.BadPFMuonFilter, process.HBHENoiseIsoFilter ) 
   process.trig_filter_seq = cms.Sequence(process.trig_filter_task)
else:
   process.trig_filter_task = cms.Task( process.HBHENoiseIsoFilter, process.triggerProducer, process.CSCTightHaloFilter, process.globalSuperTightHalo2016Filter, process.goodVerticesFilter, process.ecalBadCalibFilter, process.EcalDeadCellTriggerPrimitiveFilter, process.BadChargedCandidateFilter, process.BadPFMuonFilter ) 
   process.trig_filter_seq = cms.Sequence(process.trig_filter_task)

if options.externalFilterList:
   process.load("StopTupleMaker.SkimsAUX.EventListFilter_cfi")
   for flist in options.externalFilterList:
      tfile = tarfile.open(flist, 'r:gz')
      tfile.extractall('.')
      flist = flist.replace(".txt.tar.gz", ".txt")
      if flist.find("ecalsc") != -1:
         process.eeBadScListFilter = process.EventListFilter.clone(inputFileList=flist)
         process.trig_filter_seq += process.eeBadScListFilter
         process.stopTreeMaker.varsBool.append(cms.InputTag("eeBadScListFilter"))
      elif flist.find("csc2015") != -1:
         process.CSCTightHaloListFilter = process.EventListFilter.clone(inputFileList=flist)
         process.trig_filter_seq += process.CSCTightHaloListFilter
         process.stopTreeMaker.varsBool.append(cms.InputTag("CSCTightHaloListFilter"))
      elif flist.find("badResolutionTrack") !=-1:
         process.badResolutionTrackListFilter = process.EventListFilter.clone(inputFileList=flist)
         process.trig_filter_seq += process.badResolutionTrackListFilter
         process.stopTreeMaker.varsBool.append(cms.InputTag("badResolutionTrackListFilter"))
      elif flist.find("muonBadTrack") != -1:
         process.muonBadTrackListFilter = process.EventListFilter.clone(inputFileList=flist)
         process.trig_filter_seq += process.muonBadTrackListFilter
         process.stopTreeMaker.varsBool.append(cms.InputTag("muonBadTrackListFilter"))
      else:
         print "Do NOT support externalFilterList with name : ", flist

if options.selSMSpts == True:
   process.stopTreeMaker.vectorString.extend([cms.InputTag("smsModelFilter:fileNameStr"), cms.InputTag("smsModelFilter:smsModelStr")])
   process.stopTreeMaker.varsFloat.extend([cms.InputTag("smsModelFilter:smsMotherMass"), cms.InputTag("smsModelFilter:smsDaughterMass")])

#process.ak4Stop_Path.associate(process.myTask)


if options.mcInfo == False:

	process.comb_task = cms.Task(   process.cleanpatseq_task, process.prodMuons, process.egmGsfElectronIDTask, process.prodElectrons, process.egmPhotonIDTask, process.goodPhotons, process.QGTagger, process.QGTaggerNoLep, process.weightProducer, process.trackIsolation, process.loosetrackIsolation, process.prodIsoTrks, process.stopBJets,  process.prepareCutvars_task#, process.genHT
) #process.hltFilte process.QGAK4PFCHSr process.stopPFJets

else:
	process.comb_task = cms.Task(   process.cleanpatseq_task, process.prodMuons, process.egmGsfElectronIDTask, process.prodElectrons, process.egmPhotonIDTask, process.goodPhotons, process.QGTagger, process.QGTaggerNoLep, process.weightProducer, process.trackIsolation, process.loosetrackIsolation, process.prodIsoTrks, process.stopBJets, process.prepareCutvars_task, process.genHT, process.ISRJetProducer, process.prodGenJets #process.PDFWeights
)

# Other sequence
process.comb_seq = cms.Sequence(
  # All cleaning && all basic variables, e.g., mht, ht...     
  process.comb_task 
  # hlt requirement
  #process.QGAK4PFCHS, process.stopPFJets
)

process.dump=cms.EDAnalyzer('EventContentAnalyzer')

process.ak4Stop_Path = cms.Path( process.fullPatMetSequence * 
                                 process.comb_seq * 
                                 process.printDecayPythia8 * process.prodGenInfo * process.prodGoodVertices * 
                                 process.prodMuonsNoIso * process.prodElectronsNoIso * process.prodIsoTrks * process.prodJetIDEventFilter *
                                 process.prodJetIDEventFilterNoLep * process.METFilters *
                                 process.noBadMuonsFilter * process.badMuonsFilter * process.duplicateMuonsFilter * process.prodJetsNoLep *
                                 process.prodJets * process.prodMET * process.prodEventInfo * process.trig_filter_seq * process.prodSecondaryVertex *
                                 #process.prodBTag *
                                 #process.goodPhotons *
                                 #process.genHT *
                                 process.stopTreeMaker
                                 * process.dump
)


#process.myTask = cms.Task()
#process.myTask.add(*[getattr(process,prod) for prod in process.producers_()])
#process.myTask.add(*[getattr(process,filt) for filt in process.filters_()])
#process.ak4Stop_Path = cms.Task()
#process.ak4Stop_Path.add(*[getattr(process,prod) for prod in process.producers_()])
#process.ak4Stop_Path.add(*[getattr(process,filt) for filt in process.filters_()])
#process.ak4Stop_Path.associate(process.myTask)   
