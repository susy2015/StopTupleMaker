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

options.register('ntpVersion', "Ntp_BLAH", VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.string, "ntpVersion: to be same as the tag of the release. But can be used to produce 72X ntuple as well!")
options.register('GlobalTag', "94X_mc2017_realistic_v12", VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.string, "74X PromptReco: 74X_dataRun2_Prompt_v0")
options.register('jecDBname', "Fall17_17Nov2017_V8_MC", VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.string, "Summer15_25nsV6_DATA for data")
options.register('hltName', 'HLT', VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.string, "HLT menu to use for trigger matching")

options.register('mcInfo', True, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.bool, "process MonteCarlo data, default is data")

options.register('doPDFs', False, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.bool, "switch to enable the production of PDF weights for NNPDF3.0, CT10, MMHT2014, n.b. you need to do `scram setup lhapdf` before running (default=False)")

options.register('hltSelection', '*', VarParsing.VarParsing.multiplicity.list, VarParsing.VarParsing.varType.string, "hlTriggers (OR) used to filter events. for data: ''HLT_Mu9', 'HLT_IsoMu9', 'HLT_IsoMu13_v*''; for MC, HLT_Mu9")

options.register('debug', False, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.bool, "switch on/off debug mode")
options.register('verbose', 0, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.int, "verbose of debug")
options.register('endReport', False, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.bool, "Show report")

options.register('doPtHatWeighting', False, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.bool, "PtHat weighting for QCD flat samples, default is False")

options.register('fileslist', '', VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.string, "name of a file with input source list")

options.register('fastsim', False, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.bool, "fastsim sample or not, default is False")

options.register('selSMSpts', False, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.bool, "select model pobools")

options.parseArguments()
options._tagOrder =[]

print options

## -------------------------
## -- Check CMSSW version --
## -------------------------

procCMSSWver = os.environ['CMSSW_RELEASE_BASE'].split("/")[-1]
print "procCMSSWver : ", procCMSSWver, "\n"


## ------------------------
## -- Define the process --
## ------------------------

process = cms.Process("SUSY")

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
process.options  = cms.untracked.PSet( )
if options.endReport:
   process.options.wantSummary = cms.untracked.bool(True) 

## -- Allow unscheduled --
process.options.allowUnscheduled = cms.untracked.bool(True)

## -- Maximal Number of Events --
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )

if options.debug and options.verbose ==1:
   process.SimpleMemoryCheck = cms.Service('SimpleMemoryCheck',
                                        ignoreTotal=cms.untracked.int32(0),
                                        oncePerEventMode = cms.untracked.bool(False)
                                        )
   process.Timing = cms.Service("Timing")

## -- Input Source --
process.source = cms.Source("PoolSource",
    #Do not add files here, instead add them to "process.source.fileNames" below or using the option "options.files"
    fileNames = cms.untracked.vstring("" )
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

## -- Output source -- 

process.TFileService = cms.Service("TFileService",
   fileName = cms.string('stopFlatNtuples.root')
)

## ---------------------
## -- Calibration tag --
## ---------------------

from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag

if options.GlobalTag:
   process.GlobalTag = GlobalTag(process.GlobalTag, options.GlobalTag, '')


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
)
process.patJetsAK4PFCHSNoLep.userData.userFloats.src += ['QGTaggerNoLep:qgLikelihood','QGTaggerNoLep:ptD', 'QGTaggerNoLep:axis2', 'QGTaggerNoLep:axis1']
process.patJetsAK4PFCHSNoLep.userData.userInts.src += ['QGTaggerNoLep:mult']
process.patJetsAK4PFCHSNoLep.addGenPartonMatch = cms.bool(False)
process.patJetsAK4PFCHSNoLep.addGenJetMatch = cms.bool(False)

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


# Standard Event cleaning 
process.load("StopTupleMaker.SkimsAUX.prodFilterOutScraping_cfi")
process.load("StopTupleMaker.SkimsAUX.prodGoodVertices_cfi")
process.load("StopTupleMaker.SkimsAUX.prodSecondaryVertex_cfi")

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

process.load("StopTupleMaker.SkimsAUX.prodGenJets_cfi")
process.load("StopTupleMaker.SkimsAUX.prodGenInfo_cfi")
process.load("StopTupleMaker.SkimsAUX.prodIsoTrks_cfi")
process.load("StopTupleMaker.SkimsAUX.prodEventInfo_cfi")
process.load("StopTupleMaker.SkimsAUX.PhotonIDisoProducer_cfi")

process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

process.load("StopTupleMaker.SkimsAUX.genDecayStringMakerPythia8_cfi")
process.printDecayPythia8.src = cms.InputTag("prunedGenParticles")
process.printDecayPythia8.keyDecayStrs = cms.vstring("t", "tbar", "~chi_1+", "~chi_1-")
process.printDecayPythia8.printDecay = cms.untracked.bool(options.debug)


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
   inputStr = cms.vstring(options.ntpVersion, options.GlobalTag, options.hltName, options.jecDBname, procCMSSWver)
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

process.stopTreeMaker.vectorTLorentzVector.append(cms.InputTag("prodJets", "jetsLVec"))

process.stopTreeMaker.vectorInt.append(cms.InputTag("prodJets", "recoJetsFlavor"))

process.stopTreeMaker.vectorFloat.append(cms.InputTag("prodJets", "recoJetsJecUnc"))
process.stopTreeMaker.vectorFloat.append(cms.InputTag("prodJets", "recoJetsJecScaleRawToFull"))

process.stopTreeMaker.vectorFloat.append(cms.InputTag("prodJets", "qgLikelihood"))
process.stopTreeMaker.vectorFloat.append(cms.InputTag("prodJets", "qgPtD"))
process.stopTreeMaker.vectorFloat.append(cms.InputTag("prodJets", "qgAxis2"))
process.stopTreeMaker.vectorFloat.append(cms.InputTag("prodJets", "qgMult"))
process.stopTreeMaker.vectorFloat.append(cms.InputTag("prodJets", "qgAxis1"))

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
        
process.stopTreeMaker.vectorInt.append(cms.InputTag("prodJetsNoLep", "recoJetsFlavor"))

process.stopTreeMaker.vectorFloat.append(cms.InputTag("prodJetsNoLep", "DeepCSVb"))
process.stopTreeMaker.vectorFloat.append(cms.InputTag("prodJetsNoLep", "DeepCSVc"))
process.stopTreeMaker.vectorFloat.append(cms.InputTag("prodJetsNoLep", "DeepCSVl"))
process.stopTreeMaker.vectorFloat.append(cms.InputTag("prodJetsNoLep", "DeepCSVbb"))
process.stopTreeMaker.vectorFloat.append(cms.InputTag("prodJetsNoLep", "DeepCSVcc"))
   
    
process.stopTreeMaker.vectorFloat.append(cms.InputTag("prodJetsNoLep", "CversusL"));
process.stopTreeMaker.vectorFloat.append(cms.InputTag("prodJetsNoLep", "CversusB"));

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
   process.stopTreeMaker.varsFloat.append(cms.InputTag("genHT", "genHT"))

   process.stopTreeMaker.vectorTLorentzVector.append(cms.InputTag("prodGenJets", "genjetsLVec"))
   process.stopTreeMaker.varsFloat.extend([cms.InputTag("prodMET:genmet"), cms.InputTag("prodMET:genmetphi")])

   #NEEDS TO BE FIXED!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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


process.stopTreeMaker.varsFloat.extend([cms.InputTag("prodMET:met"), cms.InputTag("prodMET:metphi")])
process.stopTreeMaker.varsFloat.extend([cms.InputTag("prodMET:calomet"), cms.InputTag("prodMET:calometphi")])
process.stopTreeMaker.vectorFloat.extend([cms.InputTag("prodMET:metMagUp"), cms.InputTag("prodMET:metMagDown"), cms.InputTag("prodMET:metPhiUp"), cms.InputTag("prodMET:metPhiDown")])

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

if options.selSMSpts == True:
   process.stopTreeMaker.vectorString.extend([cms.InputTag("smsModelFilter:fileNameStr"), cms.InputTag("smsModelFilter:smsModelStr")])
   process.stopTreeMaker.varsFloat.extend([cms.InputTag("smsModelFilter:smsMotherMass"), cms.InputTag("smsModelFilter:smsDaughterMass")])

#process.ak4Stop_Path.associate(process.myTask)


if options.mcInfo == False:
   process.comb_task = cms.Task( process.prodMuons, process.egmGsfElectronIDTask, process.prodElectrons, process.egmPhotonIDTask, process.goodPhotons, process.QGTagger, process.QGTaggerNoLep, process.weightProducer, process.prodIsoTrks, 
                                        )
else:
   process.comb_task = cms.Task( process.prodMuons, process.egmGsfElectronIDTask, process.prodElectrons, process.egmPhotonIDTask, process.goodPhotons, process.QGTagger, process.QGTaggerNoLep, process.weightProducer, process.prodIsoTrks, process.genHT, process.ISRJetProducer, process.prodGenJets
                                        )

# Other sequence
process.comb_seq = cms.Sequence(
  # All cleaning && all basic variables, e.g., mht, ht...     
  process.comb_task 
)

process.dump=cms.EDAnalyzer('EventContentAnalyzer')

process.ak4Stop_Path = cms.Path( process.fullPatMetSequence * 
                                 process.comb_seq * 
                                 process.printDecayPythia8 * process.prodGenInfo * process.prodGoodVertices * 
                                 process.prodMuonsNoIso * process.prodElectronsNoIso * process.prodIsoTrks * process.prodJetIDEventFilter *
                                 process.prodJetIDEventFilterNoLep * process.METFilters *
                                 process.noBadMuonsFilter * process.badMuonsFilter * process.duplicateMuonsFilter * process.prodJetsNoLep *
                                 process.prodJets * process.prodMET * process.prodEventInfo * process.trig_filter_seq * process.prodSecondaryVertex *
                                 process.stopTreeMaker
#                                 * process.dump
)

