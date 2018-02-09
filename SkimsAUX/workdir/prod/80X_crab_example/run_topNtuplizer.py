import FWCore.ParameterSet.Config as cms

#Kevin's function of magic #1
def addJetInfo(process, JetTag, userFloats=[], userInts=[], btagDiscrs=cms.VInputTag(), suff=""):
    # add userfloats to jet collection
    from PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cff import updatedPatJets as patJetsUpdated

    # default suffix
    if len(suff)==0: suff = "Auxiliary"
    
    JetTagOut = cms.InputTag(JetTag.value()+suff)
    patJetsAuxiliary = patJetsUpdated.clone(
        jetSource = JetTag,
        addJetCorrFactors = cms.bool(False),
        addBTagInfo = cms.bool(False)
    )
    patJetsAuxiliary.userData.userFloats.src += userFloats
    patJetsAuxiliary.userData.userInts.src += userInts
    if len(btagDiscrs)>0:
        patJetsAuxiliary.discriminatorSources = btagDiscrs
        patJetsAuxiliary.addBTagInfo = cms.bool(True)
    setattr(process,JetTagOut.value(),patJetsAuxiliary)
    
    return (process, JetTagOut)

import FWCore.ParameterSet.VarParsing as VarParsing
### parsing job options 
import sys

options = VarParsing.VarParsing()

options.register('inputScript','',VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.string,"input Script")
options.register('outputFile','output',VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.string,"output File (w/o .root)")
options.register('maxEvents',-1,VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.int,"maximum events")
options.register('skipEvents', 0, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.int, "skip N events")
options.register('job', 0, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.int, "job number")
options.register('nJobs', 1, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.int, "total jobs")
options.register('release','8_0_1', VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.string,"release number (w/o CMSSW)")
options.register('isData', 0, VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.int,"isData flag (0 for MC, 1 for data)")

options.parseArguments()

print("Using release "+options.release)


#if hasattr(sys, "argv"):
#    options.parseArguments()


process = cms.Process("DNNFiller")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.EventContent.EventContent_cff")
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag

print "isData: ", options.isData
if options.isData:
    print "Running on Data"
    process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_data', '')
else:
    print "Running on MC"
    process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')    

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )

process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.options = cms.untracked.PSet(
   allowUnscheduled = cms.untracked.bool(True),  
   wantSummary=cms.untracked.bool(False)
)

from PhysicsTools.PatAlgos.patInputFiles_cff import filesRelValTTbarPileUpMINIAODSIM

###############################################################################################################################

process.source = cms.Source('PoolSource',
#    fileNames=cms.untracked.vstring (["/store/mc/RunIISummer16MiniAODv2/TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/110000/423685A0-BFE6-E611-B2B5-001E67DBE36D.root"]),
                            fileNames=cms.untracked.vstring (["file:/uscms_data/d3/pastika/zinv/dev/CMSSW_8_0_26_patch1/src/SusyAnaTools/SkimsAUX/workdir/prod/80X_crab_example/423685A0-BFE6-E611-B2B5-001E67DBE36D.root"]),
#            '/store/data/Run2016C/SingleMuon/MINIAOD/23Sep2016-v1/90000/109B7DBF-0C91-E611-A5EC-0CC47A4D7690.root',]),
)

#if options.inputScript != '' and options.inputScript != 'DeepNTuples.DeepNtuplizer.samples.TEST':
#    process.load(options.inputScript)

numberOfFiles = len(process.source.fileNames)
numberOfJobs = options.nJobs
jobNumber = options.job

process.source.fileNames = process.source.fileNames[jobNumber:numberOfFiles:numberOfJobs]
if options.nJobs > 1:
    print ("running over these files:")
    print (process.source.fileNames)
#process.source.fileNames = ['file:/uscms/home/verzetti/nobackup/CMSSW_8_0_25/src/DeepNTuples/copy_numEvent100.root']

process.source.skipEvents = cms.untracked.uint32(options.skipEvents)
process.maxEvents  = cms.untracked.PSet( 
    input = cms.untracked.int32 (-1) 
)

outFileName = "stopFlatNtuples.root"#options.outputFile + '_' + str(options.job) +  '.root'
print ('Using output file ' + outFileName)

process.TFileService = cms.Service("TFileService", 
                                   fileName = cms.string(outFileName))


###############################################################################################################################

process.load("StopTupleMaker.SkimsAUX.genDecayStringMakerPythia8_cfi")
process.printDecayPythia8.src = cms.InputTag("prunedGenParticles")

###############################################################################################################################

process.load("StopTupleMaker.SkimsAUX.prodFilterOutScraping_cfi")
process.load("StopTupleMaker.SkimsAUX.prodGoodVertices_cfi")

# Default is dR = 0.3, dz < 0.05, pt > 10, reliso < 0.1
process.load("StopTupleMaker.Skims.trackIsolationMaker_cfi")
process.trackIsolation = process.trackIsolationFilter.clone()
process.trackIsolation.pfCandidatesTag = cms.InputTag("packedPFCandidates")
process.trackIsolation.doTrkIsoVeto = cms.bool(False)

process.loosetrackIsolation = process.trackIsolation.clone()
process.loosetrackIsolation.minPt_PFCandidate = cms.double(5.0)
process.loosetrackIsolation.isoCut            = cms.double(0.5)

process.refalltrackIsolation = process.trackIsolation.clone()
process.refalltrackIsolation.mintPt_PFCandidate = cms.double (-1.0)
process.refalltrackIsolation.isoCut           = cms.double(9999.0)

process.load("StopTupleMaker.SkimsAUX.prodMuons_cfi")
process.prodMuonsNoIso = process.prodMuons.clone()
process.prodMuonsNoIso.DoMuonIsolation = cms.int32(0)


process.load("StopTupleMaker.SkimsAUX.prodElectrons_cfi")
process.prodElectronsNoIso = process.prodElectrons.clone()
process.prodElectronsNoIso.DoElectronIsolation = cms.int32(0)

process.load("StopTupleMaker.SkimsAUX.prodIsoTrks_cfi")

###############################################################################################################################

###### -- Add AK8 PUPPI jet collection using Jet Toolbox --
####
####from JMEAnalysis.JetToolbox.jetToolbox_cff import jetToolbox
####
##### Keep this behind the cleaned version for now, otherwise everything will be lepton cleaned
####jetToolbox( process, 'ak8', 'ak8JetSubs', 'out', 
####            runOnMC = not options.isData, 
####            PUMethod='Puppi', 
####            addSoftDropSubjets = True, 
####            addSoftDrop = True, 
####            addNsub = True, 
####            bTagDiscriminators = ['pfCombinedInclusiveSecondaryVertexV2BJetTags'], 
####            addCMSTopTagger = False)

###############################################################################################################################

#ANDRES Gamma Var  
from PhysicsTools.SelectorUtils.tools.vid_id_tools import *

switchOnVIDPhotonIdProducer(process, DataFormat.MiniAOD)

# Define which IDs we want to produce
my_photon_id_modules = ['RecoEgamma.PhotonIdentification.Identification.cutBasedPhotonID_Spring16_V2p2_cff']

# Add them to the VID producer
#if process.PhotonIDisoProducer.isFilled:
for idmod in my_photon_id_modules:
   setupAllVIDIdsInModule(process, idmod, setupVIDPhotonSelection)

process.load("StopTupleMaker.SkimsAUX.PhotonIDisoProducer_cfi")

# Set ID tags
process.goodPhotons.loosePhotonID = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-Spring16-V2p2-loose")
process.goodPhotons.mediumPhotonID = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-Spring16-V2p2-medium")
process.goodPhotons.tightPhotonID = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-Spring16-V2p2-tight")
            

###############################################################################################################################

from RecoJets.JetProducers.ak4PFJets_cfi import ak4PFJets
from RecoJets.JetProducers.ak4GenJets_cfi import ak4GenJets

## Filter out neutrinos from packed GenParticles
process.packedGenParticlesForJetsNoNu = cms.EDFilter("CandPtrSelector", 
                                                     src = cms.InputTag("prunedGenParticles"), 
                                                     cut = cms.string("abs(pdgId) != 12 && abs(pdgId) != 14 && abs(pdgId) != 16"))

## Define GenJets
process.ak4GenJetsNoNu = ak4GenJets.clone(src = 'packedGenParticlesForJetsNoNu')

## -- do projections --
process.pfCHS = cms.EDFilter("CandPtrSelector", 
                             src = cms.InputTag("packedPFCandidates"), 
                             cut = cms.string("fromPV"))

process.pfNoMuonCHSNoMu =  cms.EDProducer("CandPtrProjector", 
                                          src = cms.InputTag("pfCHS"), 
                                          veto = cms.InputTag("prodMuons", "mu2Clean"))
process.pfNoElectronCHSNoEle = cms.EDProducer("CandPtrProjector", 
                                              src = cms.InputTag("pfNoMuonCHSNoMu"), 
                                              veto = cms.InputTag("prodElectrons", "ele2Clean"))
process.ak4PFJetsCHSNoLep = ak4PFJets.clone(src = 'pfNoElectronCHSNoEle', doAreaFastjet = True) # no idea while doArea is false by default, but it's True in RECO so we have to set it

###############################################################################################################################

#if False: #int(options.release.replace("_",""))>=840 :
# bTagInfos = [
#        'pfImpactParameterTagInfos',
#        'pfInclusiveSecondaryVertexFinderTagInfos',
#        'pfDeepCSVTagInfos',
# ]
#else :
# bTagInfos = [
#        'pfImpactParameterTagInfos',
#        'pfInclusiveSecondaryVertexFinderTagInfos',
#        'deepNNTagInfos',
# ]


#if True: #int(options.release.replace("_",""))>=840 :
# bTagDiscriminators = [
#     'softPFMuonBJetTags',
#     'softPFElectronBJetTags',
#     'pfJetBProbabilityBJetTags',
#     'pfJetProbabilityBJetTags',
#     'pfCombinedInclusiveSecondaryVertexV2BJetTags',
#     'pfCombinedCvsLJetTags',     
#     'pfPositiveCombinedCvsLJetTags',
#     'pfNegativeCombinedCvsLJetTags',
#     'pfDeepCSVJetTags:probudsg', #to be fixed with new names
#     'pfDeepCSVJetTags:probb',
#     'pfDeepCSVJetTags:probc',
#     'pfDeepCSVJetTags:probbb',
#     'pfDeepCSVJetTags:probcc',
#     ]
#else :
#    bTagDiscriminators = [
#        'softPFMuonBJetTags',
#        'softPFElectronBJetTags',
#        'pfJetBProbabilityBJetTags',
#        'pfJetProbabilityBJetTags',
#        'pfCombinedInclusiveSecondaryVertexV2BJetTags',
#        'pfCombinedCvsLJetTags',     
#        'pfPositiveCombinedCvsLJetTags',
#        'pfNegativeCombinedCvsLJetTags',
#        'deepFlavourJetTags:probudsg', #to be fixed with new names
#        'deepFlavourJetTags:probb',
#        'deepFlavourJetTags:probc',
#        'deepFlavourJetTags:probbb',
#        'deepFlavourJetTags:probcc',
# ]

jetCorrectionsAK4 = ('AK4PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute'], 'None')

#from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection
#updateJetCollection(
#        process,
#        labelName = "DeepFlavour",
#        jetSource = cms.InputTag('ak4PFJetsCHSNoLep'),#'ak4Jets'
#        jetCorrections = jetCorrectionsAK4,
#        pfCandidates = cms.InputTag('pfNoElectronCHSNoEle'),
#        pvSource = cms.InputTag("offlineSlimmedPrimaryVertices"),
#
#if hasattr(process,'updatedPatJetsTransientCorrectedDeepFlavour'):
#        svSource = cms.InputTag('slimmedSecondaryVertices'),
#        muSource = cms.InputTag('slimmedMuons'),
#        elSource = cms.InputTag('slimmedElectrons'),
#        btagInfos = bTagInfos,
#        btagDiscriminators = bTagDiscriminators,
#        explicitJTA = False
#)

from PhysicsTools.PatAlgos.tools.jetTools import addJetCollection

addJetCollection(
    process,
    postfix = "",
    labelName = 'DeepFlavour',
    jetSource = cms.InputTag('ak4PFJetsCHSNoLep'),
    pvSource = cms.InputTag('offlineSlimmedPrimaryVertices'),
    pfCandidates = cms.InputTag('packedPFCandidates'),
    svSource = cms.InputTag('slimmedSecondaryVertices'),
    elSource = cms.InputTag('slimmedElectrons'),
    muSource = cms.InputTag('slimmedMuons'),
    jetCorrections = jetCorrectionsAK4,
#    btagDiscriminators = bTagDiscriminators,
    genJetCollection = cms.InputTag('ak4GenJetsNoNu'),
    genParticles = cms.InputTag('prunedGenParticles'),
    algo = 'AK', rParam = 0.4
    )


#process.patJetsDeepFlavour.trackAssociationSource = cms.InputTag("")
#process.patJetsDeepFlavour.jetChargeSource = cms.InputTag("")
#process.patJetsDeepFlavour.addAssociatedTracks = cms.bool(False)
#process.patJetsDeepFlavour.addGenJetMatch = cms.bool(False)
#process.patJetsDeepFlavour.addGenPartonMatch = cms.bool(False)
#process.patJetsDeepFlavour.addJetCharge = cms.bool(False)
#process.patJetChargeDeepFlavour.src = cms.InputTag("")
#process.patJetsDeepFlavour.userData.userFloats.src += ['QGTagger:qgLikelihood','QGTagger:ptD', 'QGTagger:axis1', 'QGTagger:axis2']
#process.patJetsDeepFlavour.userData.userInts.src += ['QGTagger:mult']

#if hasattr(process,'patJetsTransientCorrectedDeepFlavour'):
#	process.patJetsTransientCorrectedDeepFlavour.addTagInfos = cms.bool(True) 
#	process.patJetsTransientCorrectedDeepFlavour.addBTagInfo = cms.bool(True)
#else:
#	raise ValueError('I could not find updatedPatJetsTransientCorrectedDeepFlavour to embed the tagInfos, please check the cfg')

###############################################################################################################################

# QGLikelihood

qgDatabaseVersion = 'cmssw8020_v2'

#databasepath=os.environ['CMSSW_BASE']+'/src/DeepNTuples/DeepNtuplizer/data/QGL_cmssw8020_v2.db'
databasepath='QGL_cmssw8020_v2.db'

from CondCore.CondDB.CondDB_cfi import *
process.QGPoolDBESSource = cms.ESSource("PoolDBESSource",
      CondDB.DBParameters,
      toGet = cms.VPSet(),
      connect = cms.string('sqlite_file:'+databasepath),
)

for type in ['AK4PFchs','AK4PFchs_antib']:
    process.QGPoolDBESSource.toGet.extend(cms.VPSet(cms.PSet(
                record = cms.string('QGLikelihoodRcd'),
                tag    = cms.string('QGLikelihoodObject_'+qgDatabaseVersion+'_'+type),
                label  = cms.untracked.string('QGL_'+type)
                )))

process.es_prefer_jec = cms.ESPrefer("PoolDBESSource", "QGPoolDBESSource")


###############################################################################################################################

process.load('RecoJets.JetProducers.QGTagger_cfi')
process.QGTagger.srcJets   = cms.InputTag("selectedPatJetsDeepFlavour")
process.QGTagger.jetsLabel = cms.string('QGL_AK4PFchs')

process, jetTag = addJetInfo(process, cms.InputTag("selectedPatJetsDeepFlavour"), userFloats=['QGTagger:qgLikelihood','QGTagger:ptD', 'QGTagger:axis1', 'QGTagger:axis2'], userInts=['QGTagger:mult'], suff="")

###############################################################################################################################

#Deep Flavor

###process.load("RecoBTag.DeepFlavour.DeepFlavourTagInfos_cfi")
###process.load("RecoBTag.DeepFlavour.DeepFlavourJetTags_cfi")
###process.load("CommonTools.PileupAlgos.Puppi_cff")
###process.load("PhysicsTools.PatAlgos.slimming.primaryVertexAssociation_cfi")
###
###process.pfDeepFlavourTagInfos.jets = jetTag
###process.pfDeepFlavourTagInfos.vertices = cms.InputTag("offlineSlimmedPrimaryVertices")
###process.pfDeepFlavourTagInfos.secondary_vertices = cms.InputTag("slimmedSecondaryVertices")
###
###process.puppi.candName   = cms.InputTag('packedPFCandidates')
###process.puppi.vertexName = cms.InputTag('offlineSlimmedPrimaryVertices')
###
###process.primaryVertexAssociation.particles = cms.InputTag('packedPFCandidates')
###process.primaryVertexAssociation.vertices  = cms.InputTag('offlineSlimmedPrimaryVertices')
###process.primaryVertexAssociation.jets      = jetTag
###
###process, jetTag = addJetInfo(process, jetTag, userFloats=['pfDeepFlavourJetTags:probb', 'pfDeepFlavourJetTags:probbb', 'pfDeepFlavourJetTags:problepb', 'pfDeepFlavourJetTags:probc', 'pfDeepFlavourJetTags:probuds', 'pfDeepFlavourJetTags:probg'], userInts=[], suff="")

from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection

#updateJetCollection(
#        process,
#        labelName = "DeepFlavour",
#        jetSource = cms.InputTag('ak4PFJetsCHSNoLep'),#'ak4Jets'
#        jetCorrections = jetCorrectionsAK4,
#        pfCandidates = cms.InputTag('pfNoElectronCHSNoEle'),
#        pvSource = cms.InputTag("offlineSlimmedPrimaryVertices"),
#
#if hasattr(process,'updatedPatJetsTransientCorrectedDeepFlavour'):
#        svSource = cms.InputTag('slimmedSecondaryVertices'),
#        muSource = cms.InputTag('slimmedMuons'),
#        elSource = cms.InputTag('slimmedElectrons'),
#        btagInfos = bTagInfos,
#        btagDiscriminators = bTagDiscriminators,
#        explicitJTA = False
#)

updateJetCollection(
   process,
   labelName = "DeepFlavour2",
   jetSource = jetTag,
   jetCorrections = ('AK4PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute']), 'None'),
   btagDiscriminators = [
      'softPFMuonBJetTags',
      'softPFElectronBJetTags',
      'pfJetBProbabilityBJetTags',
      'pfJetProbabilityBJetTags',
      'pfCombinedInclusiveSecondaryVertexV2BJetTags',
      'pfCombinedCvsLJetTags',     
      'pfCombinedCvsBJetTags',     
      'pfCombinedSecondaryVertexV2BJetTags',
      'pfDeepCSVJetTags:probudsg', 
      'pfDeepCSVJetTags:probb', 
      'pfDeepCSVJetTags:probc', 
      'pfDeepCSVJetTags:probbb', 
      'pfDeepFlavourJetTags:probb',
      'pfDeepFlavourJetTags:probbb',
      'pfDeepFlavourJetTags:problepb',
      'pfDeepFlavourJetTags:probc',
      'pfDeepFlavourJetTags:probuds',
      'pfDeepFlavourJetTags:probg',
      ] ## to add discriminators
)

jetTag = cms.InputTag('selectedUpdatedPatJetsDeepFlavour2')

###############################################################################################################################

# update the MET to account for the new JECs
from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD
runMetCorAndUncFromMiniAOD(
    process,
    isData=bool(options.isData), # controls gen met
    #jetCollUnskimmed="updatedPatJetsUpdatedJEC",
    #jetColl="updatedPatJetsUpdatedJEC",
    #postfix="Update"
    )

###############################################################################################################################

process.load("StopTupleMaker.SkimsAUX.prodGenInfo_cfi")
process.load("StopTupleMaker.SkimsAUX.prodJets_cfi")
process.prodJets.bTagKeyString = cms.string('pfCombinedInclusiveSecondaryVertexV2BJetTags')
process.prodJets.jetPBJetTags = cms.string('pfJetBProbabilityBJetTags')
process.prodJets.jetPNegBJetTags= cms.string('pfNegativeOnlyJetBProbabilityBJetTags')
process.prodJets.jetPPosBJetTags= cms.string('pfPositiveOnlyJetBProbabilityBJetTags')
process.prodJets.jetBPBJetTags= cms.string('jetBPBJetTags')
process.prodJets.jetBPNegBJetTags= cms.string('jetBPNegBJetTags')
process.prodJets.jetBPPosBJetTags= cms.string('jetBPPosBJetTags')
process.prodJets.debug = cms.bool(False)
process.prodJets.jetSrc = jetTag
process.prodJets. ak4ptCut = cms.double(20.0)
#process.prodJets.jetOtherSrc = cms.InputTag('selectedUpdatedPatJetsDeepFlavour')

###############################################################################################################################

process.load("StopTupleMaker.SkimsAUX.prodMET_cfi")
process.prodMET.addcalomet = cms.bool(True)
process.prodMET.metSrc = cms.InputTag("slimmedMETs", "", process.name_())

###############################################################################################################################

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

process.triggerProducer.trigTagSrc = cms.InputTag("TriggerResults","","HLT")
#process.METFilters = process.filterDecisionProducer.clone( filterName  =   cms.string("Flag_METFilters") )
process.CSCTightHaloFilter = process.filterDecisionProducer.clone( filterName  =   cms.string("Flag_CSCTightHaloFilter") )
process.globalTightHalo2016Filter = process.filterDecisionProducer.clone( filterName  =   cms.string("Flag_globalTightHalo2016Filter") )
process.goodVerticesFilter = process.filterDecisionProducer.clone( filterName  =   cms.string("Flag_goodVertices") )
process.eeBadScFilter = process.filterDecisionProducer.clone( filterName  =   cms.string("Flag_eeBadScFilter") )
#process.HBHENoiseFilter = process.filterDecisionProducer.clone( filterName  =   cms.string("Flag_HBHENoiseFilter") )
process.EcalDeadCellTriggerPrimitiveFilter = process.filterDecisionProducer.clone( filterName  =   cms.string("Flag_EcalDeadCellTriggerPrimitiveFilter") )

process.filterDecisionProducerPAT = process.filterDecisionProducer.clone()
#process.filterDecisionProducerPAT.trigTagSrc = cms.InputTag("TriggerResults","","PAT")
process.filterDecisionProducerPAT.trigTagSrc = cms.InputTag("TriggerResults","")
process.noBadMuonsFilter = process.filterDecisionProducerPAT.clone( filterName  =   cms.string("Flag_noBadMuons") )
process.badMuonsFilter = process.filterDecisionProducerPAT.clone( filterName = cms.string("Flag_badMuons") )
process.duplicateMuonsFilter = process.filterDecisionProducerPAT.clone( filterName = cms.string("Flag_duplicateMuons") )

process.load('StopTupleMaker.SkimsAUX.prodJetIDEventFilter_cfi')
process.prodJetIDEventFilter.JetSource = cms.InputTag("slimmedJets")
process.prodJetIDEventFilter.MinJetPt  = cms.double(30.0)
process.prodJetIDEventFilter.MaxJetEta = cms.double(999.0)


###############################################################################################################################

process.load("StopTupleMaker.SkimsAUX.ISRJetProducer_cfi")
process.ISRJetProducer.cleanJetSrc = cms.InputTag("selectedPatJetsDeepFlavour")
process.ISRJetProducer.genParticleSrc = cms.InputTag("prunedGenParticles")
process.load("StopTupleMaker.SkimsAUX.prodEventInfo_cfi")

###############################################################################################################################

process.load("StopTupleMaker.StopTreeMaker.stopTreeMaker_cfi")
process.stopTreeMaker.debug = cms.bool(False)
process.stopTreeMaker.TreeName = cms.string("AUX")

process.stopTreeMaker.vectorTLorentzVector.append(cms.InputTag("prodJets", "jetsLVec"))
#process.stopTreeMaker.vectorTLorentzVector.append(cms.InputTag("prodGenJets", "genjetsLVec"))
process.stopTreeMaker.vectorInt.append(cms.InputTag("prodJets", "recoJetsFlavor"))
process.stopTreeMaker.vectorDouble.append(cms.InputTag("prodJets", "recoJetsJecUnc"))
process.stopTreeMaker.vectorDouble.append(cms.InputTag("prodJets", "recoJetsJecScaleRawToFull"))
process.stopTreeMaker.vectorDouble.append(cms.InputTag("prodJets", "qgLikelihood"))
process.stopTreeMaker.vectorDouble.append(cms.InputTag("prodJets", "qgPtD"))
process.stopTreeMaker.vectorDouble.append(cms.InputTag("prodJets", "qgAxis1"))
process.stopTreeMaker.vectorDouble.append(cms.InputTag("prodJets", "qgAxis2"))
process.stopTreeMaker.vectorInt.append(cms.InputTag("prodJets", "qgMult"))
process.stopTreeMaker.vectorDouble.append(cms.InputTag("prodJets", "recoJetschargedHadronEnergyFraction"))
process.stopTreeMaker.vectorDouble.append(cms.InputTag("prodJets", "recoJetschargedEmEnergyFraction"))
process.stopTreeMaker.vectorDouble.append(cms.InputTag("prodJets", "recoJetsneutralEmEnergyFraction"))
process.stopTreeMaker.vectorDouble.append(cms.InputTag("prodJets", "recoJetsmuonEnergyFraction"))
process.stopTreeMaker.vectorDouble.append(cms.InputTag("prodJets", "recoJetsHFHadronEnergyFraction"))
process.stopTreeMaker.vectorDouble.append(cms.InputTag("prodJets", "recoJetsHFEMEnergyFraction"));
#process.stopTreeMaker.vectorTLorentzVector.append(cms.InputTag("prodGenJets", "genjetsLVec"))
process.stopTreeMaker.vectorDouble.append(cms.InputTag("prodJets", "recoJetsneutralEnergyFraction"))
process.stopTreeMaker.vectorDouble.append(cms.InputTag("prodJets", "PhotonEnergyFraction"))
process.stopTreeMaker.vectorDouble.append(cms.InputTag("prodJets", "ElectronEnergyFraction"))
process.stopTreeMaker.vectorDouble.append(cms.InputTag("prodJets", "ChargedHadronMultiplicity"))
process.stopTreeMaker.vectorDouble.append(cms.InputTag("prodJets", "NeutralHadronMultiplicity"))
process.stopTreeMaker.vectorDouble.append(cms.InputTag("prodJets", "PhotonMultiplicity"))
process.stopTreeMaker.vectorDouble.append(cms.InputTag("prodJets", "ElectronMultiplicity"))
process.stopTreeMaker.vectorDouble.append(cms.InputTag("prodJets", "MuonMultiplicity"))
process.stopTreeMaker.vectorDouble.append(cms.InputTag("prodJets", "DeepCSVb"))
process.stopTreeMaker.vectorDouble.append(cms.InputTag("prodJets", "DeepCSVc"))
process.stopTreeMaker.vectorDouble.append(cms.InputTag("prodJets", "DeepCSVl"))
process.stopTreeMaker.vectorDouble.append(cms.InputTag("prodJets", "DeepCSVbb"))
process.stopTreeMaker.vectorDouble.append(cms.InputTag("prodJets", "DeepCSVcc"))
process.stopTreeMaker.vectorDouble.append(cms.InputTag("prodJets", "DeepFlavorb"))
process.stopTreeMaker.vectorDouble.append(cms.InputTag("prodJets", "DeepFlavorbb"))
process.stopTreeMaker.vectorDouble.append(cms.InputTag("prodJets", "DeepFlavorlepb"))
process.stopTreeMaker.vectorDouble.append(cms.InputTag("prodJets", "DeepFlavorc"))
process.stopTreeMaker.vectorDouble.append(cms.InputTag("prodJets", "DeepFlavoruds"))
process.stopTreeMaker.vectorDouble.append(cms.InputTag("prodJets", "DeepFlavorg"))
process.stopTreeMaker.vectorDouble.append(cms.InputTag("prodJets", "CvsL"));
process.stopTreeMaker.vectorDouble.append(cms.InputTag("prodJets", "CvsB"));
#process.stopTreeMaker.vectorDouble.append(cms.InputTag("prodJets", "DeepCSVbN"))
#process.stopTreeMaker.vectorDouble.append(cms.InputTag("prodJets", "DeepCSVcN"))
#process.stopTreeMaker.vectorDouble.append(cms.InputTag("prodJets", "DeepCSVlN"))
#process.stopTreeMaker.vectorDouble.append(cms.InputTag("prodJets", "DeepCSVbbN"))
#process.stopTreeMaker.vectorDouble.append(cms.InputTag("prodJets", "DeepCSVccN"))
#process.stopTreeMaker.vectorDouble.append(cms.InputTag("prodJets", "DeepCSVbP"))
#process.stopTreeMaker.vectorDouble.append(cms.InputTag("prodJets", "DeepCSVcP"))
#process.stopTreeMaker.vectorDouble.append(cms.InputTag("prodJets", "DeepCSVlP"))
#process.stopTreeMaker.vectorDouble.append(cms.InputTag("prodJets", "DeepCSVbbP"))
#process.stopTreeMaker.vectorDouble.append(cms.InputTag("prodJets", "DeepCSVccP"))
process.stopTreeMaker.vectorDouble.append(cms.InputTag("prodJets","CombinedSvtx"));
#process.stopTreeMaker.vectorDouble.append(cms.InputTag("prodJets","CombinedSvtxN"));
#process.stopTreeMaker.vectorDouble.append(cms.InputTag("prodJets","CombinedSvtxP"));
#process.stopTreeMaker.vectorDouble.append(cms.InputTag("prodJets","Svtx"));
#process.stopTreeMaker.vectorDouble.append(cms.InputTag("prodJets","SvtxN"));
#process.stopTreeMaker.vectorDouble.append(cms.InputTag("prodJets","SvtxHP"));
#process.stopTreeMaker.vectorDouble.append(cms.InputTag("prodJets","SvtxNHP"));
#process.stopTreeMaker.vectorDouble.append(cms.InputTag("prodJets","SoftM"));
#process.stopTreeMaker.vectorDouble.append(cms.InputTag("prodJets","SoftMN"));
#process.stopTreeMaker.vectorDouble.append(cms.InputTag("prodJets","SoftMP"));
#process.stopTreeMaker.vectorDouble.append(cms.InputTag("prodJets","SoftE"));
#process.stopTreeMaker.vectorDouble.append(cms.InputTag("prodJets","SoftEN"));
#process.stopTreeMaker.vectorDouble.append(cms.InputTag("prodJets","SoftEP"));
#process.stopTreeMaker.vectorDouble.append(cms.InputTag("prodJets","DoubleSV"));
#process.stopTreeMaker.vectorDouble.append(cms.InputTag("prodJets","cMVA"));
#process.stopTreeMaker.vectorDouble.append(cms.InputTag("prodJets","cMVAv2"));
#process.stopTreeMaker.vectorDouble.append(cms.InputTag("prodJets","cMVAv2Neg"));
#process.stopTreeMaker.vectorDouble.append(cms.InputTag("prodJets","cMVAv2Pos"));
#process.stopTreeMaker.vectorInt.append(cms.InputTag("prodJets","nTracks"));
#process.stopTreeMaker.vectorInt.append(cms.InputTag("prodJets","nSVs"));
process.stopTreeMaker.vectorDouble.append(cms.InputTag("prodJets", "JetProba"))
process.stopTreeMaker.vectorDoubleNamesInTree.append("prodJets:JetProba|JetProba_0")
#process.stopTreeMaker.vectorDouble.append(cms.InputTag("prodJets", "JetProbaN"))
#process.stopTreeMaker.vectorDoubleNamesInTree.append("prodJets:JetProbaN|JetProbaN_0")
#process.stopTreeMaker.vectorDouble.append(cms.InputTag("prodJets", "JetProbaP"))
#process.stopTreeMaker.vectorDoubleNamesInTree.append("prodJets:JetProbaP|JetProbaP_0")
process.stopTreeMaker.vectorDouble.append(cms.InputTag("prodJets", "JetBprob"))
process.stopTreeMaker.vectorDoubleNamesInTree.append("prodJets:JetProb|JetProb_0")
#process.stopTreeMaker.vectorDouble.append(cms.InputTag("prodJets", "JetBprobN"))
#process.stopTreeMaker.vectorDoubleNamesInTree.append("prodJets:JetProbN|JetProbN_0")
#process.stopTreeMaker.vectorDouble.append(cms.InputTag("prodJets", "JetBprobP"))
#process.stopTreeMaker.vectorDoubleNamesInTree.append("prodJets:JetProbP|JetProbP_0")
process.stopTreeMaker.vectorDouble.append(cms.InputTag("prodJets", "recoJetsBtag"))
process.stopTreeMaker.vectorDoubleNamesInTree.append("prodJets:recoJetsBtag|recoJetsBtag_0")
process.stopTreeMaker.vectorDouble.append(cms.InputTag("prodJets", "recoJetsCharge"))
process.stopTreeMaker.vectorDoubleNamesInTree.append("prodJets:recoJetsCharge|recoJetsCharge_0")

##CSV input variables
#process.stopTreeMaker.vectorDouble.append(cms.InputTag("prodJets", "CSVTrackJetPt"))
#process.stopTreeMaker.vectorDouble.append(cms.InputTag("prodJets", "CSVVertexCategory"))
#process.stopTreeMaker.vectorInt.append(   cms.InputTag("prodJets", "CSVJetNSecondaryVertices"))
#process.stopTreeMaker.vectorDouble.append(cms.InputTag("prodJets", "CSVTrackSumJetEtRatio"))
#process.stopTreeMaker.vectorDouble.append(cms.InputTag("prodJets", "CSVTrackSumJetDeltaR"))
#process.stopTreeMaker.vectorDouble.append(cms.InputTag("prodJets", "CSVTrackSip2dValAboveCharm"))
#process.stopTreeMaker.vectorDouble.append(cms.InputTag("prodJets", "CSVTrackSip2dSigAboveCharm"))
#process.stopTreeMaker.vectorDouble.append(cms.InputTag("prodJets", "CSVTrackSip3dValAboveCharm"))
#process.stopTreeMaker.vectorDouble.append(cms.InputTag("prodJets", "CSVTrackSip3dSigAboveCharm"))
#process.stopTreeMaker.vectorDouble.append(cms.InputTag("prodJets", "CSVVertexMass"))
#process.stopTreeMaker.vectorInt.append(   cms.InputTag("prodJets", "CSVVertexNTracks"))
#process.stopTreeMaker.vectorDouble.append(cms.InputTag("prodJets", "CSVVertexEnergyRatio"))
#process.stopTreeMaker.vectorDouble.append(cms.InputTag("prodJets", "CSVVertexJetDeltaR"))
#process.stopTreeMaker.vectorDouble.append(cms.InputTag("prodJets", "CSVFlightDistance2dVal"))
#process.stopTreeMaker.vectorDouble.append(cms.InputTag("prodJets", "CSVFlightDistance2dSig"))
#process.stopTreeMaker.vectorDouble.append(cms.InputTag("prodJets", "CSVFlightDistance3dVal"))
#process.stopTreeMaker.vectorDouble.append(cms.InputTag("prodJets", "CSVFlightDistance3dSig"))
#
##Ctagger input variables 
#process.stopTreeMaker.vectorDouble.append(cms.InputTag("prodJets", "CTagVertexCategory"))
#process.stopTreeMaker.vectorInt.append(   cms.InputTag("prodJets", "CTagJetNSecondaryVertices"))
#process.stopTreeMaker.vectorDouble.append(cms.InputTag("prodJets", "CTagTrackSumJetEtRatio"))
#process.stopTreeMaker.vectorDouble.append(cms.InputTag("prodJets", "CTagTrackSumJetDeltaR"))
#process.stopTreeMaker.vectorDouble.append(cms.InputTag("prodJets", "CTagTrackSip2dSigAboveCharm"))
#process.stopTreeMaker.vectorDouble.append(cms.InputTag("prodJets", "CTagTrackSip3dSigAboveCharm"))
#process.stopTreeMaker.vectorDouble.append(cms.InputTag("prodJets", "CTagVertexMass"))
#process.stopTreeMaker.vectorInt.append(   cms.InputTag("prodJets", "CTagVertexNTracks"))
#process.stopTreeMaker.vectorDouble.append(cms.InputTag("prodJets", "CTagVertexEnergyRatio"))
#process.stopTreeMaker.vectorDouble.append(cms.InputTag("prodJets", "CTagVertexJetDeltaR"))
#process.stopTreeMaker.vectorDouble.append(cms.InputTag("prodJets", "CTagFlightDistance2dSig"))
#process.stopTreeMaker.vectorDouble.append(cms.InputTag("prodJets", "CTagFlightDistance3dSig"))
#process.stopTreeMaker.vectorDouble.append(cms.InputTag("prodJets", "CTagMassVertexEnergyFraction"))
#process.stopTreeMaker.vectorDouble.append(cms.InputTag("prodJets", "CTagVertexBoostOverSqrtJetPt"))
#process.stopTreeMaker.vectorDouble.append(cms.InputTag("prodJets", "CTagVertexLeptonCategory"))

#pfcands
process.stopTreeMaker.vectorVectorTLorentzVector.append(cms.InputTag("prodJets", "chargedPFCandLV"))
process.stopTreeMaker.vectorVectorTLorentzVector.append(cms.InputTag("prodJets", "neutralPFCandLV"))

process.stopTreeMaker.vectorVectorDouble.append(cms.InputTag("prodJets", "chargedPFDxy"))
process.stopTreeMaker.vectorVectorDouble.append(cms.InputTag("prodJets", "chargedPFDz"))
#process.stopTreeMaker.vectorVectorDouble.append(cms.InputTag("prodJets", "chargedPFFromPV"))
process.stopTreeMaker.vectorVectorDouble.append(cms.InputTag("prodJets", "chargedPFVertexChi2"))
process.stopTreeMaker.vectorVectorDouble.append(cms.InputTag("prodJets", "chargedPFVertexNdof"))
#process.stopTreeMaker.vectorVectorDouble.append(cms.InputTag("prodJets", "chargedPFVertexMass"))
process.stopTreeMaker.vectorVectorDouble.append(cms.InputTag("prodJets", "neutralPFHCALFrac"))

####process.stopTreeMaker.vectorTLorentzVector.append(cms.InputTag("prodJets", "puppiAK8LVec"))
####process.stopTreeMaker.vectorDouble.append(cms.InputTag("prodJets", "puppiAK8Tau1"))
####process.stopTreeMaker.vectorDouble.append(cms.InputTag("prodJets", "puppiAK8Tau2"))
####process.stopTreeMaker.vectorDouble.append(cms.InputTag("prodJets", "puppiAK8Tau3"))
####process.stopTreeMaker.vectorDouble.append(cms.InputTag("prodJets", "puppiAK8SoftDropMass"))
####
####process.stopTreeMaker.vectorVectorTLorentzVector.append(cms.InputTag("prodJets", "puppiAK8SubjetLVec"))
####process.stopTreeMaker.vectorVectorDouble.append(cms.InputTag("prodJets", "puppiAK8SubjetMult"))
####process.stopTreeMaker.vectorVectorDouble.append(cms.InputTag("prodJets", "puppiAK8SubjetPtD"))
####process.stopTreeMaker.vectorVectorDouble.append(cms.InputTag("prodJets", "puppiAK8SubjetAxis1"))
####process.stopTreeMaker.vectorVectorDouble.append(cms.InputTag("prodJets", "puppiAK8SubjetAxis2"))
####process.stopTreeMaker.vectorVectorDouble.append(cms.InputTag("prodJets", "puppiAK8SubjetBDisc"))


if not options.isData:
    process.stopTreeMaker.vectorString.append(cms.InputTag("prodGenInfo", "genDecayStrVec"))
    process.stopTreeMaker.vectorInt.extend([cms.InputTag("prodGenInfo", "genDecayIdxVec"), cms.InputTag("prodGenInfo", "genDecayPdgIdVec"), cms.InputTag("prodGenInfo", "genDecayMomIdxVec"), cms.InputTag("prodGenInfo", "genDecayMomRefVec"), cms.InputTag("prodGenInfo", "WemuVec"), cms.InputTag("prodGenInfo", "WtauVec"), cms.InputTag("prodGenInfo", "WtauemuVec"), cms.InputTag("prodGenInfo", "WtauprongsVec"), cms.InputTag("prodGenInfo", "WtaunuVec"),  cms.InputTag("prodGenInfo","selPDGid")])
    process.stopTreeMaker.vectorIntNamesInTree.extend(["prodGenInfo:WemuVec|W_emuVec", "prodGenInfo:WtauVec|W_tauVec", "prodGenInfo:WtauemuVec|W_tau_emuVec", "prodGenInfo:WtauprongsVec|W_tau_prongsVec", "prodGenInfo:WtaunuVec|W_tau_nuVec"])
    process.stopTreeMaker.vectorDouble.extend([cms.InputTag("prodGenInfo", "WemupfActivityVec"), cms.InputTag("prodGenInfo", "WtauemupfActivityVec"), cms.InputTag("prodGenInfo", "WtauprongspfActivityVec")])
    process.stopTreeMaker.vectorDoubleNamesInTree.extend(["prodGenInfo:WemupfActivityVec|W_emu_pfActivityVec", "prodGenInfo:WtauemupfActivityVec|W_tau_emu_pfActivityVec", "prodGenInfo:WtauprongspfActivityVec|W_tau_prongs_pfActivityVec"])
    process.stopTreeMaker.vectorTLorentzVector.append(cms.InputTag("prodGenInfo", "genDecayLVec"))
    #process.stopTreeMaker.vectorTLorentzVector.append(cms.InputTag("prodGenInfo", "selGenParticle"))


    process.stopTreeMaker.varsDouble.extend([cms.InputTag("prodMET:genmet"), cms.InputTag("prodMET:genmetphi")])
    # Note that this default met from prodMET is both e/gamma and muon corrected which is the recommended one

process.stopTreeMaker.varsDouble.extend([cms.InputTag("prodMET", "met"), cms.InputTag("prodMET", "metphi")])


process.stopTreeMaker.varsInt.append(cms.InputTag("prodMuons", "nMuons"))
process.stopTreeMaker.varsIntNamesInTree.append("prodMuons:nMuons|nMuons_CUT")
process.stopTreeMaker.varsInt.append(cms.InputTag("prodMuonsNoIso", "nMuons"))
process.stopTreeMaker.vectorTLorentzVector.append(cms.InputTag("prodMuonsNoIso", "muonsLVec"))
process.stopTreeMaker.vectorDouble.extend([cms.InputTag("prodMuonsNoIso", "muonsCharge"), cms.InputTag("prodMuonsNoIso", "muonsMtw"), cms.InputTag("prodMuonsNoIso", "muonsRelIso"), cms.InputTag("prodMuonsNoIso", "muonsMiniIso"), cms.InputTag("prodMuonsNoIso", "muonspfActivity")])
process.stopTreeMaker.vectorInt.extend([cms.InputTag("prodMuonsNoIso", "muonsFlagMedium"), cms.InputTag("prodMuonsNoIso", "muonsFlagTight")])

process.stopTreeMaker.varsInt.append(cms.InputTag("prodElectrons", "nElectrons"))
process.stopTreeMaker.varsIntNamesInTree.append("prodElectrons:nElectrons|nElectrons_CUT")
process.stopTreeMaker.varsInt.append(cms.InputTag("prodElectronsNoIso", "nElectrons"))
process.stopTreeMaker.vectorTLorentzVector.append(cms.InputTag("prodElectronsNoIso", "elesLVec"))
process.stopTreeMaker.vectorDouble.extend([cms.InputTag("prodElectronsNoIso", "elesCharge"), cms.InputTag("prodElectronsNoIso", "elesMtw"), cms.InputTag("prodElectronsNoIso", "elesRelIso"), cms.InputTag("prodElectronsNoIso", "elesMiniIso"), cms.InputTag("prodElectronsNoIso", "elespfActivity")])
process.stopTreeMaker.vectorBool.extend([cms.InputTag("prodElectronsNoIso", "elesisEB")])
process.stopTreeMaker.vectorInt.extend([cms.InputTag("prodElectronsNoIso", "elesFlagMedium"), cms.InputTag("prodElectronsNoIso", "elesFlagVeto")])
                                   
process.stopTreeMaker.vectorTLorentzVector.append(cms.InputTag("prodIsoTrks:trksForIsoVetoLVec"))
process.stopTreeMaker.vectorDouble.extend([cms.InputTag("prodIsoTrks:trksForIsoVetocharge"), cms.InputTag("prodIsoTrks:trksForIsoVetodz"), cms.InputTag("prodIsoTrks:trksForIsoVetoiso"), cms.InputTag("prodIsoTrks:trksForIsoVetopfActivity"), cms.InputTag("prodIsoTrks:looseisoTrkscharge"), cms.InputTag("prodIsoTrks:looseisoTrksdz"), cms.InputTag("prodIsoTrks:looseisoTrksiso"), cms.InputTag("prodIsoTrks:looseisoTrksmtw"), cms.InputTag("prodIsoTrks:looseisoTrkspfActivity")])
process.stopTreeMaker.vectorDoubleNamesInTree.extend(["prodIsoTrks:trksForIsoVetocharge|trksForIsoVeto_charge", "prodIsoTrks:trksForIsoVetodz|trksForIsoVeto_dz", "prodIsoTrks:trksForIsoVetoiso|trksForIsoVeto_iso", "prodIsoTrks:trksForIsoVetopfActivity|trksForIsoVeto_pfActivity", "prodIsoTrks:looseisoTrkscharge|loose_isoTrks_charge", "prodIsoTrks:looseisoTrksdz|loose_isoTrks_dz", "prodIsoTrks:looseisoTrksiso|loose_isoTrks_iso", "prodIsoTrks:looseisoTrksmtw|loose_isoTrks_mtw", "prodIsoTrks:looseisoTrkspfActivity|loose_isoTrks_pfActivity"])
process.stopTreeMaker.vectorInt.extend([cms.InputTag("prodIsoTrks:trksForIsoVetopdgId"), cms.InputTag("prodIsoTrks:trksForIsoVetoidx"), cms.InputTag("prodIsoTrks:looseisoTrkspdgId"), cms.InputTag("prodIsoTrks:looseisoTrksidx"), cms.InputTag("prodIsoTrks:forVetoIsoTrksidx")])
process.stopTreeMaker.vectorIntNamesInTree.extend(["prodIsoTrks:trksForIsoVetopdgId|trksForIsoVeto_pdgId", "prodIsoTrks:trksForIsoVetoidx|trksForIsoVeto_idx", "prodIsoTrks:looseisoTrkspdgId|loose_isoTrks_pdgId", "prodIsoTrks:looseisoTrksidx|loose_isoTrks_idx"])
process.stopTreeMaker.vectorTLorentzVector.append(cms.InputTag("prodIsoTrks:looseisoTrksLVec"))
process.stopTreeMaker.vectorTLorentzVectorNamesInTree.append("prodIsoTrks:looseisoTrksLVec|loose_isoTrksLVec")
process.stopTreeMaker.varsInt.extend([cms.InputTag("prodIsoTrks:loosenIsoTrks"), cms.InputTag("prodIsoTrks:nIsoTrksForVeto")])
process.stopTreeMaker.varsIntNamesInTree.extend(["prodIsoTrks:loosenIsoTrks|loose_nIsoTrks", "prodIsoTrks:nIsoTrksForVeto|nIsoTrks_CUT"])

process.stopTreeMaker.vectorInt.append(cms.InputTag("triggerProducer", "PassTrigger"))
process.stopTreeMaker.vectorInt.append(cms.InputTag("triggerProducer", "TriggerPrescales"))
process.stopTreeMaker.vectorString.append(cms.InputTag("triggerProducer", "TriggerNames"))

# prodGoodVertices has the same as vtxSize in prodEventInfo...
#process.stopTreeMaker.varsInt.append(cms.InputTag("prodGoodVertices"))
#process.stopTreeMaker.varsInt.append(cms.InputTag("prodFilterOutScraping"))
process.stopTreeMaker.varsBool.append(cms.InputTag("prodJetIDEventFilter", "looseJetID"))
process.stopTreeMaker.varsBool.append(cms.InputTag("prodJetIDEventFilter", "tightJetID"))
process.stopTreeMaker.varsBool.append(cms.InputTag("prodJetIDEventFilter", "tightlepvetoJetID"))
#process.stopTreeMaker.varsInt.append(cms.InputTag("METFilters"))
#process.stopTreeMaker.varsInt.append(cms.InputTag("CSCTightHaloFilter")) # 74X txt files are ready for the 2015 working point, use this and not the flag in miniAOD 
process.stopTreeMaker.varsInt.append(cms.InputTag("globalTightHalo2016Filter"))
process.stopTreeMaker.varsInt.append(cms.InputTag("goodVerticesFilter"))
process.stopTreeMaker.varsInt.append(cms.InputTag("eeBadScFilter"))
process.stopTreeMaker.varsInt.append(cms.InputTag("EcalDeadCellTriggerPrimitiveFilter"))
process.stopTreeMaker.varsBool.append(cms.InputTag("BadChargedCandidateFilter"))
process.stopTreeMaker.varsBool.append(cms.InputTag("BadPFMuonFilter"))

process.stopTreeMaker.varsInt.append(cms.InputTag("noBadMuonsFilter"))
process.stopTreeMaker.varsInt.append(cms.InputTag("badMuonsFilter"))
process.stopTreeMaker.varsInt.append(cms.InputTag("duplicateMuonsFilter"))

process.stopTreeMaker.varsBool.append(cms.InputTag("HBHENoiseFilterResultProducer", "HBHENoiseFilterResult"))
process.stopTreeMaker.varsBoolNamesInTree.append("HBHENoiseFilterResultProducer:HBHENoiseFilterResult|HBHENoiseFilter")

process.stopTreeMaker.varsBool.append(cms.InputTag("HBHENoiseFilterResultProducer", "HBHEIsoNoiseFilterResult"))
process.stopTreeMaker.varsBoolNamesInTree.append("HBHENoiseFilterResultProducer:HBHEIsoNoiseFilterResult|HBHEIsoNoiseFilter")

process.stopTreeMaker.varsDouble.extend([cms.InputTag("prodMET:calomet"), cms.InputTag("prodMET:calometphi")])

if not options.isData:
   #isrJets
   process.stopTreeMaker.varsInt.append(cms.InputTag("ISRJetProducer", "NJetsISR"))

process.stopTreeMaker.varsInt.extend([cms.InputTag("prodEventInfo:vtxSize"), cms.InputTag("prodEventInfo:npv"), cms.InputTag("prodEventInfo:nm1"), cms.InputTag("prodEventInfo:n0"), cms.InputTag("prodEventInfo:np1")])
process.stopTreeMaker.varsDouble.extend([cms.InputTag("prodEventInfo:trunpv"), cms.InputTag("prodEventInfo:avgnpv"), cms.InputTag("prodEventInfo:storedWeight")])
process.stopTreeMaker.varsDoubleNamesInTree.extend(["prodEventInfo:trunpv|tru_npv", "prodEventInfo:avgnpv|avg_npv", "prodEventInfo:storedWeight|stored_weight"])

process.stopTreeMaker.vectorBool.append(cms.InputTag("goodPhotons", "loosePhotonID"))
process.stopTreeMaker.vectorBool.append(cms.InputTag("goodPhotons", "mediumPhotonID"))
process.stopTreeMaker.vectorBool.append(cms.InputTag("goodPhotons", "tightPhotonID"))
                                                  
process.stopTreeMaker.vectorDouble.append(cms.InputTag("goodPhotons", "pfGammaIso"))
process.stopTreeMaker.vectorDouble.append(cms.InputTag("goodPhotons", "isEB"))
process.stopTreeMaker.vectorDouble.append(cms.InputTag("goodPhotons", "genMatched"))
process.stopTreeMaker.vectorDouble.append(cms.InputTag("goodPhotons", "hadTowOverEM"))
process.stopTreeMaker.vectorDouble.append(cms.InputTag("goodPhotons", "sigmaIetaIeta"))
process.stopTreeMaker.vectorDouble.append(cms.InputTag("goodPhotons", "pfChargedIso"))
process.stopTreeMaker.vectorDouble.append(cms.InputTag("goodPhotons", "pfNeutralIso"))
process.stopTreeMaker.vectorDouble.append(cms.InputTag("goodPhotons", "pfChargedIsoRhoCorr"))
process.stopTreeMaker.vectorDouble.append(cms.InputTag("goodPhotons", "pfNeutralIsoRhoCorr"))
process.stopTreeMaker.vectorDouble.append(cms.InputTag("goodPhotons", "pfGammaIsoRhoCorr"))
process.stopTreeMaker.vectorDouble.append(cms.InputTag("goodPhotons", "hasPixelSeed"))
process.stopTreeMaker.vectorDouble.append(cms.InputTag("goodPhotons", "passElectronVeto"))
#process.stopTreeMaker.vectorBool.append(cms.InputTag("goodPhotons", "hadronization"))
process.stopTreeMaker.vectorBool.append(cms.InputTag("goodPhotons", "nonPrompt"))
#process.stopTreeMaker.vectorBool.append(cms.InputTag("goodPhotons", "fullID"))
process.stopTreeMaker.vectorDouble.append(cms.InputTag("goodPhotons", "photonPt"))
process.stopTreeMaker.vectorDouble.append(cms.InputTag("goodPhotons", "photonEta"))
process.stopTreeMaker.vectorDouble.append(cms.InputTag("goodPhotons", "photonPhi"))
process.stopTreeMaker.vectorTLorentzVector.append(cms.InputTag("goodPhotons", "gammaLVec"))
process.stopTreeMaker.vectorTLorentzVector.append(cms.InputTag("goodPhotons", "gammaLVecGen"))


###############################################################################################################################

#process.p = cms.Path(process.QGTagger * process.printDecayPythia8 
#                     * process.prodGenInfo * process.prodJets * process.prodMET 
#                     * process.prodMuons * process.prodElectrons
#                     * process.pfCHS * process.pfNoMuonCHSNoMu * process.pfNoElectronCHSNoEle * process.ak4PFJetsCHSNoLep
#                     * process.updatedPatJetsDeepFlavour
#                     * process.prodMuonsNoIso * process.prodElectronsNoIso * process.prodIsoTrks 
#                     * process.stopTreeMaker)

#process.dump=cms.EDAnalyzer('EventContentAnalyzer')

process.p = cms.Path(process.stopTreeMaker)
#process.p = cms.Path(process.pfDeepFlavourJetTags*process.stopTreeMaker*process.dump)

