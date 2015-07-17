import FWCore.ParameterSet.Config as cms

process = cms.Process("ProdFatTops")


process.load("Configuration.EventContent.EventContent_cff")
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )
process.GlobalTag.globaltag = 'PHYS14_25_V1'

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        'file:/uscms_data/d3/benwu/TTJets_PU20bx25_PHYS14_v2.root'
    )
)

process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )
process.options.allowUnscheduled = cms.untracked.bool(True)

#============================================================================#
#--------------------------     HEP v2 Producer     -------------------------#
#============================================================================#
# Select candidates that would pass CHS requirements
# This can be used as input for HTT and other jet clustering algorithms
process.chs = cms.EDFilter("CandPtrSelector", src = cms.InputTag("packedPFCandidates"), cut = cms.string("fromPV"))

from RecoJets.JetProducers.AnomalousCellParameters_cfi import *
from RecoJets.JetProducers.PFJetParameters_cfi import *
process.HTT = cms.EDProducer(
        "HTTTopJetProducer",
        PFJetParameters.clone( src = cms.InputTag("chs"),
                               doAreaFastjet = cms.bool(True),
                               doRhoFastjet = cms.bool(False),
                               jetPtMin = cms.double(100.0)
                           ),
        AnomalousCellParameters,
        optimalR = cms.bool(False),
        algorithm = cms.int32(1),
        jetAlgorithm = cms.string("CambridgeAachen"),
        rParam = cms.double(1.5),
        mode = cms.int32(4),
        minFatjetPt = cms.double(200.),
        minCandPt = cms.double(200.),
        minSubjetPt = cms.double(30.),
        writeCompound = cms.bool(True),
        minCandMass = cms.double(0.),
        maxCandMass = cms.double(1000),
        massRatioWidth = cms.double(100.),
       minM23Cut = cms.double(0.),
        minM13Cut = cms.double(0.),
        maxM13Cut = cms.double(2.))

process.ORHTT = process.HTT.clone(optimalR = cms.bool(True))

#============================================================================#
#----------------------------     Soft  Drop     ----------------------------#
#============================================================================#
from JMEAnalysis.JetToolbox.jetToolbox_cff import jetToolbox
#jetToolbox(process, 'ak8', 'ak8JetSubs', 'out', miniAOD=True, addSoftDrop=True)
jetToolbox(process, 'ak8', 'ak8JetSubs', 'out', miniAOD=True, addSoftDrop=True, addSoftDropSubjets=True, addNsub=True)
jetToolbox(process, 'ca15', 'ca15JetSubs', 'out', miniAOD=True, addSoftDrop=True, addSoftDropSubjets=True, addNsub=True)
#jetToolbox(process, 'ak8', 'ak8JetSubs', 'out', miniAOD=True, addSoftDrop=True, addSoftDropSubjets=True, addNsub=True, PUMethod='CHS')


##============================================================================#
##-----------------------     Shower Deconstruction     ----------------------#
##============================================================================#
import os
input_card = os.environ['CMSSW_BASE'] +"/src/ShowerDeconstruction/SDProducer/data/input_card.dat"
#input_card = os.environ['CMSSW_BASE'] +"/src/ShowerDeconstruction/SDProducer/data/input_card.dat"
#print '--------------------------------------------'
#print input_card
#print '/'.join(input_card.split('/')[:-1])
#print os.listdir('/'.join(input_card.split('/')[:-1]))
#print '--------------------------------------------'
sd_name = "ShowerDecon"
fj_name = "ak8PFJetsCHS"
setattr(process, sd_name+"2", cms.EDProducer("SDProducer",
                                             FatjetName = cms.string(fj_name),
                                             MicrojetCone = cms.double(0.2),
                                             InputCard = cms.string(input_card)))
setattr(process, sd_name+"3", cms.EDProducer("SDProducer",
                                             FatjetName = cms.string(fj_name),
                                             MicrojetCone = cms.double(0.3),
                                             InputCard = cms.string(input_card)))


#============================================================================#
#------------------------------     GenJet     ------------------------------#
#============================================================================#
## GenJet
from RecoJets.JetProducers.ak5GenJets_cfi import ak5GenJets
process.load('RecoJets.Configuration.GenJetParticles_cff')
process.genParticlesForJetsNoNu.src =  'packedGenParticles'
process.ak4GenJets = ak5GenJets.clone(rParam = 0.4, jetPtMin= 0.0,
                                      doAreaFastjet=True, src = 'genParticlesForJetsNoNu')
process.ak8GenJets = ak5GenJets.clone(rParam = 0.8, jetPtMin= 0.0,
                                      doAreaFastjet=True, src = 'genParticlesForJetsNoNu')

##### Prod
process.load("SusyAnaTools.SkimsAUX.prodFatTops_cfi")

process.load("SusyAnaTools.SignalScan.genDecayStringMakerPythia8_cfi")
process.printDecayPythia8.src = cms.InputTag("prunedGenParticles")
process.printDecayPythia8.keyDecayStrs = cms.vstring("t", "tbar", "~chi_1+", "~chi_1-")
process.load("SusyAnaTools.SkimsAUX.prodGenInfo_cfi")
process.load("SusyAnaTools.SkimsAUX.prodMET_cfi")

process.load("SusyAnaTools.Skims.vertex_cfi")

process.goodVertices = cms.EDFilter(
  "VertexSelector",
  filter = cms.bool(False),
  src = cms.InputTag("offlineSlimmedPrimaryVertices"),
  cut = cms.string("!isFake && ndof > 4 && abs(z) <= 24 && position.rho < 2")
)

#### Muon
process.load("SusyAnaTools.SkimsAUX.prodMuons_cfi")
process.prodMuonsNoIso = process.prodMuons.clone()
process.prodMuonsNoIso.DoMuonIsolation = cms.bool(False)
### Ele
process.load("SusyAnaTools.SkimsAUX.prodElectrons_cfi")
process.prodElectronsNoIso = process.prodElectrons.clone()
process.prodElectronsNoIso.DoElectronIsolation = cms.bool(False)

#### Iso
process.load("SusyAnaTools.Skims.trackIsolationMaker_cfi")
process.trackIsolation = process.trackIsolationFilter.clone()
process.trackIsolation.pfCandidatesTag = cms.InputTag("packedPFCandidates")
process.trackIsolation.doTrkIsoVeto = cms.bool(False)

process.load("SusyAnaTools.SkimsAUX.prodIsoTrks_cfi")


process.HEPTagger = cms.Sequence (process.chs * process.HTT )

process.p = cms.Path(process.HEPTagger * process.prodFatTops)

process.out.fileName = cms.untracked.string('myOutputFile.root')
process.out.outputCommands = cms.untracked.vstring(
    #'keep *_ak8PFJetsCHSSoftDropLinks_*_*',
    #'keep *_selectedPatJet*_*_*',
    #'keep *_Njet*_*_*',
    #'keep *_ak8PFJetsCHS_*_*',
    #'keep *_prodFatTops_*_*',
    #'keep *_prodGenInfo_*_*',
    #'keep *_ShowerDecon*_*_*',
    #'keep *_slimmed*_*_*',
    #'keep *_HTT_*_*'

    'drop *'
    #'keep *_prod*_*_*',
)

process.TFileService = cms.Service("TFileService",
   fileName = cms.string('stopFlatNtuples.root')
)

process.load("SusyAnaTools.StopTreeMaker.stopTreeMaker_cfi")
process.stopTreeMaker.TreeName = cms.string("AUX")
#process.stopTreeMaker.vectorString.append(cms.InputTag("prodGenInfo", "genDecayStrVec"))
process.stopTreeMaker.vectorInt.extend([cms.InputTag("prodGenInfo", "genDecayIdxVec"), cms.InputTag("prodGenInfo", "genDecayPdgIdVec"), cms.InputTag("prodGenInfo", "genDecayMomIdxVec"), cms.InputTag("prodGenInfo", "genDecayMomRefVec")])
#process.stopTreeMaker.vectorInt.extend([cms.InputTag("prodGenInfo", "genDecayIdxVec"), cms.InputTag("prodGenInfo", "genDecayPdgIdVec"), cms.InputTag("prodGenInfo", "genDecayMomIdxVec"), cms.InputTag("prodGenInfo", "genDecayMomRefVec"), cms.InputTag("prodGenInfo", "WemuVec"), cms.InputTag("prodGenInfo", "WtauVec"), cms.InputTag("prodGenInfo", "WtauemuVec"), cms.InputTag("prodGenInfo", "WtauprongsVec"), cms.InputTag("prodGenInfo", "WtaunuVec")])
#process.stopTreeMaker.vectorIntNamesInTree.extend(["prodGenInfo:WemuVec|W_emuVec", "prodGenInfo:WtauVec|W_tauVec", "prodGenInfo:WtauemuVec|W_tau_emuVec", "prodGenInfo:WtauprongsVec|W_tau_prongsVec", "prodGenInfo:WtaunuVec|W_tau_nuVec"])
process.stopTreeMaker.vectorTLorentzVector.append(cms.InputTag("prodGenInfo", "genDecayLVec"))

process.stopTreeMaker.vectorDouble.extend([cms.InputTag("prodFatTops", "recoJetsBtag"), cms.InputTag("prodFatTops", "AK8SD2Chi"), cms.InputTag("prodFatTops", "AK8SD3Chi"),
                                           cms.InputTag("prodFatTops", "AK8SoftDropTau32"), cms.InputTag("prodFatTops", "CA15SoftDropTau32")])
process.stopTreeMaker.varsDouble.extend([cms.InputTag("prodMET:met"), cms.InputTag("prodMET:metphi")])
#process.stopTreeMaker.vectorDoubleNamesInTree.extend(["prodFatTops:AK8SD2Chi|AK8SD2Chi_0", "prodFatTops:recoJetsBtag|recoJetsBtag_0", "prodFatTops:AK8SD3Chi|AK8SD3Chi_0"])

### Muons
process.stopTreeMaker.varsInt.append(cms.InputTag("prodMuonsNoIso", "nMuons"))
process.stopTreeMaker.vectorTLorentzVector.append(cms.InputTag("prodMuonsNoIso", "muonsLVec"))
process.stopTreeMaker.vectorDouble.extend([cms.InputTag("prodMuonsNoIso", "muonsCharge"), cms.InputTag("prodMuonsNoIso", "muonsMtw"), cms.InputTag("prodMuonsNoIso", "muonsRelIso"), cms.InputTag("prodMuonsNoIso", "muonsMiniIso")])

process.stopTreeMaker.varsInt.append(cms.InputTag("prodElectronsNoIso", "nElectrons"))
process.stopTreeMaker.vectorTLorentzVector.append(cms.InputTag("prodElectronsNoIso", "elesLVec"))
process.stopTreeMaker.vectorDouble.extend([cms.InputTag("prodElectronsNoIso", "elesCharge"), cms.InputTag("prodElectronsNoIso", "elesMtw"), cms.InputTag("prodElectronsNoIso", "elesRelIso"), cms.InputTag("prodElectronsNoIso", "elesMiniIso")])
process.stopTreeMaker.vectorBool.extend([cms.InputTag("prodElectronsNoIso", "elesisEB")])

process.stopTreeMaker.vectorTLorentzVector.append(cms.InputTag("prodIsoTrks:looseisoTrksLVec"))
process.stopTreeMaker.vectorTLorentzVectorNamesInTree.append("prodIsoTrks:looseisoTrksLVec|loose_isoTrksLVec")
process.stopTreeMaker.vectorDouble.extend([cms.InputTag("prodIsoTrks:looseisoTrksiso"), cms.InputTag("prodIsoTrks:looseisoTrksmtw")])
process.stopTreeMaker.vectorDoubleNamesInTree.extend(["prodIsoTrks:looseisoTrksiso|loose_isoTrks_iso", "prodIsoTrks:looseisoTrksmtw|loose_isoTrks_mtw"])

process.stopTreeMaker.vectorInt.extend([cms.InputTag("prodIsoTrks:looseisoTrkspdgId")])
process.stopTreeMaker.vectorIntNamesInTree.extend(["prodIsoTrks:looseisoTrkspdgId|loose_isoTrks_pdgId"])

process.stopTreeMaker.vectorTLorentzVector.append(cms.InputTag("prodFatTops", "jetsLVec"))
process.stopTreeMaker.vectorTLorentzVector.append(cms.InputTag("prodFatTops", "Gen4LVec"))
process.stopTreeMaker.vectorTLorentzVector.append(cms.InputTag("prodFatTops", "Gen8LVec"))
process.stopTreeMaker.vectorTLorentzVector.append(cms.InputTag("prodFatTops", "CMSTopLVec"))
process.stopTreeMaker.vectorTLorentzVector.append(cms.InputTag("prodFatTops", "AK8LVec"))
process.stopTreeMaker.vectorTLorentzVector.append(cms.InputTag("prodFatTops", "HEPV2TopLVec"))
process.stopTreeMaker.vectorTLorentzVector.append(cms.InputTag("prodFatTops", "ORHEPV2TopLVec"))
process.stopTreeMaker.vectorTLorentzVector.append(cms.InputTag("prodFatTops", "AK8SoftDropLVec"))
process.stopTreeMaker.vectorTLorentzVector.append(cms.InputTag("prodFatTops", "CA15SoftDropLVec"))


process.e = cms.EndPath(process.stopTreeMaker)
#process.e = cms.EndPath(process.out * process.stopTreeMaker)
