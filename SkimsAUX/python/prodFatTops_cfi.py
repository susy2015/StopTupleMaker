import FWCore.ParameterSet.Config as cms


prodFatTops = cms.EDFilter(
    "prodFatTops",
    AK4Src = cms.InputTag('slimmedJets'),
    bTagKeyString = cms.string('combinedInclusiveSecondaryVertexV2BJetTags'), # for PHYS14
    #bTagKeyString = cms.string('pfCombinedInclusiveSecondaryVertexV2BJetTags'), # For Spring15
    Gen4Src = cms.InputTag('ak4GenJets'),
    Gen8Src = cms.InputTag('ak8GenJets'),
    AK8Src = cms.InputTag('slimmedJetsAK8'),
    RecoAK8Src = cms.InputTag('ak8PFJetsCHS'),
    HTTSrc = cms.InputTag('HTT'),
    ORHTTSrc = cms.InputTag('ORHTT'),
    ak8PFJetsSoftDrop =  cms.InputTag('selectedPatJetsAK8PFCHSSoftDrop'),
    ak8PFJetsSoftDropTau =  cms.string('NjettinessAK8CHS'),
    ca15PFJetsSoftDrop =  cms.InputTag('selectedPatJetsCA15PFCHSSoftDrop'),
    ca15PFJetsSoftDropTau  = cms.string('NjettinessCA15CHS'),
    ShowerD2Src = cms.string("ShowerDecon2"),
    ShowerD3Src = cms.string("ShowerDecon3"),
    debug  = cms.bool(False)
)
