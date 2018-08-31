
import FWCore.ParameterSet.Config as cms

prodGenJets = cms.EDFilter(
  "prodGenJets",
  genjetSrc = cms.InputTag('slimmedGenJets'),
  genJetAk8Src           = cms.InputTag("slimmedGenJetsAK8"),
)
