import FWCore.ParameterSet.Config as cms
ISRJetProducer = cms.EDFilter(
  "ISRJetProducer",
  genParticleSrc = cms.InputTag("prunedGenParticles"),
  cleanJetSrc = cms.InputTag("patJetsAK4PFCHSNoLep"),
  debug  = cms.bool(False)
)
