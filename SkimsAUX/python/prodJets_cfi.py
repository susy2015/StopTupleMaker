import FWCore.ParameterSet.Config as cms

prodJets = cms.EDFilter(
  "prodJets",

  jetSrc = cms.InputTag('slimmedJets'),
  puppiJetsSrc = cms.InputTag("slimmedJetsAK8"),

  jetType = cms.string('AK4PFchs'),

  qgTaggerKey = cms.string('QGTagger'),
  deepCSVBJetTags = cms.string('pfDeepCSVJetTags'),
  deepFlavorBJetTags = cms.string('pfDeepFlavourJetTags'),
  CvsBCJetTags = cms.string('pfCombinedCvsBJetTags'),
  CvsLCJetTags = cms.string('pfCombinedCvsLJetTags'),
  bTagKeyString = cms.string('pfCombinedInclusiveSecondaryVertexV2BJetTags'),

  jetPBJetTags = cms.string("jetPBJetTags"),
  jetBPBJetTags = cms.string("jetBPBJetTags"),

  NjettinessAK8Puppi_label = cms.string('NjettinessAK8Puppi'),
  ak8PFJetsPuppi_label = cms.string('ak8PFJetsPuppi'),

  eleLVec = cms.InputTag("prodElectrons:elesLVec"), 
  muLVec = cms.InputTag("prodMuons:muonsLVec"), 
  trksForIsoVetoLVec = cms.InputTag("prodIsoTrks:trksForIsoVetoLVec"),
  looseisoTrksLVec = cms.InputTag("prodIsoTrks:looseisoTrksLVec"),

  debug  = cms.bool(False),

  AK4minPTcut_= cms.double(10.0),
  AK8minPTcut_= cms.double(170.0),

  genMatch_dR = cms.untracked.double(1.0),
  deltaRcon = cms.untracked.double(0.01),

  datapath = cms.string("NNKit/data/ak8/full")
  )

