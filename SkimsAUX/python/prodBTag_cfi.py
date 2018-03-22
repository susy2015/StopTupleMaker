import FWCore.ParameterSet.Config as cms

prodBTag = cms.EDFilter(
  "prodBTag",
  jetSrc = cms.InputTag('slimmedJets'),
  deepCSVBJetTags = cms.string('pfDeepCSVJetTags'),#'deepFlavourJetTags:probb'),#'deepFlavourJetTags'),
  deepFlavorBJetTags = cms.string('pfDeepFlavourJetTags'),
  deepCSVNegBJetTags = cms.string('negativeDeepFlavourJetTags'),#pfNegativeDeepCSVJetTags'),
  deepCSVPosBJetTags = cms.string('positiveDeepFlavourJetTags'),
  combinedSVBJetTags = cms.string('pfCombinedSecondaryVertexV2BJetTags'),
  combinedSVPosBJetTags = cms.string('pfPositiveCombinedSecondaryVertexV2BJetTags'),
  combinedSVNegBJetTags = cms.string('pfNegativeCombinedSecondaryVertexV2BJetTags'),
  combinedIVFSVBJetTags = cms.string('pfCombinedInclusiveSecondaryVertexV2BJetTags'),
  combinedIVFSVPosBJetTags = cms.string('pfPositiveCombinedInclusiveSecondaryVertexV2BJetTags'),
  combinedIVFSVNegBJetTags = cms.string('pfNegativeCombinedInclusiveSecondaryVertexV2BJetTags'),
  simpleSVHighEffBJetTags = cms.string('pfSimpleSecondaryVertexHighEffBJetTags'),
  simpleSVNegHighEffBJetTags = cms.string('pfNegativeSimpleSecondaryVertexHighEffBJetTags'),
  simpleSVHighPurBJetTags = cms.string('pfSimpleSecondaryVertexHighPurBJetTags'),
  simpleSVNegHighPurBJetTags = cms.string('pfNegativeSimpleSecondaryVertexHighPurBJetTags'),
  softPFMuonBJetTags = cms.string('softPFMuonBJetTags'),
  softPFMuonNegBJetTags = cms.string('negativeSoftPFMuonBJetTags'),
  softPFMuonPosBJetTags = cms.string('positiveSoftPFMuonBJetTags'),
  softPFElectronBJetTags = cms.string('softPFElectronBJetTags'),
  softPFElectronNegBJetTags = cms.string('negativeSoftPFElectronBJetTags'),
  softPFElectronPosBJetTags = cms.string('positiveSoftPFElectronBJetTags'),
  doubleSVBJetTags = cms.string('pfBoostedDoubleSecondaryVertexAK8BJetTags'),
  cMVABJetTags = cms.string('pfCombinedMVABJetTags'),
  cMVAv2BJetTags = cms.string('pfCombinedMVAV2BJetTags'),
  cMVAv2NegBJetTags = cms.string('pfNegativeCombinedMVAV2BJetTags'),
  cMVAv2PosBJetTags = cms.string('pfPositiveCombinedMVAV2BJetTags'),
  CvsBCJetTags = cms.string('pfCombinedCvsBJetTags'),
  CvsBNegCJetTags = cms.string('pfNegativeCombinedCvsBJetTags'),
  CvsBPosCJetTags = cms.string('pfPositiveCombinedCvsBJetTags'),
  CvsLCJetTags = cms.string('pfCombinedCvsLJetTags'),
  CvsLNegCJetTags = cms.string('pfNegativeCombinedCvsLJetTags'),
  CvsLPosCJetTags = cms.string('pfPositiveCombinedCvsLJetTags'),  
  tagInfoName = cms.string('pfDeepCSV'),#'deepNN'),
  ipTagInfos = cms.string('pfImpactParameter'),
  svTagInfos = cms.string('pfInclusiveSecondaryVertexFinder'),
  ipTagInfosCTag = cms.string('pfImpactParameter'),
  svTagInfosCTag = cms.string('pfInclusiveSecondaryVertexFinderCvsL'),
  softPFMuonTagInfosCTag = cms.string('softPFMuons'),
  softPFElectronTagInfosCTag = cms.string('softPFElectrons'),
  bDiscriminatorCSV = cms.vstring('pfCombinedInclusiveSecondaryVertexV2BJetTags'
                            ,'deepFlavourJetTags:probudsg'
                            ,'deepFlavourJetTags:probb'
                            ,'deepFlavourJetTags:probc'
                            ,'deepFlavourJetTags:probbb'
                            ,'deepFlavourJetTags:probcc'
                            ),
  bDiscriminators = cms.vstring(), 
  ak8JetSrc = cms.InputTag("selectedPatJetsAK8PFPuppi"),
  ak8ptCut = cms.double(200.0),

 )
