import FWCore.ParameterSet.Config as cms

trackIsolationFilter = cms.EDFilter("TrackIsolationFilter",
                                     pfCandidatesTag     = cms.InputTag("packedPFCandidates"),
                                     vertexInputTag      = cms.InputTag("goodVertices"),
                                     dR_ConeSize         = cms.double(0.3),
                                     dz_CutValue         = cms.double(0.1),
                                     minPt_PFCandidate   = cms.double(5.0),
                                     isoCut              = cms.double(0.2),
                                     doTrkIsoVeto        = cms.bool(True),
                                     exclPdgIdVec        = cms.vint32(),
)

trackIsolationCounter  = cms.EDFilter("PATCandViewCountFilter",
                                      minNumber = cms.uint32(0),
                                      maxNumber = cms.uint32(999999),
                                      src = cms.InputTag("trackIsolationFilter")
)

