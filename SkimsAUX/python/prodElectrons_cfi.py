
import FWCore.ParameterSet.Config as cms

prodElectrons = cms.EDFilter(
  "prodElectrons",
  VetoElectronID   =cms.InputTag('egmGsfElectronIDs:cutBasedElectronID-Fall17-94X-V2-veto'),
  LooseElectronID   =cms.InputTag('egmGsfElectronIDs:cutBasedElectronID-Fall17-94X-V2-loose'),
  MediumElectronID   =cms.InputTag('egmGsfElectronIDs:cutBasedElectronID-Fall17-94X-V2-medium'),
  TightElectronID   =cms.InputTag('egmGsfElectronIDs:cutBasedElectronID-Fall17-94X-V2-tight'),
  ElectronSource    = cms.InputTag('slimmedElectrons'),
  ConversionsSource = cms.InputTag("reducedEgamma", "reducedConversions"), #not used for WP VETO
  VertexSource      = cms.InputTag('offlineSlimmedPrimaryVertices'),#goodVertices'),
  metSource         = cms.InputTag('slimmedMETs'),
  PFCandSource = cms.InputTag('packedPFCandidates'),
  BeamSpotSource    = cms.InputTag("offlineBeamSpot"),
  RhoSource     = cms.InputTag('fixedGridRhoFastjetAll'),
  EAValues      = cms.vdouble(0.1440, 0.1562, 0.1032, 0.0859, 0.1116, 0.1321, 0.1654), #2017
  EAEtaValues   = cms.vdouble(1.0,    1.479,    2.0,    2.2,    2.3,    2.4),          #2017
  relIsoString  = cms.string("GsfEleRelPFIsoScaledCut_0"),  #2017, in 2016 or 94x V1, this was GsfEleEffAreaPFIsoCut_0
  MinElePt       = cms.double(5),
  MaxEleEta      = cms.double(2.5),
  MaxEleMiniIso  = cms.double(0.10),
  DoElectronVeto           = cms.bool(False),
  Dod0dz                   = cms.bool(True),
  DoElectronID             = cms.bool(True),
  DoElectronVtxAssociation = cms.bool(True),
  DoElectronIsolation      = cms.int32(0), # 1 for relIso; 2 for miniIso; 0 for nothing (but stores iso variables in ntuples to apply them offline)
  Debug                    = cms.bool(False)
)
