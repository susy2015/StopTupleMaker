
#from StopTupleMaker.Skims.StopJets_cff import *
#from StopTupleMaker.Skims.StopTrackIsolation_cff import *
##from StopTupleMaker.Skims.RA2Leptons_cff import *
#from StopTupleMaker.Skims.StopTauJets_cff        import *
from StopTupleMaker.Skims.StopObjects_cff        import *
from StopTupleMaker.Skims.StopDPhiSelection_cff  import *

stopFullPFchsSelectionNoMET = cms.Sequence(
  stopCountPFchsJetsPt30 *
  stopCountPFchsJetsPt70eta2p5 *
  stopCountBJets *
  jetsMETDPhiFilter *
  stopPFMuonVeto  *  
  stopElectronVeto * 
  stopTauJetVeto   
  #stopIsoTrackVeto # need to check the code
)

stopFullPFchsSelection = cms.Sequence(
  stopFullPFchsSelectionNoMET #*
  # add met filter here
)
