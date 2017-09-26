
from StopTupleMaker.Skims.RA2Cleaning_cff       import *
from StopTupleMaker.Skims.StopJets_cff           import *
from StopTupleMaker.Skims.StopLeptons_cff        import *
#from StopTupleMaker.Skims.StopPhotons_cff       import *
from StopTupleMaker.Skims.StopTrackIsolation_cff import *
from StopTupleMaker.Skims.StopTauJets_cff        import *
from StopTupleMaker.Skims.StopBTagJets_cff       import *

stopObjects = cms.Sequence(  
  stopPFJets *
  #stopMuons *     # use from RA2 for now, update soon
  #stopElectrons * # use from RA2 for now, update soon
  stopIsoTracks *
  stopTauJets *
  stopBJets
  #stopPhotons
)

stopObjectsOnPatTuples = cms.Sequence(  
  stopPFJets *
  #stopMuons *     # use from RA2 for now, update soon
  #stopElectrons * # use from RA2 for now, update soon
  stopIsoTracks *
  stopTauJets *
  stopBJets
  #stopPhotons
)

