
from StopTupleMaker.Skims.RA2Jets_cff import *
from StopTupleMaker.Skims.RA2HT_cff import *
from StopTupleMaker.Skims.RA2MHT_cff import *
from StopTupleMaker.Skims.RA2Leptons_cff import *
#from StopTupleMaker.Skims.RA2Photons_cff import *
#from StopTupleMaker.Skims.maskedECALTowers_cff import *
#from StopTupleMaker.Skims.RA2Cleaning_cff import *

ra2Objects = cms.Sequence(  
  ra2PFJets *
  htPF *
  mhtPF *
  ra2PFchsJets *
  htPFchs *
  mhtPFchs #*
#  ra2Muons *
#  ra2PFMuons *
#  ra2Electrons #*
#  ra2Photons
)
