#!/usr/bin/env python
# encoding: utf-8

# File        : MultiCrab3.py
# Author      : Ben Wu
# Contact     : benwu@fnal.gov
# Date        : 2015 Apr 01
#
# Description :


import copy, os, time

from CRABAPI.RawCommand import crabCommand
from crab3Config import config as config

workArea = 'CrabTest/TopTagTest5'
outDir =  '/store/group/lpcsusyhad/benwu/PHYS14_FatTop/TopTagTest5'
#outDir =  '/store/user/benwu/PHYS14/TopTagTest3'

jobslist = {
    #'SLHC23p1_SHCal_140PU_MET' : ['../python/METPerformance_SHCal_140PU_cfg_RECO.py',  '/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola/TP2023SHCALDR-SHCALJan23_PU140BX25_PH2_1K_FB_V6-v1/GEN-SIM-RECO',           3],
    #'SLHC25p1_SHCal_140PU_MET' : ['../python/METPerformance_SHCal_140PU_cfg_RECO.py',  '/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola/TP2023SHCALDR-SHCALMar26_PU140BX25_PH2_1K_FB_V6-v1/GEN-SIM-RECO',           2],
    #'SLHC23p1_SHCal_0PU_MET'   : ['../python/METPerformance_SHCal_NoPU_cfg_RECO.py',   '/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola/TP2023SHCALDR-SHCALJan23_NoPU_PH2_1K_FB_V6-v1/GEN-SIM-RECO',                2],
    #'SLHC23p1_Aged_140PU_MET'  : ['../python/METPerformance_Phase1_140PU_cfg_RECO.py', '/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola/GEM2019Upg14DR-Phase1age1kfixJan23_PU140BX25_PH1_1K_FB_V2-v1/GEN-SIM-RECO', 2],
    #'SLHC23p1_Aged_0PU_MET'    : ['../python/METPerformance_Phase1_140PU_cfg_RECO.py', '/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola/GEM2019Upg14DR-Phase1age1kfixJan23_NoPU_PH1_1K_FB_V2-v1/GEN-SIM-RECO',      2],
    #'SLHC23p1_NoAged_50PU_MET' : ['../python/METPerformance_Phase1_50PU_cfg_RECO.py',  '/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola/GEM2019Upg14DR-Phase1NoAgedJan23_PU50BX25_DES19_62_V8-v1/GEN-SIM-RECO',     2],
    #'SLHC23p1_NoAged_0PU_MET'  : ['../python/METPerformance_Phase1_50PU_cfg_RECO.py',  '/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola/GEM2019Upg14DR-Phase1NoAgedJan23_NoPU_DES19_62_V8-v1/GEN-SIM-RECO',         2, 2],

    #'SLHC23p1_HGCal_140PU_MET' : ['../METPerformance_HGCal_cfg_RECO.py', '/RelValZMM_14TeV/CMSSW_6_2_0_SLHC25_patch2-PU_PH2_1K_FB_V6_UPGHGCalPandPU140-v1/GEN-SIM-RECO', 1],
    #'SLHC23p1_HGCal_0PU_MET' : ['../METPerformance_HGCal_cfg_RECO.py', '/RelValZMM_14TeV/CMSSW_6_2_0_SLHC25_patch2-PH2_1K_FB_V6_UPGHGCalPand-v1/GEN-SIM-RECO', 1],

    #'SLHC25p2_HGCal_140PU_JET' : ['../python/PFAnalyzer_Phase2_cfg.py', '/QCD_Pt-15to3000_Tune4C_14TeV_pythia8/TP2023HGCALDR-HGCALMar26_PU140BX25_PH2_1K_FB_V6-v1/GEN-SIM-RECO', 3],
    #'SLHC25p2_HGCal_140PU_MET' : ['../python/METPerformance_SHCal_140PU_cfg_RECO.py', '/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola/TP2023HGCALDR-HGCALMar26_PU140BX25_PH2_1K_FB_V6-v2/GEN-SIM-RECO', 2],
    #'SLHC25p4_HGCal_140PU_MET' : ['../python/METPerformance_HGCal_140PU_cfg_RECO.py', '/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola/TP2023HGCALDR-HGCALnewsplit_PU140BX25_newsplitPU140_PH2_1K_FB_V6-v3/GEN-SIM-RECO', 3],

    #'SLHC25p6_HGCal_140PU_HF' : ['../python/PFAnalyzer_Phase2_cfg.py', '/RelValQCDForPF_14TeV/CMSSW_6_2_0_SLHC25_patch6-PU_PH2_1K_FB_V6_UPGHGCalPandPU140-v1/GEN-SIM-RECO', 1],
    #'SLHC25p4_HGCal_140PU_HF' : ['../python/PFAnalyzer_Phase2_cfg.py', '/RelValQCDForPF_14TeV/CMSSW_6_2_0_SLHC25_patch4-PU_PH2_1K_FB_V6_UPGHGCalPandPU140-v1/GEN-SIM-RECO', 1],
    #'SLHC25p6_HGCal_0PU_HF' : ['../python/PFAnalyzer_Phase2_cfg.py', '/RelValQCDForPF_14TeV/CMSSW_6_2_0_SLHC25_patch6-PH2_1K_FB_V6_UPGHGCalPand-v1/GEN-SIM-RECO', 1],
    #'SLHC25p4_HGCal_0PU_HF' : ['../python/PFAnalyzer_Phase2_cfg.py', '/RelValQCDForPF_14TeV/CMSSW_6_2_0_SLHC25_patch4-PH2_1K_FB_V6_UPGHGCalPand-v1/GEN-SIM-RECO', 1],

    #'TTbar'           : ['../treeMaker_stopRA2_PUPPI_Run2.py', '/TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v2/MINIAODSIM',                 10],
    #'WJetsHT100To200' : ['../treeMaker_stopRA2_PUPPI_Run2.py', '/WJetsToLNu_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/MINIAODSIM', 10],
    #'WJetsHT200To400' : ['../treeMaker_stopRA2_PUPPI_Run2.py', '/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/MINIAODSIM', 10],
    #'WJetsHT400To600' : ['../treeMaker_stopRA2_PUPPI_Run2.py', '/WJetsToLNu_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v3/MINIAODSIM', 10],
    #'WJetsHT600ToInf' : ['../treeMaker_stopRA2_PUPPI_Run2.py', '/WJetsToLNu_HT-600ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/MINIAODSIM', 10],

    #'T2tt425' : ['../../../python/prodFatTops_cfg.py', '/SMS-T2tt_2J_mStop-425_mLSP-325_Tune4C_13TeV-madgraph-tauola/Phys14DR-PU40bx25_PHYS14_25_V1-v1/MINIAODSIM',     5, 5],
    'TTbar'   : ['../../../python/prodFatTops_cfg.py', '/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/Phys14DR-PU20bx25_PHYS14_25_V1-v1/MINIAODSIM',         2],
    'QCD500'  : ['../../../python/prodFatTops_cfg.py', '/QCD_HT-500To1000_13TeV-madgraph/Phys14DR-PU20bx25_PHYS14_25_V1_ext1-v1/MINIAODSIM',                            2],
    'T2tt425' : ['../../../python/prodFatTops_cfg.py', '/SMS-T2tt_2J_mStop-425_mLSP-325_Tune4C_13TeV-madgraph-tauola/Phys14DR-PU40bx25_PHYS14_25_V1-v1/MINIAODSIM',     2],
    'T2tt500' : ['../../../python/prodFatTops_cfg.py', '/SMS-T2tt_2J_mStop-500_mLSP-325_Tune4C_13TeV-madgraph-tauola/Phys14DR-PU40bx25_PHYS14_25_V1-v1/MINIAODSIM',     2],
    'T2tt650' : ['../../../python/prodFatTops_cfg.py', '/SMS-T2tt_2J_mStop-650_mLSP-325_Tune4C_13TeV-madgraph-tauola/Phys14DR-PU20bx25_tsg_PHYS14_25_V1-v1/MINIAODSIM', 2],
    'T2tt850' : ['../../../python/prodFatTops_cfg.py', '/SMS-T2tt_2J_mStop-850_mLSP-100_Tune4C_13TeV-madgraph-tauola/Phys14DR-PU20bx25_tsg_PHYS14_25_V1-v1/MINIAODSIM', 2],
    'Zvv600'  : ['../../../python/prodFatTops_cfg.py', '/ZJetsToNuNu_HT-600toInf_Tune4C_13TeV-madgraph-tauola/Phys14DR-PU20bx25_PHYS14_25_V1-v1/MINIAODSIM',            2],
    'Zvv400'  : ['../../../python/prodFatTops_cfg.py', '/ZJetsToNuNu_HT-400to600_Tune4C_13TeV-madgraph-tauola/Phys14DR-PU20bx25_PHYS14_25_V1-v2/MINIAODSIM',            2],
    'Zvv200'  : ['../../../python/prodFatTops_cfg.py', '/ZJetsToNuNu_HT-200to400_Tune4C_13TeV-madgraph-tauola/Phys14DR-PU20bx25_PHYS14_25_V1-v1/MINIAODSIM',            2],
    'Zvv100'  : ['../../../python/prodFatTops_cfg.py', '/ZJetsToNuNu_HT-100to200_Tune4C_13TeV-madgraph-tauola/Phys14DR-PU20bx25_PHYS14_25_V1-v1/MINIAODSIM',            2],





}
tasklist = {}


for key, value in jobslist.items():
    tempconfig = copy.deepcopy(config)
    tempconfig.General.requestName = key
    tempconfig.General.workArea = workArea
    tempconfig.Data.publishDataName = key
    tempconfig.Data.outLFNDirBase = outDir
    if len(value) > 0:
        tempconfig.JobType.psetName = value[0]
    if len(value) > 1:
        tempconfig.Data.inputDataset = value[1]
    if len(value) > 2:
        tempconfig.Data.unitsPerJob = value[2]
    if len(value) > 3:
        tempconfig.Data.totalUnits = value[3]
    results = crabCommand('submit', config = tempconfig)
    tasklist[results['uniquerequestname']] = key
    del tempconfig

while True:
    for request, name in tasklist.items():
        dirname = './%s/crab_%s' % (workArea, name)
        fulldir = os.path.abspath(dirname)
        try:
            results = crabCommand('status', dir=fulldir)
            print "For task", request, "the job states are", results['jobsPerStatus']
        except:
            pass
        time.sleep(120)

